#ifndef _TRASERVAL_H_
#define _TRAVERSAL_H_

//#define _DEBUG_OUT


#include <set>
#include <map>
#include <list>
#include <ctime>
#include <iterator>
#include <cassert>
#include <numeric>
#include <sstream>
#include <iostream>

#include <tbb/parallel_for.h>

#include "path.h"

namespace Sibelia
{

	extern const std::string DELIMITER;

	class BlockInstance
	{
	public:
		BlockInstance() {}
		BlockInstance(int id, const size_t chr, size_t start, size_t end) : id_(id), chr_(chr), start_(start), end_(end) {}
		void Reverse();
		int GetSignedBlockId() const;
		bool GetDirection() const;
		int GetBlockId() const;
		int GetSign() const;
		size_t GetChrId() const;
		size_t GetStart() const;
		size_t GetEnd() const;
		size_t GetLength() const;
		size_t GetConventionalStart() const;
		size_t GetConventionalEnd() const;
		std::pair<size_t, size_t> CalculateOverlap(const BlockInstance & instance) const;
		bool operator < (const BlockInstance & toCompare) const;
		bool operator == (const BlockInstance & toCompare) const;
		bool operator != (const BlockInstance & toCompare) const;
	private:
		int id_;
		size_t start_;
		size_t end_;
		size_t chr_;
	};

	namespace
	{
		const bool COVERED = true;
		typedef std::vector<BlockInstance> BlockList;
		typedef std::pair<size_t, std::vector<BlockInstance> > GroupedBlock;
		typedef std::vector<GroupedBlock> GroupedBlockList;
		bool ByFirstElement(const GroupedBlock & a, const GroupedBlock & b)
		{
			return a.first < b.first;
		}

		std::string IntToStr(size_t x)
		{
			std::stringstream ss;
			ss << x;
			return ss.str();
		}

		template<class Iterator1, class Iterator2>
		void CopyN(Iterator1 it, size_t count, Iterator2 out)
		{
			for (size_t i = 0; i < count; i++)
			{
				*out++ = *it++;
			}
		}

		template<class Iterator>
		Iterator AdvanceForward(Iterator it, size_t step)
		{
			std::advance(it, step);
			return it;
		}

		template<class Iterator>
		Iterator AdvanceBackward(Iterator it, size_t step)
		{
			for (size_t i = 0; i < step; i++)
			{
				--it;
			}

			return it;
		}


		typedef std::pair<size_t, size_t> IndexPair;
		template<class T, class F, class It>
		void GroupBy(std::vector<T> & store, F pred, It out)
		{
			sort(store.begin(), store.end(), pred);
			for (size_t now = 0; now < store.size();)
			{
				size_t prev = now;
				for (; now < store.size() && !pred(store[prev], store[now]); now++);
				*out++ = std::make_pair(prev, now);
			}
		}

		template<class F>
		bool CompareBlocks(const BlockInstance & a, const BlockInstance & b, F f)
		{
			return (a.*f)() < (b.*f)();
		}

		template<class F>
		bool EqualBlocks(const BlockInstance & a, const BlockInstance & b, F f)
		{
			return f(a) == f(b);
		}

		template<class Iterator, class F, class ReturnType>
		struct FancyIterator : public std::iterator<std::forward_iterator_tag, ReturnType>
		{
		public:
			FancyIterator& operator++()
			{
				++it;
				return *this;
			}

			FancyIterator operator++(int)
			{
				FancyIterator ret(*this);
				++(*this);
				return ret;
			}

			bool operator == (FancyIterator toCompare) const
			{
				return it == toCompare.it;
			}

			bool operator != (FancyIterator toCompare) const
			{
				return !(*this == toCompare);
			}

			ReturnType operator * ()
			{
				return f(*it);
			}

			FancyIterator() {}
			FancyIterator(Iterator it, F f) : it(it), f(f) {}

		private:
			F f;
			Iterator it;
		};

		template<class Iterator, class F, class ReturnType>
		FancyIterator<Iterator, F, ReturnType> CFancyIterator(Iterator it, F f, ReturnType)
		{
			return FancyIterator<Iterator, F, ReturnType>(it, f);
		}

	}

	bool compareById(const BlockInstance & a, const BlockInstance & b);
	bool compareByChrId(const BlockInstance & a, const BlockInstance & b);
	bool compareByStart(const BlockInstance & a, const BlockInstance & b);

	void CreateOutDirectory(const std::string & path);

	class BlocksFinder
	{
	public:

		BlocksFinder(JunctionStorage & storage, size_t k) : storage_(storage), k_(k)
		{
			scoreFullChains_ = false;
		}

		struct ProcessVertexBruteForce
		{
		public:
			BlocksFinder & finder;
			std::vector<int64_t> & shuffle;

			ProcessVertexBruteForce(BlocksFinder & finder, std::vector<int64_t> & shuffle) : finder(finder), shuffle(shuffle)
			{
			}

			void operator()(tbb::blocked_range<size_t> & range) const
			{
				BestPath bestPath;
				Path currentPath(finder.storage_, finder.maxBranchSize_, finder.minBlockSize_, finder.flankingThreshold_);
				std::vector<std::vector<char > > mutexAcquired(finder.storage_.GetChrNumber(), std::vector<char>(finder.storage_.MutexNumber(), 0));
				for (size_t i = range.begin(); i != range.end(); i++)
				{
					if (finder.count_++ % 1000 == 0)
					{
						std::cerr << finder.count_ << '\t' << shuffle.size() << std::endl;
					}

					int64_t vid = shuffle[i];
					if (finder.visit_[vid + finder.storage_.GetVerticesNumber()])
					{
						continue;
					}

					for (bool explore = true; explore;)
					{
						bestPath.Init();
						currentPath.Init(vid);
						while (true)
						{
							int64_t prevBestScore = bestPath.score_;							
							finder.ExtendPathBackward(currentPath, bestPath, finder.lookingDepth_);
							bestPath.FixBackward(currentPath);
							finder.ExtendPathForward(currentPath, bestPath, finder.lookingDepth_);
							bestPath.FixForward(currentPath);
							if (bestPath.score_ <= prevBestScore)
							{
								break;
							}

						}

						if (!finder.TryFinalizeBlock(currentPath, mutexAcquired, std::cerr))
						{
							explore = false;
						}

						currentPath.Clear();
					}
				}
			}
		};		

		struct ProcessVertexSampling
		{
		public:
			BlocksFinder & finder;
			std::vector<int64_t> & shuffle;

			ProcessVertexSampling(BlocksFinder & finder, std::vector<int64_t> & shuffle) : finder(finder), shuffle(shuffle)
			{
			}

			void operator()(tbb::blocked_range<size_t> & range) const
			{	
				BestPath bestPath;
				Path currentPath(finder.storage_, finder.maxBranchSize_, finder.minBlockSize_, finder.flankingThreshold_);
				std::vector<std::vector<char > > mutexAcquired(finder.storage_.GetChrNumber(), std::vector<char>(finder.storage_.MutexNumber(), 0));				
				for (size_t i = range.begin(); i != range.end(); i++)
				{
					if (finder.count_++ % 1000 == 0)
					{
						std::cerr << finder.count_ << '\t' << shuffle.size() << std::endl;
					}

					int64_t vid = shuffle[i];
					if (finder.visit_[vid + finder.storage_.GetVerticesNumber()])
					{
						continue;
					}

					for (bool explore = true; explore;)
					{
						bestPath.Init();
						currentPath.Init(vid);
						while (true)
						{
							int64_t prevBestScore = bestPath.score_;
							finder.ExtendPathRandom(currentPath, bestPath, finder.lookingDepth_);
							if (bestPath.score_ <= prevBestScore)
							{
								break;
							}
						}

						if (!finder.TryFinalizeBlock(currentPath, mutexAcquired, std::cerr))
						{
							explore = false;
						}

						currentPath.Clear();
					}
				}
			}
		};

		void FindBlocks(int64_t minBlockSize, int64_t maxBranchSize, int64_t flankingThreshold, int64_t lookingDepth, int64_t sampleSize, int64_t threads, const std::string & debugOut)
		{
			blocksFound_ = 0;
			sampleSize_ = sampleSize;
			lookingDepth_ = lookingDepth;
			minBlockSize_ = minBlockSize;
			maxBranchSize_ = maxBranchSize;
			flankingThreshold_ = flankingThreshold;
			finishingProximity_ = 5 * k_;
			blockId_.resize(storage_.GetChrNumber());
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				blockId_[i].resize(storage_.GetChrVerticesCount(i));
			}

			std::vector<int64_t> shuffle;
			for (int64_t v = -storage_.GetVerticesNumber() + 1; v < storage_.GetVerticesNumber(); v++)
			{
				for (size_t i = 0; i < storage_.GetInstancesCount(v); i++)
				{
					if (storage_.GetJunctionInstance(v, i).IsPositiveStrand())
					{
						shuffle.push_back(v);
						break;
					}
				}
			}

			//std::random_shuffle(shuffle.begin(), shuffle.end());				
			time_t mark = time(0);
			count_ = 0;
			tbb::task_scheduler_init init(threads);
			visit_.reset(new std::atomic<bool>[storage_.GetVerticesNumber() * 2]);
			std::fill(visit_.get(), visit_.get() + storage_.GetVerticesNumber() * 2, false);
			if (sampleSize_ == 0)
			{
				tbb::parallel_for(tbb::blocked_range<size_t>(0, shuffle.size()), ProcessVertexBruteForce(*this, shuffle));
			}
			else
			{

//				tbb::serial::parallel_for(tbb::blocked_range<size_t>(0, shuffle.size()), ProcessVertexSampling(*this, shuffle));
				ProcessVertexSampling(*this, shuffle)(tbb::blocked_range<size_t>(0, shuffle.size()));
			}

			visit_.release();
			std::cerr << "Time: " << time(0) - mark << std::endl;
		}

		void Dump(std::ostream & out) const
		{
			out << "digraph G\n{\nrankdir = LR" << std::endl;
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				for (auto it = storage_.Begin(i); it != storage_.End(i) - 1; ++it)
				{
					auto jt = it + 1;
					out << it.GetVertexId(&storage_) << " -> " << jt.GetVertexId(&storage_)
						<< "[label=\"" << it.GetChar(&storage_) << ", " << it.GetChrId() << ", " << it.GetPosition(&storage_) << "\" color=blue]\n";
					out << jt.Reverse().GetVertexId(&storage_) << " -> " << it.Reverse().GetVertexId(&storage_)
						<< "[label=\"" << it.GetChar(&storage_) << ", " << it.GetChrId() << ", " << it.GetPosition(&storage_) << "\" color=red]\n";
				}
			}

			for (size_t i = 0; i < syntenyPath_.size(); i++)
			{
				for (size_t j = 0; j < syntenyPath_[i].size(); j++)
				{
					Edge e = syntenyPath_[i][j];
					out << e.GetStartVertex() << " -> " << e.GetEndVertex() <<
						"[label=\"" << e.GetChar() << ", " << i + 1 << "\" color=green]\n";
					e = e.Reverse();
					out << e.GetStartVertex() << " -> " << e.GetEndVertex() <<
						"[label=\"" << e.GetChar() << ", " << -(int64_t(i + 1)) << "\" color=green]\n";
				}
			}

			out << "}" << std::endl;
		}

		void DumpLight(std::ostream & out) const
		{
			out << "digraph G\n{\nrankdir = LR" << std::endl;
			std::vector<Edge> total;
			for (int64_t i = -storage_.GetVerticesNumber() + 1; i < storage_.GetVerticesNumber(); i++)
			{
				for (size_t j = 0; j < storage_.IngoingEdgesNumber(i); j++)
				{
					Edge e = storage_.IngoingEdge(i, j);
					total.push_back(e);
				}

				for (size_t j = 0; j < storage_.OutgoingEdgesNumber(i); j++)
				{
					Edge e = storage_.OutgoingEdge(i, j);
					total.push_back(e);
				}
			}

			std::sort(total.begin(), total.end());
			total.erase(std::unique(total.begin(), total.end()), total.end());

			for (Edge & e : total)
			{
				out << e.GetStartVertex() << " -> " << e.GetEndVertex() << "[label=\"" << e.GetChar() << ", " << e.GetCapacity() << "\"]" << std::endl;
			}

			out << "}" << std::endl;
		}

		void ListBlocksSequences(const BlockList & block, const std::string & fileName) const
		{
			std::ofstream out;
			TryOpenFile(fileName, out);
			std::vector<IndexPair> group;
			BlockList blockList = block;
			GroupBy(blockList, compareById, std::back_inserter(group));
			for (std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
			{
				for (size_t block = it->first; block < it->second; block++)
				{
					size_t length = blockList[block].GetLength();
					char strand = blockList[block].GetSignedBlockId() > 0 ? '+' : '-';
					size_t chr = blockList[block].GetChrId();
					out << ">Seq=\"" << storage_.GetChrDescription(chr) << "\",Strand='" << strand << "',";
					out << "Block_id=" << blockList[block].GetBlockId() << ",Start=";
					out << blockList[block].GetConventionalStart() << ",End=" << blockList[block].GetConventionalEnd() << std::endl;

					if (blockList[block].GetSignedBlockId() > 0)
					{
						OutputLines(storage_.GetChrSequence(chr).begin() + blockList[block].GetStart(), length, out);
					}
					else
					{
						std::string::const_reverse_iterator it(storage_.GetChrSequence(chr).begin() + blockList[block].GetEnd());
						OutputLines(CFancyIterator(it, TwoPaCo::DnaChar::ReverseChar, ' '), length, out);
					}

					out << std::endl;
				}
			}
		}


		void GenerateLegacyOutput(const std::string & outDir) const
		{
			BlockList instance;
			std::vector<std::vector<bool> > covered(storage_.GetChrNumber());
			for (size_t i = 0; i < covered.size(); i++)
			{
				covered[i].assign(storage_.GetChrSequence(i).size(), false);
			}

			for (size_t chr = 0; chr < blockId_.size(); chr++)
			{
				for (size_t i = 0; i < blockId_[chr].size();)
				{
					if (storage_.GetIterator(chr, i).IsUsed(&storage_))
					{
						int64_t bid = blockId_[chr][i].block;
						size_t j = i;
						for (; j < blockId_[chr].size() && blockId_[chr][i] == blockId_[chr][j]; j++);
						j--;
						int64_t cstart = storage_.GetIterator(chr, i, bid > 0).GetPosition(&storage_);
						int64_t cend = storage_.GetIterator(chr, j, bid > 0).GetPosition(&storage_) + (bid > 0 ? k_ : -k_);
						int64_t start = min(cstart, cend);
						int64_t end = max(cstart, cend);
						instance.push_back(BlockInstance(bid, chr, start, end));
						i = j + 1;
					}
					else
					{
						++i;
					}
				}
			}

			CreateOutDirectory(outDir);
			GenerateReport(instance, outDir + "/" + "coverage_report.txt");
			ListBlocksIndices(instance, outDir + "/" + "blocks_coords.txt");
			ListBlocksSequences(instance, outDir + "/" + "blocks_sequences.fasta");
		}


	private:

		template<class Iterator>
		void OutputLines(Iterator start, size_t length, std::ostream & out) const
		{
			for (size_t i = 1; i <= length; i++, ++start)
			{
				out << *start;
				if (i % 80 == 0 && i != length)
				{
					out << std::endl;
				}
			}
		}

		void GenerateReport(const BlockList & block, const std::string & fileName) const;
		std::vector<double> CalculateCoverage(GroupedBlockList::const_iterator start, GroupedBlockList::const_iterator end) const;
		std::string OutputIndex(const BlockInstance & block) const;
		void OutputBlocks(const std::vector<BlockInstance>& block, std::ofstream& out) const;
		void ListBlocksIndices(const BlockList & block, const std::string & fileName) const;
		void ListChromosomesAsPermutations(const BlockList & block, const std::string & fileName) const;
		void TryOpenFile(const std::string & fileName, std::ofstream & stream) const;
		void ListChrs(std::ostream & out) const;

		template<class P>
		bool TryFinalizeBlock(P & currentPath, std::vector<std::vector<char > > & mutexAcquired, std::ostream & log)
		{
			bool ret = false;
			if (currentPath.Score(true) > 0 && currentPath.MiddlePathLength() >= minBlockSize_ && currentPath.GoodInstances() > 1)
			{
				for (auto & instance : currentPath.Instances())
				{
					if (currentPath.IsGoodInstance(instance))
					{
						if (instance.Front().IsPositiveStrand())
						{
							storage_.LockRange(instance.Front(), instance.Back(), mutexAcquired);
						}
						else
						{
							storage_.LockRange(instance.Back().Reverse(), instance.Front().Reverse(), mutexAcquired);
						}
					}
				}

				std::vector<std::pair<JunctionStorage::JunctionIterator, JunctionStorage::JunctionIterator> > result;
				for (auto & instance : currentPath.Instances())
				{
					if (currentPath.IsGoodInstance(instance))
					{
						bool whole = true;
						auto start = instance.Front();
						auto end = instance.Back();
						for (; start != end && start.IsUsed(&storage_); ++start);
						for (; start != end && end.IsUsed(&storage_); --end);
						for (auto it = start; it != end + 1; it++)
						{
							if (it.IsUsed(&storage_))
							{
								whole = false;
								break;
							}
						}

						if (whole && abs(start.GetPosition(&storage_) - end.GetPosition(&storage_)) + k_ >= minBlockSize_)
						{
							result.push_back(std::make_pair(start, end));
						}
					}
				}


				if (result.size() > 1)
				{
					ret = true;
					int64_t currentBlock = ++blocksFound_;				
					int64_t instanceCount = 0;
					for (auto & instance : result)
					{
						auto end = instance.second;
						for (auto it = instance.first; it != end; ++it)
						{
							it.MarkUsed(&storage_);
							int64_t idx = it.GetIndex();
							int64_t maxidx = storage_.GetChrVerticesCount(it.GetChrId());
							blockId_[it.GetChrId()][it.GetIndex()].block = it.IsPositiveStrand() ? +currentBlock : -currentBlock;
							blockId_[it.GetChrId()][it.GetIndex()].instance = instanceCount;
						}
					}
				}

				for (auto & instance : currentPath.Instances())
				{
					if (currentPath.IsGoodInstance(instance))
					{
						if (instance.Front().IsPositiveStrand())
						{
							storage_.UnlockRange(instance.Front(), instance.Back(), mutexAcquired);
						}
						else
						{
							storage_.UnlockRange(instance.Back().Reverse(), instance.Front().Reverse(), mutexAcquired);
						}
					}
				}

				return ret;
			}

			return ret;
		}

		void ExtendPathRandom(Path & currentPath, BestPath & bestPath, int maxDepth)
		{
			int64_t startLength = currentPath.MiddlePathLength();
			for (size_t sample = 0; sample < sampleSize_; sample++)
			{
				for (size_t d = 0; ; d++)
				{
					int64_t dist = currentPath.RightDistance();
					if (dist < finishingProximity_)
					{
						visit_[currentPath.GetEndVertex() + storage_.GetVerticesNumber()] = true;
					}

					bool over = true;
					size_t attempts = storage_.OutgoingEdgesNumber(currentPath.GetEndVertex());
					for (size_t i = 0; i < attempts && (d < lookingDepth_ || currentPath.MiddlePathLength() - startLength < maxBranchSize_); i++)
					{
						Edge e = storage_.RandomForwardEdge(currentPath.GetEndVertex());
						if (e.Valid() && currentPath.PointPushBack(e))
						{
							over = false;
							break;
						}
					}

					if (over)
					{
						for (size_t i = 0; i < d; i++)
						{
							currentPath.PointPopBack();
						}

						break;
					}
					else
					{
						int64_t currentScore = currentPath.Score(scoreFullChains_);
						if (currentScore > bestPath.score_ && currentPath.Instances().size() > 1)
						{
							bestPath.UpdateForward(currentPath, currentScore);
						}
					}
				}
			}

			bestPath.FixForward(currentPath);
			startLength = currentPath.MiddlePathLength();
			for (size_t sample = 0; sample < sampleSize_; sample++)
			{
				for (size_t d = 0; ; d++)
				{
					int64_t dist = currentPath.LeftDistance();
					if (dist < finishingProximity_)
					{
						visit_[currentPath.GetStartVertex() + storage_.GetVerticesNumber()] = true;
					}

					bool over = true;
					size_t attempts = storage_.IngoingEdgesNumber(currentPath.GetStartVertex());
					for (size_t i = 0; i < attempts && (d < lookingDepth_ || currentPath.MiddlePathLength() - startLength < maxBranchSize_); i++)
					{
						Edge e = storage_.RandomBackwardEdge(currentPath.GetStartVertex());
						if (e.Valid() && currentPath.PointPushFront(e))
						{
							over = false;
							break;
						}
					}

					if (over)
					{
						for (size_t i = 0; i < d; i++)
						{
							currentPath.PointPopFront();
						}

						break;
					}
					else
					{
						int64_t currentScore = currentPath.Score(scoreFullChains_);
						if (currentScore > bestPath.score_ && currentPath.Instances().size() > 1)
						{
							bestPath.UpdateBackward(currentPath, currentScore);
						}
					}
				}
			}

			bestPath.FixBackward(currentPath);
		}				

		void ExtendPathForward(Path & currentPath, BestPath & bestPath, int maxDepth)
		{		
			if (maxDepth > 0)
			{
				int64_t prevVertex = currentPath.GetEndVertex();				
				if (currentPath.RightDistance() < finishingProximity_)
				{
					visit_[prevVertex + storage_.GetVerticesNumber()] = true;
				}

				for (int64_t idx = 0; idx < storage_.OutgoingEdgesNumber(prevVertex); idx++)
				{
					Edge e = storage_.OutgoingEdge(prevVertex, idx);
					{
						if (currentPath.PointPushBack(e))
						{
#ifdef _DEBUG_OUT
							currentPath.DebugOut(std::cerr);
#endif
							int64_t currentScore = currentPath.Score(scoreFullChains_);
							if (currentScore > bestPath.score_ && currentPath.Instances().size() > 1)
							{
								bestPath.UpdateForward(currentPath, currentScore);
							}

							ExtendPathForward(currentPath, bestPath, maxDepth - 1);
							currentPath.PointPopBack();
						}
					}
				}
			}
		}

		void ExtendPathBackward(Path & currentPath, BestPath & bestPath, int maxDepth)
		{	
			if (maxDepth > 0)
			{
				int64_t prevVertex = currentPath.GetStartVertex();
				if (currentPath.LeftDistance() < finishingProximity_)
				{
					visit_[prevVertex + storage_.GetVerticesNumber()] = true;
				}

				for (int64_t idx = 0; idx < storage_.IngoingEdgesNumber(prevVertex); idx++)
				{
					Edge e = storage_.IngoingEdge(prevVertex, idx);
					{
						if (currentPath.PointPushFront(e))
						{
#ifdef _DEBUG_OUT
							currentPath.DebugOut(std::cerr);
#endif
							int64_t currentScore = currentPath.Score(scoreFullChains_);
							if (currentScore > bestPath.score_ && currentPath.Instances().size() > 1)
							{
								bestPath.UpdateBackward(currentPath, currentScore);
							}

							ExtendPathBackward(currentPath, bestPath, maxDepth - 1);
							currentPath.PointPopFront();
						}
					}
				}
			}
		}

		int64_t k_;
		std::atomic<int64_t> count_;
		std::atomic<int64_t> blocksFound_;
		int64_t sampleSize_;
		int64_t scalingFactor_;
		bool scoreFullChains_;
		int64_t lookingDepth_;
		int64_t minBlockSize_;
		int64_t maxBranchSize_;
		int64_t finishingProximity_;
		int64_t flankingThreshold_;
		JunctionStorage & storage_;
		std::unique_ptr<std::atomic<bool>[] > visit_;
		std::vector<std::vector<Edge> > syntenyPath_;
		std::vector<std::vector<Assignment> > blockId_;
	};
}

#endif