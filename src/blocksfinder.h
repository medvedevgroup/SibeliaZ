#ifndef _TRASERVAL_H_
#define _TRAVERSAL_H_

//#define _DEBUG_OUT


#include <set>
#include <map>
#include <list>
#include <ctime>
#include <queue>
#include <iterator>
#include <cassert>
#include <numeric>
#include <sstream>
#include <iostream>
#include <unordered_map>

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

		struct ProcessVertexDijkstra
		{
		public:
			BlocksFinder & finder;
			std::vector<int64_t> & shuffle;

			ProcessVertexDijkstra(BlocksFinder & finder, std::vector<int64_t> & shuffle) : finder(finder), shuffle(shuffle)
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
						std::cout << finder.count_ << '\t' << shuffle.size() << std::endl;
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
							finder.ExtendPathDijkstra(currentPath, bestPath, finder.lookingDepth_);
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
			tbb::parallel_for(tbb::blocked_range<size_t>(0, shuffle.size()), ProcessVertexDijkstra(*this, shuffle));
			visit_.release();
			std::cout << "Time: " << time(0) - mark << std::endl;
		}

		void Dump(std::ostream & out) const
		{
			out << "digraph G\n{\nrankdir = LR" << std::endl;
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				auto end = storage_.End(i).Prev();
				for (auto it = storage_.Begin(i); it != end; ++it)
				{
					auto jt = it.Next();
					out << it.GetVertexId() << " -> " << jt.GetVertexId() << "[label=\"" << it.GetChar() << ", " << it.GetChrId() << ", " << it.GetPosition() << "\" color=blue]\n";
					out << jt.Reverse().GetVertexId() << " -> " << it.Reverse().GetVertexId() << "[label=\"" << it.GetChar() << ", " << it.GetChrId() << ", " << it.GetPosition() << "\" color=red]\n";
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
					if (storage_.GetIterator(chr, i).IsUsed())
					{
						int64_t bid = blockId_[chr][i].block;
						size_t j = i;
						for (; j < blockId_[chr].size() && blockId_[chr][i] == blockId_[chr][j]; j++);
						j--;
						int64_t cstart = storage_.GetIterator(chr, i, bid > 0).GetPosition();
						int64_t cend = storage_.GetIterator(chr, j, bid > 0).GetPosition() + (bid > 0 ? k_ : -k_);
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
						for (; start != end && start.IsUsed(); ++start);
						for (; start != end && end.IsUsed(); --end);
						for (auto it = start; it != end.Next(); it++)
						{
							if (it.IsUsed())
							{
								whole = false;
								break;
							}
						}

						if (whole && abs(start.GetPosition() - end.GetPosition()) + k_ >= minBlockSize_)
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
							it.MarkUsed();
							int64_t idx = it.GetIndex();
							int64_t maxidx = storage_.GetChrVerticesCount(it.GetChrId());
							blockId_[it.GetChrId()][it.GetIndex()].block = it.IsPositiveStrand() ? +currentBlock : -currentBlock;
							blockId_[it.GetChrId()][it.GetIndex()].instance = instanceCount;
						}
					}
				}
				
				/*
				std::vector<Edge> epath;
				currentPath.DumpPath(epath);				
				std::cerr << "Winner: " << std::endl;
				for (auto e : epath)
				{
					std::cerr << "(" << e.GetEndVertex() << ", " << storage_.GetInstancesCount(e.GetEndVertex()) << ") ";
				}

				std::cerr << std::endl;
				*/

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
			}

			return ret;
		}

		struct DijkstraState
		{
			int64_t prevIdx;
			int64_t vertex;
			int64_t distance;
			Edge e;
			DijkstraState(int64_t vertex, int64_t distance, int64_t prevIdx, Edge e) : vertex(vertex), distance(distance), prevIdx(prevIdx), e(e)
			{
			}

			bool operator < (const DijkstraState & state) const
			{
				return distance > state.distance;
			}
		};

		void Dijkstra(int64_t source, int64_t sink, std::vector<Edge> & path, bool forward)
		{
			path.clear();
			std::priority_queue<DijkstraState> q;
			std::unordered_map<int64_t, int64_t> vertexIdx;
			std::vector<int64_t> distance;
			std::vector<int64_t> prevIdx;
			std::vector<Edge> nextEdge;
			int64_t prevVertex = source;
			q.push(DijkstraState(prevVertex, 0, -1, Edge()));

			while (q.size() > 0)
			{
				DijkstraState state = q.top();
				q.pop();
				int64_t nowIdx;
				auto it = vertexIdx.find(state.vertex);
				if (it == vertexIdx.end())
				{
					nowIdx = vertexIdx.size();
					vertexIdx[state.vertex] = nowIdx;
					distance.push_back(INT64_MAX);
					nextEdge.push_back(Edge());
					prevIdx.push_back(-1);
				}
				else
				{
					nowIdx = it->second;
				}

				if (state.distance < distance[nowIdx])
				{
					prevIdx[nowIdx] = state.prevIdx;
					distance[nowIdx] = state.distance;
					if (state.prevIdx >= 0)
					{
						nextEdge[state.prevIdx] = state.e;
						if (state.vertex == sink)
						{
							break;
						}
					}
					if (forward)
					{
						for (int64_t idx = 0; idx < storage_.OutgoingEdgesNumber(state.vertex); idx++)
						{
							Edge e = storage_.OutgoingEdge(state.vertex, idx);
							int64_t nextVertex = e.GetEndVertex();
							DijkstraState next(nextVertex, state.vertex == prevVertex ? 0 : state.distance + e.GetLength(), nowIdx, e);
							if (next.distance <= maxBranchSize_)
							{
								q.push(next);
							}
						}
					}
					else
					{
						for (int64_t idx = 0; idx < storage_.IngoingEdgesNumber(state.vertex); idx++)
						{
							Edge e = storage_.IngoingEdge(state.vertex, idx);
							int64_t nextVertex = e.GetStartVertex();
							DijkstraState next(nextVertex, state.vertex == prevVertex ? 0 : state.distance + e.GetLength(), nowIdx, e);
							if (next.distance <= maxBranchSize_)
							{
								q.push(next);
							}
						}
					}
				}
			}

			for (int64_t idx = prevIdx.back(); idx >= 0; idx = prevIdx[idx])
			{
				path.push_back(nextEdge[idx]);
			}
		}

		int64_t MostPopularVertex(const Path & currentPath, bool forward)
		{
			std::vector<int64_t> popularVid;
			for (auto inst : currentPath.Instances())
			{
				auto it = forward ? inst.Back().Next() : inst.Front().Prev();
				for (size_t d = 0; it.Valid() && (d < lookingDepth_ || abs(it.GetPosition() - inst.Back().GetPosition()) < maxBranchSize_); d++)
				{
					if (!currentPath.IsInPath(it.GetVertexId()))
					{
						popularVid.push_back(it.GetVertexId());
					}

					if (forward)
					{
						++it;
					}
					else
					{
						--it;
					}
				}
			}

			int64_t bestVid = 0;
			int64_t bestCount = 0;
			std::sort(popularVid.begin(), popularVid.end());
			for (size_t i = 0; i < popularVid.size(); )
			{
				size_t j = i + 1;
				for (; j < popularVid.size() && popularVid[j] == popularVid[i]; ++j);
				if (j - i > bestCount)
				{
					bestCount = j - i;
					bestVid = popularVid[i];
				}

				i = j;
			}

			return bestVid;
		}


		void ExtendPathDijkstra(Path & currentPath, BestPath & bestPath, int maxDepth)
		{			
			int64_t nextForwardVid = MostPopularVertex(currentPath, true);
			if (nextForwardVid != 0)
			{
				std::vector<Edge> newPath;
				Dijkstra(currentPath.GetEndVertex(), nextForwardVid, newPath, true);
				for (auto it = newPath.rbegin(); it != newPath.rend(); ++it)
				{
					if (currentPath.PointPushBack(*it))
					{
						int64_t currentScore = currentPath.Score(scoreFullChains_);
						if (currentScore > bestPath.score_ && currentPath.Instances().size() > 1)
						{
							bestPath.UpdateForward(currentPath, currentScore);
						}
					}
					else
					{
						break;
					}
				}
			}

			int64_t nextBackwardVid = MostPopularVertex(currentPath, false);
			if (nextBackwardVid != 0)
			{
				std::vector<Edge> newPath;
				Dijkstra(currentPath.GetStartVertex(), nextForwardVid, newPath, false);
				for (auto it = newPath.rbegin(); it != newPath.rend(); ++it)
				{
					if (currentPath.PointPushFront(*it))
					{
						int64_t currentScore = currentPath.Score(scoreFullChains_);
						if (currentScore > bestPath.score_ && currentPath.Instances().size() > 1)
						{
							bestPath.UpdateBackward(currentPath, currentScore);
						}
					}
					else
					{
						break;
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