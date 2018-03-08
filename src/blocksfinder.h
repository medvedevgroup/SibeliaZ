#ifndef _TRASERVAL_H_
#define _TRAVERSAL_H_

//#define _DEBUG_OUT_

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
			scoreFullChains_ = true;
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
				std::vector<uint32_t> data;
				std::vector<uint32_t> count(finder.storage_.GetVerticesNumber() * 2 + 1, 0);
				std::pair<int64_t, std::vector<Path::Instance> > goodInstance;
				Path currentPath(finder.storage_, finder.maxBranchSize_, finder.minBlockSize_, finder.k_);
				for (size_t i = range.begin(); i != range.end(); i++)
				{
					if (finder.count_++ % 1000 == 0)
					{
						std::cout << finder.count_ << '\t' << shuffle.size() << std::endl;
					}
					
					int64_t score;
					int64_t vid = shuffle[i];
#ifdef _DEBUG_OUT_
					std::cerr << "Vid: " << vid << std::endl;
#endif
					for (bool explore = true; explore;)
					{
						currentPath.Init(vid);
						int64_t bestScore = 0;
						size_t bestRightSize = currentPath.RightSize();
						size_t bestLeftSize = currentPath.LeftSize();
#ifdef _DEBUG_OUT_
						std::cerr << "Going forward:" << std::endl;
#endif
						int64_t minRun = max(finder.minBlockSize_, finder.maxBranchSize_) * 3;
						while (true)
						{
							int64_t prevBestScore = currentPath.Score(finder.scoreFullChains_);
							bool ret = finder.ExtendPathDijkstraForward(currentPath, count, data, bestRightSize, bestScore, score);
							if (!ret || (score <= 0 && currentPath.MiddlePathLength() >= minRun))
							{
								break;
							}
						}
						
						while (currentPath.RightSize() > bestRightSize)
						{
							currentPath.PointPopBack();
						}

#ifdef _DEBUG_OUT_
						std::cerr << "Going backward:" << std::endl;
#endif
						while (true)
						{
							int64_t prevBestScore = currentPath.Score(finder.scoreFullChains_);
							bool ret = finder.ExtendPathDijkstraBackward(currentPath, count, data, bestLeftSize, bestScore, score);
							if (!ret || (score <= 0 && currentPath.MiddlePathLength() >= minRun))
							{
								break;
							}
						}

						while (currentPath.LeftSize() > bestLeftSize)
						{
							currentPath.PointPopFront();
						}
						
						if (bestScore > 0)
						{							
							goodInstance.second.clear();
							goodInstance.first = bestScore;
							for (auto it : currentPath.AllInstances())
							{
								if (currentPath.IsGoodInstance(*it))
								{
									goodInstance.second.push_back(*it);
								}
							}

							if (!finder.TryFinalizeBlock(currentPath, goodInstance, std::cerr))
							{
								explore = false;
							}
						}
						else
						{
							explore = false;
						}

						currentPath.Clear();
					}
				}
			}
		};


		void FindBlocks(int64_t minBlockSize, int64_t maxBranchSize, int64_t lookingDepth, int64_t sampleSize, int64_t threads, const std::string & debugOut)
		{
			blocksFound_ = 0;
			sampleSize_ = sampleSize;
			lookingDepth_ = lookingDepth;
			minBlockSize_ = minBlockSize;
			maxBranchSize_ = maxBranchSize;		
			blockId_.resize(storage_.GetChrNumber());
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				blockId_[i].resize(storage_.GetChrVerticesCount(i));
			}

			std::vector<int64_t> shuffle;			
			for (int64_t v = -storage_.GetVerticesNumber() + 1; v < storage_.GetVerticesNumber(); v++)
			{
				for (JunctionStorage::JunctionIterator it(v); it.Valid(); ++it)
				{
					if(it.IsPositiveStrand())
					{
						shuffle.push_back(v);
						break;
					}					
				}
			}
			
			srand(time(0));
			//std::random_shuffle(shuffle.begin(), shuffle.end());			
			time_t mark = time(0);
			count_ = 0;
			tbb::task_scheduler_init init(threads);
			tbb::parallel_for(tbb::blocked_range<size_t>(0, shuffle.size()), ProcessVertexDijkstra(*this, shuffle));
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
			ListBlocksIndices(instance, outDir + "/" + "blocks_coords.txt");
			ListBlocksSequences(instance, outDir + "/" + "blocks_sequences.fasta");
			GenerateReport(instance, outDir + "/" + "coverage_report.txt");

			std::string buf;
			std::set<int64_t> vertex;
			std::ifstream oldCoordsIn(outDir + "/missing.maf");
			if (oldCoordsIn)
			{
				std::ofstream missingDot(outDir + "/missing.dot");
				while (std::getline(oldCoordsIn, buf) && buf[0] != '-')
				{
					std::string seq;
					std::stringstream ss(buf);
					char sign;
					int seqId, start, length, end, seqSize;
					ss >> seq >> seq >> start >> length >> sign >> seqSize;
					seqId = atoi(seq.substr(2).c_str()) - 1;
					end = start + length;
					if (sign == '-')
					{
						start = seqSize - start;
						end = seqSize - end;
						std::swap(start, end);
						assert(start < end);
					}

					for (auto it = storage_.Begin(seqId); it.Valid(); ++it)
					{
						int64_t pos = it.GetPosition();
						if (pos >= start && pos < end)
						{
							vertex.insert(it.GetVertexId());
						}
					}

				}

				missingDot << "digraph G\n{\nrankdir = LR" << std::endl;
				std::vector<std::pair<JunctionStorage::JunctionSequentialIterator, JunctionStorage::JunctionSequentialIterator> > vvisit;
				for (auto vid : vertex)
				{
					DumpVertex(vid, missingDot, vvisit);
					missingDot << vid << "[shape=square]" << std::endl;
				}

				missingDot << "}" << std::endl;
			}			
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
		void CalculateCoverage(GroupedBlockList::const_iterator start, GroupedBlockList::const_iterator end, std::vector<bool> & cover, std::vector<double> & ret) const;
		std::string OutputIndex(const BlockInstance & block) const;
		void OutputBlocks(const std::vector<BlockInstance>& block, std::ofstream& out) const;
		void ListBlocksIndices(const BlockList & block, const std::string & fileName) const;
		void ListChromosomesAsPermutations(const BlockList & block, const std::string & fileName) const;
		void TryOpenFile(const std::string & fileName, std::ofstream & stream) const;
		void ListChrs(std::ostream & out) const;
		
		template<class T>
		void DumpVertex(int64_t id, std::ostream & out, T & visit, int64_t cnt = 5) const
		{
			for (auto kt = JunctionStorage::JunctionIterator(id); kt.Valid(); ++kt)
			{
				auto jt = kt.SequentialIterator();
				for (int64_t i = 0; i < cnt; i++)
				{
					auto it = jt - 1;
					auto pr = std::make_pair(it, jt);
					if (it.Valid() && std::find(visit.begin(), visit.end(), pr) == visit.end())
					{
						int64_t length = it.GetPosition() - jt.GetPosition();
						out << it.GetVertexId() << " -> " << jt.GetVertexId()
							<< "[label=\"" << it.GetChar() << ", " << it.GetChrId() << ", " << it.GetPosition() << "," << length << "\""
							<< (it.IsPositiveStrand() ? "color=blue" : "color=red") << "]\n";
						visit.push_back(pr);
					}

					jt = it;
				}
			}

			for (auto kt = JunctionStorage::JunctionIterator(id); kt.Valid(); ++kt)
			{
				auto it = kt.SequentialIterator();
				for (int64_t i = 0; i < cnt; i++)
				{
					auto jt = it + 1;
					auto pr = std::make_pair(it, jt);
					if (jt.Valid() && std::find(visit.begin(), visit.end(), pr) == visit.end())
					{
						int64_t length = it.GetPosition() - jt.GetPosition();
						out << it.GetVertexId() << " -> " << jt.GetVertexId()
							<< "[label=\"" << it.GetChar() << ", " << it.GetChrId() << ", " << it.GetPosition() << "," << length << "\""
							<< (it.IsPositiveStrand() ? "color=blue" : "color=red") << "]\n";
						visit.push_back(pr);
					}

					it = jt;
				}
			}
		}


		template<class P>
		bool TryFinalizeBlock(P & currentPath, std::pair<int64_t, std::vector<Path::Instance> > & goodInstance, std::ostream & log)
		{
			bool ret = false;
			if (goodInstance.first > 0 && goodInstance.second.size() > 1)
			{
				std::sort(goodInstance.second.begin(), goodInstance.second.end(), Path::Instance::OldComparator);
				{
					std::pair<size_t, size_t> idx(SIZE_MAX, SIZE_MAX);
					for (auto & instIt : goodInstance.second)
					{
						auto & instance = instIt;
						{
							if (instance.Front().IsPositiveStrand())
							{
								storage_.LockRange(instance.Front(), instance.Back(), idx);
							}
							else
							{
								storage_.LockRange(instance.Back().Reverse(), instance.Front().Reverse(), idx);
							}
						}						
					}
				}

				std::vector<std::pair<JunctionStorage::JunctionSequentialIterator, JunctionStorage::JunctionSequentialIterator> > result;
				for (auto & instIt : goodInstance.second)
				{
					auto & instance = instIt;
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
					int64_t instanceCount = 0;
					int64_t currentBlock = ++blocksFound_;
					for (auto & instance : result)
					{
						auto it = instance.first;
						do
						{
							it.MarkUsed();
							int64_t idx = it.GetIndex();
							int64_t maxidx = storage_.GetChrVerticesCount(it.GetChrId());
							blockId_[it.GetChrId()][it.GetIndex()].block = it.IsPositiveStrand() ? +currentBlock : -currentBlock;
							blockId_[it.GetChrId()][it.GetIndex()].instance = instanceCount;

						} while (it++ != instance.second);	

						instanceCount++;
					}
				}

				{
					std::pair<size_t, size_t> idx(SIZE_MAX, SIZE_MAX);
					for (auto & instIt : goodInstance.second)
					{
						auto & instance = instIt;
						{
							if (instance.Front().IsPositiveStrand())
							{
								storage_.UnlockRange(instance.Front(), instance.Back(), idx);
							}
							else
							{
								storage_.UnlockRange(instance.Back().Reverse(), instance.Front().Reverse(), idx);
							}
						}
					}
				} 
				
			}

			return ret;
		}

		struct NextVertex
		{
			int32_t diff;
			int32_t count;
			JunctionStorage::JunctionSequentialIterator origin;
			NextVertex(): count(0)
			{

			}

			NextVertex(int64_t diff, JunctionStorage::JunctionSequentialIterator origin) : origin(origin), diff(diff), count(1)
			{

			}
		};


		std::pair<int32_t, NextVertex> MostPopularVertex(const Path & currentPath, bool forward, std::vector<uint32_t> & count, std::vector<uint32_t> & data)
		{
			NextVertex ret;
			int32_t bestVid = 0;
			int64_t startVid = forward ? currentPath.RightVertex() : currentPath.LeftVertex();
			for (auto & inst : currentPath.AllInstances())
			{
				int64_t nowVid = forward ?  inst->Back().GetVertexId() : inst->Front().GetVertexId();
				if (nowVid == startVid)
				{
					int64_t weight = abs(inst->Front().GetPosition() - inst->Back().GetPosition()) + 1;
					auto origin = forward ? inst->Back() : inst->Front();
					auto it = forward ? origin.Next() : origin.Prev();
					for (size_t d = 1; it.Valid() && (d < lookingDepth_); d++)
					{
						int32_t vid = it.GetVertexId();
						if (!currentPath.IsInPath(vid) && !it.IsUsed())
						{
							auto adjVid = vid + storage_.GetVerticesNumber();
							if (count[adjVid] == 0)
							{
								data.push_back(adjVid);
							}

							count[adjVid] += weight;
							auto diff = abs(it.GetAbsolutePosition() - origin.GetAbsolutePosition());
							if (count[adjVid] > ret.count || (count[adjVid] == ret.count && diff < ret.diff))
							{
								ret.diff = diff;
								ret.origin = origin;
								ret.count = count[adjVid];
								bestVid = vid;
							}
						}
						else
						{
							break;
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
			}
			

			for (auto vid : data)
			{
				count[vid] = 0;
			}

			data.clear();
			return std::make_pair(bestVid, ret);
		}				

		std::pair<int32_t, NextVertex> MostPopularVertexSampling(const Path & currentPath, bool forward, std::vector<uint32_t> & count, std::vector<uint32_t> & data)
		{
			NextVertex ret;
			int32_t bestVid = 0;
			return std::make_pair(bestVid, ret);
		}

		bool ExtendPathDijkstraForward(Path & currentPath,
			std::vector<uint32_t> & count,
			std::vector<uint32_t> & data,
			size_t & bestRightSize,
			int64_t & bestScore,
			int64_t & nowScore)
		{	
			bool success = false;
			int64_t origin = currentPath.Origin();
			std::pair<int32_t, NextVertex> nextForwardVid;
			if (sampleSize_ == 0 || storage_.GetInstancesCount(currentPath.RightVertex()) <= sampleSize_)
			{
				nextForwardVid = MostPopularVertex(currentPath, true, count, data);
			}
			else
			{
				nextForwardVid = MostPopularVertexSampling(currentPath, true, count, data);
			}
			
			if (nextForwardVid.first != 0)
			{				
				for(auto it = nextForwardVid.second.origin; it.GetVertexId() != nextForwardVid.first; ++it)
				{
#ifdef _DEBUG_OUT_
					std::cerr << "Attempting to push back the vertex:" << it.GetVertexId() << std::endl;
#endif
					success = currentPath.PointPushBack(it.OutgoingEdge());
					if (success)
					{
						nowScore = currentPath.Score(scoreFullChains_);
#ifdef _DEBUG_OUT_
						std::cerr << "Success! New score:" << nowScore << std::endl;
						currentPath.DumpPath(std::cerr);
						currentPath.DumpInstances(std::cerr);
#endif												
						if (nowScore > bestScore)
						{
							bestScore = nowScore;
							bestRightSize = currentPath.RightSize();
						}
					}					
				}
			}

			return success;
		}

		bool ExtendPathDijkstraBackward(Path & currentPath,
			std::vector<uint32_t> & count,
			std::vector<uint32_t> & data,
			size_t & bestLeftSize,
			int64_t & bestScore,
			int64_t & nowScore)
		{
			bool success = false;
			std::pair<int32_t, NextVertex> nextBackwardVid;
			if (sampleSize_ == 0 || storage_.GetInstancesCount(currentPath.RightVertex()) <= sampleSize_)
			{
				nextBackwardVid = MostPopularVertex(currentPath, false, count, data);
			}
			else
			{
				nextBackwardVid = MostPopularVertexSampling(currentPath, false, count, data);
			}

			if (nextBackwardVid.first != 0)
			{
				for (auto it = nextBackwardVid.second.origin; it.GetVertexId() != nextBackwardVid.first; --it)
				{
#ifdef _DEBUG_OUT_
					std::cerr << "Attempting to push front the vertex:" << it.GetVertexId() << std::endl;
#endif
					success = currentPath.PointPushFront(it.IngoingEdge());
					if (success)
					{
						nowScore = currentPath.Score(scoreFullChains_);
#ifdef _DEBUG_OUT_
						std::cerr << "Success! New score:" << nowScore << std::endl;
						currentPath.DumpPath(std::cerr);
						currentPath.DumpInstances(std::cerr);
#endif		
						if (nowScore > bestScore)
						{
							bestScore = nowScore;
							bestLeftSize = currentPath.LeftSize();
						}
					}
				}
			}
						
			return success;
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
		JunctionStorage & storage_;
		std::vector<std::vector<Edge> > syntenyPath_;
		std::vector<std::vector<Assignment> > blockId_;
	};
}

#endif