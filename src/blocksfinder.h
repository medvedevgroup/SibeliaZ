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
#include <functional>
#include <unordered_map>


#include <tbb/mutex.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>

#include "path.h"

namespace Sibelia
{
	extern const std::string DELIMITER;
	extern const std::string VERSION;

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

		struct ProcessVertexIsFork
		{
		public:
			BlocksFinder & finder;
			std::vector<int64_t> & shuffle;

			ProcessVertexIsFork(BlocksFinder & finder, std::vector<int64_t> & shuffle) : finder(finder), shuffle(shuffle)
			{
			}

			void operator()(tbb::blocked_range<size_t> & range) const
			{
				BubbledBranches forwardBubble;
				BubbledBranches backwardBubble;
				std::vector<JunctionStorage::JunctionSequentialIterator> instance;				
				for (size_t i = range.begin(); i != range.end(); i++)
				{
					if (finder.count_++ % 1000 == 0)
					{
						std::cout << finder.count_ << '\t' << shuffle.size() << std::endl;
					}

					instance.clear();
					int64_t minStart = INT64_MAX;
					int64_t vertex = shuffle[i];
					for (auto it = JunctionStorage::JunctionIterator(vertex); it.Valid(); ++it)
					{
						instance.push_back(it.SequentialIterator());
						auto pos = it.GetPosition();
						if (pos < minStart)
						{
							minStart = pos;
						}
					}

					finder.BubbledBranchesForward(vertex, instance, forwardBubble);
					finder.BubbledBranchesBackward(vertex, instance, backwardBubble);
					bool isFork = false;
					for (size_t i = 0; i < forwardBubble.size() && !isFork; i++)
					{
						for (size_t j = 0; j < forwardBubble[i].size() && !isFork; j++)
						{
							size_t k = forwardBubble[i][j];
							if (std::find(backwardBubble[i].begin(), backwardBubble[i].end(), k) == backwardBubble[i].end())
							{
								finder.source_.push_back(std::make_pair(minStart, vertex));
								isFork = true;
							}
						}
					}
				}
			}
		};

		struct ProcessVertexDijkstra
		{
		public:
			BlocksFinder & finder;

			ProcessVertexDijkstra(BlocksFinder & finder) : finder(finder)
			{
			}

			void operator()(tbb::blocked_range<size_t> & range) const
			{
				std::vector<uint32_t> data;
				std::pair<int64_t, std::vector<Path::Instance> > goodInstance;
				Path currentPath(finder.storage_, finder.maxBranchSize_, finder.minBlockSize_, finder.minBlockSize_, finder.maxFlankingSize_);
				for (size_t i = range.begin(); i != range.end(); i++)
				{
					if (finder.count_++ % 1000 == 0)
					{
						std::cout << finder.count_ << '\t' << finder.source_.size() << std::endl;
					}
					
					int64_t score;
					int64_t vid = finder.source_[i].second;
#ifdef _DEBUG_OUT_
					finder.debug_ = finder.missingVertex_.count(vid);
					if (finder.debug_)
					{
						std::cerr << "Vid: " << vid << std::endl;
					}
#endif
					for (bool explore = true; explore;)
					{
						currentPath.Init(vid);
						int64_t bestScore = 0;
						size_t bestRightSize = currentPath.RightSize();
#ifdef _DEBUG_OUT_
						if (finder.debug_)
						{
							std::cerr << "Going forward:" << std::endl;
						}
#endif
						goodInstance.first = 0;
						goodInstance.second.clear();
						int64_t minRun = max(finder.minBlockSize_, finder.maxBranchSize_) * 2;
						while (true)
						{
							bool ret = true;
							bool positive = false;
							int64_t prevLength = currentPath.MiddlePathLength();
							while ((ret = finder.ExtendPathDijkstraForward(currentPath, bestRightSize, bestScore, score, goodInstance)) && currentPath.MiddlePathLength() - prevLength <= minRun)
							{
								positive = positive || (score > 0);
							}

							if (!ret || !positive)
							{
								break;
							}
						}
						
						if (bestScore > 0)
						{			

#ifdef _DEBUG_OUT_
							if (finder.debug_)
							{
								std::cerr << "Setting a new block. Best score:" << bestScore << std::endl;
								currentPath.DumpPath(std::cerr);
								currentPath.DumpInstances(std::cerr);
							}
#endif
				
							if (!finder.TryFinalizeBlock(currentPath, goodInstance, std::cerr))
							{
								explore = false;
							}
							else
							{
								int x = 1;
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

		void MissingSet(const std::string & fileName, std::set<int64_t> & result) const
		{
			std::string buf;
			std::ifstream oldCoordsIn(fileName);
			if (oldCoordsIn)
			{
				std::ofstream missingDot("missing.dot");
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
							result.insert(it.GetVertexId());
							result.insert(-it.GetVertexId());
						}
					}

				}
			}
		}

		static bool DegreeCompare(const JunctionStorage & storage, int64_t v1, int64_t v2)
		{
			return storage.GetInstancesCount(v1) > storage.GetInstancesCount(v2);
		}

		void FindBlocks(int64_t minBlockSize, int64_t maxBranchSize, int64_t maxFlankingSize, int64_t lookingDepth, int64_t sampleSize, int64_t threads, const std::string & debugOut)
		{
			blocksFound_ = 0;
			sampleSize_ = sampleSize;
			lookingDepth_ = lookingDepth;
			minBlockSize_ = minBlockSize;
			maxBranchSize_ = maxBranchSize;	
			maxFlankingSize_ = maxFlankingSize;
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

			using namespace std::placeholders;
			std::sort(shuffle.begin(), shuffle.end(), std::bind(DegreeCompare, std::cref(storage_), _1, _2));
			srand(time(0));
			time_t mark = time(0);
			count_ = 0;
			tbb::task_scheduler_init init(threads);
			tbb::parallel_for(tbb::blocked_range<size_t>(0, shuffle.size()), ProcessVertexIsFork(*this, shuffle));
			count_ = 0;			
			tbb::parallel_for(tbb::blocked_range<size_t>(0, source_.size()), ProcessVertexDijkstra(*this));
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

		void ListBlocksSequences(const BlockList & block, const std::string & directory) const
		{
			std::vector<IndexPair> group;
			BlockList blockList = block;
			GroupBy(blockList, compareById, std::back_inserter(group));
			for (std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
			{
				std::ofstream out;
				std::stringstream ss;
				ss << directory << "/" << blockList[it->first].GetBlockId() << ".fa";
				TryOpenFile(ss.str(), out);
				for (size_t block = it->first; block < it->second; block++)
				{
					size_t length = blockList[block].GetLength();
					size_t chr = blockList[block].GetChrId();
					size_t chrSize = storage_.GetChrSequence(chr).size();
					out << ">" << blockList[block].GetBlockId() << "_" << block - it->first << " ";
					out << storage_.GetChrDescription(chr) << ";";
					if (blockList[block].GetSignedBlockId() > 0)
					{
						out << blockList[block].GetStart() << ";" << length << ";" << "1;" << chrSize << std::endl;
						OutputLines(storage_.GetChrSequence(chr).begin() + blockList[block].GetStart(), length, out);
					}
					else
					{
						size_t start = chrSize - blockList[block].GetEnd();
						out << start << ";" << length << ";" << "-1;" << chrSize << std::endl;
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
			std::string blocksDir = outDir + "/blocks";
			CreateOutDirectory(blocksDir);
			ListBlocksIndicesGFF(instance, outDir + "/" + "blocks_coords.gff");
			ListBlocksIndices(instance, outDir + "/" + "blocks_coords.txt");
			ListBlocksSequences(instance, blocksDir);
			GenerateReport(instance, outDir + "/" + "coverage_report.txt");			
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
		void CalculateCoverage(GroupedBlockList::const_iterator start, GroupedBlockList::const_iterator end, std::vector<double> & ret) const;
		std::string OutputIndex(const BlockInstance & block) const;
		void OutputBlocks(const std::vector<BlockInstance>& block, std::ofstream& out) const;
		void ListBlocksIndices(const BlockList & block, const std::string & fileName) const;
		void ListBlocksIndicesGFF(const BlockList & blockList, const std::string & fileName) const;
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
					currentPath.Score();
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


		std::pair<int32_t, NextVertex> MostPopularVertex(const Path & currentPath, bool forward, std::unordered_map<int32_t, int32_t> & count)
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
			
			return std::make_pair(bestVid, ret);
		}				

		bool ExtendPathDijkstraForward(Path & currentPath,
			size_t & bestRightSize,
			int64_t & bestScore,
			int64_t & nowScore,
			std::pair<int64_t, std::vector<Path::Instance> > & goodInstance)
		{	
			bool success = false;
			int64_t origin = currentPath.Origin();
			std::pair<int32_t, NextVertex> nextForwardVid;
			std::unordered_map<int32_t, int32_t> count;
			nextForwardVid = MostPopularVertex(currentPath, true, count);
			if (nextForwardVid.first != 0)
			{				
				for(auto it = nextForwardVid.second.origin; it.GetVertexId() != nextForwardVid.first; ++it)
				{
#ifdef _DEBUG_OUT_
					if (debug_)
					{
						std::cerr << "Attempting to push back the vertex:" << it.GetVertexId() << std::endl;
					}
					
					if (missingVertex_.count(it.GetVertexId()))
					{
						std::cerr << "Alert: " << it.GetVertexId() << ", origin: " << currentPath.Origin() << std::endl;
					}
#endif
					success = currentPath.PointPushBack(it.OutgoingEdge());
					if (success)
					{
						nowScore = currentPath.Score(scoreFullChains_);
#ifdef _DEBUG_OUT_
						if (debug_)
						{
							std::cerr << "Success! New score:" << nowScore << std::endl;
							currentPath.DumpPath(std::cerr);
							currentPath.DumpInstances(std::cerr);
						}						
#endif												
						if (nowScore > bestScore)
						{
							bestScore = nowScore;
							bestRightSize = currentPath.RightSize();

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
							}							
						}
					}					
				}
			}

			return success;
		}

		struct BranchData
		{
			std::vector<size_t> branchId;
		};

		typedef std::vector<std::vector<size_t> > BubbledBranches;

		struct Fork
		{
			Fork(JunctionStorage::JunctionSequentialIterator it, JunctionStorage::JunctionSequentialIterator jt)
			{
				if (it < jt)
				{
					branch[0] = it;
					branch[1] = jt;
				}
				else
				{
					branch[0] = jt;
					branch[1] = it;
				}
			}

			JunctionStorage::JunctionSequentialIterator branch[2];
		};

		int64_t ChainLength(const Fork & now, const Fork & next) const
		{
			return min(abs(now.branch[0].GetPosition() - next.branch[0].GetPosition()), abs(now.branch[1].GetPosition() - next.branch[1].GetPosition()));
		}

		Fork ExpandSourceFork(const Fork & source) const
		{
			for (auto now = source; ; )
			{
				auto next = TakeBubbleStep(now);
				if (next.branch[0].Valid())
				{
					int64_t vid0 = now.branch[0].GetVertexId();
					int64_t vid1 = now.branch[1].GetVertexId();
					assert(vid0 == vid1 && abs(now.branch[0].GetPosition() - next.branch[0].GetPosition()) < maxBranchSize_ &&
						abs(now.branch[1].GetPosition() - next.branch[1].GetPosition()) < maxBranchSize_);
					now = next;
				}
				else
				{
					return now;
				}
			}

			return source;
		}

		Fork TakeBubbleStep(const Fork & source) const
		{
			auto it = source.branch[0];
			std::map<int64_t, int64_t> firstBranch;
			for (int64_t i = 1; abs(it.GetPosition() - source.branch[0].GetPosition()) < maxBranchSize_ && (++it).Valid(); i++)
			{
				int64_t d = abs(it.GetPosition() - source.branch[0].GetPosition());
				firstBranch[it.GetVertexId()] = i;
			}

			it = source.branch[1];
			for (int64_t i = 1; abs(it.GetPosition() - source.branch[1].GetPosition()) < maxBranchSize_ && (++it).Valid(); i++)
			{
				auto kt = firstBranch.find(it.GetVertexId());
				if (kt != firstBranch.end())
				{
					return Fork(source.branch[0] + kt->second, it);
				}
			}

			return Fork(JunctionStorage::JunctionSequentialIterator(), JunctionStorage::JunctionSequentialIterator());
		}

		void BubbledBranchesForward(int64_t vertexId, const std::vector<JunctionStorage::JunctionSequentialIterator> & instance, BubbledBranches & bulges) const
		{
			std::vector<size_t> parallelEdge[5];
			std::map<int64_t, BranchData> visit;
			bulges.assign(instance.size(), std::vector<size_t>());
			for (size_t i = 0; i < instance.size(); i++)
			{
				auto vertex = instance[i];
				if ((vertex + 1).Valid())
				{
					parallelEdge[TwoPaCo::DnaChar::MakeUpChar(vertex.GetChar())].push_back(i);
				}

				for (int64_t startPosition = vertex++.GetPosition(); vertex.Valid() && abs(startPosition - vertex.GetPosition()) <= maxBranchSize_; ++vertex)
				{
					int64_t nowVertexId = vertex.GetVertexId();
					auto point = visit.find(nowVertexId);
					if (point == visit.end())
					{
						BranchData bData;
						bData.branchId.push_back(i);
						visit[nowVertexId] = bData;
					}
					else
					{
						point->second.branchId.push_back(i);
					}
				}
			}

			for (size_t i = 0; i < 5; i++)
			{
				for (size_t j = 0; j < parallelEdge[i].size(); j++)
				{
					for (size_t k = j + 1; k < parallelEdge[i].size(); k++)
					{
						size_t smallBranch = parallelEdge[i][j];
						size_t largeBranch = parallelEdge[i][k];
						bulges[smallBranch].push_back(largeBranch);
					}
				}
			}

			for (auto point = visit.begin(); point != visit.end(); ++point)
			{
				std::sort(point->second.branchId.begin(), point->second.branchId.end());
				for (size_t j = 0; j < point->second.branchId.size(); j++)
				{
					for (size_t k = j + 1; k < point->second.branchId.size(); k++)
					{
						size_t smallBranch = point->second.branchId[j];
						size_t largeBranch = point->second.branchId[k];
						if (smallBranch != largeBranch && std::find(bulges[smallBranch].begin(), bulges[smallBranch].end(), largeBranch) == bulges[smallBranch].end())
						{
							bulges[smallBranch].push_back(largeBranch);
						}
					}
				}
			}
		}

		void BubbledBranchesBackward(int64_t vertexId, const std::vector<JunctionStorage::JunctionSequentialIterator> & instance, BubbledBranches & bulges) const
		{
			std::vector<size_t> parallelEdge[5];
			std::map<int64_t, BranchData> visit;
			bulges.assign(instance.size(), std::vector<size_t>());
			for (size_t i = 0; i < instance.size(); i++)
			{
				auto vertex = instance[i];
				auto prev = vertex - 1;
				if (prev.Valid())
				{
					parallelEdge[TwoPaCo::DnaChar::MakeUpChar(prev.GetChar())].push_back(i);
				}

				for (int64_t startPosition = vertex--.GetPosition(); vertex.Valid() && abs(startPosition - vertex.GetPosition()) <= maxBranchSize_; --vertex)
				{
					int64_t nowVertexId = vertex.GetVertexId();
					auto point = visit.find(nowVertexId);
					if (point == visit.end())
					{
						BranchData bData;
						bData.branchId.push_back(i);
						visit[nowVertexId] = bData;
					}
					else
					{
						point->second.branchId.push_back(i);
					}
				}
			}

			for (size_t i = 0; i < 5; i++)
			{
				for (size_t j = 0; j < parallelEdge[i].size(); j++)
				{
					for (size_t k = j + 1; k < parallelEdge[i].size(); k++)
					{
						size_t smallBranch = parallelEdge[i][j];
						size_t largeBranch = parallelEdge[i][k];
						bulges[smallBranch].push_back(largeBranch);
					}
				}
			}

			for (auto point = visit.begin(); point != visit.end(); ++point)
			{
				std::sort(point->second.branchId.begin(), point->second.branchId.end());
				for (size_t j = 0; j < point->second.branchId.size(); j++)
				{
					for (size_t k = j + 1; k < point->second.branchId.size(); k++)
					{
						size_t smallBranch = point->second.branchId[j];
						size_t largeBranch = point->second.branchId[k];
						if (smallBranch != largeBranch && std::find(bulges[smallBranch].begin(), bulges[smallBranch].end(), largeBranch) == bulges[smallBranch].end())
						{
							bulges[smallBranch].push_back(largeBranch);
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
		int64_t maxFlankingSize_;
		JunctionStorage & storage_;
		std::vector<std::vector<Edge> > syntenyPath_;
		std::vector<std::vector<Assignment> > blockId_;
		tbb::concurrent_vector<std::pair<int64_t, int64_t> > source_;
#ifdef _DEBUG_OUT_
		bool debug_;
		std::set<int64_t> missingVertex_;
#endif
	};
}

#endif
