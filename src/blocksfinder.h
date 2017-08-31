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

		BlocksFinder(JunctionStorage & storage, size_t k) : storage_(storage), k_(k), forbidden_(storage)
		{
			scoreFullChains_ = false;
		}

		void FindBlocks(int64_t minBlockSize, int64_t maxBranchSize, int64_t flankingThreshold, int64_t lookingDepth, int64_t sampleSize, const std::string & debugOut)
		{			
			blocksFound_ = 0;
			sampleSize_ = sampleSize;
			lookingDepth_ = lookingDepth;
			minBlockSize_ = minBlockSize;
			maxBranchSize_ = maxBranchSize;
			flankingThreshold_ = flankingThreshold;
			std::vector<std::vector<bool> > junctionInWork;
			std::vector<std::pair<int64_t, int64_t> > bubbleCountVector;
			blockId_.resize(storage_.GetChrNumber());
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				blockId_[i].resize(storage_.GetChrVerticesCount(i));
			}

			int64_t count = 0;
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

			std::random_shuffle(shuffle.begin(), shuffle.end());
			std::ofstream debugStream(debugOut.c_str());
			BestPath bestPath;
			Path currentPath(storage_, maxBranchSize_, minBlockSize_, flankingThreshold_, blockId_);
			time_t mark = time(0);
			for(auto vid : shuffle)
			{
				if (count++ % 1000 == 0)
				{
					std::cerr << count << '\t' << shuffle.size() << std::endl;
				}

				ExtendSeed(vid, currentPath, bestPath, debugStream);
			}

			std::cout << "Time: " << time(0) - mark << std::endl;
		}

		void Dump(std::ostream & out) const
		{
			return;
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
					if (blockId_[chr][i].block != Assignment::UNKNOWN_BLOCK)
					{
						int64_t bid = blockId_[chr][i].block;
						size_t j = i;
						for (; j < blockId_[chr].size() && blockId_[chr][i] == blockId_[chr][j]; j++);
						j--;
						int64_t cstart = storage_.GetIterator(chr, i, bid > 0).GetPosition(&storage_);
						int64_t cend = storage_.GetIterator(chr, j, bid > 0).GetPosition(&storage_) + (bid > 0 ? k_ : -k_);
						int64_t start = std::min(cstart, cend);
						int64_t end = std::max(cstart, cend);
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

		struct BranchData
		{
			BranchData() {}
			std::vector<size_t> branchId;
		};

		typedef std::vector< std::vector<size_t> > BubbledBranches;
		
	
		void ExtendSeed(int64_t vid,
			Path & currentPath,
			BestPath & bestPath,
			std::ostream & debugOut)
		{			
			bestPath.Init();
			currentPath.Init(vid);
			while (true)
			{				
				int64_t prevBestScore = bestPath.score_;
				if (sampleSize_ > 0)
				{
					ExtendPathRandom(currentPath, bestPath, lookingDepth_);
					if (bestPath.score_ <= prevBestScore)
					{
						break;
					}
				}
				else
				{
					ExtendPathBackward(currentPath, bestPath, lookingDepth_);
					bestPath.FixBackward(currentPath);					
					ExtendPathForward(currentPath, bestPath, lookingDepth_);
					bestPath.FixForward(currentPath);
					if (bestPath.score_ <= prevBestScore)
					{
						break;
					}
				}							
			}
						
			if (currentPath.Score(true) > 0 && currentPath.MiddlePathLength() >= minBlockSize_ && currentPath.GoodInstances() > 1)
			{
				++blocksFound_;
//				debugOut << "Block No." << blocksFound_ << ":"  << std::endl;
//				currentPath.DebugOut(debugOut, false);
				syntenyPath_.push_back(std::vector<Edge>());
				std::vector<Edge> nowPathBody;
				currentPath.DumpPath(nowPathBody);
				for (auto pt : nowPathBody)
				{					
					syntenyPath_.back().push_back(pt);
				}

				for (size_t i = 0; i < syntenyPath_.back().size() - 1; i++)
				{
					forbidden_.Add(syntenyPath_.back()[i]);
				}

				int64_t instanceCount = 0;
				for (auto & instance : currentPath.Instances())
				{
					if (currentPath.IsGoodInstance(instance))
					{
						auto end = instance.Back() + 1;
						for (auto it = instance.Front(); it != end; ++it)
						{
							int64_t idx = it.GetIndex();
							int64_t maxidx = storage_.GetChrVerticesCount(it.GetChrId());							
							blockId_[it.GetChrId()][it.GetIndex()].block = it.IsPositiveStrand() ? blocksFound_ : -blocksFound_;
							blockId_[it.GetChrId()][it.GetIndex()].instance = instanceCount;
						}

						instanceCount++;
					}					
				}
			}			
			
			currentPath.Clean();
		}

		void ExtendPathRandom(Path & currentPath, BestPath & bestPath, int maxDepth)
		{
			for (size_t sample = 0; sample < sampleSize_; sample++)
			{
				for (size_t d = 0; ; d++)
				{
					bool over = true;
					for (size_t i = 0; i < 4 && d < lookingDepth_; i++)
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
			
			for (size_t sample = 0; sample < sampleSize_; sample++)
			{
				for (size_t d = 0; ; d++)
				{
					bool over = true;
					for (size_t i = 0; i < 4 && d < lookingDepth_; i++)
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
				for (int64_t idx = 0; idx < storage_.OutgoingEdgesNumber(prevVertex); idx++)
				{
					Edge e = storage_.OutgoingEdge(prevVertex, idx);
					if (forbidden_.Notin(e))
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
				for (int64_t idx = 0; idx < storage_.IngoingEdgesNumber(prevVertex); idx++)
				{
					Edge e = storage_.IngoingEdge(prevVertex, idx);
					if (forbidden_.Notin(e))
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

		void CountBubbles(int64_t vertexId, std::map<int64_t, int64_t> & bubbleCount)
		{
			BubbledBranches bulges;
			std::map<int64_t, BranchData> visit;
			std::vector<JunctionStorage::JunctionIterator> instance;
			for (size_t i = 0; i < storage_.GetInstancesCount(vertexId); i++)
			{				
				JunctionStorage::JunctionIterator vertex = storage_.GetJunctionInstance(vertexId, i);
				instance.push_back(vertex);				
				for (int64_t startPosition = vertex++.GetPosition(&storage_); vertex.Valid(&storage_) && abs(startPosition - vertex.GetPosition(&storage_)) < maxBranchSize_; ++vertex)
				{					
					int64_t nowVertexId = vertex.GetVertexId(&storage_);
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
						break;
					}
				}
			}

			for (auto point = visit.begin(); point != visit.end(); ++point)
			{
				if (point->second.branchId.size() > 1)
				{
					size_t n = point->second.branchId.size();
					for (size_t i = 0; i < point->second.branchId.size(); ++i)
					{
						auto & edgeIt = instance[point->second.branchId[i]];
						if (edgeIt.IsPositiveStrand())
						{
							bubbleCount[edgeIt.GetVertexId(&storage_)] += n * (n - 1) / 2;
						}
					}
				}
			}
		}
	
		int64_t k_;
		int64_t sampleSize_;
		int64_t blocksFound_;
		Forbidden forbidden_;
		bool scoreFullChains_;		
		int64_t lookingDepth_;		
		int64_t minBlockSize_;
		int64_t maxBranchSize_;
		int64_t flankingThreshold_;
		JunctionStorage & storage_;
		std::vector<std::vector<Edge> > syntenyPath_;
		std::vector<std::vector<Assignment> > blockId_;	
	};
}

#endif