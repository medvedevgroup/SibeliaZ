#ifndef _TRASERVAL_H_
#define _TRAVERSAL_H_

#include "junctionstorage.h"

#include <set>
#include <map>
#include <list>
#include <deque>
#include <numeric>
#include <sstream>
#include <iostream>

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
		const char COVERED = 1;
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

		BlocksFinder(const JunctionStorage & storage, size_t k) : storage_(storage), k_(k)
		{
			
		}

		void FindBlocks(int64_t minBlockSize, int64_t maxBranchSize, int64_t flankingThreshold)
		{
			blocksFound_ = 0;
			minBlockSize_ = minBlockSize;
			maxBranchSize_ = maxBranchSize;
			flankingThreshold_ = flankingThreshold;
			std::vector<std::pair<int64_t, int64_t> > bubbleCountVector;

			blockId_.resize(storage_.GetChrNumber());
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				blockId_[i].assign(storage_.GetChrVerticesCount(i), UNKNOWN_BLOCK);
			}
			
			std::map<int64_t, int64_t> bubbleCount;
			for (int64_t vid = -storage_.GetVerticesNumber() + 1; vid < storage_.GetVerticesNumber(); vid++)
			{
				CountBubbles(vid, bubbleCount);
			}

			std::map<Edge, int64_t> edgeLength;			
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				for (auto it = storage_.Begin(i); it != storage_.End(i) - 1; )
				{
					auto jt = it + 1;
					Edge e(it.GetVertexId(), jt.GetVertexId());
					edgeLength[e] = abs(jt.GetPosition() - it.GetPosition());
				}
			}

			for (auto it = bubbleCount.rbegin(); it != bubbleCount.rend(); ++it)
			{
				bubbleCountVector.push_back(std::make_pair(it->second, it->first));
			}

			std::sort(bubbleCountVector.begin(), bubbleCountVector.end());
			int count = 0;
			for (auto it = bubbleCountVector.rbegin(); it != bubbleCountVector.rend(); ++it)
			{
				if (count++ % 100 == 0)
				{
					std::cerr << count << '\t' << bubbleCountVector.size() << std::endl;
				}

				ExtendSeed(it->second, edgeLength, bubbleCount);
			}

			edgeLength.clear();
		}

		void Dump(std::ostream & out) const
		{
			out << "digraph G\n{\nrankdir = LR" << std::endl;
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				for (auto it = storage_.Begin(i); it != storage_.End(i) - 1; ++it)
				{
					auto jt = it + 1;
					out << it.GetVertexId() << " -> " << jt.GetVertexId() << "[label=\"" << it.GetChrId() << ", " << it.GetPosition() << "\" color=blue]" << std::endl;
					out << -it.GetVertexId() << " -> " << -it.GetVertexId() << "[label=\"" << it.GetChrId() << ", " << it.GetPosition() << "\" color=red]" << std::endl;
				}
			}

			for (size_t i = 0; i < syntenyPath_.size(); i++)
			{
				for (size_t j = 0; j < syntenyPath_[i].size() - 1; j++)
				{
					out << syntenyPath_[i][j] << " -> " << syntenyPath_[i][j + 1] << "[label=\"" << i + 1 << "\" color=green]" << std::endl;
					out << -syntenyPath_[i][j + 1] << " -> " << -syntenyPath_[i][j] << "[label=\"" << -int64_t(i + 1) << "\" color=green]" << std::endl;
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
			for (size_t chr = 0; chr < blockId_.size(); chr++)
			{
				for (size_t i = 0; i < blockId_[chr].size();)
				{
					if (blockId_[chr][i] != UNKNOWN_BLOCK)
					{
						int64_t bid = blockId_[chr][i];
						size_t j = i;
						for (; j < blockId_[chr].size() && blockId_[chr][i] == blockId_[chr][j]; j++);
						j--;
						int64_t start = storage_.GetIterator(chr, i, bid > 0).GetPosition();
						int64_t end = storage_.GetIterator(chr, j, bid > 0).GetPosition();
						instance.push_back(BlockInstance(blockId_[chr][i], chr, std::min(start, end), std::max(start, end)));
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

		void GenerateReport(const BlockList & block, const std::string & fileName) const
		{
			std::ofstream out;
			TryOpenFile(fileName, out);
			GroupedBlockList sepBlock;
			std::vector<IndexPair> group;
			BlockList blockList = block;
			GroupBy(blockList, compareById, std::back_inserter(group));
			for (std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
			{
				sepBlock.push_back(std::make_pair(it->second - it->first, std::vector<BlockInstance>(blockList.begin() + it->first, blockList.begin() + it->second)));
			}

			ListChrs(out);
			out << "Degree\tCount\tTotal";
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				out << "\tSeq " << i + 1;
			}

			out << std::endl;
			group.clear();
			GroupBy(sepBlock, ByFirstElement, std::back_inserter(group));
			group.push_back(IndexPair(0, sepBlock.size()));
			for (std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
			{
				if (it != group.end() - 1)
				{
					out << sepBlock[it->first].first << '\t' << it->second - it->first << '\t';
				}
				else
				{
					out << "All\t" << it->second - it->first << "\t";
				}

				out.precision(2);
				out.setf(std::ostream::fixed);
				std::vector<double> coverage = CalculateCoverage(sepBlock.begin() + it->first, sepBlock.begin() + it->second);
				std::copy(coverage.begin(), coverage.end(), std::ostream_iterator<double>(out, "%\t"));
				out << std::endl;
			}

			out << DELIMITER << std::endl;
		}

		std::vector<double> CalculateCoverage(GroupedBlockList::const_iterator start, GroupedBlockList::const_iterator end) const
		{
			std::vector<double> ret;
			std::vector<char> cover;
			double totalBp = 0;
			double totalCoveredBp = 0;
			for (size_t chr = 0; chr < storage_.GetChrNumber(); chr++)
			{
				totalBp += storage_.GetChrSequence(chr).size();
				cover.assign(storage_.GetChrSequence(chr).size(), 0);
				for (GroupedBlockList::const_iterator it = start; it != end; ++it)
				{
					for (size_t i = 0; i < it->second.size(); i++)
					{
						if (it->second[i].GetChrId() == chr)
						{
							std::fill(cover.begin() + it->second[i].GetStart(), cover.begin() + it->second[i].GetEnd(), COVERED);
						}
					}
				}

				double nowCoveredBp = static_cast<double>(std::count(cover.begin(), cover.end(), COVERED));
				ret.push_back(nowCoveredBp / cover.size() * 100);
				totalCoveredBp += nowCoveredBp;
			}

			ret.insert(ret.begin(), totalCoveredBp / totalBp * 100);
			return ret;
		}


		std::string OutputIndex(const BlockInstance & block) const
		{
			std::stringstream out;
			out << block.GetChrId() + 1 << '\t' << (block.GetSignedBlockId() < 0 ? '-' : '+') << '\t';
			out << block.GetConventionalStart() << '\t' << block.GetConventionalEnd() << '\t' << block.GetEnd() - block.GetStart();
			return out.str();
		}


		void OutputBlocks(const std::vector<BlockInstance>& block, std::ofstream& out) const
		{
			std::vector<IndexPair> group;
			std::vector<BlockInstance> blockList = block;
			GroupBy(blockList, compareById, std::back_inserter(group));
			for (std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
			{
				size_t length = it->second - it->first;
				std::sort(blockList.begin() + it->first, blockList.begin() + it->second, compareByChrId);
				out << "Block #" << blockList[it->first].GetBlockId() << std::endl;
				out << "Seq_id\tStrand\tStart\tEnd\tLength" << std::endl;
				for (auto jt = blockList.begin() + it->first; jt < blockList.begin() + it->first + length; ++jt)
				{
					out << OutputIndex(*jt) << std::endl;
				}

				out << DELIMITER << std::endl;
			}
		}


		void ListBlocksIndices(const BlockList & block, const std::string & fileName) const
		{
			std::ofstream out;
			TryOpenFile(fileName, out);
			ListChrs(out);
			OutputBlocks(block, out);
		}

		void ListChromosomesAsPermutations(const BlockList & block, const std::string & fileName) const
		{
			std::ofstream out;
			TryOpenFile(fileName, out);
			std::vector<IndexPair> group;
			BlockList blockList = block;
			GroupBy(blockList, compareByChrId, std::back_inserter(group));
			for (std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
			{
				out.setf(std::ios_base::showpos);
				size_t length = it->second - it->first;
				size_t chr = blockList[it->first].GetChrId();
				out << '>' << storage_.GetChrDescription(chr) << std::endl;
				std::sort(blockList.begin() + it->first, blockList.begin() + it->second);
				for (auto jt = blockList.begin() + it->first; jt < blockList.begin() + it->first + length; ++jt)
				{
					out << jt->GetSignedBlockId() << " ";
				}

				out << "$" << std::endl;
			}
		}


		void TryOpenFile(const std::string & fileName, std::ofstream & stream) const
		{
			stream.open(fileName.c_str());
			if (!stream)
			{
				throw std::runtime_error(("Cannot open file " + fileName).c_str());
			}
		}


		void ListChrs(std::ostream & out) const
		{
			out << "Seq_id\tSize\tDescription" << std::endl;
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				out << i + 1 << '\t' << storage_.GetChrSequence(i).size() << '\t' << storage_.GetChrDescription(i) << std::endl;
			}

			out << DELIMITER << std::endl;
		}

		struct BranchData
		{
			BranchData() {}
			std::vector<size_t> branchId;
		};

		typedef std::vector< std::vector<size_t> > BubbledBranches;

		struct Path
		{
			Path() : score(0) {}
			Path(int64_t start, const JunctionStorage & storage) : score(0)
			{
				body.push_back(Point(start, 0));
				CheckPushBack(storage, 0);
			}

			static bool Compatible(EdgeStorage::EdgeIterator start, EdgeStorage::EdgeIterator end, int64_t maxBranchSize)
			{
				if (!(start.GetChrId() == end.GetChrId() && start.IsPositiveStrand() == end.IsPositiveStrand()))
				{
					return false;
				}

				if ((start.IsPositiveStrand() && end.GetEndPosition() - start.GetStartPosition() > maxBranchSize) ||
					(!start.IsPositiveStrand() && -(end.GetEndPosition() - start.GetStartPosition()) > maxBranchSize))
				{
					return false;
				}

				return true;
			}

			bool Score(EdgeStorage::EdgeIterator start, EdgeStorage::EdgeIterator end)
			{
				return abs(end.GetEndPosition() - start.GetStartPosition());
			}

			void CheckPushBack(const EdgeStorage & storage, int64_t maxBranchSize)
			{
				int64_t v = body.back().vertex;
				for (size_t i = 0; i < storage.GetOutgoingEdgesCount(v); i++)
				{
					EdgeStorage::EdgeIterator it = storage.GetOutgoingEdge(v, i);
					for (auto & p : instance)
					{
						if (Compatible(p.back(), it, maxBranchSize))
						{
							score += Score(p.back(), it);
							p.push_back(it);
							break;
						}
					}
				}				
			}

			int64_t score;

			struct Point
			{
				int64_t vertex;
				int64_t distance;
				Point() {}
				Point(int64_t vertex, int64_t distance) : vertex(vertex), distance(distance) {}
			};

			typedef std::deque<JunctionStorage::JunctionIterator> Instance;

			std::deque<Point> body;
			std::list<Instance> instance;

			bool IsInPath(int64_t vertex) const
			{
				for (auto & it : body)
				{
					if (it.vertex == vertex)
					{
						return true;
					}
				}

				return false;
			}

			void PointPushBack(int64_t vertex, int64_t length, const JunctionStorage & storage, int64_t maxBranchSize)
			{
				body.push_back(Point(vertex, body.back().distance + length));
				CheckPushBack(storage, maxBranchSize);
			}

			void PointPushFront(int64_t vertex, int64_t length, const JunctionStorage & storage, int64_t maxBranchSize)
			{
				body.push_front(Point(vertex, body.front().distance - length));
			}

			void PointPopBack()
			{
				body.pop_back();
			}

			void PointPopFront()
			{
				body.pop_front();
			}

			void ExportPath(std::vector<std::vector<int64_t> > & result) const
			{
				result.push_back(std::vector<int64_t>());
				for (auto & p : body)
				{
					result.back().push_back(p.vertex);
				}
			}

			int64_t Length() const
			{
				return body.back().distance - body.front().distance;
			}
		};

		void ExtendSeed(int64_t vertex, std::map<Edge, int64_t> & edgeLength, const std::map<int64_t, int64_t> & bubbleCount)
		{
			Path bestPath(vertex, storage_);
			while (true)
			{
				Path currentPath = bestPath;
				int64_t prevBestScore = bestPath.score;
				ExtendPathBackward(currentPath, bestPath, edgeLength, 8);
				currentPath = bestPath;
				ExtendPathForward(currentPath, bestPath, edgeLength, 8);
				if (bestPath.score <= prevBestScore)
				{
					break;
				}
			}

			if (bestPath.score > 0 && bestPath.Length() >= minBlockSize_ && bestPath.instance.size() > 1)
			{
				blocksFound_++;
				bestPath.ExportPath(syntenyPath_);
				for (size_t i = 0; i < syntenyPath_.back().size(); i++)
				{
					forbidden_.insert(syntenyPath_.back()[i]);
				}

				for (auto & instance : bestPath.instance)
				{
					for (auto it : instance)
					{						
						blockId_[it.GetChrId()][it.GetIdx()] = it.IsPositiveStrand() ? blocksFound_ : -blocksFound_;
					}
				}				
			}
		}
		/*
		void RescorePath(Path & path, std::map<std::pair<uint64_t, uint64_t>, bool> & assignment, bool finalScoring = true)
		{
			path.score = 0;
			path.instances = 0;
			assignment.clear();
			std::set<std::pair<uint64_t, uint64_t> >  seen;
			for (size_t startVertexIdx = 0; startVertexIdx < path.vertex.size(); startVertexIdx++)
			{				
				int64_t startVertex = path.vertex[startVertexIdx];
				for (size_t e = 0; e < storage_.GetOutgoingEdgesCount(startVertex); e++)
				{
					bool anyHit = false;
					int64_t localScore = 0;
					size_t lastHitVertexIdx = startVertexIdx;
					EdgeStorage::EdgeIterator edge = storage_.GetOutgoingEdge(startVertex, e);					
					int64_t lastHitPosition = edge.GetStartPosition();
					std::vector<std::pair<uint64_t, uint64_t> > history;
					for (; edge.Valid() && abs(edge.GetEndPosition() - lastHitPosition) < maxBranchSize_; ++edge)
					{
						auto place = std::make_pair(edge.GetChrId(), edge.GetIdx());						
						if (seen.count(place) > 0 || blockId_[edge.GetChrId()][edge.GetIdx()] != UNKNOWN_BLOCK)
						{
							break;
						}

						seen.insert(place);
						history.push_back(place);
						int64_t endVertex = edge.GetEndVertexId();
						size_t it = std::find(path.vertex.begin() + lastHitVertexIdx + 1, path.vertex.end(), endVertex) - path.vertex.begin();

						if (it != path.vertex.size())
						{
							int64_t distance = path.distance[it] - path.distance[lastHitVertexIdx];
							if (distance <= maxBranchSize_)
							{						
								anyHit = true;
								localScore += distance;
								lastHitVertexIdx = it;
								lastHitPosition = edge.GetEndPosition();								
							}
							else
							{
								break;
							}
						}
					}
						
					if (anyHit)
					{
						int64_t leftFlank = path.distance[startVertexIdx] - path.distance.front();
						int64_t chainLength = path.distance[lastHitVertexIdx] - path.distance[startVertexIdx];
						int64_t rightFlank = path.distance.back() - path.distance[lastHitVertexIdx];
						if (leftFlank > flankingThreshold_ || rightFlank > flankingThreshold_)
						{
							if (chainLength >= minBlockSize_ - flankingThreshold_ * 2)
							{
								path.score = -INT_MAX;
								return;
							}
						}
						else
						{
							path.instances++;
							path.score += localScore - leftFlank - rightFlank;							
							for (auto place : history)
							
								assignment[place] = edge.IsPositiveStrand();
							}
						}
					}								
				}
			}
			
		}*/

		void ExtendPathForward(Path & currentPath, Path & bestPath, std::map<Edge, int64_t> & edgeLength, int maxDepth)
		{
			if (maxDepth > 0)
			{
				std::vector<int64_t> adjList;
				int64_t prevVertex = currentPath.body.back().vertex;
				storage_.SuccessorsList(prevVertex, adjList);
				for (auto nextVertex : adjList)
				{
					if (currentPath.IsInPath(nextVertex) && forbidden_.count(Edge(prevVertex, nextVertex)) == 0)
					{
						currentPath.PointPushBack(nextVertex, edgeLength[Edge(currentPath.body.back().vertex, nextVertex)], storage_, maxBranchSize_);
						ExtendPathForward(currentPath, bestPath, edgeLength, maxDepth - 1);
						if (currentPath.score > bestPath.score && currentPath.instance.size() > 1)
						{
							bestPath = currentPath;
						}

						currentPath.PointPopBack();
					}
				}
			}
		}

		void ExtendPathBackward(Path & currentPath, Path & bestPath, std::map<Edge, int64_t> & edgeLength, int maxDepth)
		{
			if (maxDepth > 0)
			{
				std::vector<int64_t> adjList;
				int64_t prevVertex = currentPath.body.front().vertex;
				storage_.PredecessorsList(prevVertex, adjList);
				for (auto nextVertex : adjList)
				{
					if (currentPath.IsInPath(nextVertex) && forbidden_.count(Edge(nextVertex, prevVertex)) == 0)
					{							
						currentPath.PointPushFront(nextVertex, edgeLength[Edge(nextVertex, prevVertex)], storage_, maxBranchSize_);
						if (currentPath.score > bestPath.score && currentPath.instance.size() > 1)
						{
							bestPath = currentPath;
						}

						ExtendPathBackward(currentPath, bestPath, edgeLength, maxDepth - 1);
						currentPath.PointPopFront();
					}
				}
			}
		}
		
		/*
		double EdgeWeight(const std::map<EdgeStorage::Edge, uint32_t> & bubbleCount, EdgeStorage::Edge e) const
		{
			if (bubbleCount.count(e) > 0)
			{
				return bubbleCount.find(e)->second;
			}

			if (bubbleCount.count(e.Reverse()) > 0)
			{
				return bubbleCount.find(e.Reverse())->second;
			}

			return .5;
		}
		
		bool RandomChoice(Path & currentPath, std::map<EdgeStorage::Edge, int64_t> & edgeLength, const std::map<EdgeStorage::Edge, uint32_t> & bubbleCount)
		{
			std::vector<bool> isForward;
			std::vector<double> probability;
			std::vector<double> weight;
			std::vector<int64_t> nextVertex;
			std::vector<int64_t> adjForwardList;
			std::vector<int64_t> adjBackwardList;
			storage_.SuccessorsList(currentPath.vertex.back(), adjForwardList);
			storage_.PredecessorsList(currentPath.vertex.front(), adjBackwardList);
			for (int64_t v : adjForwardList)
			{
				EdgeStorage::Edge e(currentPath.vertex.back(), v);
				if (!forbidden_.count(e) && !forbidden_.count(e.Reverse()) && 
					std::find(currentPath.vertex.begin(), currentPath.vertex.end(), v) == currentPath.vertex.end())
				{
					nextVertex.push_back(v);
					isForward.push_back(true);
					weight.push_back(EdgeWeigth(bubbleCount, e));
				}
			}

			for (int64_t v : adjBackwardList)
			{
				EdgeStorage::Edge e(v, currentPath.vertex.front());
				if (!forbidden_.count(e) && !forbidden_.count(e.Reverse()) && 
					std::find(currentPath.vertex.begin(), currentPath.vertex.end(), v) == currentPath.vertex.end())
				{
					nextVertex.push_back(v);
					isForward.push_back(false);
					weight.push_back(EdgeWeigth(bubbleCount, e));
				}
			}

			if (weight.size() > 0)
			{
				double total = std::accumulate(weight.begin(), weight.end(), 0.0);
				for (size_t i = 0; i < weight.size(); i++)
				{
					probability.push_back(weight[i] / total + (i > 0 ? probability[i - 1] : 0));
				}

				probability.back() = 1.01;
				double coin = double(rand()) / RAND_MAX;
				size_t it = std::lower_bound(probability.begin(), probability.end(), coin) - probability.begin();
				if (isForward[it])
				{
					currentPath.distance.push_back(currentPath.distance.back() + edgeLength[EdgeStorage::Edge(currentPath.vertex.back(), nextVertex[it])]);
					currentPath.vertex.push_back(nextVertex[it]);
				}
				else
				{
					currentPath.distance.push_front(currentPath.distance.front() -edgeLength[EdgeStorage::Edge(nextVertex[it], currentPath.vertex.front())]);
					currentPath.vertex.push_front(nextVertex[it]);
				}

				return true;
			}			

			return false;
		}

		void ExtendPathRandom(Path & startPath, Path & bestPath, std::map<EdgeStorage::Edge, int64_t> & edgeLength, int sampleSize, int maxDepth, const std::map<EdgeStorage::Edge, uint32_t> & bubbleCount)		
		{
			std::map<std::pair<uint64_t, uint64_t>, bool> seen;
			for (size_t it = 0; it < sampleSize; it++)
			{
				Path currentPath = startPath;
				for (size_t i = 0; i < maxDepth; i++)
				{
					if (RandomChoice(currentPath, edgeLength, bubbleCount))
					{
						seen.clear();
						RescorePath(currentPath, seen);
						if (currentPath.score > bestPath.score)
						{
							bestPath = currentPath;
						}
					}
					else
					{
						break;
					}
				}
			}						
		}*/		

		void CountBubbles(int64_t vertexId, std::map<int64_t, int64_t> & bubbleCount)
		{
			//THERE IS A BUG HERE: NOT ACCOUNTING FOR BUBBLE'S MINIMALITY
			BubbledBranches bulges;
			std::map<int64_t, BranchData> visit;
			std::vector<JunctionStorage::JunctionIterator> instance;
			for (size_t i = 0; i < storage_.GetInstancesCount(vertexId); i++)
			{				
				JunctionStorage::JunctionIterator vertex = storage_.GetJunctionInstance(vertexId, i);
				instance.push_back(vertex++);
				for (int64_t startPosition = vertex.GetPosition(); abs(startPosition - vertex.GetPosition()) < maxBranchSize_ && vertex.Valid(); ++vertex)
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
							bubbleCount[edgeIt.GetVertexId()] += n * (n - 1) / 2;
						}
					}
				}
			}
		}
		
		size_t k_;
		int64_t blocksFound_;
		int64_t minBlockSize_;
		int64_t maxBranchSize_;
		int64_t flankingThreshold_;
		const JunctionStorage & storage_;
		std::set<Edge> forbidden_;
		std::vector<std::vector<int64_t> > blockId_;
		std::vector<std::vector<int64_t> > syntenyPath_;		

		static const int32_t UNKNOWN_BLOCK;
	};
}

#endif