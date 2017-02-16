#ifndef _TRASERVAL_H_
#define _TRAVERSAL_H_

#include "edgestorage.h"

#include <set>
#include <map>
#include <deque>
#include <iostream>

namespace Sibelia
{

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

	class BlocksFinder
	{
	public:

		BlocksFinder(const EdgeStorage & storage, size_t k) : storage_(storage), k_(k)
		{
			
		}

		void FindBlocks(size_t minBlockSize, size_t maxBranchSize)
		{
			blocksFound_ = 0;
			minBlockSize_ = minBlockSize;
			maxBranchSize_ = maxBranchSize;
			std::vector<std::pair<uint32_t, EdgeStorage::Edge> > bubbleCountVector;

			blockId_.resize(storage_.GetChrNumber());
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				blockId_[i].assign(storage_.GetChrEdgeCount(i), UNKNOWN_BLOCK);
			}
			
			std::map<EdgeStorage::Edge, uint32_t> bubbleCount;
			for (int64_t vid = -storage_.GetVerticesNumber() + 1; vid < storage_.GetVerticesNumber(); vid++)
			{
				if (!storage_.EdgesOnlyOnNegativeStrand(vid))
				{
					CountBubbles(vid, bubbleCount);
				}
			}

			std::map<EdgeStorage::Edge, uint64_t> edgeLength;			
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				for (auto it = storage_.GetIterator(i, 0); it.CanInc(); ++it)
				{
					edgeLength[*it] = it.GetLength();
				}
			}

			for (auto it = bubbleCount.rbegin(); it != bubbleCount.rend(); ++it)
			{
				bubbleCountVector.push_back(std::make_pair(it->second, it->first));
			}

			std::sort(bubbleCountVector.begin(), bubbleCountVector.end());
			for (auto it = bubbleCountVector.rbegin(); it != bubbleCountVector.rend(); ++it)
			{
				if (forbidden_.count(it->second) == 0)
				{
					ExtendSeedEdge(it->second, edgeLength);
				}				
			}

			edgeLength.clear();
		}

		void Dump(std::ostream & out) const
		{
			out << "digraph G\n{\nrankdir = LR" << std::endl;
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				for (EdgeStorage::EdgeIterator it = storage_.GetIterator(i, 0); it.CanInc(); ++it)
				{
					out << it.GetStartVertexId() << " -> " << it.GetEndVertexId() << "[label=\"" << it.GetChrId() << ", " << it.GetStartPosition() << "\" color=blue]" << std::endl;
				}
			}

			for (size_t i = 0; i < syntenyPath_.size(); i++)
			{
				for (size_t j = 0; j < syntenyPath_[i].size() - 1; j++)
				{
					out << syntenyPath_[i][j] << " -> " << syntenyPath_[i][j + 1] << "[label=\"" << i << "\" color=red]" << std::endl;
				}
			}

			out << "}" << std::endl;
		}

		void GenerateOutput(std::ostream & out) const
		{
			std::vector<BlockInstance> instance;
			for (size_t chr = 0; chr < blockId_.size(); chr++)
			{
				for (size_t i = 0; i < blockId_[chr].size(); )
				{
					if (blockId_[chr][i] != UNKNOWN_BLOCK)
					{
						size_t j = i;
						for (; j < blockId_[chr].size() && blockId_[chr][i] == blockId_[chr][j]; j++);
						uint64_t start = storage_.GetIterator(chr, i).GetStartPosition();
						uint64_t end = storage_.GetIterator(chr, j).GetStartPosition();
						if (blockId_[chr][i] > 0)
						{
							end += k_;
						}
						else
						{
							start -= k_;
						}

						instance.push_back(BlockInstance(blockId_[chr][i], chr, start, end));
						i = j;
					}
					else
					{
						++i;
					}
				}
			}

			std::sort(instance.begin(), instance.end());
			for (size_t i = 0; i < instance.size(); )
			{
				size_t j = i;
				for (; j < instance.size() && instance[i].GetBlockId() == instance[j].GetBlockId(); ++j)
				{
					out << instance[j].GetSignedBlockId() << '\t' << instance[j].GetChrId() << '\t' << instance[j].GetStart() << '\t' << instance[j].GetEnd() << std::endl;
				}

				i = j;
				out << std::endl;
			}
		}

	private:

		struct BranchData
		{
			BranchData() {}
			std::vector<size_t> branchId;
		};

		typedef std::vector< std::vector<size_t> > BubbledBranches;

		struct Path
		{
			Path() : score(0) {}
			int64_t score;
			std::deque<int64_t> vertex;
			std::deque<int64_t> distance;
		};

		void ExtendSeedEdge(EdgeStorage::Edge edge, std::map<EdgeStorage::Edge, uint64_t> & edgeLength)
		{			
			Path bestPath;
			bestPath.vertex.push_back(edge.GetStartVertex());
			bestPath.vertex.push_back(edge.GetEndVertex());
			bestPath.distance.push_back(0);
			for (size_t i = 0; i < bestPath.vertex.size() - 1; i++)
			{
				bestPath.distance.push_back(edgeLength[Sibelia::EdgeStorage::Edge(bestPath.vertex[i], bestPath.vertex[i + 1])]);
			}

			std::map<std::pair<uint64_t, uint64_t>, bool> seen;
			bool x = edge.GetStartVertex() == 28;
			while (true)
			{
				Path currentPath = bestPath;
				int64_t prevBestScore = bestPath.score;			
				ExtendPathBackward(currentPath, bestPath, edgeLength, 12, seen);
				currentPath = bestPath;				
				ExtendPathForward(currentPath, bestPath, edgeLength, 12, seen);
				if (bestPath.score <= prevBestScore)
				{
					break;
				}
			}

			if (bestPath.score > 0)
			{
				blocksFound_++;
				RescorePath(bestPath, seen);
				syntenyPath_.push_back(bestPath.vertex);
				for (size_t i = 0; i < syntenyPath_.back().size() - 1; i++)
				{
					forbidden_.insert(EdgeStorage::Edge(syntenyPath_.back()[i], syntenyPath_.back()[i + 1]));
				}

				for (auto it : seen)
				{					
					blockId_[it.first.first][it.first.second] = it.second ? blocksFound_ : -blocksFound_;
				}
			}			
		}

		void RescorePath(Path & path, std::map<std::pair<uint64_t, uint64_t>, bool> & seen)
		{
			seen.clear();
			path.score = 0;
			for (size_t startVertexIdx = 0; startVertexIdx < path.vertex.size(); startVertexIdx++)
			{				
				int64_t startVertex = path.vertex[startVertexIdx];
				for (size_t e = 0; e < storage_.GetOutgoingEdgesCount(startVertex); e++)
				{
					bool anyHit = false;
					size_t lastHitVertexIdx = startVertexIdx;
					EdgeStorage::EdgeIterator edge = storage_.GetOutgoingEdge(startVertex, e);					
					int64_t lastHitPosition = edge.GetStartPosition();					
					for (; edge.CanInc() && abs(edge.GetEndPosition() - lastHitPosition) < maxBranchSize_; ++edge)
					{
						auto place = std::make_pair(edge.GetChrId(), edge.GetIdx());
						if (seen.count(place) > 0 || blockId_[edge.GetChrId()][edge.GetIdx()] != UNKNOWN_BLOCK)
						{
							break;
						}

						seen[place] = edge.IsPositiveStrand();
						size_t it = std::find(path.vertex.begin() + lastHitVertexIdx, path.vertex.end(), edge.GetEndVertexId()) - path.vertex.begin();
						if (it != path.vertex.size())
						{
							if (path.distance[it] <= maxBranchSize_)
							{								
								for (size_t j = lastHitVertexIdx + 1; j <= it; ++j)
								{
									path.score += path.distance[j];
									anyHit = true;
								}

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
						for (size_t i = 0; i <= startVertexIdx; i++)
						{
							path.score -= path.distance[i];
						}

						for (size_t i = lastHitVertexIdx + 1; i < path.distance.size(); i++)
						{
							path.score -= path.distance[i];
						}
					}				
				}
			}
			
		}

		void ExtendPathBackward(Path & currentPath, Path & bestPath, std::map<EdgeStorage::Edge, uint64_t> & edgeLength, int maxDepth, std::map<std::pair<uint64_t, uint64_t>, bool> & seen)
		{
			seen.clear();
			RescorePath(currentPath, seen);
			if (currentPath.score > bestPath.score)
			{
				bestPath = currentPath;
			}

			if (maxDepth > 0)
			{
				std::vector<int64_t> adjList;
				storage_.PredecessorsList(currentPath.vertex.front(), adjList);
				for (auto nextVertex : adjList)
				{
					if (std::find(currentPath.vertex.begin(), currentPath.vertex.end(), nextVertex) == currentPath.vertex.end() &&
						forbidden_.count(EdgeStorage::Edge(nextVertex, currentPath.vertex.front())) == 0)
					{						
						currentPath.vertex.push_front(nextVertex);
						currentPath.distance.clear();
						currentPath.distance.push_back(0);
						for (size_t i = 0; i < currentPath.vertex.size() - 1; i++)
						{
							currentPath.distance.push_back(edgeLength[Sibelia::EdgeStorage::Edge(currentPath.vertex[i], currentPath.vertex[i + 1])]);
						}						

						ExtendPathBackward(currentPath, bestPath, edgeLength, maxDepth - 1, seen);
						currentPath.distance.pop_front();
						currentPath.vertex.pop_front();
					}
				}
			}
		}

		void ExtendPathForward(Path & currentPath, Path & bestPath, std::map<EdgeStorage::Edge, uint64_t> & edgeLength, int maxDepth, std::map<std::pair<uint64_t, uint64_t>, bool> & seen)
		{
			seen.clear();
			RescorePath(currentPath, seen);
			if (currentPath.score > bestPath.score)
			{
				bestPath = currentPath;
			}

			if (maxDepth > 0)
			{
				std::vector<int64_t> adjList;
				storage_.SuccessorsList(currentPath.vertex.back(), adjList);
				for (auto nextVertex : adjList)
				{
					if (std::find(currentPath.vertex.begin(), currentPath.vertex.end(), nextVertex) == currentPath.vertex.end() &&
						forbidden_.count(EdgeStorage::Edge(currentPath.vertex.back(), nextVertex)) == 0)
					{
						currentPath.vertex.push_back(nextVertex);
						currentPath.distance.clear();
						currentPath.distance.push_back(0);
						for (size_t i = 0; i < currentPath.vertex.size() - 1; i++)
						{
							currentPath.distance.push_back(edgeLength[EdgeStorage::Edge(currentPath.vertex[i], currentPath.vertex[i + 1])]);
						}

						ExtendPathForward(currentPath, bestPath, edgeLength, maxDepth - 1, seen);
						currentPath.distance.pop_back();
						currentPath.vertex.pop_back();
					}
				}
			}
		}

		void CountBubbles(int64_t vertexId, std::map<EdgeStorage::Edge, uint32_t> & bubbleCount)
		{			
			BubbledBranches bulges;
			std::map<int64_t, BranchData> visit;
			std::vector<EdgeStorage::EdgeIterator> outEdges;
			for (size_t i = 0; i < storage_.GetOutgoingEdgesCount(vertexId); i++)
			{				
				EdgeStorage::EdgeIterator edge = storage_.GetOutgoingEdge(vertexId, i);
				outEdges.push_back(edge);
				for (int64_t startPosition = edge.GetStartPosition(); abs(startPosition - edge.GetEndPosition()) < maxBranchSize_ && edge.CanInc(); ++edge)
				{					
					int64_t nowVertexId = edge.GetEndVertexId();
					if (nowVertexId == vertexId)
					{
						break;
					}

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
						auto & edgeIt = outEdges[point->second.branchId[i]];
						if (edgeIt.IsPositiveStrand())
						{
							bubbleCount[*edgeIt] += n * (n - 1) / 2;
						}
					}
				}
			}
		}
		
		size_t k_;
		int64_t blocksFound_;
		size_t minBlockSize_;
		size_t maxBranchSize_;		
		const EdgeStorage & storage_;
		std::set<EdgeStorage::Edge> forbidden_;
		std::vector<std::vector<int64_t> > blockId_;
		std::vector<std::deque<int64_t> > syntenyPath_;		

		static const int32_t UNKNOWN_BLOCK;
	};
}

#endif