#ifndef _TRASERVAL_H_
#define _TRAVERSAL_H_

#include "edgestorage.h"

#include <set>
#include <map>

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
			blockId_.resize(storage_.GetChrNumber());
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				blockId_[i].assign(storage_.GetChrEdgeCount(i), UNKNOWN_BLOCK);
			}
			
			std::map<EdgeStorage::Edge, uint32_t> bubbleCount;
			for (size_t vid = 0; vid < storage_.GetVerticesNumber(); vid++)
			{
				CountBubbles(vid, bubbleCount);
			}

			std::map<EdgeStorage::Edge, uint32_t> edgeLength;
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				for (auto it = storage_.GetIterator(i, 0); it.CanInc(); ++it)
				{
					edgeLength[*it] = it.GetLength();
				}
			}

			for (auto it = bubbleCount.rbegin(); it != bubbleCount.rend(); ++it)
			{
				ExtendSeedEdge(it->first, edgeLength);
			}

			edgeLength.clear();
		}

		void GenerateOutput(std::ostream & out) const
		{
			std::vector<BlockInstance> instance;
			for (size_t chr = 0; chr < blockId_.size(); chr++)
			{
				for (size_t i = 0; i < blockId_[i].size(); )
				{
					if (blockId_[chr][i] != UNKNOWN_BLOCK)
					{
						size_t j = i;
						for (; j < blockId_[chr].size() && blockId_[chr][i] == blockId_[chr][j]; j++);
						uint64_t start = storage_.GetIterator(chr, i).GetPosition();
						uint64_t end = storage_.GetIterator(chr, j).GetPosition();
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
			int score;
			std::vector<uint32_t> vertex;
			std::vector<uint32_t> distance;
		};

		void ExtendSeedEdge(EdgeStorage::Edge edge, std::map<EdgeStorage::Edge, uint32_t> & edgeLength)
		{
			Path bestPath;			
			bestPath.vertex.push_back(edge.GetStartVertex());
			bestPath.vertex.push_back(edge.GetEndVertex());
			bestPath.distance.push_back(edgeLength[edge]);
			std::map<std::pair<uint64_t, uint64_t>, bool> seen;
			while (true)
			{
				Path currentPath = bestPath;
				int prevBestScore = bestPath.score;				
				ExtendPath(currentPath, bestPath, edgeLength, 10, seen);
				if (bestPath.score <= prevBestScore)
				{
					break;
				}
			}

			blocksFound_++;
			ReScorePath(bestPath, seen);
			for (auto it : seen)
			{
				blockId_[it.first.first][it.first.second] = blocksFound_ ? it.second : -blocksFound_;
			}
		}

		uint64_t TraceForward(const Path & path, uint64_t lastHitIdxPath, EdgeStorage::EdgeIterator e, std::map<std::pair<uint64_t, uint64_t>, bool> & seen)
		{
			uint64_t score = 0;
			EdgeStorage::EdgeIterator lastHitEdge = e;
			for (; e.CanInc() && abs(e.GetPosition() - lastHitEdge.GetPosition()) <= maxBranchSize_; ++e)
			{
				auto place = std::make_pair(e.GetChrId(), e.GetIdx());
				if (seen.count(place) > 0 || blockId_[e.GetChrId()][e.GetIdx()] != UNKNOWN_BLOCK)
				{
					break;
				}

				seen[place] = e.IsPositiveStrand();
				auto it = std::find(path.vertex.begin() + lastHitIdxPath, path.vertex.end(), e.GetEndVertexId()) - path.vertex.begin();
				if (it != path.vertex.size())
				{
					if (path.distance[it] <= maxBranchSize_)
					{
						lastHitEdge = e;
						lastHitIdxPath = it;
						score += path.distance[it];
					}
					else
					{
						break;
					}
				}
			}

			return score;
		}

		void ReScorePath(Path & path, std::map<std::pair<uint64_t, uint64_t>, bool> & seen)
		{
			seen.clear();
			path.score = 0;
			for (uint64_t i = 0; i < path.vertex.size(); i++)
			{
				uint64_t v = path.vertex[i];
				for (size_t e = 0; e < storage_.GetOutgoingEdgesCount(path.vertex[i]); e++)
				{
					path.score += TraceForward(path, i, storage_.GetOutgoingEdge(v, e), seen);
				}
			}
		}

		void ExtendPath(Path & currentPath, Path & bestPath, std::map<EdgeStorage::Edge, uint32_t> & edgeLength, int maxDepth, std::map<std::pair<uint64_t, uint64_t>, bool> & seen)
		{
			seen.clear();
			ReScorePath(currentPath, seen);
			if (currentPath.score > bestPath.score)
			{
				bestPath = currentPath;
			}

			if (maxDepth > 0)
			{
				std::vector<uint64_t> adjList;
				storage_.AdjacencyList(currentPath.vertex.back(), adjList);
				for (auto nextVertex : adjList)
				{
					if (std::find(currentPath.vertex.begin(), currentPath.vertex.end(), nextVertex) == currentPath.vertex.end())
					{
						currentPath.distance.push_back(edgeLength[EdgeStorage::Edge(currentPath.vertex.back(), nextVertex)]);
						currentPath.vertex.push_back(nextVertex);
						ExtendPath(currentPath, bestPath, edgeLength, maxDepth - 1, seen);
					}
				}
			}			
		}

		void CountBubbles(uint64_t vertexId, std::map<EdgeStorage::Edge, uint32_t> & bubbleCount)
		{			
			BubbledBranches bulges;
			std::map<size_t, BranchData> visit;
			std::vector<EdgeStorage::EdgeIterator> outEdges;
			for (size_t i = 0; i < storage_.GetOutgoingEdgesCount(vertexId); i++)
			{				
				EdgeStorage::EdgeIterator edge = storage_.GetOutgoingEdge(vertexId, i);
				outEdges.push_back(edge);
				size_t startVertex = edge.GetStartVertexId();				
				for (size_t step = 1; step < maxBranchSize_ && edge.CanInc(); ++edge)
				{
					size_t nowVertexId = edge.GetEndVertexId();
					if (nowVertexId == startVertex)
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
					uint64_t n = point->second.branchId.size();
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
		size_t blocksFound_;
		size_t minBlockSize_;
		size_t maxBranchSize_;		
		const EdgeStorage & storage_;
		std::vector<std::vector<int32_t> > blockId_;

		static const int32_t UNKNOWN_BLOCK;
	};
}

#endif