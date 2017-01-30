#ifndef _TRASERVAL_H_
#define _TRAVERSAL_H_

#include "edgestorage.h"

#include <unordered_map>

namespace Sibelia
{
	class BlocksFinder
	{
	public:

		BlocksFinder(const EdgeStorage storage) : storage_(storage)
		{
			blockId_.resize(storage.GetChrNumber());
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				blockId_[i].assign(storage_.GetChrEdgeCount(i), false);
			}
		}

		void FindBlocks(size_t minBlockSize, size_t maxBranchSize)
		{
			std::vector<std::vector<uint32_t> > bubbleCountPerEdge(storage_.GetChrNumber());
			for (size_t i = 0; i < bubbleCountPerEdge.size(); i++)
			{
				bubbleCountPerEdge[i].assign(storage_.GetChrEdgeCount(i), 0);
			}
		}

	private:

		struct BranchData
		{
			BranchData() {}
			std::vector<size_t> branchId;
		};

		typedef std::vector< std::vector<size_t> > BubbledBranches;

		void AnyBubbles(uint64_t vertexId, std::vector<std::vector<uint32_t> > & bubbleCount)
		{			
			BubbledBranches bulges;
			std::unordered_map<size_t, BranchData> visit;
			std::vector<EdgeStorage::EdgeIterator> outEdges;
			for (size_t i = 0; i < storage_.GetOutgoingEdgesCount(i); i++)
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
						auto & edge = outEdges[point->second.branchId[i]];
						if (edge.IsPositiveStrand())
						{
							bubbleCount[edge.GetChrId()][edge.GetIdx()] += n * (n - 1) / 2;
						}
					}
				}
			}
		}

		void GoThroughBubbles(size_t minBlockSize, size_t maxBranchSize, std::vector<std::vector<uint32_t> > & bubbleCount)
		{
			for (uint64_t vertexId = 0; vertexId < storage_.GetVerticesNumber(); vertexId++)
			{
				AnyBubbles(vertexId, bubbleCount);
			}
		}

		size_t minBlockSize_;
		size_t maxBranchSize_;		
		const EdgeStorage & storage_;
		std::vector<std::vector<uint32_t> > blockId_;
	};
}

#endif