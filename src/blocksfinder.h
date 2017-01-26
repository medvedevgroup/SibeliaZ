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
			mark_.resize(storage.GetChrNumber());
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				mark_[i].assign(storage_.GetChrEdgeCount(i), false);
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
			BranchData(char ch) : endChar(ch) {}
			char endChar;
			std::vector<size_t> branchId;
		};

		typedef std::vector< std::vector<size_t> > BubbledBranches;

		bool AnyBubbles(EdgeStorage & edgeStorage,
			uint64_t vertexId,
			BubbledBranches & bulges,
			size_t minBranchSize)
		{
			bulges.clear();			
			std::unordered_map<size_t, BranchData> visit;
			for (size_t i = 0; i < edgeStorage.GetOutgoingEdgesCount(i); i++)
			{				
				EdgeStorage::EdgeIterator edge = edgeStorage.GetOutgoingEdge(vertexId, i);
				size_t startVertex = edge.GetStartVertexId();				
				char endChar = edge.GetChar();
				for (size_t step = 1; step < minBranchSize && edge.CanInc(); ++edge)
				{
					size_t nowVertexId = edge.GetEndVertexId();
					if (nowVertexId == startVertex)
					{
						break;
					}

					auto point = visit.find(nowVertexId);
					if (point == visit.end())
					{
						BranchData bData(endChar);
						bData.branchId.push_back(i);
						visit[nowVertexId] = bData;
					}
					else if (point->second.endChar != endChar)
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
					/*std::cerr << "[";
					for (size_t i = 0; i < kt->second.branchIds.size(); ++i)
					{
					std::cerr << kt->second.branchIds[i] << ",";
					}
					std::cerr << "]";*/
				}
			}
			//std::cerr << std::endl;
			return !bulges.empty();
		}

		void GoThroughBubbles(size_t minBlockSize, size_t maxBranchSize, std::vector<std::vector<uint32_t> > & bubbleCount)
		{
			for (uint64_t vertexId = 0; vertexId < storage_.GetVerticesNumber(); vertexId++)
			{
				
			}
		}

		size_t minBlockSize_;
		size_t maxBranchSize_;		
		const EdgeStorage & storage_;
		std::vector<std::vector<bool> > mark_;		
	};
}

#endif