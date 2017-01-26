#ifndef _TRASERVAL_H_
#define _TRAVERSAL_H_

#include "edgestorage.h"

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

		void CountBubbles(size_t minBlockSize, size_t maxBranchSize, std::vector<std::vector<uint32_t> > & bubbleCount)
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