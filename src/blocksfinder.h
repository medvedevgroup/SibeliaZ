#ifndef _TRASERVAL_H_
#define _TRAVERSAL_H_

#include "vertexstorage.h"

namespace Sibelia
{
	class BlocksFinder
	{
	public:

		BlocksFinder(const VertexStorage storage, size_t minBlockSize, size_t maxBranchSize) : storage_(storage),
			minBlockSize_(minBlockSize), maxBranchSize_(maxBranchSize)
		{

		}

	private:



		size_t minBlockSize_;
		size_t maxBranchSize_;		
		const VertexStorage & storage_;
	};
}

#endif