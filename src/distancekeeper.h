#ifndef _DISTANCE_KEEPER_H_
#define _DISTANCE_KEEPER_H_

#include "junctionstorage.h"

namespace Sibelia
{
	class DistanceKeeper
	{
	public:
		DistanceKeeper(int64_t vertices) : vertices_(vertices)
		{
			distance_.assign(vertices_ * 2, false);
		}

		bool IsSet(int64_t v) const
		{
			return distance_[v + vertices_];
		}

		void Set(int64_t v)
		{
			distance_[v + vertices_] = true;
		}

		void Unset(int64_t v)
		{
			distance_[v + vertices_] = false;
		}

	private:
		int64_t vertices_;
		std::vector<bool> distance_;
	};
}

#endif
