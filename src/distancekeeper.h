#ifndef _DISTANCE_KEEPER_H_
#define _DISTANCE_KEEPER_H_

#include "junctionstorage.h"

namespace Sibelia
{
	class DistanceKeeper
	{
	public:
		DistanceKeeper(int64_t vertices) : vertices_(vertices), NO_DISTANCE(INT32_MAX)
		{
			distance_.assign(vertices_ * 2, NO_DISTANCE);
		}

		bool IsSet(int64_t v) const
		{
			return distance_[v + vertices_] != NO_DISTANCE;
		}

		void Set(int64_t v, int64_t d)
		{
			distance_[v + vertices_] = d;
		}

		int64_t Get(int64_t v) const
		{
			assert(distance_[v + vertices_] != NO_DISTANCE);
			return distance_[v + vertices_];
		}

		void Unset(int64_t v)
		{
			distance_[v + vertices_] = NO_DISTANCE;
		}

	private:
		int64_t vertices_;
		int32_t NO_DISTANCE;
		std::vector<int32_t> distance_;
	};
}

#endif
