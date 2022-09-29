#ifndef _DISTANCE_KEEPER_H_
#define _DISTANCE_KEEPER_H_

#include "junctionstorage.h"
#include <climits>

namespace Sibelia
{
	class DistanceKeeper
	{
	public:
		DistanceKeeper(int64_t vertices) : vertices_(vertices), NOT_SET(INT_MAX)
		{
			distance_.assign(vertices_ * 2, NOT_SET);
		}

		bool IsSet(int64_t v) const
		{
			return distance_[v + vertices_] != NOT_SET;
		}

		void Set(int64_t v, int distance)
		{
			distance_[v + vertices_] = distance;
		}

		int Get(int64_t v) const
		{
			return distance_[v + vertices_];
		}

		void Unset(int64_t v)
		{
			distance_[v + vertices_] = NOT_SET;
		}

	private:
		int64_t vertices_;
		const int NOT_SET;
		std::vector<int> distance_;
	};
}

#endif
