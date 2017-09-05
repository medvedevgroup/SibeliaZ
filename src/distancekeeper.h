#ifndef _DISTANCE_KEEPER_H_
#define _DISTANCE_KEEPER_H_

#include "forbidden.h"
#include <unordered_map>

namespace Sibelia
{
	class DistanceKeeper
	{
	public:
		DistanceKeeper(int64_t vertices) : vertices_(vertices)
		{
		}

		bool IsSet(int64_t v) const
		{
			return distance_.count(v) > 0;
		}

		void Set(int64_t v, int64_t d)
		{
			distance_[v] = d;
		}

		int64_t Get(int64_t v) const
		{
			assert(IsSet(v));
			return distance_.find(v)->second;
		}

		void Unset(int64_t v)
		{
			distance_.erase(v);
		}

	private:
		int64_t vertices_;
		std::unordered_map<int64_t, int32_t> distance_;		
	};
}

#endif
