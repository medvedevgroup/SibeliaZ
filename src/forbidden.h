#ifndef _FORBIDDEN_H_

#include "junctionstorage.h"

namespace Sibelia
{
	class Forbidden
	{
	public:
		Forbidden(const JunctionStorage & storage): vertices_(storage.GetVerticesNumber()), forbidden_(vertices_ * 2)
		{

		}

		void Add(const Edge & e)
		{
			int64_t v = e.GetStartVertex() + vertices_;
			if (std::find(forbidden_[v].begin(), forbidden_[v].end(), e.GetChar()) == forbidden_[v].end())
			{
				forbidden_[v].push_back(e.GetChar());
				forbidden_[-e.GetEndVertex() + vertices_].push_back(e.GetChar());
			}
		}

		bool Notin(const Edge & e) const
		{
			int64_t v = e.GetStartVertex() + vertices_;
			return std::find(forbidden_[v].begin(), forbidden_[v].end(), e.GetChar()) == forbidden_[v].end();
		}

	private:
		int64_t vertices_;
		std::vector<std::string> forbidden_;
	};
}
#endif // !_FORBIDDEN_H_

