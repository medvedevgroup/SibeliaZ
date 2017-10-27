#ifndef _FORBIDDEN_H_
#define _FORBIDDEN_H_

#include "junctionstorage.h"

namespace Sibelia
{
	class Forbidden
	{
	public:
		Forbidden(const JunctionStorage & storage): vertices_(storage.GetVerticesNumber()), forbidden_(vertices_ * 2)
		{
			for (auto it = forbidden_.begin(); it != forbidden_.end(); it++)
			{
				std::fill(it->value, it->value + 4, 0);
			}
		}

		void Add(const Edge & e)
		{			
			if (e.GetChar() != 'N')
			{
				Edge er = e.Reverse();
				forbidden_[e.GetStartVertex() + vertices_].value[TwoPaCo::DnaChar::MakeUpChar(e.GetChar())] = true;
				forbidden_[er.GetStartVertex() + vertices_].value[TwoPaCo::DnaChar::MakeUpChar(er.GetChar())] = true;
			}			
		}

		bool IsForbidden(const Edge & e)
		{
			if (e.GetChar() == 'N')
			{
				return false;
			}

	
			bool ret = forbidden_[e.GetStartVertex() + vertices_].value[TwoPaCo::DnaChar::MakeUpChar(e.GetChar())];
			return ret;
		}

	private:
		int64_t vertices_;

		struct BoolArray
		{
			bool value[4];
		};

		std::vector<BoolArray> forbidden_;
		
	};
}
#endif // !_FORBIDDEN_H_

