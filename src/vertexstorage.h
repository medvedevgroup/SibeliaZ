#ifndef _VERTEX_STORAGE_H_
#define _VERTEX_STORAGE_H_

#include <vector>
#include <cstdint>

#include "junctionapi.h"

namespace Sibelia
{
	class VertexStorage
	{
	private:
		struct Vertex
		{
			int64_t id;
			uint64_t pos;
		};

		typedef std::vector<Vertex> VertexVector;

	public:		

		class StrandIterator
		{
		public:
			StrandIterator() : idx_(0), chr_(0)
			{

			}

			bool IsPositiveStrand() const
			{
				return (*chr_)[idx_].id > 0;
			}

			int64_t GetId() const
			{
				return IsPositiveStrand() ? (*chr_)[idx_].id : -(*chr_)[idx_].id;
			}

			int64_t GetPosition() const
			{
				return (*chr_)[idx_].pos;
			}			

			bool CanInc() const
			{
				if (IsPositiveStrand())
				{
					return idx_ < chr_->size();
				}
				else
				{
					return idx_ > 0;
				}
			}

			bool CanDec() const
			{
				if (IsPositiveStrand())
				{
					return idx_ > 0;
				}
				else
				{
					return idx_ < chr_->size();
				}
			}

			StrandIterator& operator++ ()
			{
				Inc();
				return *this;
			}
			
			StrandIterator operator++ (int)
			{
				StrandIterator ret(*this);
				Inc();
				return ret;
			}

			StrandIterator& operator-- ()
			{
				Dec();
				return *this;
			}

			StrandIterator operator-- (int)
			{				
				StrandIterator ret(*this);
				Dec();
				return ret;
			}

		private:

			void Inc()
			{
				if (CanInc())
				{
					idx_ += IsPositiveStrand() ? +1 : -1;
				}
			}

			void Dec()
			{
				if (CanDec())
				{
					idx_ += IsPositiveStrand() ? -1 : +1;
				}
			}

			StrandIterator(int64_t idx, const VertexVector * chr) : idx_(idx),  chr_(chr)
			{

			}

			friend class VertexStorage;
			int64_t idx_;
			const VertexVector * chr_;
		};

		uint64_t VerticesNumber() const
		{
			return coordinate_.size();
		}

		void Positions(uint64_t vertexId, std::vector<StrandIterator> & pos) const
		{
			pos.clear();
			for (auto coord : coordinate_)
			{				
				pos.push_back(StrandIterator(coord.idx, &posChr_[coord.chr]));
			}
		}

		void Init()
		{

		}

	private:

		struct Coordinate
		{
			uint32_t chr;
			uint32_t idx;
		};

		std::vector<VertexVector> posChr_;
		std::vector<Coordinate> coordinate_;
	};
}

#endif