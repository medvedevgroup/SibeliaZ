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

			Vertex(const TwoPaCo::JunctionPosition & junction) : id(junction.GetId()), pos(junction.GetPos())
			{

			}
		};
		
		struct Coordinate
		{
			uint32_t chr;
			uint32_t idx;

			Coordinate() {}
			Coordinate(uint32_t chr, uint32_t idx) : chr(chr), idx(idx)
			{

			}
		};

		typedef std::vector<Vertex> VertexVector;
		typedef std::vector<Coordinate> CoordinateVector;

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
			for (auto coord : coordinate_[vertexId])
			{				
				pos.push_back(StrandIterator(coord.idx, &posChr_[coord.chr]));
			}
		}

		void Init(const std::string & inFileName)
		{
			TwoPaCo::JunctionPositionReader reader(inFileName);
			for (TwoPaCo::JunctionPosition junction; reader.NextJunctionPosition(junction);)
			{
				if (junction.GetChr() >= posChr_.size())
				{
					posChr_.push_back(VertexVector());
				}

				posChr_[junction.GetChr()].push_back(Vertex(junction));
				int64_t absId = abs(junction.GetId());

				while (absId >= coordinate_.size())
				{
					coordinate_.push_back(CoordinateVector());
				}

				coordinate_[absId].push_back(Coordinate(junction.GetChr(), posChr_[junction.GetChr()].size() - 1));
			}
		}

		VertexStorage() {}
		VertexStorage(const std::string & fileName)
		{
			Init(fileName);
		}

	private:
		
		std::vector<VertexVector> posChr_;
		std::vector<CoordinateVector> coordinate_;
	};
}

#endif