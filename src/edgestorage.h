#ifndef _EDGE_STORAGE_H_
#define _EDGE_STORAGE_H_

#include <vector>
#include <cstdint>

#include "junctionapi.h"

namespace Sibelia
{
	class EdgeStorage
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

		class EdgeIterator
		{
		public:
			EdgeIterator() : idx_(0), chr_(0)
			{

			}

			bool IsPositiveStrand() const
			{
				return (*chr_)[idx_].id > 0;
			}

			int64_t GetStartVertexId() const
			{
				return IsPositiveStrand() ? (*chr_)[idx_].id : -(*chr_)[idx_].id;
			}

			int64_t GetEndVertexId() const
			{
				return IsPositiveStrand() ? (*chr_)[idx_ + 1].id : -(*chr_)[idx_ - 1].id;
			}

			char GetChar() const;

			int64_t GetPosition() const
			{
				return (*chr_)[idx_].pos;
			}

			uint64_t GetIdx() const
			{
				return idx_;
			}

			uint64_t GetChrId() const
			{
				return chrId_;
			}			

			bool CanInc() const
			{
				if (IsPositiveStrand())
				{
					return idx_ < chr_->size() - 1;
				}
				else
				{
					return idx_ > 1;
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

			EdgeIterator& operator++ ()
			{
				Inc();
				return *this;
			}
			
			EdgeIterator operator++ (int)
			{
				EdgeIterator ret(*this);
				Inc();
				return ret;
			}

			EdgeIterator& operator-- ()
			{
				Dec();
				return *this;
			}

			EdgeIterator operator-- (int)
			{				
				EdgeIterator ret(*this);
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

			EdgeIterator(int64_t idx, const VertexVector * chr, uint64_t chrId) : idx_(idx), chr_(chr), chrId_(chrId)
			{

			}

			friend class EdgeStorage;
			int64_t idx_;			
			const VertexVector * chr_;
			uint64_t chrId_;
		};

		uint64_t GetChrNumber() const
		{
			return posChr_.size();
		}

		uint64_t GetChrEdgeCount(uint64_t chrId) const
		{
			return posChr_[chrId].size() - 1;
		}

		EdgeIterator GetIterator(uint64_t chrId, uint64_t idx) const
		{
			return EdgeIterator(idx, &posChr_[chrId], chrId);
		}

		uint64_t GetVerticesNumber() const
		{
			return coordinate_.size();
		}

		uint64_t GetOutgoingEdgesCount(uint64_t vertexId) const
		{
			return coordinate_[vertexId].size();
		}

		EdgeIterator GetOutgoingEdge(uint64_t vertexId, uint64_t idx) const
		{
			auto coord = coordinate_[vertexId][idx];
			return EdgeIterator(coord.idx, &posChr_[coord.chr], coord.chr);
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

		EdgeStorage() {}
		EdgeStorage(const std::string & fileName)
		{
			Init(fileName);
		}

	private:
		
		std::vector<VertexVector> posChr_;
		std::vector<CoordinateVector> coordinate_;
	};
}

#endif