#ifndef _EDGE_STORAGE_H_
#define _EDGE_STORAGE_H_

#include <string>
#include <vector>
#include <cstdint>

#include "junctionapi.h"
#include "streamfastaparser.h"

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
			EdgeIterator() : idx_(0), storage_(0)
			{

			}

			bool IsPositiveStrand() const
			{
				return storage_->posChr_[chrId_][idx_].id > 0;
			}

			int64_t GetStartVertexId() const
			{
				return IsPositiveStrand() ? storage_->posChr_[chrId_][idx_].id : -storage_->posChr_[chrId_][idx_].id;
			}

			int64_t GetEndVertexId() const
			{
				return IsPositiveStrand() ? storage_->posChr_[chrId_][idx_ + 1].id : storage_->posChr_[chrId_][idx_ - 1].id;
			}

			char GetChar() const
			{
				if (IsPositiveStrand())
				{
					return storage_->seq_[chrId_][idx_ + storage_->k_];
				}

				return TwoPaCo::DnaChar::ReverseChar(storage_->seq_[chrId_][idx_ - storage_->k_]);
			}

			int64_t GetPosition() const
			{
				return storage_->posChr_[chrId_][idx_].pos;
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
					return idx_ < storage_->posChr_[chrId_].size() - 1;
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
					return idx_ < storage_->posChr_[chrId_].size();
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

			EdgeIterator(int64_t idx, const EdgeStorage * storage, uint64_t chrId) : idx_(idx), storage_(storage), chrId_(chrId)
			{

			}

			friend class EdgeStorage;
			int64_t idx_;			
			const EdgeStorage * storage_;
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
			return EdgeIterator(idx, this, chrId);
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
			return EdgeIterator(coord.idx, this, coord.chr);
		}

		void Init(const std::string & inFileName, const std::string & genomesFileName)
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

			size_t record = 0;
			seq_.resize(posChr_.size());
			for (TwoPaCo::StreamFastaParser parser(genomesFileName); parser.ReadRecord(); record++)
			{
				for (char ch; parser.GetChar(ch); )
				{
					seq_[record].push_back(ch);
				}
			}
		}

		EdgeStorage() {}
		EdgeStorage(const std::string & fileName, const std::string & genomesFileName, uint64_t k) : k_(k)
		{
			Init(fileName, genomesFileName);
		}

	private:
		
		uint64_t k_;
		std::vector<std::string> seq_;
		std::vector<VertexVector> posChr_;
		std::vector<CoordinateVector> coordinate_;
	};
}

#endif