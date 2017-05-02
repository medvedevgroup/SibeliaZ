#ifndef _JUNCTION_STORAGE_H_
#define _JUNCTION_STORAGE_H_

#include <string>
#include <vector>
#include <cstdint>

#include "junctionapi.h"
#include "streamfastaparser.h"

namespace Sibelia
{
	class Edge
	{
	public:
		Edge() {}

		Edge(int64_t startVertex, int64_t endVertex) : startVertex_(startVertex), endVertex_(endVertex)
		{

		}

		int64_t GetStartVertex() const
		{
			return startVertex_;
		}

		int64_t GetEndVertex() const
		{
			return endVertex_;
		}

		bool operator < (const Edge & e) const
		{
			return std::make_pair(startVertex_, endVertex_) < std::make_pair(e.startVertex_, e.endVertex_);
		}

		Edge Reverse() const
		{
			return Edge(-endVertex_, -startVertex_);
		}

	private:
		int64_t startVertex_;
		int64_t endVertex_;
	};

	class JunctionStorage
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

		class JunctionIterator
		{
		public:
			JunctionIterator() : idx_(0), storage_(0)
			{

			}

			bool IsPositiveStrand() const
			{
				return isPositiveStrand_;
			}
			
			int64_t operator * () const
			{
				return GetVertexId();
			}

			/*
			Edge operator -> () const
			{
				return Edge(GetStartVertexId(), GetEndVertexId());
			}

			int64_t GetEndVertexId() const
			{
				return IsPositiveStrand() ? storage_->posChr_[chrId_][idx_ + 1].id : -storage_->posChr_[chrId_][idx_ - 1].id;
			}

			char GetChar() const
			{
				if (IsPositiveStrand())
				{
					return storage_->sequence_[chrId_][idx_ + storage_->k_];
				}

				return TwoPaCo::DnaChar::ReverseChar(storage_->sequence_[chrId_][idx_ - storage_->k_]);
			}			
		
			int64_t GetEndPosition() const
			{
				return (++EdgeIterator(*this)).GetStartPosition();
			}*/


			int64_t GetVertexId() const
			{
				return IsPositiveStrand() ? storage_->posChr_[chrId_][idx_].id : -storage_->posChr_[chrId_][idx_].id;
			}

			int64_t GetPosition() const
			{
				if (IsPositiveStrand())
				{
					return storage_->posChr_[chrId_][idx_].pos;
				}

				return storage_->posChr_[chrId_][idx_].pos + storage_->k_;
			}

			uint64_t GetIndex() const
			{
				return idx_;
			}

			uint64_t GetChrId() const
			{
				return chrId_;
			}						

			bool Valid() const
			{
				return idx_ >= 0 && idx_ < storage_->posChr_[chrId_].size();
			}

			JunctionIterator& operator++ ()
			{
				Inc();
				return *this;
			}
			
			JunctionIterator operator++ (int)
			{
				JunctionIterator ret(*this);
				Inc();
				return ret;
			}

			JunctionIterator operator + (size_t step) const
			{
				JunctionIterator ret(*this);
				ret.Inc(step);
				return ret;
			}

			JunctionIterator operator - (size_t step) const
			{
				JunctionIterator ret(*this);
				ret.Dec(step);
				return ret;
			}

			JunctionIterator& operator-- ()
			{
				Dec();
				return *this;
			}

			JunctionIterator operator-- (int)
			{				
				JunctionIterator ret(*this);
				Dec();
				return ret;
			}

			bool operator == (const JunctionIterator & arg) const
			{
				return this->chrId_ == arg.chrId_ && this->idx_ == arg.idx_;
			}

			bool operator != (const JunctionIterator & arg) const
			{
				return !(*this == arg);
			}

		private:

			void Inc(int64_t step = 1)
			{
				if (Valid())
				{
					idx_ += IsPositiveStrand() ? +step : -step;
				}
			}

			void Dec(int64_t step = 1)
			{
				if (Valid())
				{
					idx_ += IsPositiveStrand() ? -step : +step;
				}
			}

			JunctionIterator(const JunctionStorage * storage, uint64_t chrId, int64_t idx, bool isPositiveStrand) : storage_(storage), idx_(idx), chrId_(chrId), isPositiveStrand_(isPositiveStrand)
			{

			}

			friend class JunctionStorage;
			const JunctionStorage * storage_;
			uint64_t chrId_;
			int64_t idx_;
			bool isPositiveStrand_;
			
		};

		uint64_t GetChrNumber() const
		{
			return posChr_.size();
		}

		const std::string& GetChrSequence(uint64_t idx) const
		{
			return sequence_[idx];
		}

		const std::string& GetChrDescription(uint64_t idx) const
		{
			return sequenceDescription_[idx];
		}

		uint64_t GetChrVerticesCount(uint64_t chrId) const
		{
			return posChr_[chrId].size();
		}

		JunctionIterator GetIterator(uint64_t chrId, uint64_t idx, bool isPositiveStrand = true) const
		{
			if (isPositiveStrand)
			{
				return JunctionIterator(this, chrId, idx, isPositiveStrand);
			}

			return JunctionIterator(this, chrId, posChr_[chrId].size() - idx - 1, isPositiveStrand);
		}

		JunctionIterator Begin(uint64_t chrId, bool isPositiveStrand = true) const
		{
			return JunctionIterator(this, chrId, 0, isPositiveStrand);
		}

		JunctionIterator End(uint64_t chrId, bool isPositiveStrand = true) const
		{
			return JunctionIterator(this, chrId, posChr_[chrId].size() - 1, isPositiveStrand);
		}

		int64_t GetVerticesNumber() const
		{
			return coordinate_.size();
		}
				
		uint64_t GetInstancesCount(int64_t vertexId) const
		{
			return coordinate_[abs(vertexId)].size();
		}

		JunctionIterator GetJunctionInstance(int64_t vertexId, uint64_t n) const
		{
			auto coord = coordinate_[abs(vertexId)][n];			
			return JunctionIterator(this, coord.chr, coord.idx, posChr_[coord.chr][coord.idx].id == vertexId);
		}

		void PredecessorsList(int64_t vertexId, std::vector<int64_t> & list) const
		{
			list.clear();
			for (auto coord : coordinate_[abs(vertexId)])
			{
				if (posChr_[coord.chr][coord.idx].id == vertexId)
				{
					if (coord.idx > 0)
					{
						list.push_back(posChr_[coord.chr][coord.idx - 1].id);
					}
				}
				else
				{
					if (coord.idx + 1 < posChr_[coord.chr].size())
					{
						list.push_back(-posChr_[coord.chr][coord.idx + 1].id);
					}					
				}
			}

			std::sort(list.begin(), list.end());
			list.erase(std::unique(list.begin(), list.end()), list.end());
		}

		void SuccessorsList(int64_t vertexId, std::vector<int64_t> & list) const
		{
			list.clear();
			for (auto coord : coordinate_[abs(vertexId)])
			{
				if (posChr_[coord.chr][coord.idx].id == vertexId)
				{
					if (coord.idx + 1 < posChr_[coord.chr].size())
					{
						list.push_back(posChr_[coord.chr][coord.idx + 1].id);
					}
				}
				else
				{
					if (coord.idx > 0)
					{
						list.push_back(-posChr_[coord.chr][coord.idx - 1].id);
					}
				}
			}

			std::sort(list.begin(), list.end());
			list.erase(std::unique(list.begin(), list.end()), list.end());
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
				size_t absId = abs(junction.GetId());

				while (absId >= coordinate_.size())
				{
					coordinate_.push_back(CoordinateVector());
				}

				coordinate_[absId].push_back(Coordinate(junction.GetChr(), posChr_[junction.GetChr()].size() - 1));
			}

			size_t record = 0;
			sequence_.resize(posChr_.size());
			for (TwoPaCo::StreamFastaParser parser(genomesFileName); parser.ReadRecord(); record++)
			{
				sequenceDescription_.push_back(parser.GetCurrentHeader());
				for (char ch; parser.GetChar(ch); )
				{
					sequence_[record].push_back(ch);
				}
			}
		}

		JunctionStorage() {}
		JunctionStorage(const std::string & fileName, const std::string & genomesFileName, uint64_t k) : k_(k)
		{
			Init(fileName, genomesFileName);
		}

	private:
		
		uint64_t k_;
		std::vector<std::string> sequence_;
		std::vector<std::string> sequenceDescription_;
		std::vector<VertexVector> posChr_;
		std::vector<CoordinateVector> coordinate_;
	};
}

#endif