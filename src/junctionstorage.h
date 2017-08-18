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
		Edge() : startVertex_(INT64_MAX), endVertex_(INT64_MAX) {}

		Edge(int64_t startVertex, int64_t endVertex, char ch, char revCh, int64_t length) :
			startVertex_(startVertex), endVertex_(endVertex), ch_(ch), revCh_(revCh), length_(length)
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

		char GetChar() const
		{
			return ch_;
		}

		int64_t GetLength() const
		{
			return length_;
		}

		Edge Reverse() const
		{
			return Edge(-endVertex_, -startVertex_, revCh_, ch_, length_);
		}

		char GetRevChar() const
		{
			return revCh_;
		}

		bool operator < (const Edge & e) const
		{
			if (startVertex_ != e.startVertex_)
			{
				return startVertex_ < e.startVertex_;
			}

			if (endVertex_ != e.endVertex_)
			{
				return endVertex_ < e.endVertex_;
			}

			if (ch_ != e.ch_)
			{
				return ch_ < e.ch_;
			}

			return false;
		}

		bool Valid() const
		{
			return startVertex_ != INT64_MAX;
		}

		bool operator == (const Edge & e) const
		{
			return startVertex_ == e.startVertex_ && endVertex_ == e.endVertex_ && ch_ == e.ch_;
		}

		bool operator != (const Edge & e) const
		{
			return !(*this == e);
		}

	private:
		int64_t startVertex_;
		int64_t endVertex_;
		char ch_;
		char revCh_;
		int64_t length_;
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
			JunctionIterator() : idx_(0)
			{

			}

			bool IsPositiveStrand() const
			{
				return chrId_ > 0;
			}

			int64_t GetVertexId(const JunctionStorage * storage_) const
			{
				return IsPositiveStrand() ? storage_->posChr_[GetChrId()][idx_].id : -storage_->posChr_[GetChrId()][idx_].id;
			}

			int64_t GetPosition(const JunctionStorage * storage_) const
			{
				if (IsPositiveStrand())
				{
					return storage_->posChr_[GetChrId()][idx_].pos;
				}

				return storage_->posChr_[GetChrId()][idx_].pos + storage_->k_;
			}

			JunctionIterator Reverse()
			{
				return JunctionIterator(GetChrId(), idx_, !IsPositiveStrand());
			}

			char GetChar(const JunctionStorage * storage_) const
			{
				int64_t pos = storage_->posChr_[GetChrId()][idx_].pos;
				if (IsPositiveStrand())
				{
					return storage_->sequence_[GetChrId()][pos + storage_->k_];
				}

				return TwoPaCo::DnaChar::ReverseChar(storage_->sequence_[GetChrId()][pos - 1]);
			}

			uint64_t GetIndex() const
			{
				return idx_;
			}

			uint64_t GetRelativeIndex(const JunctionStorage * storage_) const
			{
				if (IsPositiveStrand())
				{
					return idx_;
				}

				return storage_->posChr_[GetChrId()].size() - idx_ - 1;
			}

			uint64_t GetChrId() const
			{
				return abs(chrId_) - 1;
			}

			bool Valid(const JunctionStorage * storage_) const
			{
				return idx_ >= 0 && idx_ < storage_->posChr_[GetChrId()].size();
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

				idx_ += IsPositiveStrand() ? +step : -step;
			}

			void Dec(int64_t step = 1)
			{

				idx_ += IsPositiveStrand() ? -step : +step;
			}

			JunctionIterator(uint64_t chrId, int64_t idx, bool isPositiveStrand) : idx_(idx), chrId_(isPositiveStrand ? chrId + 1 : -(chrId + 1))
			{

			}

			friend class JunctionStorage;
			int64_t chrId_;
			int64_t idx_;
		};

		int64_t GetChrNumber() const
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

		int64_t GetChrVerticesCount(uint64_t chrId) const
		{
			return posChr_[chrId].size();
		}

		JunctionIterator GetIterator(uint64_t chrId, uint64_t idx, bool isPositiveStrand = true) const
		{
			return JunctionIterator(chrId, idx, isPositiveStrand);
		}

		JunctionIterator Begin(uint64_t chrId, bool isPositiveStrand = true) const
		{
			return JunctionIterator(chrId, 0, isPositiveStrand);
		}

		JunctionIterator End(uint64_t chrId, bool isPositiveStrand = true) const
		{
			return JunctionIterator(chrId, posChr_[chrId].size(), isPositiveStrand);
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
			return JunctionIterator(coord.chr, coord.idx, posChr_[coord.chr][coord.idx].id == vertexId);
		}

		int64_t IngoingEdgesNumber(int64_t vertexId) const
		{
			return ingoingEdge_[vertexId + GetVerticesNumber()].size();
		}

		int64_t OutgoingEdgesNumber(int64_t vertexId) const
		{
			return outgoingEdge_[vertexId + GetVerticesNumber()].size();
		}
		
		Edge IngoingEdge(int64_t vertexId, int64_t idx) const
		{
			return ingoingEdge_[vertexId + GetVerticesNumber()][idx];
		}

		Edge OutgoingEdge(int64_t vertexId, int64_t idx) const
		{
			return outgoingEdge_[vertexId + GetVerticesNumber()][idx];
		}

		void IngoingEdges(int64_t vertexId, std::vector<Edge> & list) const
		{
			list.clear();
			for (auto coord : coordinate_[abs(vertexId)])
			{
				const Vertex & now = posChr_[coord.chr][coord.idx];
				if (now.id == vertexId)
				{
					if (coord.idx > 0)
					{
						const Vertex & prev = posChr_[coord.chr][coord.idx - 1];
						char ch = sequence_[coord.chr][prev.pos + k_];
						char revCh = TwoPaCo::DnaChar::ReverseChar(sequence_[coord.chr][now.pos - 1]);
						list.push_back(Edge(prev.id, now.id, ch, revCh, now.pos - prev.pos));
					}
				}
				else
				{					
					if (coord.idx + 1 < posChr_[coord.chr].size())
					{
						const Vertex & prev = posChr_[coord.chr][coord.idx + 1];
						char ch = TwoPaCo::DnaChar::ReverseChar(sequence_[coord.chr][prev.pos - 1]);
						char revCh = sequence_[coord.chr][now.pos + k_];
						list.push_back(Edge(-prev.id, -now.id, ch, revCh, prev.pos - now.pos));
					}					
				}
			}

			std::sort(list.begin(), list.end());
			list.erase(std::unique(list.begin(), list.end()), list.end());
		}

		void OutgoingEdges(int64_t vertexId, std::vector<Edge> & list) const
		{
			list.clear();
			for (auto coord : coordinate_[abs(vertexId)])
			{
				const Vertex & now = posChr_[coord.chr][coord.idx];
				if (now.id == vertexId)
				{
					if (coord.idx + 1 < posChr_[coord.chr].size())
					{
						const Vertex & next = posChr_[coord.chr][coord.idx + 1];
						char ch = sequence_[coord.chr][now.pos + k_];
						char revCh = TwoPaCo::DnaChar::ReverseChar(sequence_[coord.chr][next.pos - 1]);
						list.push_back(Edge(now.id, next.id, ch, revCh, next.pos - now.pos));
					}
				}
				else
				{
					if (coord.idx > 0)
					{
						const Vertex & next = posChr_[coord.chr][coord.idx - 1];
						char ch = TwoPaCo::DnaChar::ReverseChar(sequence_[coord.chr][now.pos - 1]);
						char revCh = sequence_[coord.chr][now.pos + k_];
						list.push_back(Edge(-now.id, -next.id, ch, revCh, now.pos - next.pos));
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

			int64_t vertices = GetVerticesNumber();
			ingoingEdge_.resize(vertices * 2 + 1);
			outgoingEdge_.resize(vertices * 2 + 1);
			for (int64_t vertexId = -vertices + 1; vertexId < vertices; vertexId++)
			{
				IngoingEdges(vertexId, ingoingEdge_[vertexId + vertices]);
				OutgoingEdges(vertexId, outgoingEdge_[vertexId + vertices]);
			}
		}

			
		Edge RandomForwardEdge(int64_t vid) const
		{
			int64_t adjVid = vid + GetVerticesNumber();
			if (outgoingEdge_[adjVid].size() > 0)
			{
				return outgoingEdge_[adjVid][rand() % outgoingEdge_[adjVid].size()];
			}

			return Edge();
		}
		
		Edge RandomBackwardEdge(int64_t vid) const
		{
			int64_t adjVid = vid + GetVerticesNumber();
			if (ingoingEdge_[adjVid].size() > 0)
			{
				return ingoingEdge_[adjVid][rand() % ingoingEdge_[adjVid].size()];
			}

			return Edge();
		}

		JunctionStorage() {}
		JunctionStorage(const std::string & fileName, const std::string & genomesFileName, uint64_t k) : k_(k)
		{
			Init(fileName, genomesFileName);
		}


	private:
		
		struct LightEdge
		{
			int64_t vertex;
			char ch;
		};

		int64_t k_;
		std::vector<std::vector<Edge> > ingoingEdge_;
		std::vector<std::vector<Edge> > outgoingEdge_;
		std::vector<std::string> sequence_;
		std::vector<std::string> sequenceDescription_;
		std::vector<VertexVector> posChr_;
		std::vector<CoordinateVector> coordinate_;
	};
}

#endif