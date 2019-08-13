#ifndef _JUNCTION_STORAGE_H_
#define _JUNCTION_STORAGE_H_

#include <set>
#include <atomic>
#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include <stdexcept>
#include <algorithm>

#include <tbb/mutex.h>

#include <streamfastaparser.h>
#include <junctionapi.h>

namespace Sibelia
{	
	using std::min;
	using std::max;

	class Edge
	{
	public:
		Edge() : startVertex_(INT64_MAX), endVertex_(INT64_MAX) {}

		Edge(int64_t startVertex, int64_t endVertex, char ch, char revCh, int64_t length, int64_t capacity) :
			startVertex_(startVertex), endVertex_(endVertex), ch_(ch), revCh_(revCh), length_(length), capacity_(capacity)
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

		int64_t GetCapacity() const
		{
			return capacity_;
		}

		Edge Reverse() const
		{
			return Edge(-endVertex_, -startVertex_, revCh_, ch_, length_, capacity_);
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
		
		void Inc()
		{
			capacity_++;
		}

	private:
		int64_t startVertex_;
		int64_t endVertex_;		
		int64_t length_;
		int64_t capacity_;
		char ch_;
		char revCh_;
	};

	class JunctionStorage
	{
	private:

		struct Vertex
		{
			int32_t id;
			uint32_t chr;
			uint32_t idx;
			uint32_t pos;
			char ch;
			char revCh;

			Vertex(const TwoPaCo::JunctionPosition & junction) : id(static_cast<int32_t>(junction.GetId())), chr(junction.GetChr()), pos(junction.GetPos())
			{

			}

			bool operator < (const Vertex & a) const
			{
				return std::make_pair(chr, idx) < std::make_pair(a.chr, a.idx);
			}
		};

		struct Position
		{
			int32_t id;
			uint32_t pos;
			bool used;

			Position(const TwoPaCo::JunctionPosition & junction) : used(false)
			{
				id = static_cast<int32_t>(junction.GetId());
				pos = junction.GetPos();
			}
		};

		typedef std::vector<Vertex> VertexVector;
		typedef std::vector<Position> PositionVector;

	public:

		class JunctionSequentialIterator
		{
		public:
			JunctionSequentialIterator() : idx_(0)
			{

			}

			bool IsPositiveStrand() const
			{
				return chrId_ > 0;
			}

			int64_t GetVertexId() const
			{
				return IsPositiveStrand() ? JunctionStorage::this_->position_[GetChrId()][idx_].id : -JunctionStorage::this_->position_[GetChrId()][idx_].id;
			}

			int64_t GetPosition() const
			{
				if (IsPositiveStrand())
				{
					return JunctionStorage::this_->position_[GetChrId()][idx_].pos;
				}

				return JunctionStorage::this_->position_[GetChrId()][idx_].pos + JunctionStorage::this_->k_;
			}

			int64_t GetAbsolutePosition() const
			{
				return JunctionStorage::this_->position_[GetChrId()][idx_].pos;
			}

			Edge OutgoingEdge() const
			{
				const Position & now = JunctionStorage::this_->position_[GetChrId()][idx_];
				if (IsPositiveStrand())
				{
					const Position & next = JunctionStorage::this_->position_[GetChrId()][idx_ + 1];
					char ch = JunctionStorage::this_->sequence_[GetChrId()][now.pos + JunctionStorage::this_->k_];
					char revCh = TwoPaCo::DnaChar::ReverseChar(JunctionStorage::this_->sequence_[GetChrId()][next.pos - 1]);
					return Edge(now.id, next.id, ch, revCh, next.pos - now.pos, 1);
				}
				else
				{
					const Position & next = JunctionStorage::this_->position_[GetChrId()][idx_ - 1];
					char ch = TwoPaCo::DnaChar::ReverseChar(JunctionStorage::this_->sequence_[GetChrId()][now.pos - 1]);
					char revCh = JunctionStorage::this_->sequence_[GetChrId()][now.pos + JunctionStorage::this_->k_];
					return Edge(-now.id, -next.id, ch, revCh, now.pos - next.pos, 1);
				}
			}

			Edge IngoingEdge() const
			{
				const Position & now = JunctionStorage::this_->position_[GetChrId()][idx_];
				if (IsPositiveStrand())
				{
					const Position & prev = JunctionStorage::this_->position_[GetChrId()][idx_ - 1];
					char ch = JunctionStorage::this_->sequence_[GetChrId()][prev.pos + JunctionStorage::this_->k_];
					char revCh = TwoPaCo::DnaChar::ReverseChar(JunctionStorage::this_->sequence_[GetChrId()][now.pos - 1]);
					return Edge(prev.id, now.id, ch, revCh, now.pos - prev.pos, 1);
				}
				else
				{
					const Position & prev = JunctionStorage::this_->position_[GetChrId()][idx_ + 1];
					char ch = TwoPaCo::DnaChar::ReverseChar(JunctionStorage::this_->sequence_[GetChrId()][prev.pos - 1]);
					char revCh = JunctionStorage::this_->sequence_[GetChrId()][now.pos + JunctionStorage::this_->k_];
					return Edge(-prev.id, -now.id, ch, revCh, prev.pos - now.pos, 1);
				}
			}

			JunctionSequentialIterator Reverse()
			{
				return JunctionSequentialIterator(GetChrId(), idx_, !IsPositiveStrand());
			}

			char GetChar() const
			{
				int64_t pos = JunctionStorage::this_->position_[GetChrId()][idx_].pos;
				if (IsPositiveStrand())
				{
					return JunctionStorage::this_->sequence_[GetChrId()][pos + JunctionStorage::this_->k_];
				}

				return TwoPaCo::DnaChar::ReverseChar(JunctionStorage::this_->sequence_[GetChrId()][pos - 1]);
			}

			uint64_t GetIndex() const
			{
				return idx_;
			}

			uint64_t GetRelativeIndex() const
			{
				if (IsPositiveStrand())
				{
					return idx_;
				}

				return JunctionStorage::this_->position_[GetChrId()].size() - idx_ - 1;
			}

			uint64_t GetChrId() const
			{
				return abs(chrId_) - 1;
			}

			bool Valid() const
			{
				return idx_ >= 0 && size_t(idx_) < JunctionStorage::this_->position_[GetChrId()].size();
			}

			bool IsUsed() const
			{
				if (IsPositiveStrand())
				{
					return JunctionStorage::this_->position_[GetChrId()][idx_].used;
				}
				
				if (idx_ > 0)
				{
					return JunctionStorage::this_->position_[GetChrId()][idx_ - 1].used;
				}

				return false;
			}

			void MarkUsed() const
			{
				if (IsPositiveStrand())
				{
					JunctionStorage::this_->position_[GetChrId()][idx_].used = true;
				}
				else if (idx_ > 0)
				{
					JunctionStorage::this_->position_[GetChrId()][idx_ - 1].used = true;
				}
			}

			JunctionSequentialIterator& operator++ ()
			{
				Inc();
				return *this;
			}

			JunctionSequentialIterator operator++ (int)
			{
				JunctionSequentialIterator ret(*this);
				Inc();
				return ret;
			}

			JunctionSequentialIterator operator + (size_t step) const
			{
				JunctionSequentialIterator ret(*this);
				ret.Inc(step);
				return ret;
			}

			JunctionSequentialIterator operator - (size_t step) const
			{
				JunctionSequentialIterator ret(*this);
				ret.Dec(step);
				return ret;
			}

			JunctionSequentialIterator& operator-- ()
			{
				Dec();
				return *this;
			}

			JunctionSequentialIterator operator-- (int)
			{
				JunctionSequentialIterator ret(*this);
				Dec();
				return ret;
			}

			JunctionSequentialIterator Next() const
			{
				JunctionSequentialIterator ret(*this);
				return ++ret;
			}

			JunctionSequentialIterator Prev() const
			{
				JunctionSequentialIterator ret(*this);
				return --ret;
			}

			bool operator < (const JunctionSequentialIterator & arg) const
			{
				if (GetChrId() != arg.GetChrId())
				{
					return GetChrId() < arg.GetChrId();
				}

				return GetIndex() < arg.GetIndex();
			}

			bool operator == (const JunctionSequentialIterator & arg) const
			{
				return this->chrId_ == arg.chrId_ && this->idx_ == arg.idx_;
			}

			bool operator != (const JunctionSequentialIterator & arg) const
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

			JunctionSequentialIterator(int64_t chrId, int64_t idx, bool isPositiveStrand) : idx_(idx), chrId_(isPositiveStrand ? chrId + 1 : -(chrId + 1))
			{

			}

			friend class JunctionStorage;
			int64_t chrId_;
			int64_t idx_;
		};



		class JunctionIterator
		{
		public:
			JunctionIterator() : vid_(0)
			{

			}

			bool IsPositiveStrand() const
			{
				return JunctionStorage::this_->vertex_[abs(vid_)][iidx_].id == vid_;
			}

			int64_t GetVertexId() const
			{
				return vid_;
			}

			int64_t GetPosition() const
			{
				return JunctionStorage::this_->vertex_[abs(vid_)][iidx_].pos;
			}

			char GetChar() const
			{
				if (IsPositiveStrand())
				{
					return JunctionStorage::this_->vertex_[abs(vid_)][iidx_].ch;
				}

				return JunctionStorage::this_->vertex_[abs(vid_)][iidx_].revCh;
			}

			JunctionSequentialIterator SequentialIterator() const
			{
				return JunctionSequentialIterator(GetChrId(), GetIndex(), IsPositiveStrand());
			}

			uint64_t GetIndex() const
			{
				return JunctionStorage::this_->vertex_[abs(vid_)][iidx_].idx;
			}

			uint64_t GetRelativeIndex() const
			{
				if (IsPositiveStrand())
				{
					return JunctionStorage::this_->vertex_[abs(vid_)][iidx_].idx;;
				}

				return JunctionStorage::this_->position_[GetChrId()].size() - JunctionStorage::this_->vertex_[abs(vid_)][iidx_].idx; -1;
			}

			uint64_t GetChrId() const
			{
				return JunctionStorage::this_->vertex_[abs(vid_)][iidx_].chr;
			}

			bool Valid() const
			{
				return iidx_ < JunctionStorage::this_->vertex_[abs(vid_)].size();
			}

			size_t InstancesCount() const
			{
				return JunctionStorage::this_->vertex_[abs(vid_)].size();
			}

			JunctionIterator operator + (size_t inc) const
			{
				return JunctionIterator(vid_, iidx_ + inc);
			}

			JunctionIterator& operator++ ()
			{
				++iidx_;
				return *this;
			}

			JunctionIterator operator++ (int)
			{
				JunctionIterator ret(*this);
				++iidx_;
				return ret;
			}

			bool operator < (const JunctionIterator & arg) const
			{
				if (GetChrId() != arg.GetChrId())
				{
					return GetChrId() < arg.GetChrId();
				}

				return GetIndex() < arg.GetIndex();
			}

			bool operator == (const JunctionIterator & arg) const
			{
				return this->vid_ == arg.vid_ && this->iidx_ == arg.iidx_;
			}

			bool operator != (const JunctionIterator & arg) const
			{
				return !(*this == arg);
			}

			JunctionIterator(int64_t vid) : iidx_(0), vid_(vid)
			{
			}

		private:

			JunctionIterator(int64_t vid, size_t iidx) : iidx_(iidx), vid_(vid)
			{
			}

			friend class JunctionStorage;
			size_t iidx_;
			int64_t vid_;

		};

		int64_t GetChrNumber() const
		{
			return position_.size();
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
			return position_[chrId].size();
		}

		JunctionSequentialIterator GetIterator(uint64_t chrId, uint64_t idx, bool isPositiveStrand = true) const
		{
			return JunctionSequentialIterator(chrId, idx, isPositiveStrand);
		}

		JunctionSequentialIterator Begin(uint64_t chrId, bool isPositiveStrand = true) const
		{
			return JunctionSequentialIterator(chrId, 0, isPositiveStrand);
		}

		JunctionSequentialIterator End(uint64_t chrId, bool isPositiveStrand = true) const
		{
			return JunctionSequentialIterator(chrId, position_[chrId].size(), isPositiveStrand);
		}

		int64_t GetVerticesNumber() const
		{
			return vertex_.size();
		}

		uint64_t GetInstancesCount(int64_t vertexId) const
		{
			return vertex_[abs(vertexId)].size();
		}

		const std::string& GetSequence(size_t idx) const
		{
			return sequence_[idx];
		}

		void Init(const std::string & inFileName, const std::vector<std::string> & genomesFileName, int64_t threads, int64_t abundanceThreshold, int64_t loopThreshold)
		{
			this_ = this;
			std::vector<size_t> abundance;
			{
				TwoPaCo::JunctionPositionReader reader(inFileName);
				for (TwoPaCo::JunctionPosition junction; reader.NextJunctionPosition(junction);)
				{
					if (junction.GetChr() >= position_.size())
					{
						position_.push_back(std::vector<Position>());
					}

					size_t absId = abs(junction.GetId());
					while (absId >= vertex_.size())
					{
						vertex_.push_back(VertexVector());
						abundance.push_back(0);
					}

					++abundance[absId];
				}
			}


			{
				size_t chr = 0;
				uint32_t idx = 0;
				TwoPaCo::JunctionPositionReader reader(inFileName);
				for (TwoPaCo::JunctionPosition junction; reader.NextJunctionPosition(junction);)
				{
					if (junction.GetChr() > chr)
					{
						chr++;
						idx = 0;
					}

					size_t absId = abs(junction.GetId());
					if (abundance[absId] < size_t(abundanceThreshold) )
					{
						position_[junction.GetChr()].push_back(Position(junction));
						vertex_[absId].push_back(Vertex(junction));
						vertex_[absId].back().idx = idx++;
					}
				}
			}


			size_t record = 0;
			sequence_.resize(position_.size());
			for (const auto & fastaFileName : genomesFileName)
			{
				for (TwoPaCo::StreamFastaParser parser(fastaFileName); parser.ReadRecord(); record++)
				{
					sequenceDescription_.push_back(parser.GetCurrentHeader());
					sequenceId_[parser.GetCurrentHeader()] = sequenceDescription_.size() - 1;
					for (char ch; parser.GetChar(ch); )
					{
						sequence_[record].push_back(ch);
					}
				}
			}

			for (size_t i = 0; i < vertex_.size(); i++)
			{
				for (size_t j = 0; j < vertex_[i].size(); j++)
				{
					int64_t chr = vertex_[i][j].chr;
					int64_t pos_ = vertex_[i][j].pos;
					vertex_[i][j].ch = sequence_[chr][pos_ + JunctionStorage::this_->k_];
					vertex_[i][j].revCh = pos_ > 0 ? TwoPaCo::DnaChar::ReverseChar(sequence_[chr][pos_ - 1]) : 'N';
				}
			}

			for (auto & v : vertex_)
			{
				std::sort(v.begin(), v.end());
			}
		}

		JunctionStorage() {}
		JunctionStorage(const std::string & fileName, const std::vector<std::string> & genomesFileName, uint64_t k, int64_t threads, int64_t abundanceThreshold, int64_t loopThreshold) : k_(k)
		{
			Init(fileName, genomesFileName, threads, abundanceThreshold, loopThreshold);
		}

		bool IsSequencePresent(const std::string & str) const
		{
			return sequenceId_.count(str) > 0;
		}

		size_t GetSequenceId(const std::string & str) const
		{
			return sequenceId_.find(str)->second;
		}

		void DebugUsed() const
		{
			for (size_t i = 0; i < position_.size(); i++)
			{
				for (size_t j = 0; j < position_[i].size(); j++)
				{
					std::cout << (position_[i][j].used ? 1 : 0);
				}

				std::cout << std::endl;
			}
		}

	private:

		struct LightEdge
		{
			int64_t vertex;
			char ch;
		};


		int64_t k_;
		int64_t mutexBits_;
		std::map<std::string, size_t> sequenceId_;
		std::vector<std::string> sequence_;
		std::vector<std::string> sequenceDescription_;		
		std::vector<VertexVector> vertex_;
		std::vector<std::vector<Position> > position_;
		static JunctionStorage * this_;
	};
}

#endif
