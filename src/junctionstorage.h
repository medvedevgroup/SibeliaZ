#ifndef _JUNCTION_STORAGE_H_
#define _JUNCTION_STORAGE_H_

#include <atomic>
#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include <stdexcept>
#include <algorithm>

#include <tbb/mutex.h>

#include "junctionapi.h"
#include "streamfastaparser.h"


namespace Sibelia
{	
	using std::min;
	using std::max;

	class Edge
	{
	public:
		Edge() : startVertex_(INT64_MAX), endVertex_(INT64_MAX) {}

		Edge(int64_t startVertex, int64_t endVertex, char ch, char revCh, int32_t length, int32_t capacity) :
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
		char ch_;
		char revCh_;		
		int32_t length_;
		int32_t capacity_;
	};

	class JunctionStorage
	{
	private:
		struct Vertex
		{
			int64_t id;
			int64_t pos;
			char ch_;
			char revCh_;

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

			int64_t GetVertexId() const
			{
				return vid_;
			}

			int64_t GetPosition() const
			{
				return pos_;
			}

			JunctionIterator Reverse()
			{
				return JunctionIterator(GetChrId(), idx_, !IsPositiveStrand());
			}

			char GetChar() const
			{
				return ch_;
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

				return JunctionStorage::this_->posChr_[GetChrId()].size() - idx_ - 1;
			}

			uint64_t GetChrId() const
			{
				return abs(chrId_) - 1;
			}

			JunctionIterator Next() const
			{
				JunctionIterator ret(*this);
				return ++ret;
			}

			JunctionIterator Prev() const
			{
				JunctionIterator ret(*this);
				return --ret;
			}

			bool Valid() const
			{
				return idx_ >= 0 && idx_ < JunctionStorage::this_->posChr_[GetChrId()].size();
			}

			bool IsUsed() const
			{
				bool ret = JunctionStorage::this_->used_[GetChrId()][idx_];
				return ret;
			}			

			void MarkUsed() const
			{
				JunctionStorage::this_->used_[GetChrId()][idx_] = true;
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
				if (Valid())
				{
					Init();
				}
			}

			void Dec(int64_t step = 1)
			{
				idx_ += IsPositiveStrand() ? -step : +step;
				if (Valid())
				{
					Init();
				}				
			}

			JunctionIterator(int64_t chrId, int64_t idx, bool isPositiveStrand) : idx_(idx), chrId_(isPositiveStrand ? chrId + 1 : -(chrId + 1))
			{
				Init();
			}
			
			void Init()
			{
				if (IsPositiveStrand()) 
				{
					vid_ = JunctionStorage::this_->posChr_[GetChrId()][idx_].id;
					pos_ = JunctionStorage::this_->posChr_[GetChrId()][idx_].pos;
					ch_ = JunctionStorage::this_->posChr_[GetChrId()][idx_].ch_;
				}
				else
				{
					vid_ = -JunctionStorage::this_->posChr_[GetChrId()][idx_].id;
					pos_ = JunctionStorage::this_->posChr_[GetChrId()][idx_].pos + JunctionStorage::this_->k_;
					ch_ = JunctionStorage::this_->posChr_[GetChrId()][idx_].revCh_;
				}
			}

			friend class JunctionStorage;
			int64_t chrId_;
			int64_t idx_;
			int64_t vid_;
			int64_t pos_;
			char ch_;

		};

		void LockRange(JunctionIterator start, JunctionIterator end, std::vector<std::vector<char > > & mutexAcquired)
		{
			do
			{
				size_t idx = MutexIdx(start.GetChrId(), start.GetIndex());
				if (!mutexAcquired[start.GetChrId()][idx])
				{
					mutex_[start.GetChrId()][idx].lock();
					mutexAcquired[start.GetChrId()][idx] = true;
				}
				
				
			} while (start++ != end);
		}

		void UnlockRange(JunctionIterator start, JunctionIterator end, std::vector<std::vector<char > > & mutexAcquired)
		{
			do
			{
				size_t idx = MutexIdx(start.GetChrId(), start.GetIndex());
				if (mutexAcquired[start.GetChrId()][idx])
				{
					mutex_[start.GetChrId()][idx].unlock();
					mutexAcquired[start.GetChrId()][idx] = false;
				}
				

			} while (start++ != end);
		}

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
						Edge newEdge(prev.id, now.id, ch, revCh, now.pos - prev.pos, 1);
						auto it = std::find(list.begin(), list.end(), newEdge);
						if (it == list.end())
						{
							list.push_back(newEdge);
						}
						else
						{
							it->Inc();
						}
					}
				}
				else
				{					
					if (coord.idx + 1 < posChr_[coord.chr].size())
					{
						const Vertex & prev = posChr_[coord.chr][coord.idx + 1];
						char ch = TwoPaCo::DnaChar::ReverseChar(sequence_[coord.chr][prev.pos - 1]);
						char revCh = sequence_[coord.chr][now.pos + k_];
						Edge newEdge(-prev.id, -now.id, ch, revCh, prev.pos - now.pos, 1);
						auto it = std::find(list.begin(), list.end(), newEdge);
						if (it == list.end())
						{
							list.push_back(newEdge);
						}
						else
						{
							it->Inc();
						}
					}					
				}
			}

			std::sort(list.begin(), list.end());
			list.erase(std::unique(list.begin(), list.end()), list.end());
		}

		size_t MutexNumber() const
		{
			return 1 << mutexBits_;
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
						Edge newEdge = Edge(now.id, next.id, ch, revCh, next.pos - now.pos, 1);
						auto it = std::find(list.begin(), list.end(), newEdge);
						if (it == list.end())
						{
							list.push_back(newEdge);
						}
						else
						{
							it->Inc();
						}
						
					}
				}
				else
				{
					if (coord.idx > 0)
					{
						const Vertex & next = posChr_[coord.chr][coord.idx - 1];
						char ch = TwoPaCo::DnaChar::ReverseChar(sequence_[coord.chr][now.pos - 1]);
						char revCh = sequence_[coord.chr][now.pos + k_];
						Edge newEdge(-now.id, -next.id, ch, revCh, now.pos - next.pos, 1);
						auto it = std::find(list.begin(), list.end(), newEdge);
						if (it == list.end())
						{
							list.push_back(newEdge);
						}
						else
						{
							it->Inc();
						}
					}
				}
			}

			std::sort(list.begin(), list.end());
			list.erase(std::unique(list.begin(), list.end()), list.end());
		}

		void Init(const std::string & inFileName, const std::string & genomesFileName, int64_t threads)
		{
			this_ = this;
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
			used_.resize(posChr_.size());
			for (size_t i = 0; i < posChr_.size(); i++)
			{
				used_[i].reset(new std::atomic<bool>[posChr_[i].size()]);
				std::fill(used_[i].get(), used_[i].get() + posChr_[i].size(), false);
			}

			for (TwoPaCo::StreamFastaParser parser(genomesFileName); parser.ReadRecord(); record++)
			{
				sequenceDescription_.push_back(parser.GetCurrentHeader());
				for (char ch; parser.GetChar(ch); )
				{
					sequence_[record].push_back(ch);
				}				
			}

			for (size_t i = 0; i < posChr_.size(); i++)
			{
				for (size_t j = 0; j < posChr_[i].size(); j++)
				{
					int64_t pos_ = posChr_[i][j].pos;
					posChr_[i][j].ch_ = sequence_[i][pos_ + JunctionStorage::this_->k_];
					posChr_[i][j].revCh_ = pos_ > 0 ? TwoPaCo::DnaChar::ReverseChar(sequence_[i][pos_ - 1]) : 'N';
				}
			}
			
			mutex_.resize(GetChrNumber());
			chrSizeBits_.resize(GetChrNumber(), 1);
			for (mutexBits_ = 3; (1 << mutexBits_) < threads * 8; mutexBits_++);
			for (size_t i = 0; i < mutex_.size(); i++) 
			{
				mutex_[i].reset(new tbb::mutex[1 << mutexBits_]);
				for (; (int64_t(1) << chrSizeBits_[i]) <= posChr_[i].size(); chrSizeBits_[i]++);
				chrSizeBits_[i] = max(int64_t(0), chrSizeBits_[i] - mutexBits_);
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
		JunctionStorage(const std::string & fileName, const std::string & genomesFileName, uint64_t k, int64_t threads) : k_(k)
		{
			Init(fileName, genomesFileName, threads);
		}


	private:
		
		struct LightEdge
		{
			int64_t vertex;
			char ch;
		};

		size_t MutexIdx(size_t chrId, size_t idx) const
		{
			size_t ret = idx >> chrSizeBits_[chrId];
			assert(ret < MutexNumber());
			return ret;
		}

		int64_t k_;
		int64_t mutexBits_;
		std::vector<std::vector<Edge> > ingoingEdge_;
		std::vector<std::vector<Edge> > outgoingEdge_;
		std::vector<std::string> sequence_;
		std::vector<std::string> sequenceDescription_;
		std::vector<VertexVector> posChr_;
		std::vector<int64_t> chrSizeBits_;
		std::vector<CoordinateVector> coordinate_;
		std::vector<std::unique_ptr<std::atomic<bool>[] > > used_;
		std::vector<std::unique_ptr<tbb::mutex[]> > mutex_;
		static JunctionStorage * this_;
	};
}

#endif
