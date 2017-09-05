#ifndef _PATH_H_
#define _PATH_H_

#include <set>
#include <cassert>
#include <algorithm>
#include "distancekeeper.h"


namespace Sibelia
{
	struct Assignment
	{
		static const int64_t UNKNOWN_BLOCK;
		int32_t block;
		int32_t instance;
		Assignment() : block(UNKNOWN_BLOCK), instance(UNKNOWN_BLOCK)
		{

		}

		bool operator == (const Assignment & assignment) const
		{
			return block == assignment.block && instance == assignment.instance;
		}
	};

	struct BestPath;

	struct Path
	{
	public:
		Path(const JunctionStorage & storage,			
			int64_t maxBranchSize,
			int64_t minBlockSize,
			int64_t maxFlankingSize,
			std::vector<std::vector<Assignment> > & blockId,
			bool checkConsistency = false) :
			maxBranchSize_(maxBranchSize), minBlockSize_(minBlockSize), maxFlankingSize_(maxFlankingSize), storage_(&storage),
			distanceKeeper_(storage.GetVerticesNumber()), minChainSize_(minBlockSize - 2 * maxFlankingSize), blockId_(blockId),
			checkConsistency_(checkConsistency)
		{
			
		}

		void Init(int64_t vid)
		{
			leftBody_.clear();
			rightBody_.clear();
			instance_.clear();
			distanceKeeper_.Clear();

			origin_ = vid;			
			distanceKeeper_.Set(vid, 0);
			for (size_t i = 0; i < storage_->GetInstancesCount(vid); i++)
			{
				JunctionStorage::JunctionIterator it = storage_->GetJunctionInstance(vid, i);
				if (blockId_[it.GetChrId()][it.GetIndex()].block == Assignment::UNKNOWN_BLOCK)
				{
					instance_.insert(Instance(it));
				}
			}
		}		

		struct Instance
		{	
		private:
			JunctionStorage::JunctionIterator front_;
			JunctionStorage::JunctionIterator back_;
		public:			
		
			Instance()
			{

			}

			Instance(JunctionStorage::JunctionIterator & it) : front_(it), back_(it)
			{

			}

			void ChangeFront(const JunctionStorage::JunctionIterator & it)
			{
				front_ = it;
			}

			void ChangeBack(const JunctionStorage::JunctionIterator & it)
			{
				back_ = it;
			}			

			bool SinglePoint() const
			{
				return front_ == back_;
			}

			JunctionStorage::JunctionIterator Front() const
			{
				return front_;
			}

			JunctionStorage::JunctionIterator Back() const
			{
				return back_;
			}

			int64_t LeftFlankDistance(const DistanceKeeper & keeper, const JunctionStorage * storage_) const
			{
				return keeper.Get(front_.GetVertexId(storage_));
			}

			int64_t RightFlankDistance(const DistanceKeeper & keeper, const JunctionStorage * storage_) const
			{
				return keeper.Get(back_.GetVertexId(storage_));
			}

			bool Within(const JunctionStorage::JunctionIterator it) const
			{
				if (it.GetChrId() == front_.GetChrId())
				{
					int64_t left = min(front_.GetIndex(), back_.GetIndex());
					int64_t right = max(front_.GetIndex(), back_.GetIndex());
					return it.GetIndex() >= left && it.GetIndex() <= right;
				}

				return false;
			}

			bool operator < (const Instance & inst) const
			{
				if (front_.GetChrId() != inst.front_.GetChrId())
				{
					return front_.GetChrId() < inst.front_.GetChrId();
				}

				int64_t idx1 = back_.IsPositiveStrand() ? back_.GetIndex() : front_.GetIndex();
				int64_t idx2 = inst.back_.IsPositiveStrand() ? inst.back_.GetIndex() : inst.front_.GetIndex();
				return idx1 < idx2;
			}
		};

		struct Point
		{
		private:
			Edge edge;
			int64_t startDistance;
		public:
			Point() {}
			Point(Edge edge, int64_t startDistance) : edge(edge), startDistance(startDistance) {}

			Edge Edge() const
			{
				return edge;
			}

			int64_t StartDistance() const
			{
				return startDistance;
			}

			int64_t EndDistance() const
			{
				return startDistance + edge.GetLength();
			}

			bool operator == (const Point & p) const
			{
				return startDistance == p.startDistance && edge == p.edge;
			}

			bool operator != (const Point & p) const
			{
				return p != *this;
			}
		};

		int64_t Origin() const
		{
			return origin_;
		}		

		const std::multiset<Instance> & Instances() const
		{
			return instance_;
		}

		int64_t MiddlePathLength() const
		{
			return (rightBody_.size() > 0 ? rightBody_.back().EndDistance() : 0) - (leftBody_.size() > 0 ? leftBody_.back().StartDistance() : 0);
		}

		int64_t GetEndVertex() const
		{
			if (rightBody_.size() > 0)
			{
				return rightBody_.back().Edge().GetEndVertex();
			}

			return origin_;
		}

		int64_t GetStartVertex() const
		{
			if (leftBody_.size() > 0)
			{
				return leftBody_.back().Edge().GetStartVertex();
			}

			return origin_;
		}

		void DumpPath(std::vector<Edge> & ret) const
		{
			ret.clear();
			for (auto it = leftBody_.rbegin(); it != leftBody_.rend(); ++it)
			{
				ret.push_back(it->Edge());
			}

			for (auto it = rightBody_.rbegin(); it != rightBody_.rend(); ++it)
			{
				ret.push_back(it->Edge());
			}
		}

		bool Compatible(const JunctionStorage::JunctionIterator & start, const JunctionStorage::JunctionIterator & end, const Edge & e) const
		{
			if (start.GetChrId() != end.GetChrId() || start.IsPositiveStrand() != end.IsPositiveStrand())
			{
				return false;
			}

			int64_t diff = end.GetPosition(storage_) - start.GetPosition(storage_);
			if (start.IsPositiveStrand())
			{
				if (diff < 0)
				{
					return false;
				}

				auto start1 = start + 1;
				if (diff > maxBranchSize_ && (start.GetChar(storage_) != e.GetChar() || end != start1 || start1.GetVertexId(storage_) != e.GetEndVertex()))
				{
					return false;
				}
			}
			else
			{
				if (-diff < 0)
				{
					return false;
				}

				auto start1 = start + 1;
				if (-diff > maxBranchSize_ && (start.GetChar(storage_) != e.GetChar() || end != start1 || start1.GetVertexId(storage_) != e.GetEndVertex()))
				{
					return false;
				}
			}

			return true;
		}

		bool PointPushBack(const Edge & e)
		{
			int64_t vertex = e.GetEndVertex();
			if (distanceKeeper_.IsSet(vertex))
			{
				return false;
			}

			bool fail = false;
			int64_t startVertexDistance = rightBody_.empty() ? 0 : rightBody_.back().EndDistance();
			int64_t endVertexDistance = startVertexDistance + e.GetLength();

			for (size_t i = 0; i < storage_->GetInstancesCount(vertex); i++)
			{
				bool newInstance = true;
				JunctionStorage::JunctionIterator nowIt = storage_->GetJunctionInstance(vertex, i);
				if (blockId_[nowIt.GetChrId()][nowIt.GetIndex()].block == Assignment::UNKNOWN_BLOCK)
				{
					auto inst = instance_.upper_bound(Instance(nowIt));
					if (inst != instance_.end() && inst->Within(nowIt))
					{
						continue;
					}

					if (nowIt.IsPositiveStrand())
					{
						if (inst != instance_.begin() && Compatible((--inst)->Back(), nowIt, e))
						{
							newInstance = false;
						}
					}
					else
					{
						if (inst != instance_.end() && Compatible(inst->Back(), nowIt, e))
						{
							newInstance = false;
						}
					}
					
					if (!newInstance && inst->Back().GetVertexId(storage_) != vertex)
					{
						int64_t nextLength = abs(nowIt.GetPosition(storage_) - inst->Front().GetPosition(storage_));
						int64_t leftFlankSize = abs(inst->LeftFlankDistance(distanceKeeper_, storage_) - (leftBody_.empty() ? 0 : leftBody_.back().StartDistance()));
						if (nextLength >= minChainSize_ && leftFlankSize > maxFlankingSize_)
						{
							fail = true;
							break;
						}
						
						const_cast<Instance&>(*inst).ChangeBack(nowIt);
					}
					else
					{
						instance_.insert(Instance(nowIt));
					}
				}
			}

			rightBody_.push_back(Point(e, startVertexDistance));
			distanceKeeper_.Set(e.GetEndVertex(), endVertexDistance);

			if (fail)
			{
				PointPopBack();
				return false;
			}
				
			return true;
		}

		void PointPopBack()
		{
			int64_t lastVertex = rightBody_.back().Edge().GetEndVertex();
			rightBody_.pop_back();
			distanceKeeper_.Unset(lastVertex);
			for (auto it = instance_.begin(); it != instance_.end(); )
			{
				if (it->Back().GetVertexId(storage_) == lastVertex)
				{
					if (it->Front() == it->Back())
					{						
						it = instance_.erase(it);
					}
					else
					{						
						auto jt = it->Back();
						while (true)
						{
							if (distanceKeeper_.IsSet(jt.GetVertexId(storage_)))
							{
								const_cast<Instance&>(*it).ChangeBack(jt);
								break;
							}
							else
							{
								--jt;
							}
						}

						it++;
					}				
				}
				else
				{
					++it;
				}
			}			
		}

		bool PointPushFront(const Edge & e)
		{
			int64_t vertex = e.GetStartVertex();
			if (distanceKeeper_.IsSet(vertex))
			{
				return false;
			}

			bool fail = false;
			int64_t endVertexDistance = leftBody_.empty() ? 0 : leftBody_.back().StartDistance();
			int64_t startVertexDistance = endVertexDistance - e.GetLength();

			for (size_t i = 0; i < storage_->GetInstancesCount(vertex); i++)
			{
				bool newInstance = true;
				JunctionStorage::JunctionIterator nowIt = storage_->GetJunctionInstance(vertex, i);
				if (blockId_[nowIt.GetChrId()][nowIt.GetIndex()].block == Assignment::UNKNOWN_BLOCK)
				{
					auto inst = instance_.upper_bound(Instance(nowIt));
					if (inst != instance_.end() && inst->Within(nowIt))
					{
						continue;
					}

					if (nowIt.IsPositiveStrand())
					{
						if (inst != instance_.end() && Compatible(nowIt, inst->Front(), e))
						{
							newInstance = false;
						}
					}
					else
					{
						if (inst != instance_.begin() && Compatible(nowIt, (--inst)->Front(), e))
						{
							newInstance = false;
						}
					}

					if (!newInstance && inst->Front().GetVertexId(storage_) != vertex)
					{
						int64_t nextLength = abs(nowIt.GetPosition(storage_) - inst->Back().GetPosition(storage_));
						int64_t rightFlankSize = abs(inst->RightFlankDistance(distanceKeeper_, storage_) - (rightBody_.empty() ? 0 : rightBody_.back().EndDistance()));
						if (nextLength >= minChainSize_ && rightFlankSize > maxFlankingSize_)
						{
							fail = true;
							break;
						}

						const_cast<Instance&>(*inst).ChangeFront(nowIt);
					}
					else
					{
						instance_.insert(Instance(nowIt));
					}
				}
			}

			leftBody_.push_back(Point(e, startVertexDistance));
			distanceKeeper_.Set(e.GetStartVertex(), startVertexDistance);
			
			if (fail)
			{
				PointPopFront();
				return false;
			}

			return true;
		}

		void PointPopFront()
		{
			int64_t lastVertex = leftBody_.back().Edge().GetStartVertex();
			leftBody_.pop_back();
			distanceKeeper_.Unset(lastVertex);
			for (auto it = instance_.begin(); it != instance_.end(); )
			{
				if (it->Front().GetVertexId(storage_) == lastVertex)
				{	
					if (it->Front() == it->Back())
					{
						it = instance_.erase(it);
					}
					else
					{
						auto jt = it->Front();
						while (true)
						{
							if (distanceKeeper_.IsSet(jt.GetVertexId(storage_)))
							{
								const_cast<Instance&>(*it).ChangeFront(jt);
								break;
							}
							else
							{
								++jt;
							}
						}

						it++;
					}
				}
				else
				{
					++it;
				}
			}
		}

		int64_t Score(bool final = false) const
		{
			int64_t score;
			int64_t length;
			int64_t ret = 0;
			for (auto & inst : instance_)
			{
				InstanceScore(inst, length, score);
				if (!final || length >= minChainSize_)
				{
					ret += score;
				}
			}

			return ret;
		}

		int64_t GoodInstances() const
		{
			int64_t ret = 0;
			for (auto & it : instance_)
			{
				if (IsGoodInstance(it))
				{
					ret++;
				}
			}

			return ret;
		}

		bool IsGoodInstance(const Instance & it) const
		{
			int64_t score;
			int64_t length;
			InstanceScore(it, length, score);
			return length >= minChainSize_;
		}

		void InstanceScore(const Instance & inst, int64_t & length, int64_t & score) const
		{
			int64_t leftFlank = abs(inst.LeftFlankDistance(distanceKeeper_, storage_) - (leftBody_.size() > 0 ? leftBody_.back().StartDistance() : 0));
			int64_t rightFlank = abs(inst.RightFlankDistance(distanceKeeper_, storage_) - (rightBody_.size() > 0 ? rightBody_.back().EndDistance() : 0));
			length = abs(inst.Front().GetPosition(storage_) - inst.Back().GetPosition(storage_));
			score = length - leftFlank - rightFlank;
		}

	private:


		friend class BestPath;

		std::vector<Point> leftBody_;
		std::vector<Point> rightBody_;
		std::multiset<Instance> instance_;

		int64_t origin_;
		int64_t minChainSize_;
		int64_t minBlockSize_;
		int64_t maxBranchSize_;
		int64_t maxFlankingSize_;
		bool checkConsistency_;
		DistanceKeeper distanceKeeper_;
		const JunctionStorage * storage_;		
		std::vector<std::vector<Assignment> > & blockId_;		
	};

	struct BestPath
	{
		int64_t score_;
		int64_t leftFlank_;
		int64_t rightFlank_;
		std::vector<Path::Point> newLeftBody_;
		std::vector<Path::Point> newRightBody_;

		void FixForward(Path & path)
		{
			for (auto & pt : newRightBody_)
			{
				bool ret = path.PointPushBack(pt.Edge());
				assert(ret);
			}

			newRightBody_.clear();
			rightFlank_ = path.rightBody_.size();
		}

		void FixBackward(Path & path)
		{
			for (auto & pt : newLeftBody_)
			{
				bool ret = path.PointPushFront(pt.Edge());
				assert(ret);
			}

			newLeftBody_.clear();
			leftFlank_ = path.leftBody_.size();
		}

		void UpdateForward(const Path & path, int64_t newScore)
		{
			score_ = newScore;
			newRightBody_.clear();
			std::copy(path.rightBody_.begin() + rightFlank_, path.rightBody_.end(), std::back_inserter(newRightBody_));
		}

		void UpdateBackward(const Path & path, int64_t newScore)
		{
			score_ = newScore;
			newLeftBody_.clear();
			std::copy(path.leftBody_.begin() + leftFlank_, path.leftBody_.end(), std::back_inserter(newLeftBody_));
		}

		BestPath() 
		{
			Init();
		}

		void Init()
		{
			score_ = leftFlank_ = rightFlank_ = 0;
		}
	};
}

#endif

