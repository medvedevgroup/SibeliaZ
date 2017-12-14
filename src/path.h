#ifndef _PATH_H_
#define _PATH_H_

#include <set>
#include <cassert>
#include <algorithm>
#include "distancekeeper.h"


#include <tbb/mutex.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_scheduler_init.h>

namespace Sibelia
{
	struct Assignment
	{
		int32_t block;
		int32_t instance;
		Assignment()
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
			bool checkConsistency = false) :
			maxBranchSize_(maxBranchSize),
			minBlockSize_(minBlockSize),
			maxFlankingSize_(maxFlankingSize),
			storage_(&storage),
			distanceKeeper_(storage.GetVerticesNumber()),
			instance_(storage.GetChrNumber())
		{
			
		}

		void Init(int64_t vid)
		{
			origin_ = vid;			
			distanceKeeper_.Set(vid, 0);
			leftBodyFlank_ = rightBodyFlank_ = 0;
			for (auto & instanceSet : instance_)
			{
				instanceSet.clear();
			}

			for (JunctionStorage::JunctionIterator it(vid); it.Valid(); ++it)
			{				
				if (!it.IsUsed())
				{
					instance_[it.GetChrId()].insert(Instance(it, 0));
				}
			}
		}

		bool IsInPath(int64_t vertex) const
		{
			return distanceKeeper_.IsSet(vertex);
		}

		struct Instance
		{	
		private:
			int64_t frontDistance_;
			int64_t backDistance_;
			JunctionStorage::JunctionIterator front_;
			JunctionStorage::JunctionIterator back_;
		public:			
		
			Instance()
			{

			}

			Instance(JunctionStorage::JunctionIterator & it, int64_t distance) : front_(it), back_(it), frontDistance_(distance), backDistance_(distance)
			{

			}

			void ChangeFront(const JunctionStorage::JunctionIterator & it, int64_t distance)
			{
				front_ = it;
				frontDistance_ = distance;
			}

			void ChangeBack(const JunctionStorage::JunctionIterator & it, int64_t distance)
			{
				back_ = it;
				backDistance_ = distance;
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

			int64_t LeftFlankDistance() const
			{
				return frontDistance_;
			}

			int64_t RightFlankDistance() const
			{
				return backDistance_;
			}

			bool Within(const JunctionStorage::JunctionIterator it) const
			{
				int64_t left = min(front_.GetIndex(), back_.GetIndex());
				int64_t right = max(front_.GetIndex(), back_.GetIndex());
				return it.GetIndex() >= left && it.GetIndex() <= right;
			}

			bool operator < (const Instance & inst) const
			{				
				int64_t idx1 = back_.IsPositiveStrand() ? back_.GetIndex() : front_.GetIndex();
				int64_t idx2 = inst.back_.IsPositiveStrand() ? inst.back_.GetIndex() : inst.front_.GetIndex();
				return idx1 < idx2;
			}
		};

		typedef std::multiset<Instance> InstanceSet;

		struct Point
		{
		private:
			Edge edge;
			int64_t startDistance;
		public:
			Point() {}
			Point(Edge edge, int64_t startDistance) : edge(edge), startDistance(startDistance) {}

			Edge GetEdge() const
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

		const std::vector<InstanceSet> & Instances() const
		{
			return instance_;
		}

		int64_t LeftDistance() const
		{
			return -leftBodyFlank_;
		}

		int64_t RightDistance() const
		{
			return rightBodyFlank_;
		}

		int64_t MiddlePathLength() const
		{
			return LeftDistance() + RightDistance();
		}

		int64_t GetEndVertex() const
		{
			if (rightBody_.size() > 0)
			{
				return rightBody_.back().GetEdge().GetEndVertex();
			}

			return origin_;
		}

		int64_t GetStartVertex() const
		{
			if (leftBody_.size() > 0)
			{
				return leftBody_.back().GetEdge().GetStartVertex();
			}

			return origin_;
		}

		size_t RightSize() const
		{
			return rightBody_.size() + 1;
		}

		int64_t RightVertex(size_t idx) const
		{
			if (idx == 0)
			{
				return origin_;
			}
			
			return rightBody_[idx - 1].GetEdge().GetEndVertex();
		}

		size_t LeftSize() const
		{
			return leftBody_.size() + 1;
		}

		int64_t LeftVertex(size_t idx) const
		{
			if (idx == 0)
			{
				return origin_;
			}

			return leftBody_[idx - 1].GetEdge().GetStartVertex();
		}

		void DumpInstances(std::ostream & out) const
		{
			for (auto & instanceSet : instance_)
			{
				for (auto inst : instanceSet)
				{
					out << "(" << (inst.Front().IsPositiveStrand() ? '+' : '-') << inst.Front().GetChrId() << ' ' << inst.Front().GetIndex() << ' ' << inst.Back().GetIndex() << ')' << std::endl;
				}
			}			
		}

		void DumpPath(std::vector<Edge> & ret) const
		{
			ret.clear();
			for (auto it = leftBody_.rbegin(); it != leftBody_.rend(); ++it)
			{
				ret.push_back(it->GetEdge());
			}

			for (auto it = rightBody_.rbegin(); it != rightBody_.rend(); ++it)
			{
				ret.push_back(it->GetEdge());
			}
		}

		bool Compatible(const JunctionStorage::JunctionIterator & start, const JunctionStorage::JunctionIterator & end, const Edge & e) const
		{
			if (start.IsPositiveStrand() != end.IsPositiveStrand())
			{
				return false;
			}

			int64_t diff = end.GetPosition() - start.GetPosition();
			if (start.IsPositiveStrand())
			{
				if (diff < 0)
				{
					return false;
				}

				auto start1 = start.SequentialIterator().Next();
				if (diff > maxBranchSize_ && (!start1.Valid() || start.GetChar() != e.GetChar() || start1.GetVertexId() != e.GetEndVertex()))
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

				auto start1 = start.SequentialIterator().Next();
				if (-diff > maxBranchSize_ && (!start1.Valid() || start.GetChar() != e.GetChar() || start1.GetVertexId() != e.GetEndVertex()))
				{
					return false;
				}
			}

			return true;
		}		

		class PointPushFrontWorker
		{
		public:
			Edge e;
			Path * path;
			int64_t vertex;
			int64_t distance;
			bool & failFlag;

			PointPushFrontWorker(Path * path, int64_t vertex, int64_t distance, Edge e, bool & failFlag) : path(path), vertex(vertex), e(e), failFlag(failFlag), distance(distance)
			{

			}

			void operator()() const
			{
				for (JunctionStorage::JunctionIterator nowIt(vertex); nowIt.Valid() && !failFlag; nowIt++)
				{
					bool newInstance = true;
					if (!nowIt.IsUsed())
					{
						auto & instanceSet = path->instance_[nowIt.GetChrId()];
						auto inst = instanceSet.upper_bound(Instance(nowIt, 0));
						if (inst != instanceSet.end() && inst->Within(nowIt))
						{
							continue;
						}

						if (nowIt.IsPositiveStrand())
						{
							if (inst != instanceSet.end() && path->Compatible(nowIt, inst->Front(), e))
							{
								newInstance = false;
							}
						}
						else
						{
							if (inst != instanceSet.begin() && path->Compatible(nowIt, (--inst)->Front(), e))
							{
								newInstance = false;
							}
						}

						if (!newInstance && inst->Front().GetVertexId() != vertex)
						{
							int64_t nextLength = abs(nowIt.GetPosition() - inst->Back().GetPosition());
							int64_t rightFlankSize = path->rightBodyFlank_ - inst->RightFlankDistance();
							assert(rightFlankSize >= 0);
							if (nextLength >= path->minBlockSize_ && rightFlankSize > path->maxFlankingSize_)
							{
								failFlag = true;
								break;
							}

							const_cast<Instance&>(*inst).ChangeFront(nowIt, distance);
						}
						else
						{
							instanceSet.insert(Instance(nowIt, distance));
						}
					}
				}
			}

		};

		class PointPushBackWorker
		{
		public:
			Edge e;
			Path * path;
			int64_t vertex;
			int64_t distance;
			bool & failFlag;

			PointPushBackWorker(Path * path, int64_t vertex, int64_t distance, Edge e, bool & failFlag) : path(path), vertex(vertex), e(e), failFlag(failFlag), distance(distance)
			{

			}

			void operator()() const
			{
				for (JunctionStorage::JunctionIterator nowIt(vertex); nowIt.Valid() && !failFlag; nowIt++)
				{
					bool newInstance = true;
					if (!nowIt.IsUsed())
					{						
						auto & instanceSet = path->instance_[nowIt.GetChrId()];
						auto inst = instanceSet.upper_bound(Instance(nowIt, 0));
						if (inst != instanceSet.end() && inst->Within(nowIt))
						{
							continue;
						}

						if (nowIt.IsPositiveStrand())
						{
							if (inst != instanceSet.begin() && path->Compatible((--inst)->Back(), nowIt, e))
							{
								newInstance = false;
							}
						}
						else
						{
							if (inst != instanceSet.end() && path->Compatible(inst->Back(), nowIt, e))
							{
								newInstance = false;
							}
						}

						if (!newInstance && inst->Back().GetVertexId() != vertex)
						{
							int64_t nextLength = abs(nowIt.GetPosition() - inst->Front().GetPosition());
							int64_t leftFlankSize = -(path->leftBodyFlank_ - inst->LeftFlankDistance());
							assert(leftFlankSize >= 0);
							if (nextLength >= path->minBlockSize_ && leftFlankSize > path->maxFlankingSize_)
							{
								failFlag = true;
								break;
							}

							const_cast<Instance&>(*inst).ChangeBack(nowIt, distance);
						}
						else
						{
							instanceSet.insert(Instance(nowIt, distance));
						}
					}
				}
			}

		};

		bool PointPushBack(const Edge & e)
		{
			int64_t vertex = e.GetEndVertex();
			if (distanceKeeper_.IsSet(vertex))
			{
				return false;
			}

			bool failFlag = false;
			int64_t startVertexDistance = rightBodyFlank_;
			int64_t endVertexDistance = startVertexDistance + e.GetLength();
			PointPushBackWorker(this, vertex, endVertexDistance, e, failFlag)();
			rightBody_.push_back(Point(e, startVertexDistance));
			distanceKeeper_.Set(e.GetEndVertex(), endVertexDistance);
			rightBodyFlank_ = rightBody_.back().EndDistance();			
			return !failFlag;
		}

		
		bool PointPushFront(const Edge & e)
		{
			int64_t vertex = e.GetStartVertex();
			if (distanceKeeper_.IsSet(vertex))
			{
				return false;
			}

			bool failFlag = false;
			int64_t endVertexDistance = leftBodyFlank_;
			int64_t startVertexDistance = endVertexDistance - e.GetLength();
			PointPushFrontWorker(this, vertex, startVertexDistance, e, failFlag)();
			leftBody_.push_back(Point(e, startVertexDistance));
			distanceKeeper_.Set(e.GetStartVertex(), startVertexDistance);
			leftBodyFlank_ = leftBody_.back().StartDistance();		
			return !failFlag;
		}

		int64_t Score(bool final = false) const
		{
			int64_t score;
			int64_t length;
			int64_t ret = 0;
			int64_t middlePath = MiddlePathLength();
			for(auto & instanceSet : instance_)
			{
				for (auto & inst : instanceSet)
				{
					InstanceScore(inst, length, score, middlePath);
					if (!final || length >= minBlockSize_)
					{
						ret += score;
					}
				}
			}

			return ret;
		}

		int64_t GoodInstances() const
		{
			int64_t ret = 0;
			for (auto & instanceSet : instance_)
			{
				for (auto & inst : instanceSet)
				{
					if (IsGoodInstance(inst))
					{
						ret++;
					}
				}
			}

			return ret;
		}

		bool IsGoodInstance(const Instance & it) const
		{
			int64_t score;
			int64_t length;
			InstanceScore(it, length, score, MiddlePathLength());
			return length >= minBlockSize_;
		}

		void InstanceScore(const Instance & inst, int64_t & length, int64_t & score, int64_t middlePath) const
		{			
			length = abs(inst.Front().GetPosition() - inst.Back().GetPosition());
			score = length - (middlePath - length);
		}		

		void Clear()
		{
			for (auto pt : leftBody_)
			{
				distanceKeeper_.Unset(pt.GetEdge().GetStartVertex());
			}

			for (auto pt : rightBody_)
			{
				distanceKeeper_.Unset(pt.GetEdge().GetEndVertex());
			}

			leftBody_.clear();
			rightBody_.clear();
			distanceKeeper_.Unset(origin_);			
			for (int64_t v1 = -storage_->GetVerticesNumber() + 1; v1 < storage_->GetVerticesNumber(); v1++)
			{
				assert(!distanceKeeper_.IsSet(v1));
			}

			for (auto & instanceSet : instance_)
			{
				instanceSet.clear();
			}
		}

	private:


		std::vector<Point> leftBody_;
		std::vector<Point> rightBody_;
		std::vector<InstanceSet> instance_;

		int64_t origin_;
		int64_t minBlockSize_;
		int64_t maxBranchSize_;
		int64_t maxFlankingSize_;
		int64_t leftBodyFlank_;
		int64_t rightBodyFlank_;		
		DistanceKeeper distanceKeeper_;
		const JunctionStorage * storage_;
	};
}

#endif

