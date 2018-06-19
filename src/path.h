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
			int64_t minScoringUnit,
			int64_t maxFlankingSize) :
			maxBranchSize_(maxBranchSize),
			minBlockSize_(minBlockSize),
			minScoringUnit_(minScoringUnit),
			maxFlankingSize_(maxFlankingSize),
			storage_(&storage),
			distanceKeeper_(storage.GetVerticesNumber()),
			instance_(storage.GetChrNumber())
		{

		}

		void Init(int64_t vid)
		{
			origin_ = vid;
			distanceKeeper_.Set(vid);
			leftBodyFlank_ = rightBodyFlank_ = 0;
			for (JunctionStorage::JunctionIterator it(vid); it.Valid(); ++it)
			{
				if (!it.IsUsed())
				{
					allInstance_.push_back(instance_[it.GetChrId()].insert(Instance(it.SequentialIterator(), 0)));
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
			int32_t compareIdx_;
			int32_t frontDistance_;
			int32_t backDistance_;
			JunctionStorage::JunctionSequentialIterator front_;
			JunctionStorage::JunctionSequentialIterator back_;
		public:

			static bool OldComparator(const Instance & a, const Instance & b)
			{
				if (a.front_.GetChrId() != b.front_.GetChrId())
				{
					return a.front_.GetChrId() < b.front_.GetChrId();
				}

				int64_t idx1 = a.back_.IsPositiveStrand() ? a.back_.GetIndex() : a.front_.GetIndex();
				int64_t idx2 = b.back_.IsPositiveStrand() ? b.back_.GetIndex() : b.front_.GetIndex();
				return idx1 < idx2;
			}

			Instance()
			{

			}

			Instance(const JunctionStorage::JunctionSequentialIterator & it, int64_t distance) : front_(it),
				back_(it),
				frontDistance_(distance),
				backDistance_(distance),
				compareIdx_(it.GetIndex())
			{

			}

			void ChangeFront(const JunctionStorage::JunctionSequentialIterator & it, int64_t distance)
			{
				front_ = it;
				frontDistance_ = distance;
				assert(backDistance_ >= frontDistance_);
				if (!back_.IsPositiveStrand())
				{
					compareIdx_ = front_.GetIndex();
				}
			}

			void ChangeBack(const JunctionStorage::JunctionSequentialIterator & it, int64_t distance)
			{
				back_ = it;
				backDistance_ = distance;
				assert(backDistance_ >= frontDistance_);
				if (back_.IsPositiveStrand())
				{
					compareIdx_ = back_.GetIndex();
				}
			}

			bool SinglePoint() const
			{
				return front_ == back_;
			}

			JunctionStorage::JunctionSequentialIterator Front() const
			{
				return front_;
			}

			JunctionStorage::JunctionSequentialIterator Back() const
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

			int64_t UtilityLength() const
			{
				return backDistance_ - frontDistance_;
			}

			int64_t RealLength() const
			{
				return abs(front_.GetPosition() - back_.GetPosition());
			}

			bool Within(const JunctionStorage::JunctionIterator it) const
			{
				int64_t left = min(front_.GetIndex(), back_.GetIndex());
				int64_t right = max(front_.GetIndex(), back_.GetIndex());
				return it.GetIndex() >= left && it.GetIndex() <= right;
			}

			bool operator < (const Instance & inst) const
			{
				return compareIdx_ < inst.compareIdx_;
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

		const std::vector<InstanceSet::iterator> & AllInstances() const
		{
			return allInstance_;
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

		Point RightPoint(size_t idx) const
		{
			return rightBody_[idx];
		}

		Point LeftPoint(size_t idx) const
		{
			return leftBody_[idx];
		}

		int64_t RightVertex() const
		{
			if (rightBody_.size() == 0)
			{
				return origin_;
			}

			return rightBody_.back().GetEdge().GetEndVertex();
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

		int64_t LeftVertex() const
		{
			if (leftBody_.size() == 0)
			{
				return origin_;
			}

			return leftBody_.back().GetEdge().GetStartVertex();
		}

		void DumpPath(std::ostream & out) const
		{
			out << "Left:" << std::endl;
			for (auto it = leftBody_.rbegin(); it != leftBody_.rend(); ++it)
			{
				out << it->GetEdge().GetStartVertex() << " -> " << it->GetEdge().GetEndVertex() << ", " << it->GetEdge().GetChar() << ", " << it->StartDistance() << ", " << it->EndDistance() << std::endl;
			}

			out << "Right:" << std::endl;
			for (auto it = rightBody_.begin(); it != rightBody_.end(); ++it)
			{
				out << it->GetEdge().GetStartVertex() << " -> " << it->GetEdge().GetEndVertex() << ", " << it->GetEdge().GetChar() << ", " << it->StartDistance() << ", " << it->EndDistance() << std::endl;
			}
		}

		void DumpInstances(std::ostream & out) const
		{
			size_t total = 0;
			for (auto & instanceSet : instance_)
			{
				total += instanceSet.size();
				for (auto inst : instanceSet)
				{
					int64_t middlePath = MiddlePathLength();
					int64_t length = inst.UtilityLength();

					int64_t start = inst.Front().GetPosition();
					int64_t end = inst.Back().GetPosition();
					out << "(" << (inst.Front().IsPositiveStrand() ? '+' : '-') <<
						inst.Front().GetChrId() << ' ' << start << ' ' << end << ' ' << end - start << ';' <<
						inst.LeftFlankDistance() << ' ' << inst.RightFlankDistance() << ')' << std::endl;
				}
			}

			out << "Total: " << total << std::endl;
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

		bool Compatible(const JunctionStorage::JunctionSequentialIterator & start, const JunctionStorage::JunctionSequentialIterator & end, const Edge & e) const
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

				auto start1 = start.Next();
				if (diff > maxBranchSize_ && (!start1.Valid() || start.GetChar() != e.GetChar() || end != start1 || start1.GetVertexId() != e.GetEndVertex()))
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

				auto start1 = start.Next();
				if (-diff > maxBranchSize_ && (!start1.Valid() || start.GetChar() != e.GetChar() || end != start1 || start1.GetVertexId() != e.GetEndVertex()))
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
						auto inst = instanceSet.upper_bound(Instance(nowIt.SequentialIterator(), 0));
						if (inst != instanceSet.end() && inst->Within(nowIt))
						{
							continue;
						}

						if (nowIt.IsPositiveStrand())
						{
							if (inst != instanceSet.end() && path->Compatible(nowIt.SequentialIterator(), inst->Front(), e))
							{
								newInstance = false;
							}
						}
						else
						{
							if (inst != instanceSet.begin() && path->Compatible(nowIt.SequentialIterator(), (--inst)->Front(), e))
							{
								newInstance = false;
							}
						}

						if (!newInstance && inst->Front().GetVertexId() != vertex)
						{
							bool prevGoodInstance = path->IsGoodInstance(*inst);
							const_cast<Instance&>(*inst).ChangeFront(nowIt.SequentialIterator(), distance);
							if (!prevGoodInstance && path->IsGoodInstance(*inst))
							{
								path->goodInstance_.push_back(inst);
							}
						}
						else
						{
							path->allInstance_.push_back(instanceSet.insert(Instance(nowIt.SequentialIterator(), distance)));
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
						auto inst = instanceSet.upper_bound(Instance(nowIt.SequentialIterator(), 0));
						if (inst != instanceSet.end() && inst->Within(nowIt))
						{
							continue;
						}

						if (nowIt.IsPositiveStrand())
						{
							if (inst != instanceSet.begin() && path->Compatible((--inst)->Back(), nowIt.SequentialIterator(), e))
							{
								newInstance = false;
							}
						}
						else
						{
							if (inst != instanceSet.end() && path->Compatible(inst->Back(), nowIt.SequentialIterator(), e))
							{
								newInstance = false;
							}
						}

						if (!newInstance && inst->Back().GetVertexId() != vertex)
						{
							bool prevGoodInstance = path->IsGoodInstance(*inst);
							const_cast<Instance&>(*inst).ChangeBack(nowIt.SequentialIterator(), distance);
							if (!prevGoodInstance && path->IsGoodInstance(*inst))
							{
								path->goodInstance_.push_back(inst);
							}
						}
						else
						{
							path->allInstance_.push_back(instanceSet.insert(Instance(nowIt.SequentialIterator(), distance)));
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
			distanceKeeper_.Set(e.GetEndVertex());
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
			distanceKeeper_.Set(e.GetStartVertex());
			leftBodyFlank_ = leftBody_.back().StartDistance();
			return !failFlag;
		}

		int64_t Score(bool final = false) const
		{
			int64_t ret = 0;
			int64_t multiplier = goodInstance_.size();
			for (auto & instanceIt : goodInstance_)
			{
				int64_t score = instanceIt->RealLength();
				int64_t penalty = MiddlePathLength() - instanceIt->UtilityLength();
				assert(penalty >= 0);
				if (penalty >= maxFlankingSize_)
				{
					score -= penalty * multiplier;
				}
				else
				{
					score -= penalty;
				}

				ret += score;
			}

			return ret;
		}

		int64_t GoodInstances() const
		{
			return goodInstance_.size();
		}

		static bool CmpInstance(const InstanceSet::iterator & a, const InstanceSet::iterator & b)
		{
			return Path::Instance::OldComparator(*a, *b);
		}

		const std::vector<InstanceSet::iterator> & GoodInstancesList() const
		{
			return goodInstance_;
		}

		bool IsGoodInstance(const Instance & inst) const
		{
			return inst.RealLength() >= minBlockSize_;
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

			for (auto it : allInstance_)
			{
				instance_[it->Front().GetChrId()].erase(it);
			}

			allInstance_.clear();
			goodInstance_.clear();
		}

	private:

		std::vector<Point> leftBody_;
		std::vector<Point> rightBody_;
		std::vector<InstanceSet> instance_;
		std::vector<InstanceSet::iterator> allInstance_;
		std::vector<InstanceSet::iterator> goodInstance_;

		int64_t origin_;
		int64_t minBlockSize_;
		int64_t minScoringUnit_;
		int64_t maxBranchSize_;
		int64_t leftBodyFlank_;
		int64_t rightBodyFlank_;
		int64_t maxFlankingSize_;
		DistanceKeeper distanceKeeper_;
		const JunctionStorage * storage_;
		friend class BestPath;
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
				bool ret = path.PointPushBack(pt.GetEdge());
				assert(ret);
			}

			newRightBody_.clear();
			rightFlank_ = path.rightBody_.size();
		}

		void FixBackward(Path & path)
		{
			for (auto & pt : newLeftBody_)
			{
				bool ret = path.PointPushFront(pt.GetEdge());
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

