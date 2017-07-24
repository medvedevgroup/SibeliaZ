#ifndef _PATH_H_
#define _PATH_H_

#include "forbidden.h"

namespace Sibelia
{
	struct Assignment
	{
		static const int64_t UNKNOWN_BLOCK;
		int64_t block;
		int64_t instance;
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
		Path() {}
		Path(int64_t vid,
			const JunctionStorage & storage,
			int64_t maxBranchSize,
			int64_t minBlockSize,
			int64_t maxFlankingSize,
			const std::vector<std::vector<Assignment> > & blockId,
			bool checkConsistency = false) :
			maxBranchSize_(maxBranchSize), minBlockSize_(minBlockSize), maxFlankingSize_(maxFlankingSize), storage_(&storage),
			minChainSize_(minBlockSize - 2 * maxFlankingSize), blockId_(&blockId), checkConsistency_(checkConsistency), origin_(vid)
		{
			FillBuffer(vid);
			for (auto & it : junctionBuffer_)
			{
				instance_.push_back(Instance());
				instance_.back().seq.push_back(it);
				instance_.back().leftFlankDistance = instance_.back().rightFlankDistance = 0;
			}
		}

		struct Instance
		{			
			int64_t nextFlankDistance;
			int64_t leftFlankDistance;
			int64_t rightFlankDistance;
			JunctionStorage::JunctionIterator nextJunction;
			std::deque<JunctionStorage::JunctionIterator> seq;

			bool operator == (const Instance & inst) const
			{
				return leftFlankDistance == inst.leftFlankDistance && rightFlankDistance == inst.rightFlankDistance && seq == inst.seq;
			}

			bool operator != (const Instance & inst) const
			{
				return inst != *this;
			}
		};

		struct Point
		{
			Edge edge;
			int64_t startDistance;
			int64_t endDistance;
			Point() {}
			Point(Edge edge, int64_t startDistance) : edge(edge), startDistance(startDistance), endDistance(startDistance + edge.GetLength()) {}
			bool operator == (const Point & p) const
			{
				return startDistance == p.startDistance && edge == p.edge;
			}

			bool operator != (const Point & p) const
			{
				return p != *this;
			}
		};

		void PrintInstance(const Instance & inst, std::ostream & out) const
		{
			for (auto & it : inst.seq)
			{
				out << "{vid:" << it.GetVertexId()
					<< " str:" << inst.seq.front().IsPositiveStrand()
					<< " chr:" << inst.seq.front().GetChrId()
					<< " pos:" << it.GetPosition() << "} ";
			}

			out << std::endl;
		}

		void DebugOut(std::ostream & out, bool all = true) const
		{
			out << "Path: ";
			for (auto & pt : body_)
			{
				out << "(" << pt.edge.GetStartVertex() << "," << pt.edge.GetEndVertex() << "," << pt.edge.GetChar() << ") ";
			}

			out << std::endl << "Instances: " << std::endl;
			for (auto & inst : instance_)
			{
				if (all || IsGoodInstance(inst))
				{
					PrintInstance(inst, out);
				}
			}

			out << std::endl;		
		}

		size_t PathLength() const
		{
			return body_.size();
		}

		const std::list<Instance> & Instances() const
		{
			return instance_;
		}

		const std::deque<Point> & PathBody() const
		{
			return body_;
		}

		int64_t MiddlePathLength() const
		{
			return body_.back().endDistance - body_.front().startDistance;
		}

		Edge GetEdge(size_t index) const
		{
			return body_[index].edge;
		}		

		int64_t GetEndVertex() const
		{
			if (body_.size() > 0)
			{
				return body_.back().edge.GetEndVertex();
			}

			return origin_;
		}

		int64_t GetStartVertex() const
		{
			if (body_.size() > 0)
			{
				return body_.back().edge.GetStartVertex();
			}

			return origin_;
		}
		
		bool PointPushBack(const Edge & e)
		{
			int64_t vertex = e.GetEndVertex();
			if (FindVertexInPath(vertex) != body_.end())
			{
				return false;
			}

			int64_t startVertexDistance = body_.empty() ? 0 : body_.back().endDistance;
			int64_t endVertexDistance = startVertexDistance + e.GetLength();
			for (auto & inst : instance_)
			{
				inst.nextFlankDistance = inst.rightFlankDistance;
				inst.nextJunction = storage_->End(inst.seq.back().GetChrId(), inst.seq.back().IsPositiveStrand());				
			}

			FillBuffer(vertex);
			std::vector<JunctionStorage::JunctionIterator> newInstance;
			for (auto it : junctionBuffer_)
			{
				bool nowNewInstance = true;
				for (auto & inst : instance_)
				{
					if (Compatible(inst.seq.back(), it, e))
					{
						nowNewInstance = false;						
						int64_t leftFlank = abs(inst.leftFlankDistance - body_.front().startDistance);
						int64_t nextLength = abs(it.GetPosition() - inst.seq.front().GetPosition()) >= minChainSize_;
						if (nextLength >= minChainSize_ && leftFlank > maxFlankingSize_)
						{
							return false;
						}
						
						if (it.GetRelativeIndex() < inst.nextJunction.GetRelativeIndex())
						{
							inst.nextJunction = it;
							inst.nextFlankDistance = endVertexDistance;
							break;
						}												
					}
				}

				if (nowNewInstance)
				{
					newInstance.push_back(it);
				}
			}

			for (auto & it : instance_)
			{
				int64_t nextRightFlank = abs(endVertexDistance - it.nextFlankDistance);
				int64_t length = abs(it.seq.front().GetPosition() - it.seq.back().GetPosition());				
				if (length >= minChainSize_ && nextRightFlank > maxFlankingSize_)
				{
					return false;
				}
			}

			for (auto & inst : instance_)
			{
				if (inst.nextJunction.GetIndex() != storage_->End(inst.seq.front().GetChrId()).GetIndex())
				{					
					inst.seq.push_back(inst.nextJunction);
					inst.rightFlankDistance = endVertexDistance;
				}
			}
				
			for (auto & it : newInstance)
			{
				instance_.push_back(Instance());
				instance_.back().seq.push_back(it);
				instance_.back().leftFlankDistance = instance_.back().rightFlankDistance = endVertexDistance;
			}

			body_.push_back(Point(e, startVertexDistance));
			assert(IsConsistent());
			return true;
		}

		bool PointPushFront(const Edge & e)
		{
			int64_t vertex = e.GetStartVertex();
			if (FindVertexInPath(vertex) != body_.end())
			{
				return false;
			}

			int64_t endVertexDistance = body_.empty() ? 0 : body_.front().startDistance;
			int64_t startVertexDistance = endVertexDistance - e.GetLength();
			for (auto & inst : instance_)
			{
				inst.nextFlankDistance = inst.leftFlankDistance;
				inst.nextJunction = storage_->End(inst.seq.back().GetChrId(), inst.seq.back().IsPositiveStrand());
			}

			FillBuffer(vertex);
			std::vector<JunctionStorage::JunctionIterator> newInstance;
			for (auto it : junctionBuffer_)
			{
				bool nowNewInstance = true;
				for (auto & inst : instance_)
				{
					if (Compatible(inst.seq.back(), it, e))
					{
						nowNewInstance = false;
						int64_t rightFlank = abs(inst.rightFlankDistance - body_.back().endDistance);
						int64_t nextLength = abs(it.GetPosition() - inst.seq.back().GetPosition());
						if (nextLength >= minChainSize_ && rightFlank > maxFlankingSize_)
						{
							return false;
						}
						
						if (it.GetRelativeIndex() > inst.nextJunction.GetRelativeIndex() || inst.nextJunction.GetIndex() == storage_->End(inst.seq.front().GetChrId()).GetIndex())
						{
							inst.nextJunction = it;
							inst.nextFlankDistance = startVertexDistance;
							break;
						}
					}
				}

				if (nowNewInstance)
				{
					newInstance.push_back(it);
				}
			}

			for (auto & it : instance_)
			{
				int64_t length = abs(it.seq.front().GetPosition() - it.seq.back().GetPosition());
				int64_t nextLeftFlank = abs(startVertexDistance - it.nextFlankDistance);
				if (length >= minChainSize_ && nextLeftFlank > maxFlankingSize_)
				{
					return false;
				}
			}

			for (auto & inst : instance_)
			{
				if (inst.nextJunction.GetIndex() != storage_->End(inst.seq.front().GetChrId()).GetIndex())
				{
					inst.seq.push_front(inst.nextJunction);
					inst.leftFlankDistance = startVertexDistance;
				}
			}

			for (auto & it : newInstance)
			{
				instance_.push_back(Instance());
				instance_.back().seq.push_back(it);
				instance_.back().leftFlankDistance = instance_.back().rightFlankDistance = startVertexDistance;
			}

			body_.push_front(Point(e, startVertexDistance));
			assert(IsConsistent());
			return true;
		}		

		void PointPopFront()
		{
			int64_t lastVertex = body_.front().edge.GetStartVertex();
			body_.pop_front();
			for (auto it = instance_.begin(); it != instance_.end();)
			{
				if (it->seq.front().GetVertexId() == lastVertex)
				{
					JunctionStorage::JunctionIterator & j = it->seq.front();					
					it->seq.pop_front();
					if (it->seq.empty())
					{
						it = instance_.erase(it);
					}
					else
					{
						for (auto jt = body_.begin(); ; ++jt)
						{
							if (jt->edge.GetStartVertex() == it->seq.front().GetVertexId())
							{
								it++->leftFlankDistance = jt->startDistance;
								break;
							}
						}						
					}
				}
				else
				{
					++it;
				}
			}

			assert(IsConsistent());
		}	

		void PointPopBack()
		{
			int64_t lastVertex = body_.back().edge.GetEndVertex();
			body_.pop_back();
			for (auto it = instance_.begin(); it != instance_.end();)
			{
				if (it->seq.back().GetVertexId() == lastVertex)
				{
					JunctionStorage::JunctionIterator & j = it->seq.back();					
					it->seq.pop_back();
					if (it->seq.empty())
					{
						it = instance_.erase(it);						
					}
					else
					{
						for (auto jt = body_.rbegin(); ; ++jt)
						{
							if (jt->edge.GetEndVertex() == it->seq.back().GetVertexId())
							{
								it++->rightFlankDistance = jt->endDistance;
								break;
							}
						}			
					}
				}
				else
				{
					++it;
				}
			}

			assert(IsConsistent());
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
			int64_t leftFlank = abs(inst.leftFlankDistance - body_.front().startDistance);
			int64_t rightFlank = abs(inst.rightFlankDistance - body_.back().endDistance);
			length = abs(inst.seq.front().GetPosition() - inst.seq.back().GetPosition());
			score = length - leftFlank - rightFlank;
		}

		bool IsConsistent() const
		{
			if (checkConsistency_)
			{
				Path testPath(body_.front().edge.GetStartVertex(), *storage_, maxBranchSize_, minBlockSize_, maxFlankingSize_, *blockId_, false);
				for (auto it = body_.begin() + 1; it != body_.end(); ++it)
				{
					if (!testPath.PointPushBack(it->edge))
					{
						testPath.PointPushBack(it->edge);
						assert(false);
					}
				}

				assert(body_ == testPath.body_);				
				assert(instance_ == testPath.instance_);				
			}			

			return true;
		}

	private:

		bool Compatible(const JunctionStorage::JunctionIterator & start, const JunctionStorage::JunctionIterator & end, const Edge & e) const
		{
			if (start.GetChrId() != end.GetChrId() || start.IsPositiveStrand() != end.IsPositiveStrand())
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

				auto start1 = start + 1;
				if (diff > maxBranchSize_ && (start.GetChar() != e.GetChar() || end != start1 || start1.GetVertexId() != e.GetEndVertex()))
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
				if (-diff > maxBranchSize_ && (start.GetChar() != e.GetChar() || end != start1 || start1.GetVertexId() != e.GetEndVertex()))
				{
					return false;
				}
			}

			return true;
		}

		friend class BestPath;
		
		std::deque<Point> body_;
		std::list<Instance> instance_;
		int64_t origin_;
		int64_t minChainSize_;
		int64_t minBlockSize_;
		int64_t maxBranchSize_;
		int64_t maxFlankingSize_;
		bool checkConsistency_;
		const JunctionStorage * storage_;		
		const std::vector<std::vector<Assignment> > * blockId_;
		std::vector<JunctionStorage::JunctionIterator> junctionBuffer_;
		
		void FillBuffer(int64_t vertex)
		{
			junctionBuffer_.clear();
			for (size_t i = 0; i < storage_->GetInstancesCount(vertex); i++)
			{
				JunctionStorage::JunctionIterator it = storage_->GetJunctionInstance(vertex, i);
				if ((*blockId_)[it.GetChrId()][it.GetIndex()].block != Assignment::UNKNOWN_BLOCK)
				{
					continue;
				}

				bool within = false;
				for (auto & instance : instance_)
				{
					int64_t chrId = instance.seq.front().GetChrId();
					int64_t leftIdx = instance.seq.front().GetIndex();
					int64_t rightIdx = instance.seq.back().GetIndex();
					if (it.GetChrId() == chrId && it.GetIndex() >= leftIdx && it.GetIndex() <= rightIdx)
					{
						within = true;
						break;
					}
				}

				if (!within)
				{
					junctionBuffer_.push_back(it);
				}
			}
		}		

		std::deque<Point>::const_iterator FindVertexInPath(int64_t vertex) const
		{
			for (auto it = body_.begin(); it != body_.end(); ++it)
			{
				if (it->edge.GetEndVertex() == vertex || it->edge.GetStartVertex() == vertex)
				{
					return it;
				}
			}

			return body_.end();
		}
	};

	struct BestPath
	{
		int64_t score;		
		int64_t leftFlank;
		int64_t rightFlank;
		std::deque<Edge> left;
		std::deque<Edge> right;

		void FixForward(Path & path)
		{			
			auto it = right.begin() + rightFlank;
			for (; it != right.end(); ++it)
			{
				assert(path.PointPushBack(*it));
				/*
				if (!path.PointPushBack(*it))
				{
					path.PointPushBack(*it);
					Path testPath(right.front(), *path.storage_, path.maxBranchSize_, path.minBlockSize_, path.maxFlankingSize_, *path.blockId_);
					for (auto it = right.begin() + 1; it != right.end(); ++it)
					{
						if (!testPath.PointPushBack(*it))
						{
							assert(false);
						}
					}
				}
				*/
			}

			rightFlank = right.size();
			assert(left.front() == right.front());
		}

		void FixBackward(Path & path)
		{
			auto it = left.begin() + leftFlank;
			for (; it != left.end(); ++it)
			{
				assert(path.PointPushFront(*it));				
			}

			leftFlank = left.size();
			assert(left.front() == right.front());
		}

		void UpdateForward(const Path & path, int64_t newScore)
		{
			score = newScore;
			right.erase(right.begin() + rightFlank, right.end());
			auto it = path.PathBody().end();
			for (--it; it->edge != right[rightFlank - 1]; --it);
			for (++it; it != path.PathBody().end(); ++it)
			{
				right.push_back(it->edge);
			}

			assert(left.front() == right.front());
		}

		void UpdateBackward(const Path & path, int64_t newScore)
		{
			score = newScore;
			left.erase(left.begin() + leftFlank, left.end());
			auto it = path.PathBody().rend();			

			for (--it; it->edge != left[leftFlank - 1]; --it);
			for (++it; it != path.PathBody().rend(); ++it)
			{
				left.push_back(it->edge);
			}

			assert(left.front() == right.front());
		}

		BestPath(int64_t vid) : score(0)
		{			
			rightFlank = right.size();
			leftFlank = left.size();
		}
	};
}

#endif

