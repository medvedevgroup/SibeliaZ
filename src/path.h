#ifndef _PATH_H_
#define _PATH_H_

#include "junctionstorage.h"

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


	struct Path
	{
	public:
		Path() {}
		Path(int64_t start,
			const JunctionStorage & storage,
			int64_t maxBranchSize,
			int64_t minBlockSize,
			int64_t maxFlankingSize,
			const std::vector<std::vector<Assignment> > & blockId) :
			maxBranchSize_(maxBranchSize), minBlockSize_(minBlockSize), maxFlankingSize_(maxFlankingSize), storage_(&storage),
			minChainSize_(minBlockSize - 2 * maxFlankingSize), blockId_(&blockId)
		{
			PointPushBack(Edge(0, start, 0, 0));
		}

		struct Instance
		{
			int64_t leftFlankDistance;
			int64_t rightFlankDistance;
			std::deque<JunctionStorage::JunctionIterator> seq;
		};

		struct Point
		{
			int64_t vertex;
			int64_t distance;
			Point() {}
			Point(int64_t vertex, int64_t distance) : vertex(vertex), distance(distance) {}
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
				out << pt.vertex << ' ';
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

		bool PointPushBack(const Edge & e)
		{
			int64_t vertex = e.GetEndVertex();
			if (FindVertexInPath(vertex) != body_.end())
			{
				return false;
			}

			size_t j = 0;
			int64_t vertexDistance = body_.empty() ? 0 : body_.back().distance + e.GetLength();
			for (auto & inst : instance_)
			{
				nextFlank_[j++] = inst.rightFlankDistance;
			}

			for (size_t i = 0; i < storage_->GetInstancesCount(vertex); i++)
			{
				j = 0;
				JunctionStorage::JunctionIterator it = storage_->GetJunctionInstance(vertex, i);
				if ((*blockId_)[it.GetChrId()][it.GetIndex()].block != Assignment::UNKNOWN_BLOCK ||
					inside_.count(std::make_pair(it.GetChrId(), it.GetIndex())))
				{
					continue;
				}

				for (auto & inst : instance_)
				{
					if (Compatible(inst.seq.back(), it, e))
					{
						int64_t leftFlank = abs(inst.leftFlankDistance - body_.front().distance);
						if (abs(it.GetPosition() - inst.seq.front().GetPosition()) >= minChainSize_ && leftFlank > maxFlankingSize_)
						{
							return false;
						}

						nextFlank_[j] = vertexDistance;
						break;
					}

					j++;
				}
			}

			auto it = instance_.begin();
			for (int64_t & rightFlank : nextFlank_)
			{
				if (abs(it->seq.front().GetPosition() - it->seq.back().GetPosition()) >= minChainSize_ && abs(vertexDistance - rightFlank) > maxFlankingSize_)
				{
					return false;
				}

				++it;
			}

			for (size_t i = 0; i < storage_->GetInstancesCount(vertex); i++)
			{
				bool newInstance = true;
				JunctionStorage::JunctionIterator it = storage_->GetJunctionInstance(vertex, i);
				if ((*blockId_)[it.GetChrId()][it.GetIndex()].block != Assignment::UNKNOWN_BLOCK ||
					inside_.count(std::make_pair(it.GetChrId(), it.GetIndex())))
				{
					continue;
				}

				for (auto & inst : instance_)
				{
					if (Compatible(inst.seq.back(), it, e))
					{
						newInstance = false;
						inst.seq.push_back(it);
						inst.rightFlankDistance = vertexDistance;
						inside_.insert(std::make_pair(it.GetChrId(), it.GetIndex()));
						break;
					}
				}

				if (newInstance)
				{
					nextFlank_.push_back(0);
					instance_.push_back(Instance());
					instance_.back().seq.push_back(it);
					instance_.back().leftFlankDistance = instance_.back().rightFlankDistance = vertexDistance;
				}
			}

			body_.push_back(Point(vertex, vertexDistance));
			return true;
		}

		bool PointPushFront(const Edge & e)
		{
			int64_t vertex = e.GetStartVertex();
			if (FindVertexInPath(vertex) != body_.end())
			{
				return false;
			}

			size_t j = 0;
			int64_t vertexDistance = body_.empty() ? 0 : body_.front().distance - e.GetLength();
			for (auto & inst : instance_)
			{
				nextFlank_[j++] = inst.leftFlankDistance;
			}

			for (size_t i = 0; i < storage_->GetInstancesCount(vertex); i++)
			{
				j = 0;
				JunctionStorage::JunctionIterator it = storage_->GetJunctionInstance(vertex, i);
				if ((*blockId_)[it.GetChrId()][it.GetIndex()].block != Assignment::UNKNOWN_BLOCK ||
					inside_.count(std::make_pair(it.GetChrId(), it.GetIndex())))
				{
					continue;
				}

				for (auto & inst : instance_)
				{
					if (Compatible(it, inst.seq.front(), e))
					{
						int64_t rightFlank = abs(inst.rightFlankDistance - body_.back().distance);
						if (abs(it.GetPosition() - inst.seq.back().GetPosition()) >= minChainSize_ && rightFlank > maxFlankingSize_)
						{
							return false;
						}

						nextFlank_[j] = vertexDistance;
						break;
					}

					j++;
				}
			}

			auto it = instance_.begin();
			for (int64_t leftFlank : nextFlank_)
			{
				if (abs(it->seq.front().GetPosition() - it->seq.back().GetPosition()) >= minChainSize_ && abs(vertexDistance - leftFlank) > maxFlankingSize_)
				{
					return false;
				}

				++it;
			}

			for (size_t i = 0; i < storage_->GetInstancesCount(vertex); i++)
			{
				bool newInstance = true;
				JunctionStorage::JunctionIterator it = storage_->GetJunctionInstance(vertex, i);
				if ((*blockId_)[it.GetChrId()][it.GetIndex()].block != Assignment::UNKNOWN_BLOCK ||
					inside_.count(std::make_pair(it.GetChrId(), it.GetIndex())))
				{
					continue;
				}

				for (auto & inst : instance_)
				{
					if (Compatible(it, inst.seq.front(), e))
					{
						newInstance = false;
						inst.seq.push_front(it);
						inst.leftFlankDistance = vertexDistance;
						inside_.insert(std::make_pair(it.GetChrId(), it.GetIndex()));
						break;
					}
				}

				if (newInstance)
				{
					nextFlank_.push_back(0);
					instance_.push_back(Instance());
					instance_.back().seq.push_back(it);
					instance_.back().leftFlankDistance = instance_.back().rightFlankDistance = vertexDistance;
				}
			}

			body_.push_front(Point(vertex, vertexDistance));
			return true;
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
			return body_.back().distance - body_.front().distance;
		}

		int64_t GetVertex(size_t index) const
		{
			return body_[index].vertex;
		}

		void PointPopFront()
		{
			int64_t newLeftFlankDistance = (++body_.begin())->distance;
			for (auto it = instance_.begin(); it != instance_.end();)
			{
				if (it->seq.front().GetVertexId() == body_.front().vertex)
				{
					JunctionStorage::JunctionIterator & j = it->seq.front();
					inside_.erase(std::make_pair(j.GetChrId(), j.GetIndex()));
					it->seq.pop_front();
					if (it->seq.empty())
					{
						it = instance_.erase(it);
						nextFlank_.pop_back();
					}
					else
					{
						++it->leftFlankDistance = newLeftFlankDistance;
					}
				}
				else
				{
					++it;
				}
			}

			body_.pop_front();
		}

		void PointPopBack()
		{
			int64_t newRightFlankDistance = (--(--body_.end()))->distance;
			for (auto it = instance_.begin(); it != instance_.end();)
			{
				if (it->seq.back().GetVertexId() == body_.back().vertex)
				{
					JunctionStorage::JunctionIterator & j = it->seq.back();
					inside_.erase(std::make_pair(j.GetChrId(), j.GetIndex()));
					it->seq.pop_back();
					if (it->seq.empty())
					{
						it = instance_.erase(it);
						nextFlank_.pop_back();
					}
					else
					{
						it++->rightFlankDistance = newRightFlankDistance;
					}
				}
				else
				{
					++it;
				}
			}

			body_.pop_back();
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
			int64_t leftFlank = abs(inst.leftFlankDistance - body_.front().distance);
			int64_t rightFlank = abs(inst.rightFlankDistance - body_.back().distance);
			length = abs(inst.seq.front().GetPosition() - inst.seq.back().GetPosition());
			score = length - leftFlank - rightFlank;
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

		std::set<std::pair<int64_t, int64_t> > inside_;
		std::deque<Point> body_;
		std::list<Instance> instance_;
		int64_t minChainSize_;
		int64_t minBlockSize_;
		int64_t maxBranchSize_;
		int64_t maxFlankingSize_;
		std::vector<int64_t> nextFlank_;
		const JunctionStorage * storage_;
		const std::vector<std::vector<Assignment> > * blockId_;

		std::deque<Point>::const_iterator FindVertexInPath(int64_t vertex) const
		{
			for (auto it = body_.begin(); it != body_.end(); ++it)
			{
				if (it->vertex == vertex)
				{
					return it;
				}
			}

			return body_.end();
		}
	};
}

#endif

