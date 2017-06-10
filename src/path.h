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
		Path(Edge start,
			const JunctionStorage & storage,
			int64_t maxBranchSize,
			int64_t minBlockSize,
			int64_t maxFlankingSize,
			const std::vector<std::vector<Assignment> > & blockId) :
			maxBranchSize_(maxBranchSize), minBlockSize_(minBlockSize), maxFlankingSize_(maxFlankingSize), storage_(&storage),
			minChainSize_(minBlockSize - 2 * maxFlankingSize), blockId_(&blockId)
		{
			PointPushBack(start);
		}

		struct Instance
		{
			int64_t leftFlankDistance;
			int64_t rightFlankDistance;
			std::deque<JunctionStorage::JunctionIterator> seq;
		};

		struct Point
		{
			Edge edge;
			int64_t distance;
			Point() {}
			Point(Edge edge, int64_t distance) : edge(edge), distance(distance) {}
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
				//out << pt.vertex << ' ';
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

			FillBuffer(vertex);
			for (auto it : junctionBuffer_)
			{
				j = 0;			
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

			for (auto it : junctionBuffer_)
			{
				bool newInstance = true;				
				for (auto & inst : instance_)
				{
					if (Compatible(inst.seq.back(), it, e))
					{
						newInstance = false;
						inst.seq.push_back(it);
						inst.rightFlankDistance = vertexDistance;						
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

			body_.push_back(Point(e, vertexDistance));
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

			FillBuffer(vertex);
			for (auto it : junctionBuffer_)
			{
				j = 0;				
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

			for (auto it : junctionBuffer_)
			{
				bool newInstance = true;			
				for (auto & inst : instance_)
				{
					if (Compatible(it, inst.seq.front(), e))
					{
						newInstance = false;
						inst.seq.push_front(it);
						inst.leftFlankDistance = vertexDistance;						
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

			body_.push_front(Point(e, vertexDistance));
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

		Edge GetEdge(size_t index) const
		{
			return body_[index].edge;
		}

		void PointPopFront()
		{
			int64_t newLeftFlankDistance = (++body_.begin())->distance;
			for (auto it = instance_.begin(); it != instance_.end();)
			{
				if (it->seq.front().GetVertexId() == body_.front().edge.GetStartVertex())
				{
					JunctionStorage::JunctionIterator & j = it->seq.front();					
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
				if (it->seq.back().GetVertexId() == body_.back().edge.GetEndVertex())
				{
					JunctionStorage::JunctionIterator & j = it->seq.back();					
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
		
		std::deque<Point> body_;
		std::list<Instance> instance_;
		int64_t minChainSize_;
		int64_t minBlockSize_;
		int64_t maxBranchSize_;
		int64_t maxFlankingSize_;
		std::vector<int64_t> nextFlank_;
		std::vector<JunctionStorage::JunctionIterator> junctionBuffer_;
		const JunctionStorage * storage_;
		const std::vector<std::vector<Assignment> > * blockId_;

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
		int64_t leftFlankVertex;
		int64_t rightFlankVertex;
		std::deque<Edge> body;
		void Fix()
		{
			if (body.size() > 0)
			{
				leftFlankVertex = body.front().GetStartVertex();
				rightFlankVertex = body.back().GetEndVertex();
			}
		}

		void UpdateForward(const Path & path, int64_t newScore)
		{
			score = newScore;
			while (body.size() > 0 && body.back().GetEndVertex() != rightFlankVertex)
			{
				body.pop_back();
			}

			if (body.size() > 0)
			{
				auto it = --path.PathBody().end();
				for (; it->edge.GetEndVertex() != rightFlankVertex; --it);
				for (++it; it != path.PathBody().end(); ++it)
				{
					body.push_back(it->edge);
				}
			}
			else
			{
				for (auto it = path.PathBody().begin(); it != path.PathBody().end(); ++it)
				{
					body.push_back(it->edge);
				}
			}			
		}

		void UpdateBackward(const Path & path, int64_t newScore)
		{
			score = newScore;
			while (body.size() > 0 && body.front().GetStartVertex() != leftFlankVertex)
			{
				body.pop_front();
			}

			if (body.size() > 0)
			{
				auto it = path.PathBody().begin();
				for (; it->edge.GetStartVertex() != leftFlankVertex; ++it);
				for (++it; it != path.PathBody().end(); ++it)
				{
					body.push_front(it->edge);
				}
			}
			else
			{
				for (auto it = path.PathBody().begin(); it != path.PathBody().end(); ++it)
				{
					body.push_back(it->edge);
				}
			}			
		}

		BestPath() : score(0), leftFlankVertex(0), rightFlankVertex(0) {}
	};
}

#endif

