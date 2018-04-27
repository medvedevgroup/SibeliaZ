#ifndef _TRASERVAL_H_
#define _TRAVERSAL_H_

//#define _DEBUG_OUT_

#include <set>
#include <map>
#include <list>
#include <ctime>
#include <queue>
#include <iterator>
#include <cassert>
#include <numeric>
#include <sstream>
#include <iostream>
#include <functional>
#include <unordered_map>

#include <tbb/mutex.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>

#include "path.h"

namespace Sibelia
{
	extern const std::string DELIMITER;

	class BlockInstance
	{
	public:
		BlockInstance() {}
		BlockInstance(int id, const size_t chr, size_t start, size_t end) : id_(id), chr_(chr), start_(start), end_(end) {}
		void Reverse();
		int GetSignedBlockId() const;
		bool GetDirection() const;
		int GetBlockId() const;
		int GetSign() const;
		size_t GetChrId() const;
		size_t GetStart() const;
		size_t GetEnd() const;
		size_t GetLength() const;
		size_t GetConventionalStart() const;
		size_t GetConventionalEnd() const;
		std::pair<size_t, size_t> CalculateOverlap(const BlockInstance & instance) const;
		bool operator < (const BlockInstance & toCompare) const;
		bool operator == (const BlockInstance & toCompare) const;
		bool operator != (const BlockInstance & toCompare) const;
	private:
		int id_;
		size_t start_;
		size_t end_;
		size_t chr_;
	};

	namespace
	{
		const bool COVERED = true;
		typedef std::vector<BlockInstance> BlockList;
		typedef std::pair<size_t, std::vector<BlockInstance> > GroupedBlock;
		typedef std::vector<GroupedBlock> GroupedBlockList;
		bool ByFirstElement(const GroupedBlock & a, const GroupedBlock & b)
		{
			return a.first < b.first;
		}

		std::string IntToStr(size_t x)
		{
			std::stringstream ss;
			ss << x;
			return ss.str();
		}

		template<class Iterator1, class Iterator2>
		void CopyN(Iterator1 it, size_t count, Iterator2 out)
		{
			for (size_t i = 0; i < count; i++)
			{
				*out++ = *it++;
			}
		}

		template<class Iterator>
		Iterator AdvanceForward(Iterator it, size_t step)
		{
			std::advance(it, step);
			return it;
		}

		template<class Iterator>
		Iterator AdvanceBackward(Iterator it, size_t step)
		{
			for (size_t i = 0; i < step; i++)
			{
				--it;
			}

			return it;
		}


		typedef std::pair<size_t, size_t> IndexPair;
		template<class T, class F, class It>
		void GroupBy(std::vector<T> & store, F pred, It out)
		{
			sort(store.begin(), store.end(), pred);
			for (size_t now = 0; now < store.size();)
			{
				size_t prev = now;
				for (; now < store.size() && !pred(store[prev], store[now]); now++);
				*out++ = std::make_pair(prev, now);
			}
		}

		template<class F>
		bool CompareBlocks(const BlockInstance & a, const BlockInstance & b, F f)
		{
			return (a.*f)() < (b.*f)();
		}

		template<class F>
		bool EqualBlocks(const BlockInstance & a, const BlockInstance & b, F f)
		{
			return f(a) == f(b);
		}

		template<class Iterator, class F, class ReturnType>
		struct FancyIterator : public std::iterator<std::forward_iterator_tag, ReturnType>
		{
		public:
			FancyIterator& operator++()
			{
				++it;
				return *this;
			}

			FancyIterator operator++(int)
			{
				FancyIterator ret(*this);
				++(*this);
				return ret;
			}

			bool operator == (FancyIterator toCompare) const
			{
				return it == toCompare.it;
			}

			bool operator != (FancyIterator toCompare) const
			{
				return !(*this == toCompare);
			}

			ReturnType operator * ()
			{
				return f(*it);
			}

			FancyIterator() {}
			FancyIterator(Iterator it, F f) : it(it), f(f) {}

		private:
			F f;
			Iterator it;
		};

		template<class Iterator, class F, class ReturnType>
		FancyIterator<Iterator, F, ReturnType> CFancyIterator(Iterator it, F f, ReturnType)
		{
			return FancyIterator<Iterator, F, ReturnType>(it, f);
		}

	}

	bool compareById(const BlockInstance & a, const BlockInstance & b);
	bool compareByChrId(const BlockInstance & a, const BlockInstance & b);
	bool compareByStart(const BlockInstance & a, const BlockInstance & b);

	void CreateOutDirectory(const std::string & path);

	class BlocksFinder
	{
	public:

		BlocksFinder(JunctionStorage & storage, size_t k) : storage_(storage), k_(k)
		{
			scoreFullChains_ = true;
		}		


		void MissingSet(const std::string & fileName, std::set<int64_t> & result) const
		{
			std::string buf;
			std::ifstream oldCoordsIn(fileName);
			if (oldCoordsIn)
			{
				std::ofstream missingDot("missing.dot");
				while (std::getline(oldCoordsIn, buf) && buf[0] != '-')
				{
					std::string seq;
					std::stringstream ss(buf);
					char sign;
					int seqId, start, length, end, seqSize;
					ss >> seq >> seq >> start >> length >> sign >> seqSize;
					seqId = atoi(seq.substr(2).c_str()) - 1;
					end = start + length;
					if (sign == '-')
					{
						start = seqSize - start;
						end = seqSize - end;
						std::swap(start, end);
						assert(start < end);
					}

					for (auto it = storage_.Begin(seqId); it.Valid(); ++it)
					{
						int64_t pos = it.GetPosition();
						if (pos >= start && pos < end)
						{
							result.insert(it.GetVertexId());
							result.insert(-it.GetVertexId());
						}
					}

				}
			}
		}

		static bool DegreeCompare(const JunctionStorage & storage, int64_t v1, int64_t v2)
		{
			return storage.GetInstancesCount(v1) > storage.GetInstancesCount(v2);
		}

		struct ProcessVertexDijkstra
		{
		public:
			BlocksFinder & finder;
			std::vector<int64_t> & shuffle;

			ProcessVertexDijkstra(BlocksFinder & finder, std::vector<int64_t> & shuffle) : finder(finder), shuffle(shuffle)
			{
			}

			void operator()(tbb::blocked_range<size_t> & range) const
			{
				BubbledBranches forwardBubble;
				BubbledBranches backwardBubble;
				std::vector<JunctionStorage::JunctionSequentialIterator> instance;
				for (size_t i = range.begin(); i != range.end(); i++)
				{
					if (finder.count_++ % 1000 == 0)
					{
						tbb::mutex::scoped_lock lock(finder.globalMutex_);
						std::cout << finder.count_ << '\t' << shuffle.size() << std::endl;
					}

					instance.clear();
					int64_t vertex = shuffle[i];
					for (auto it = JunctionStorage::JunctionIterator(vertex); it.Valid(); ++it)
					{
						instance.push_back(it.SequentialIterator());
					}

					finder.BubbledBranchesForward(vertex, instance, forwardBubble);
					finder.BubbledBranchesBackward(vertex, instance, backwardBubble);
					bool isFork = false;
					for (size_t i = 0; i < forwardBubble.size() && !isFork; i++)
					{
						for (size_t j = 0; j < forwardBubble[i].size() && !isFork; j++)
						{
							size_t k = forwardBubble[i][j];
							if (std::find(backwardBubble[i].begin(), backwardBubble[i].end(), k) == backwardBubble[i].end())
							{
								finder.source_.push_back(vertex);
								isFork = true;
							}
						}
					}
				}
			}
		};

		static void AddIfNotExists(std::vector<int64_t> & adj, int64_t value)
		{
			if (std::find(adj.begin(), adj.end(), value) == adj.end())
			{
				adj.push_back(value);
			}
		}

		void FindBlocks(int64_t minBlockSize, int64_t maxBranchSize, int64_t lookingDepth, int64_t sampleSize, int64_t threads, const std::string & debugOut)
		{
			blocksFound_ = 0;
			sampleSize_ = sampleSize;
			lookingDepth_ = lookingDepth;
			minBlockSize_ = minBlockSize;
			maxBranchSize_ = maxBranchSize;
			blockId_.resize(storage_.GetChrNumber());
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				blockId_[i].resize(storage_.GetChrVerticesCount(i));
			}

			std::vector<int64_t> shuffle;
			for (int64_t v = -storage_.GetVerticesNumber() + 1; v < storage_.GetVerticesNumber(); v++)
			{
				for (JunctionStorage::JunctionIterator it(v); it.Valid(); ++it)
				{
					if (it.IsPositiveStrand())
					{
						shuffle.push_back(v);
						break;
					}
				}
			}
			
#ifdef _DEBUG_OUT_
			MissingSet("missing.maf", missingVertex_);
			std::ofstream missingDot("missing.dot");
			missingDot << "digraph G\n{\nrankdir = LR" << std::endl;
			std::vector<std::pair<JunctionStorage::JunctionSequentialIterator, JunctionStorage::JunctionSequentialIterator> > vvisit;
			for (auto vid : missingVertex_)
			{
				DumpVertex(vid, missingDot, vvisit, 10);
				missingDot << vid << "[shape=square]" << std::endl;
			}

			missingDot << "}" << std::endl;
#endif
			count_ = 0;
			time_t mark = time(0);
			tbb::task_scheduler_init init(threads);
			tbb::parallel_for(tbb::blocked_range<size_t>(0, shuffle.size()), ProcessVertexDijkstra(*this, shuffle));
			std::cout << source_.size() << ' ' << shuffle.size() << std::endl;
			
			std::vector<std::vector<int64_t> > adj(source_.size());
			{
				std::vector<std::vector<int64_t> > visit(source_.size());
				std::vector<int64_t> sourceLabel(storage_.GetVerticesNumber() * 2 + 1, -1);
				for (size_t i = 0; i < source_.size(); i++)
				{
					sourceLabel[source_[i] + storage_.GetVerticesNumber()] = i;
				}
				
				for (size_t chr = 0; chr < storage_.GetChrNumber(); chr++)
				{
					for (size_t strand = 0; strand < 2; strand++)
					{
						for (JunctionStorage::JunctionSequentialIterator it = storage_.Begin(chr, strand == 0); it.Valid(); ++it)
						{
							int64_t itSrcId = sourceLabel[it.GetVertexId() + storage_.GetVerticesNumber()];
							if (itSrcId == -1)
							{
								continue;
							}

							for (auto jt = it + 1; jt.Valid() && abs(it.GetPosition() - jt.GetPosition()) <= maxBranchSize_; ++jt)
							{
								int64_t jtSrcId = sourceLabel[jt.GetVertexId() + storage_.GetVerticesNumber()];
								if (jtSrcId == -1)
								{
									continue;
								}

								AddIfNotExists(adj[itSrcId], jtSrcId);
								AddIfNotExists(adj[jtSrcId], itSrcId);
								AddIfNotExists(visit[jtSrcId], itSrcId);
							}
						}
					}
				}

				for (auto & v : visit)
				{
					for (size_t i = 0; i < v.size(); i++)
					{
						for (size_t j = 0; j < v.size(); j++)
						{
							AddIfNotExists(adj[v[i]], v[j]);
							AddIfNotExists(adj[v[i]], v[j]);
						}
					}
				}
			}

			int components = 0;
			std::vector<int> visit(source_.size(), -1);
			{
				for (int v = 0; v < source_.size(); v++)
				{
					if (visit[v] == -1)
					{
						auto nowComponent = components++;
						std::vector<int64_t> stack;
						stack.push_back(v);
						while (stack.size() > 0)
						{
							auto now = stack.back();
							stack.pop_back();
							if (visit[now] == -1)
							{
								visit[now] = nowComponent;
								for (auto next : adj[now])
								{
									if (visit[next] == -1)
									{
										stack.push_back(next);
									}
								}
							}
						}
					}
				}
			}

			std::vector<std::vector<int64_t> > componentContent(components);
			std::cout << "Components: " << components << std::endl;
			for (size_t i = 0; i < visit.size(); i++)
			{
				componentContent[visit[i]].push_back(i);
			}


			for(int i = 0; i < components; i++)
			{
				std::stringstream fn;
				fn << "comps/" << i << ".dot";
				std::ofstream out(fn.str().c_str());
				out << "strict graph G\n{\nrankdir = LR" << std::endl;
				fn.str("");
				fn << "comps/" << i << "_graph.dot";
				std::ofstream grout(fn.str().c_str());
				grout << "digraph G\n{\nrankdir = LR" << std::endl;

				std::vector<std::pair<JunctionStorage::JunctionSequentialIterator, JunctionStorage::JunctionSequentialIterator> > vvisit;
				for (auto v : componentContent[i])
				{
					auto vLabel = source_[v];
					out << vLabel << std::endl;
					for (auto u : adj[v])
					{
						auto uLabel = source_[u];
						out << vLabel << " -- " << uLabel << std::endl;
					}

					DumpVertex(vLabel, grout, vvisit, 10);
					grout << vLabel << "[shape=square]" << std::endl;
				}

				out << "}" << std::endl;
				grout << "}" << std::endl;
			}
		}

		void Dump(std::ostream & out) const
		{
			out << "digraph G\n{\nrankdir = LR" << std::endl;
			for (size_t i = 0; i < storage_.GetChrNumber(); i++)
			{
				auto end = storage_.End(i).Prev();
				for (auto it = storage_.Begin(i); it != end; ++it)
				{
					auto jt = it.Next();
					out << it.GetVertexId() << " -> " << jt.GetVertexId() << "[label=\"" << it.GetChar() << ", " << it.GetChrId() << ", " << it.GetPosition() << "\" color=blue]\n";
					out << jt.Reverse().GetVertexId() << " -> " << it.Reverse().GetVertexId() << "[label=\"" << it.GetChar() << ", " << it.GetChrId() << ", " << it.GetPosition() << "\" color=red]\n";
				}
			}

			for (size_t i = 0; i < syntenyPath_.size(); i++)
			{
				for (size_t j = 0; j < syntenyPath_[i].size(); j++)
				{
					Edge e = syntenyPath_[i][j];
					out << e.GetStartVertex() << " -> " << e.GetEndVertex() <<
						"[label=\"" << e.GetChar() << ", " << i + 1 << "\" color=green]\n";
					e = e.Reverse();
					out << e.GetStartVertex() << " -> " << e.GetEndVertex() <<
						"[label=\"" << e.GetChar() << ", " << -(int64_t(i + 1)) << "\" color=green]\n";
				}
			}

			out << "}" << std::endl;
		}

		void GenerateLegacyOutput(const std::string & outDir) const;

	private:

		template<class Iterator>
		void OutputLines(Iterator start, size_t length, std::ostream & out) const
		{
			for (size_t i = 1; i <= length; i++, ++start)
			{
				out << *start;
				if (i % 80 == 0 && i != length)
				{
					out << std::endl;
				}
			}
		}

		void GenerateReport(const BlockList & block, const std::string & fileName) const;
		void CalculateCoverage(GroupedBlockList::const_iterator start, GroupedBlockList::const_iterator end, std::vector<bool> & cover, std::vector<double> & ret) const;
		std::string OutputIndex(const BlockInstance & block) const;
		void OutputBlocks(const std::vector<BlockInstance>& block, std::ofstream& out) const;
		void ListBlocksIndices(const BlockList & block, const std::string & fileName) const;
		void ListBlocksSequences(const BlockList & block, const std::string & directory) const;
		void ListChromosomesAsPermutations(const BlockList & block, const std::string & fileName) const;
		void TryOpenFile(const std::string & fileName, std::ofstream & stream) const;
		void ListChrs(std::ostream & out) const;
		
		template<class T>
		void DumpVertex(int64_t id, std::ostream & out, T & visit, int64_t cnt = 5) const
		{
			for (auto kt = JunctionStorage::JunctionIterator(id); kt.Valid(); ++kt)
			{
				auto jt = kt.SequentialIterator();
				for (int64_t i = 0; i < cnt; i++)
				{
					auto it = jt - 1;
					auto pr = std::make_pair(it, jt);
					if (it.Valid() && std::find(visit.begin(), visit.end(), pr) == visit.end())
					{
						int64_t length = it.GetPosition() - jt.GetPosition();
						out << it.GetVertexId() << " -> " << jt.GetVertexId()
							<< "[label=\"" << it.GetChar() << ", " << it.GetChrId() << ", " << it.GetPosition() << "," << length << "\""
							<< (it.IsPositiveStrand() ? "color=blue" : "color=red") << "]\n";
						visit.push_back(pr);
					}

					jt = it;
				}
			}

			for (auto kt = JunctionStorage::JunctionIterator(id); kt.Valid(); ++kt)
			{
				auto it = kt.SequentialIterator();
				for (int64_t i = 0; i < cnt; i++)
				{
					auto jt = it + 1;
					auto pr = std::make_pair(it, jt);
					if (jt.Valid() && std::find(visit.begin(), visit.end(), pr) == visit.end())
					{
						int64_t length = it.GetPosition() - jt.GetPosition();
						out << it.GetVertexId() << " -> " << jt.GetVertexId()
							<< "[label=\"" << it.GetChar() << ", " << it.GetChrId() << ", " << it.GetPosition() << "," << length << "\""
							<< (it.IsPositiveStrand() ? "color=blue" : "color=red") << "]\n";
						visit.push_back(pr);
					}

					it = jt;
				}
			}
		}

		struct BranchData
		{
			std::vector<size_t> branchId;
		};

		typedef std::vector<std::vector<size_t> > BubbledBranches;

		struct Fork
		{
			Fork(JunctionStorage::JunctionSequentialIterator it, JunctionStorage::JunctionSequentialIterator jt)
			{
				if (it < jt)
				{
					branch[0] = it;
					branch[1] = jt;
				}
				else
				{
					branch[0] = jt;
					branch[1] = it;
				}
			}

			JunctionStorage::JunctionSequentialIterator branch[2];
		};

		int64_t ChainLength(const Fork & now, const Fork & next) const
		{
			return min(abs(now.branch[0].GetPosition() - next.branch[0].GetPosition()), abs(now.branch[1].GetPosition() - next.branch[1].GetPosition()));
		}

		Fork ExpandSourceFork(const Fork & source) const
		{
			for (auto now = source; ; )
			{
				auto next = TakeBubbleStep(now);
				if (next.branch[0].Valid())
				{
					int64_t vid0 = now.branch[0].GetVertexId();
					int64_t vid1 = now.branch[1].GetVertexId();
					assert(vid0 == vid1 && abs(now.branch[0].GetPosition() - next.branch[0].GetPosition()) < maxBranchSize_ && 
						abs(now.branch[1].GetPosition() - next.branch[1].GetPosition()) < maxBranchSize_);
					now = next;
				}
				else
				{
					return now;
				}
			}

			return source;
		}

		Fork TakeBubbleStep(const Fork & source) const
		{
			auto it = source.branch[0];
			std::map<int64_t, int64_t> firstBranch;
			for (int64_t i = 1; abs(it.GetPosition() - source.branch[0].GetPosition()) < maxBranchSize_ && (++it).Valid(); i++)
			{
				int64_t d = abs(it.GetPosition() - source.branch[0].GetPosition());
				firstBranch[it.GetVertexId()] = i;
			}

			it = source.branch[1];
			for (int64_t i = 1; abs(it.GetPosition() - source.branch[1].GetPosition()) < maxBranchSize_ && (++it).Valid(); i++)
			{
				auto kt = firstBranch.find(it.GetVertexId());
				if (kt != firstBranch.end())
				{
					return Fork(source.branch[0] + kt->second, it);
				}
			}

			return Fork(JunctionStorage::JunctionSequentialIterator(), JunctionStorage::JunctionSequentialIterator());
		}

		void BubbledBranchesForward(int64_t vertexId, const std::vector<JunctionStorage::JunctionSequentialIterator> & instance, BubbledBranches & bulges) const
		{
			std::vector<size_t> parallelEdge[5];
			std::map<int64_t, BranchData> visit;
			bulges.assign(instance.size(), std::vector<size_t>());
			for (size_t i = 0; i < instance.size(); i++)
			{
				auto vertex = instance[i];
				if ((vertex + 1).Valid())
				{
					parallelEdge[TwoPaCo::DnaChar::MakeUpChar(vertex.GetChar())].push_back(i);
				}

				for (int64_t startPosition = vertex++.GetPosition(); vertex.Valid() && abs(startPosition - vertex.GetPosition()) <= maxBranchSize_; ++vertex)
				{
					int64_t nowVertexId = vertex.GetVertexId();
					auto point = visit.find(nowVertexId);
					if (point == visit.end())
					{
						BranchData bData;
						bData.branchId.push_back(i);
						visit[nowVertexId] = bData;
					}
					else
					{
						point->second.branchId.push_back(i);
					}
				}
			}

			for (size_t i = 0; i < 5; i++)
			{
				for (size_t j = 0; j < parallelEdge[i].size(); j++)
				{
					for (size_t k = j + 1; k < parallelEdge[i].size(); k++)
					{
						size_t smallBranch = parallelEdge[i][j];
						size_t largeBranch = parallelEdge[i][k];
						bulges[smallBranch].push_back(largeBranch);
					}
				}
			}

			for (auto point = visit.begin(); point != visit.end(); ++point)
			{
				std::sort(point->second.branchId.begin(), point->second.branchId.end());
				for (size_t j = 0; j < point->second.branchId.size(); j++)
				{
					for (size_t k = j + 1; k < point->second.branchId.size(); k++)
					{
						size_t smallBranch = point->second.branchId[j];
						size_t largeBranch = point->second.branchId[k];
						if (smallBranch != largeBranch && std::find(bulges[smallBranch].begin(), bulges[smallBranch].end(), largeBranch) == bulges[smallBranch].end())
						{
							bulges[smallBranch].push_back(largeBranch);
						}
					}
				}
			}
		}

		void BubbledBranchesBackward(int64_t vertexId, const std::vector<JunctionStorage::JunctionSequentialIterator> & instance, BubbledBranches & bulges) const
		{
			std::vector<size_t> parallelEdge[5];
			std::map<int64_t, BranchData> visit;
			bulges.assign(instance.size(), std::vector<size_t>());
			for (size_t i = 0; i < instance.size(); i++)
			{
				auto vertex = instance[i];
				auto prev = vertex - 1;
				if (prev.Valid())
				{
					parallelEdge[TwoPaCo::DnaChar::MakeUpChar(prev.GetChar())].push_back(i);
				}

				for (int64_t startPosition = vertex--.GetPosition(); vertex.Valid() && abs(startPosition - vertex.GetPosition()) <= maxBranchSize_; --vertex)
				{
					int64_t nowVertexId = vertex.GetVertexId();
					auto point = visit.find(nowVertexId);
					if (point == visit.end())
					{
						BranchData bData;
						bData.branchId.push_back(i);
						visit[nowVertexId] = bData;
					}
					else
					{
						point->second.branchId.push_back(i);
					}
				}
			}

			for (size_t i = 0; i < 5; i++)
			{
				for (size_t j = 0; j < parallelEdge[i].size(); j++)
				{
					for (size_t k = j + 1; k < parallelEdge[i].size(); k++)
					{
						size_t smallBranch = parallelEdge[i][j];
						size_t largeBranch = parallelEdge[i][k];
						bulges[smallBranch].push_back(largeBranch);
					}
				}
			}

			for (auto point = visit.begin(); point != visit.end(); ++point)
			{
				std::sort(point->second.branchId.begin(), point->second.branchId.end());
				for (size_t j = 0; j < point->second.branchId.size(); j++)
				{
					for (size_t k = j + 1; k < point->second.branchId.size(); k++)
					{
						size_t smallBranch = point->second.branchId[j];
						size_t largeBranch = point->second.branchId[k];
						if (smallBranch != largeBranch && std::find(bulges[smallBranch].begin(), bulges[smallBranch].end(), largeBranch) == bulges[smallBranch].end())
						{
							bulges[smallBranch].push_back(largeBranch);
						}
					}
				}
			}
		}


		int64_t k_;
		std::atomic<int64_t> count_;
		std::atomic<int64_t> blocksFound_;
		int64_t sampleSize_;
		int64_t scalingFactor_;
		bool scoreFullChains_;
		int64_t lookingDepth_;
		int64_t minBlockSize_;
		int64_t maxBranchSize_;
		JunctionStorage & storage_;
		std::vector<std::vector<Edge> > syntenyPath_;
		std::vector<std::vector<Assignment> > blockId_;
		tbb::concurrent_vector<int64_t> source_;
		tbb::mutex globalMutex_;
#ifdef _DEBUG_OUT_
		bool debug_;
		std::set<int64_t> missingVertex_;
#endif
	};
}

#endif