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

#include <tbb/tbb.h>
#include <tbb/parallel_for.h>

#include "path.h"

namespace Sibelia
{
	extern const std::string DELIMITER;
	extern const std::string VERSION;
	extern const size_t GAME_OVER;

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
			FancyIterator & operator++()
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

		struct Bundle
		{
			int64_t vid;
			char ch;
			size_t count = 0;

			bool operator < (const Bundle & a) const
			{
				return count > a.count;
			}
		};

		typedef std::vector<Path::Instance> * InstanceVectorPtr;

		BlocksFinder(JunctionStorage & storage, size_t k) : storage_(storage), k_(k)
		{
			progressCount_ = 50;
			scoreFullChains_ = true;
		}

		struct ProcessVertex
		{
		public:
			BlocksFinder & finder;

			ProcessVertex(BlocksFinder & finder) : finder(finder)
			{
			}

			void operator()()
			{
				size_t i;
				std::vector<size_t> data;
				std::vector<uint32_t> count(finder.storage_.GetVerticesNumber() * 2 + 1, 0);
				InstanceVectorPtr bestInstance;
				Path finalizer(finder.storage_, finder.maxBranchSize_, finder.minBlockSize_, finder.minBlockSize_, finder.maxFlankingSize_, true);
				Path currentPath(finder.storage_, finder.maxBranchSize_, finder.minBlockSize_, finder.minBlockSize_, finder.maxFlankingSize_, true);
				for(bool run = true; run; )
				{
					bool go = false;
					finder.taskMutex_.lock();
					if (finder.taskQueue_.size() > 0)
					{
						i = finder.taskQueue_.top();
						if (i != GAME_OVER)
						{
							finder.taskQueue_.pop();
							go = true;
						}
						else
						{
							run = false;
						}
					}

					finder.taskMutex_.unlock();

					if (!go)
					{
						continue;
					}

					bestInstance = new std::vector<Path::Instance>();

					int64_t score;
					int64_t vid = finder.bundle_[i].vid;
					char initChar = finder.bundle_[i].ch;
#ifdef _DEBUG_OUT_
					finder.debug_ = finder.missingVertex_.count(vid);
					if (finder.debug_)
					{
						std::cerr << "Vid: " << vid << std::endl;
					}
#endif
					{
						currentPath.Init(vid, initChar);

						int64_t bestScore = 0;
						size_t bestRightSize = currentPath.RightSize();
						size_t bestLeftSize = currentPath.LeftSize();
#ifdef _DEBUG_OUT_
						if (finder.debug_)
						{
							std::cerr << "Going forward:" << std::endl;
						}
#endif
						int64_t minRun = max(finder.minBlockSize_, finder.maxBranchSize_) * 2;
						while (true)
						{
							bool ret = true;
							bool positive = false;
							int64_t prevLength = currentPath.MiddlePathLength();
							while ((ret = finder.ExtendPathForward(currentPath, count, data, bestRightSize, bestScore, score, *bestInstance)) && currentPath.MiddlePathLength() - prevLength <= minRun)
							{
								positive = positive || (score > 0);
							}

							if (!ret || !positive)
							{
								break;
							}
						}

						{
							std::vector<Edge> bestEdge;
							for (size_t i = 0; i < bestRightSize - 1; i++)
							{
								bestEdge.push_back(currentPath.RightPoint(i).GetEdge());
							}

							currentPath.Clear();
							currentPath.Init(vid, initChar);
							for (auto & e : bestEdge)
							{
								currentPath.PointPushBack(e);
							}
						}
#ifdef _DEBUG_OUT_
						if (finder.debug_)
						{
							std::cerr << "Going backward:" << std::endl;
						}
#endif
						while (true)
						{
							bool ret = true;
							bool positive = false;
							int64_t prevLength = currentPath.MiddlePathLength();
							while ((ret = finder.ExtendPathBackward(currentPath, count, data, bestLeftSize, bestScore, score, *bestInstance)) && currentPath.MiddlePathLength() - prevLength <= minRun);
							{
								positive = positive || (score > 0);
							}

							if (!ret || !positive)
							{
								break;
							}
						}

						finder.resultMutex_.lock();
						finder.resultQueue_.push(std::make_pair(i, bestInstance));
						finder.resultMutex_.unlock();
						
						currentPath.Clear();
					}
				}
			}

			/*
			void FinalizeBlock(const Path & currentPath, Path & finalizer, size_t bestRightSize, size_t bestLeftSize, char initChar)
			{
				bool ret = false;

				finalizer.Init(currentPath.Origin(), initChar);
				for (size_t i = 0; i < bestRightSize - 1 && finalizer.PointPushBack(currentPath.RightPoint(i).GetEdge()); i++);
				for (size_t i = 0; i < bestLeftSize - 1 && finalizer.PointPushFront(currentPath.LeftPoint(i).GetEdge()); i++);

				int64_t finalScore = finalizer.Score();
				int64_t finalInstances = finalizer.GoodInstances();
				if (finalScore > 0 && finalInstances > 1)
				{
					ret = true;
					int64_t currentBlock = ++blocksFound_;
					for (auto jt : finalizer.AllInstances())
					{
						if (finalizer.IsGoodInstance(*jt))
						{
							blocksMutex_.lock();
							if (jt->Front().IsPositiveStrand())
							{
								blocksInstance_.push_back(BlockInstance(+currentBlock, jt->Front().GetChrId(), jt->Front().GetPosition(), jt->Back().GetPosition() + k_));
							}
							else
							{
								blocksInstance_.push_back(BlockInstance(-currentBlock, jt->Front().GetChrId(), jt->Back().GetPosition() - k_, jt->Front().GetPosition()));
							}

							blocksMutex_.unlock();
							for (auto it = jt->Front(); it != jt->Back(); ++it)
							{
								it.MarkUsed();
							}
						}
					}
				}

				finalizer.Clear();
			}*/
		};

		static bool DegreeCompare(const JunctionStorage & storage, int64_t v1, int64_t v2)
		{
			return storage.GetInstancesCount(v1) > storage.GetInstancesCount(v2);
		}

		void Split(std::string & source, std::vector<std::string> & result)
		{
			std::stringstream ss;
			ss << source;
			result.clear();
			while (ss >> source)
			{
				result.push_back(source);
			}
		}

		void FindBlocks(int64_t minBlockSize, int64_t maxBranchSize, int64_t maxFlankingSize, int64_t lookingDepth, int64_t sampleSize, int64_t threads, const std::string & debugOut)
		{
			sampleSize_ = sampleSize;
			lookingDepth_ = lookingDepth;
			minBlockSize_ = minBlockSize;
			maxBranchSize_ = maxBranchSize;
			maxFlankingSize_ = maxFlankingSize;

			for (int64_t v = -storage_.GetVerticesNumber() + 1; v < storage_.GetVerticesNumber(); v++)
			{
				std::set<char> good;
				std::map<char, size_t> count;
				for (JunctionStorage::JunctionIterator it(v); it.Valid(); ++it)
				{
					if (it.IsPositiveStrand())
					{
						good.insert(it.GetChar());
					}

					count[it.GetChar()] += 1;
				}

				for (auto p : count)
				{
					Bundle b;
					b.vid = v;
					b.ch = p.first;
					b.count = p.second;
					if (b.count > 1 && good.count(b.ch))
					{
						bundle_.push_back(b);
					}
				}
			}

			count_ = 0;
			std::cout << '[' << std::flush;
			progressPortion_ = bundle_.size() / progressCount_;
			if (progressPortion_ == 0)
			{
				progressPortion_ = 1;
			}

			std::sort(bundle_.begin(), bundle_.end());
			tbb::task_scheduler_init init(static_cast<int>(threads));

			currentBundle_ = 0;
			clock_t mark = clock();
			std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);

			size_t step = threads * 8;
			for (size_t task = 0; task < (step < bundle_.size() ? step : bundle_.size()); task++)
			{
				taskQueue_.push(task);
			}

			for (size_t i = 0; i < threads; i++)
			{
				ProcessVertex process(*this);
				workerThread[i].reset(new tbb::tbb_thread(process));
			}

			size_t top;
			blocksFound_ = 0;
			for (size_t task = 0; task < bundle_.size(); )
			{
				bool skip = true;
				InstanceVectorPtr instance = 0;
				resultMutex_.lock();
				if(resultQueue_.size() > 0 && (top = resultQueue_.top().first) == task)
				{
					skip = false;
					instance = resultQueue_.top().second;
					resultQueue_.pop();
				}

				resultMutex_.unlock();

				if (skip)
				{
					continue;
				}

				bool isGood = true;
				if (instance != 0)
				{
					for (auto inst : *instance)
					{
						for (auto it = inst.Front(); it != inst.Back(); ++it)
						{
							if (it.IsUsed())
							{
								isGood = false;
								break;
							}
						}

						if (!isGood)
						{
							break;
						}
					}
				}

				if (isGood)
				{
					if (count_++ % progressPortion_ == 0)
					{
						std::cout << '.' << std::flush;
					}

					if (instance->size() > 1)
					{
						int64_t currentBlock = ++blocksFound_;
						for (auto jt : *instance)
						{
							if (jt.Front().IsPositiveStrand())
							{
								blocksInstance_.push_back(BlockInstance(+currentBlock, jt.Front().GetChrId(), jt.Front().GetPosition(), jt.Back().GetPosition() + k_));
							}
							else
							{
								blocksInstance_.push_back(BlockInstance(-currentBlock, jt.Front().GetChrId(), jt.Back().GetPosition() - k_, jt.Front().GetPosition()));
							}

							for (auto it = jt.Front(); it != jt.Back(); ++it)
							{
								it.MarkUsed();
							}
						}
					}

					size_t nextTask = task + step;
					if (nextTask < bundle_.size())
					{
						taskMutex_.lock();
						taskQueue_.push(nextTask);
						taskMutex_.unlock();
					}

					++task;
				}
				else
				{
					taskMutex_.lock();
					taskQueue_.push(task);
					taskMutex_.unlock();
				}

				delete instance;
			}

			taskMutex_.lock();
			taskQueue_.push(GAME_OVER);
			taskMutex_.unlock();

			for (auto & thread : workerThread)
			{
				thread->join();
			}

			std::cout << ']' << std::endl;
			std::cout << "Time: " << clock() - mark << std::endl;
		}


		void ListBlocksSequences(const BlockList & block, const std::string & directory) const
		{
			std::vector<IndexPair> group;
			BlockList blockList = block;
			GroupBy(blockList, compareById, std::back_inserter(group));
			for (std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
			{
				std::ofstream out;
				std::stringstream ss;
				ss << directory << "/" << blockList[it->first].GetBlockId() << ".fa";
				TryOpenFile(ss.str(), out);
				for (size_t block = it->first; block < it->second; block++)
				{
					size_t length = blockList[block].GetLength();
					size_t chr = blockList[block].GetChrId();
					size_t chrSize = storage_.GetChrSequence(chr).size();
					out << ">" << blockList[block].GetBlockId() << "_" << block - it->first << " ";
					out << storage_.GetChrDescription(chr) << ";";
					if (blockList[block].GetSignedBlockId() > 0)
					{
						out << blockList[block].GetStart() << ";" << length << ";" << "+;" << chrSize << std::endl;
						OutputLines(storage_.GetChrSequence(chr).begin() + blockList[block].GetStart(), length, out);
					}
					else
					{
						size_t start = chrSize - blockList[block].GetEnd();
						out << start << ";" << length << ";" << "-;" << chrSize << std::endl;
						std::string::const_reverse_iterator it(storage_.GetChrSequence(chr).begin() + blockList[block].GetEnd());
						OutputLines(CFancyIterator(it, TwoPaCo::DnaChar::ReverseChar, ' '), length, out);
					}

					out << std::endl;
				}
			}
		}

		struct SortByMultiplicity
		{
			SortByMultiplicity(const std::vector<int> & multiplicityOrigin) : multiplicity(multiplicityOrigin)
			{
			}

			bool operator () (const BlockInstance & a, const BlockInstance & b) const
			{
				auto mlp1 = multiplicity[a.GetBlockId()];
				auto mlp2 = multiplicity[b.GetBlockId()];
				if (mlp1 != mlp2)
				{
					return mlp1 > mlp2;
				}

				return a.GetBlockId() < b.GetBlockId();
			}

			const std::vector<int> & multiplicity;
		};

		void GenerateOutput(const std::string & outDir, bool genSeq)
		{
			std::vector<std::vector<bool> > covered(storage_.GetChrNumber());
			for (size_t i = 0; i < covered.size(); i++)
			{
				covered[i].assign(storage_.GetChrSequence(i).size() + 1, false);
			}

			int64_t trimmedId = 1;
			std::vector<IndexPair> group;
			std::vector<BlockInstance> buffer;
			std::vector<BlockInstance> trimmedBlocks;
			std::vector<int> copiesCount_(blocksFound_ + 1, 0);
			for (auto b : blocksInstance_)
			{
				copiesCount_[b.GetBlockId()]++;
			}

			GroupBy(blocksInstance_, SortByMultiplicity(copiesCount_), std::back_inserter(group));
			for (auto g : group)
			{
				buffer.clear();
				for (size_t i = g.first; i < g.second; i++)
				{
					size_t chr = blocksInstance_[i].GetChrId();
					size_t start = blocksInstance_[i].GetStart();
					size_t end = blocksInstance_[i].GetEnd();
					for (; covered[chr][start] && start < end; start++);
					for (; covered[chr][end] && end > start; end--);
					if (end - start >= minBlockSize_)
					{
						buffer.push_back(BlockInstance(blocksInstance_[i].GetSign() * trimmedId, chr, start, end));
						std::fill(covered[chr].begin() + start, covered[chr].begin() + end, true);
					}
				}

				if (buffer.size() > 1)
				{
					trimmedId++;
					for (const auto & it : buffer)
					{
						trimmedBlocks.push_back(it);
					}
				}
				else
				{
					for (const auto & it : buffer)
					{
						std::fill(covered[it.GetChrId()].begin() + it.GetStart(), covered[it.GetChrId()].begin() + it.GetEnd(), false);
					}
				}
			}

			std::cout.setf(std::cout.fixed);
			std::cout.precision(2);
			std::cout << "Blocks found: " << blocksFound_ << std::endl;
			std::cout << "Coverage: " << CalculateCoverage(trimmedBlocks) << std::endl;

			CreateOutDirectory(outDir);
			std::string blocksDir = outDir + "/blocks";
			ListBlocksIndicesGFF(trimmedBlocks, outDir + "/" + "blocks_coords.gff");
			if (genSeq)
			{
				CreateOutDirectory(blocksDir);
				ListBlocksSequences(trimmedBlocks, blocksDir);
			}
		}


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

		double CalculateCoverage(const BlockList & block) const;
		void ListBlocksIndicesGFF(BlockList & blockList, const std::string & fileName);
		void TryOpenFile(const std::string & fileName, std::ofstream & stream) const;

		struct NextVertex
		{
			int64_t diff;
			int64_t count;
			JunctionStorage::JunctionSequentialIterator origin;
			NextVertex() : count(0)
			{

			}

			NextVertex(int64_t diff, JunctionStorage::JunctionSequentialIterator origin) : origin(origin), diff(diff), count(1)
			{

			}
		};

		std::pair<int64_t, NextVertex> MostPopularVertex(const Path & currentPath, bool forward, std::vector<uint32_t> & count, std::vector<size_t> & data)
		{
			NextVertex ret;
			int64_t bestVid = 0;
			int64_t startVid = forward ? currentPath.RightVertex() : currentPath.LeftVertex();
			const auto & instList = currentPath.GoodInstancesList().size() >= 2 ? currentPath.GoodInstancesList() : currentPath.AllInstances();
			for (auto & inst : instList)
			{
				int64_t nowVid = forward ? inst->Back().GetVertexId() : inst->Front().GetVertexId();
				if (nowVid == startVid)
				{
					int64_t weight = abs(inst->Front().GetPosition() - inst->Back().GetPosition()) + 1;
					auto origin = forward ? inst->Back() : inst->Front();
					auto it = forward ? origin.Next() : origin.Prev();
					for (size_t d = 1; it.Valid() && (d < size_t(lookingDepth_) || abs(it.GetPosition() - origin.GetPosition()) <= maxBranchSize_); d++)
					{
						int64_t vid = it.GetVertexId();
						if (!currentPath.IsInPath(vid) && !it.IsUsed())
						{
							auto adjVid = vid + storage_.GetVerticesNumber();
							if (count[adjVid] == 0)
							{
								data.push_back(adjVid);
							}

							count[adjVid] += static_cast<uint32_t>(weight);
							auto diff = abs(it.GetAbsolutePosition() - origin.GetAbsolutePosition());
							if (count[adjVid] > ret.count || (count[adjVid] == ret.count && diff < ret.diff))
							{
								ret.diff = diff;
								ret.origin = origin;
								ret.count = count[adjVid];
								bestVid = vid;
							}
						}
						else
						{
							break;
						}

						if (forward)
						{
							++it;
						}
						else
						{
							--it;
						}
					}
				}
			}


			for (auto vid : data)
			{
				count[vid] = 0;
			}

			data.clear();
			return std::make_pair(bestVid, ret);
		}

		bool ExtendPathForward(Path & currentPath,
			std::vector<uint32_t> & count,
			std::vector<size_t> & data,
			size_t & bestRightSize,
			int64_t & bestScore,
			int64_t & nowScore,
			std::vector<Path::Instance> & bestInstance)
		{
			bool success = false;
			int64_t origin = currentPath.Origin();
			std::pair<int64_t, NextVertex> nextForwardVid;
			nextForwardVid = MostPopularVertex(currentPath, true, count, data);
			if (nextForwardVid.first != 0)
			{
				for (auto it = nextForwardVid.second.origin; it.GetVertexId() != nextForwardVid.first; ++it)
				{
#ifdef _DEBUG_OUT_
					if (debug_)
					{
						std::cerr << "Attempting to push back the vertex:" << it.GetVertexId() << std::endl;
					}

					if (missingVertex_.count(it.GetVertexId()))
					{
						std::cerr << "Alert: " << it.GetVertexId() << ", origin: " << currentPath.Origin() << std::endl;
					}
#endif
					success = currentPath.PointPushBack(it.OutgoingEdge());
					if (success)
					{
						nowScore = currentPath.Score(scoreFullChains_);
#ifdef _DEBUG_OUT_
						if (debug_)
						{
							std::cerr << "Success! New score:" << nowScore << std::endl;
							currentPath.DumpPath(std::cerr);
							currentPath.DumpInstances(std::cerr);
						}
#endif												
						if (nowScore > bestScore)
						{
							bestScore = nowScore;
							bestRightSize = currentPath.RightSize();
							if (nowScore > 0)
							{
								bestInstance.clear();
								for (auto it : currentPath.GoodInstancesList())
								{
									bestInstance.push_back(*it);
								}
							}
						}
					}
				}
			}

			return success;
		}

		bool ExtendPathBackward(Path & currentPath,
			std::vector<uint32_t> & count,
			std::vector<size_t> & data,
			size_t & bestLeftSize,
			int64_t & bestScore,
			int64_t & nowScore,
			std::vector<Path::Instance> & bestInstance)
		{
			bool success = false;
			std::pair<int64_t, NextVertex> nextBackwardVid;
			nextBackwardVid = MostPopularVertex(currentPath, false, count, data);
			if (nextBackwardVid.first != 0)
			{
				for (auto it = nextBackwardVid.second.origin; it.GetVertexId() != nextBackwardVid.first; --it)
				{
#ifdef _DEBUG_OUT_
					if (debug_)
					{
						std::cerr << "Attempting to push front the vertex:" << it.GetVertexId() << std::endl;
					}

					if (missingVertex_.count(it.GetVertexId()))
					{
						std::cerr << "Alert: " << it.GetVertexId() << ", origin: " << currentPath.Origin() << std::endl;
					}
#endif
					success = currentPath.PointPushFront(it.IngoingEdge());
					if (success)
					{
						nowScore = currentPath.Score(scoreFullChains_);
#ifdef _DEBUG_OUT_
						if (debug_)
						{
							std::cerr << "Success! New score:" << nowScore << std::endl;
							currentPath.DumpPath(std::cerr);
							currentPath.DumpInstances(std::cerr);
						}
#endif		
						if (nowScore > bestScore)
						{
							bestScore = nowScore;
							bestLeftSize = currentPath.LeftSize();
							if (nowScore > 0)
							{
								bestInstance.clear();
								for (auto it : currentPath.GoodInstancesList())
								{
									bestInstance.push_back(*it);
								}
							}
						}
					}
				}
			}

			return success;
		}

		std::vector<Bundle> bundle_;
		tbb::spin_mutex taskMutex_;
		std::priority_queue<size_t, std::vector<size_t>, std::greater<size_t> > taskQueue_;
		tbb::spin_mutex resultMutex_;
		typedef std::pair<size_t, InstanceVectorPtr> Pair;
		std::priority_queue<Pair, std::vector<Pair>, std::greater<Pair> > resultQueue_;


		int64_t k_;
		size_t progressCount_;
		size_t progressPortion_;
		std::atomic<int64_t> count_;
		std::atomic<size_t> currentBundle_;
		std::atomic<int64_t> blocksFound_;
		int64_t sampleSize_;
		int64_t scalingFactor_;
		bool scoreFullChains_;
		int64_t lookingDepth_;
		int64_t minBlockSize_;
		int64_t maxBranchSize_;
		int64_t maxFlankingSize_;
		JunctionStorage & storage_;
		tbb::mutex progressMutex_;
		tbb::mutex blocksMutex_;
		std::ofstream debugOut_;
		std::vector<BlockInstance> blocksInstance_;
		std::vector<std::vector<Edge> > syntenyPath_;
#ifdef _DEBUG_OUT_
		bool debug_;
		std::set<int64_t> missingVertex_;
#endif
	};
}

#endif