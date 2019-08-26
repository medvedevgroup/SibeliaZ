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

#include <tbb/parallel_for.h>

#include "sweeper.h"

namespace Sibelia
{
	extern const std::string DELIMITER;
	extern const std::string VERSION;

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
			progressCount_ = 50;
			scoreFullChains_ = true;
		}

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
			blocksFound_ = 0;
			minBlockSize_ = minBlockSize;
			maxBranchSize_ = maxBranchSize;
			maxFlankingSize_ = maxFlankingSize;

			using namespace std::placeholders;
			tbb::task_scheduler_init init(static_cast<int>(threads));

			time_t mark = time(0);
			count_ = 0;
			
			time_t start = clock();
			starter_ = 0;
			progressCount_ = 0;
			tbb::parallel_for(tbb::blocked_range<size_t>(0, storage_.GetChrNumber()), ChrSweep(*this));
			std::sort(template_.begin(), template_.end());
			RunTemplate runner(*this);
			runner(tbb::blocked_range<size_t>(0, template_.size()));
			std::cout << double(clock() - start) / CLOCKS_PER_SEC << std::endl;
		}

		struct ChrSweep
		{
		public:
			BlocksFinder & finder;

			ChrSweep(BlocksFinder & finder) : finder(finder)
			{

			}

			void operator()(tbb::blocked_range<size_t> & range) const
			{
				for (size_t r = range.begin(); r != range.end(); r++)
				{
					Sweeper sweeper(finder.storage_.Begin(r));
					sweeper.Sweep(finder.minBlockSize_, finder.maxBranchSize_, finder.k_, finder.blocksFound_, finder.template_, finder.globalMutex_);
					{
						tbb::mutex::scoped_lock lock(finder.globalMutex_);
						std::cout << ++finder.progressCount_ << ' ' << finder.storage_.GetChrNumber() << std::endl;
					}
				}
			}
		};

		struct RunTemplate
		{
		public:
			BlocksFinder & finder;

			RunTemplate(BlocksFinder & finder) : finder(finder)
			{

			}

			void Finalize(Path & currentPath) const
			{
				int64_t currentBlock = ++finder.blocksFound_;
				//finder.log_ << '*' << currentBlock << std::endl;
				if (currentPath.GoodInstancesList().size() > 1)
				{
					std::cout << "g" << std::endl;
					for (auto it : currentPath.GoodInstancesList())
					{
						auto & jt = *it;
						if (jt.Front().IsPositiveStrand())
						{
							finder.blocksInstance_.push_back(BlockInstance(+currentBlock, jt.Front().GetChrId(), jt.Front().GetPosition(), jt.Back().GetPosition() + finder.k_));
						}
						else
						{
							finder.blocksInstance_.push_back(BlockInstance(-currentBlock, jt.Front().GetChrId(), jt.Back().GetPosition() - finder.k_, jt.Front().GetPosition()));
						}

						//finder.log_ << jt.Front().GetChrId() << ' ' << jt.Front().GetPosition() << ' ' << jt.Back().GetPosition() << std::endl;
						for (auto it = jt.Front(); it != jt.Back(); ++it)
						{
							it.MarkUsed();
						}
					}
				}
				
			}

			void operator()(tbb::blocked_range<size_t> & range) const
			{
				std::vector<size_t> data;
				std::vector<uint32_t> count(finder.storage_.GetVerticesNumber() * 2 + 1, 0);
				Path currentPath(finder.storage_, finder.maxBranchSize_, finder.minBlockSize_, finder.minBlockSize_, finder.maxFlankingSize_, true);
				for (size_t i = range.begin(); i < range.end(); i++)
				{
					std::cout << i << ' ' << range.end() << std::endl;
					auto run = finder.template_[i];
					while (run.first < run.second)
					{
						for (; run.first.IsUsed() && run.first < run.second; ++run.first);
						auto end = run.first;
						for (; !end.IsUsed() && end < run.second; ++end);
						if (end.GetPosition() - run.first.GetPosition() >= finder.minBlockSize_)
						{
							currentPath.Init(run.first.GetVertexId(), 'N');
							for (; run.first != end; ++run.first)
							{
								currentPath.PointPushBack(run.first.OutgoingEdge());
							}
							
							Finalize(currentPath);
							currentPath.Clear();
						}

						run.first = end;
					}
				}
			}
		};
		
		

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
					int64_t length = blockList[block].GetLength();
					size_t chr = blockList[block].GetChrId();
					int64_t chrSize = storage_.GetChrSequence(chr).size();
					out << ">" << blockList[block].GetBlockId() << "_" << block - it->first << " ";
					out << storage_.GetChrDescription(chr) << ";";
					if (blockList[block].GetSignedBlockId() > 0)
					{
						out << blockList[block].GetStart() << ";" << length << ";" << "+;" << chrSize << std::endl;
						OutputLines(storage_.GetChrSequence(chr).begin() + blockList[block].GetStart(), length, out);
					}
					else
					{
						int64_t start = chrSize - blockList[block].GetEnd();
						out << start << ";" << length << ";" << "-;" << chrSize << std::endl;
						std::string::const_reverse_iterator it(storage_.GetChrSequence(chr).begin() + blockList[block].GetEnd());
						OutputLines(CFancyIterator(it, TwoPaCo::DnaChar::ReverseChar, ' '), length, out);
					}

					out << std::endl;
				}

				//std::cout << std::endl;
			}
		}

		void GenerateOutput(const std::string & outDir, bool genSeq)
		{
			const auto & trimmedBlocks = blocksInstance_;
			std::vector<std::vector<bool> > covered(storage_.GetChrNumber());
			for (size_t i = 0; i < covered.size(); i++)
			{
				covered[i].assign(storage_.GetChrSequence(i).size() + 1, false);
			}
			
			for (auto & b : blocksInstance_)
			{
				for (size_t i = b.GetStart(); i < b.GetEnd(); i++)
				{
					covered[b.GetChrId()][i] = true;
				}
			}
			
			size_t total = 0;
			size_t totalCovered = 0;
			for (auto & chr : covered)
			{
				total += chr.size();
				totalCovered += std::count(chr.begin(), chr.end(), true);
			}

			std::cout.setf(std::cout.fixed);
			std::cout.precision(2);
			std::cout << "Blocks found: " << blocksFound_ << std::endl;
			std::cout << "Coverage: " << double(totalCovered) / total << std::endl;

			CreateOutDirectory(outDir);
			std::string blocksDir = outDir + "/blocks";
			ListBlocksIndicesGFF(trimmedBlocks, outDir + "/" + "blocks_coords.gff");
			if (genSeq)
			{
				CreateOutDirectory(blocksDir);
				ListBlocksSequences(trimmedBlocks, blocksDir);
			}

			ListBlocksIndices(trimmedBlocks, outDir + "/" + "blocks_coords.txt");
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

		void ListChrs(std::ostream & out) const;
		std::string OutputIndex(const BlockInstance & block) const;
		void OutputBlocks(const std::vector<BlockInstance>& block, std::ofstream& out) const;
		void ListBlocksIndices(const BlockList & block, const std::string & fileName) const;
		void ListBlocksIndicesGFF(const BlockList & blockList, const std::string & fileName) const;
		void TryOpenFile(const std::string & fileName, std::ofstream & stream) const;


		int64_t k_;
		size_t progressCount_;
		size_t progressPortion_;
		std::atomic<int64_t> count_;
		std::atomic<int64_t> starter_;
		std::atomic<int64_t> blocksFound_;
		std::vector<size_t> pointComponent_;

		bool scoreFullChains_;
		int64_t scalingFactor_;
		int64_t minBlockSize_;
		int64_t maxBranchSize_;
		int64_t maxFlankingSize_;
		JunctionStorage & storage_;
		tbb::mutex progressMutex_;
		tbb::mutex globalMutex_;
		std::ofstream debugOut_;
		std::vector<Template> template_;
		std::vector<BlockInstance> blocksInstance_;

		//std::ofstream forkLog;


#ifdef _DEBUG_OUT_
		bool debug_;
		std::set<int64_t> missingVertex_;
#endif
	};
}

#endif
