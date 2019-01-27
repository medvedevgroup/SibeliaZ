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

#include "path.h"

namespace Sibelia
{
	extern const std::string DELIMITER;
	extern const std::string VERSION;

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
			progressCount_ = 50;
			scoreFullChains_ = true;			
		}

		struct ProcessVertex
		{
		public:
			BlocksFinder & finder;
			std::vector<int64_t> & shuffle;

			ProcessVertex(BlocksFinder & finder, std::vector<int64_t> & shuffle) : finder(finder), shuffle(shuffle)
			{
			}

			void operator()(tbb::blocked_range<size_t> & range) const
			{
				std::vector<uint32_t> data;
				std::vector<uint32_t> count(finder.storage_.GetVerticesNumber() * 2 + 1, 0);
				std::pair<int64_t, std::vector<Path::Instance> > goodInstance;
				Path finalizer(finder.storage_, finder.maxBranchSize_, finder.minBlockSize_, finder.minBlockSize_, finder.maxFlankingSize_);
				Path currentPath(finder.storage_, finder.maxBranchSize_, finder.minBlockSize_, finder.minBlockSize_, finder.maxFlankingSize_);
				for (size_t i = range.begin(); i != range.end(); i++)
				{
					if (finder.count_++ % finder.progressPortion_ == 0)
					{
						finder.progressMutex_.lock();
						std::cout << '.' << std::flush;
						finder.progressMutex_.unlock();
					}

					int64_t score;
					int64_t vid = shuffle[i];
#ifdef _DEBUG_OUT_
					finder.debug_ = finder.missingVertex_.count(vid);
					if (finder.debug_)
					{
						std::cerr << "Vid: " << vid << std::endl;
					}
#endif
					for (bool explore = true; explore;)
					{
						currentPath.Init(vid);
						if (currentPath.AllInstances().size() < 2)
						{
							currentPath.Clear();
							break;
						}

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
							while ((ret = finder.ExtendPathForward(currentPath, count, data, bestRightSize, bestScore, score)) && currentPath.MiddlePathLength() - prevLength <= minRun)
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
							currentPath.Init(vid);
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
							while ((ret = finder.ExtendPathBackward(currentPath, count, data, bestLeftSize, bestScore, score)) && currentPath.MiddlePathLength() - prevLength <= minRun);
							{
								positive = positive || (score > 0);
							}

							if (!ret || !positive)
							{
								break;
							}
						}

						if (bestScore > 0)
						{
#ifdef _DEBUG_OUT_
							if (finder.debug_)
							{
								std::cerr << "Setting a new block. Best score:" << bestScore << std::endl;
								currentPath.DumpPath(std::cerr);
								currentPath.DumpInstances(std::cerr);
							}
#endif
							if (!finder.TryFinalizeBlock(currentPath, finalizer, bestRightSize, bestLeftSize))
							{
								explore = false;
							}
						}
						else
						{
							explore = false;
						}

						currentPath.Clear();
					}
				}
			}
		};

		static bool DegreeCompare(const JunctionStorage & storage, int64_t v1, int64_t v2)
		{
			return storage.GetInstancesCount(v1) > storage.GetInstancesCount(v2);
		}

		void OutputMissing(const std::string & missingFile, const std::string & missingOutDir)
		{
			std::string buf;
			std::string id, sign, seq;
			std::set<int64_t> result;
			CreateOutDirectory(missingOutDir);
			std::ifstream oldCoordsIn(missingFile);
			while (std::getline(oldCoordsIn, buf))
			{
				if (buf == "END")
				{
					break;
				}

				if (buf.empty())
				{
					std::ofstream missingDot(missingOutDir + std::string("/missing") + id + ".dot");
					missingDot << "digraph G\n{\nrankdir = LR" << std::endl;
					std::vector<std::pair<JunctionStorage::JunctionSequentialIterator, JunctionStorage::JunctionSequentialIterator> > vvisit;
					for (auto vid : result)
					{
						DumpVertex(vid, missingDot, vvisit, 2);
						missingDot << vid << "[shape=square]" << std::endl;
					}

					result.clear();
					missingDot << "}" << std::endl;
				}

				std::stringstream ss(buf);
				int seqId, start, end, seqSize;
				ss >> id >> seq >> start >> end >> sign;
				for (auto it = storage_.Begin(storage_.GetSequenceId(seq)); it.Valid(); ++it)
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

		struct MafRecord
		{
			std::string seq;
			int64_t start;
			int64_t bodySize;
			int64_t seqSize;
			bool isPositive;
			std::string body;

			int64_t PositiveStart() const
			{
				if (isPositive)
				{
					return start;
				}

				return (seqSize - (start + bodySize)) + bodySize - 1;
			}
		};

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


		void OutputSubgraph(const std::string & mafFile, const std::string & outFile, const std::string & pairFiles)
		{
			std::string buf;
			std::ifstream mafIn(mafFile.c_str());
			std::multimap<std::string, std::string> neighbour;
			std::ofstream outStream(outFile.c_str());
			{
				std::vector<std::string> edgeBuf;
				std::ifstream pairFilesIn(pairFiles.c_str());
				while (std::getline(pairFilesIn, buf))
				{
					std::ifstream pairIn(buf.c_str());
					while (std::getline(pairIn, buf))
					{
						Split(buf, edgeBuf);
						if (edgeBuf.size() > 1)
						{
							neighbour.insert(std::make_pair(edgeBuf[0], edgeBuf[1]));
						}
					}
				}
			}

			std::cerr << "Total pairs: " << neighbour.size() << std::endl;

			if (!mafIn)
			{
				throw std::runtime_error("Can't open the MAF file");
			}

			std::string buffer;
			bool checkBlock = false;
			std::vector<std::string> data;
			std::vector<std::string> geneId;
			std::vector<std::string> newGeneId;
			std::vector<MafRecord> record;
			for (bool over = false; !over;)
			{
				if (!std::getline(mafIn, buffer))
				{
					over = checkBlock = true;
				}
				else
				{
					if (buffer.size() > 0 && buffer[0] != '#')
					{
						if (buffer[0] == 'a')
						{
							newGeneId.clear();
							size_t pos = buffer.find('=');
							if (pos != buffer.npos)
							{
								std::stringstream ss(buffer.substr(pos + 1));
								while (std::getline(ss, buf, ';'))
								{
									newGeneId.push_back(buf);
								}
							}
							
							checkBlock = true;
						}
						else
						{
							Split(buffer, data);
							if (data.size() != 7)
							{
								std::cerr << buffer << std::endl;
								for (auto & str : data)
								{
									std::cerr << str << std::endl;
								}

								throw std::runtime_error("A wrong line in the MAF file");
							}

							MafRecord nowRecord;
							nowRecord.seq = data[1];
							if (storage_.IsSequencePresent(nowRecord.seq))
							{
								nowRecord.start = std::atol(data[2].c_str());
								nowRecord.bodySize = std::atol(data[3].c_str());
								nowRecord.isPositive = data[4] == "+";
								nowRecord.seqSize = std::atol(data[5].c_str());
								nowRecord.body = data[6];
								record.push_back(nowRecord);
							}
						}
					}
				}

				if (checkBlock)
				{
					checkBlock = false;
					if (record.size() > 1 && geneId.size() == record.size())
					{
						int total = 0;
						int exceed = 0;
						for(size_t k = 0; k < record.size(); ++k)
						{
							const auto & rec = record[k];
							int64_t start = rec.PositiveStart();
							int64_t end = rec.PositiveStart();
							if (rec.isPositive)
							{
								end += rec.bodySize;
							}
							else
							{
								start -= rec.bodySize;
							}

							int nowNeighbours = neighbour.count(rec.seq);
							for (auto it = storage_.Begin(storage_.GetSequenceId(rec.seq)); it.Valid(); ++it)
							{
								int64_t pos = it.GetPosition();
								if (pos >= start && pos < end)
								{
									total++;
									if (nowNeighbours < storage_.GetInstancesCount(it.GetVertexId()))
									{
										exceed++;
									}
							//		outStream << it.GetPosition() << ' ' << it.GetVertexId() << ' ' << storage_.GetInstancesCount(it.GetVertexId()) << ';';
								}
							}

							outStream << geneId[0] << ';' << geneId[1] << '\t' <<  double(exceed) / total << std::endl;
						}
					}

					record.clear();
					newGeneId.swap(geneId);
				}
			}
		}

		void FindBlocks(int64_t minBlockSize, int64_t maxBranchSize, int64_t maxFlankingSize, int64_t lookingDepth, int64_t sampleSize, int64_t threads, const std::string & debugOut)
		{
			blocksFound_ = 0;
			sampleSize_ = sampleSize;
			lookingDepth_ = lookingDepth;
			minBlockSize_ = minBlockSize;
			maxBranchSize_ = maxBranchSize;
			maxFlankingSize_ = maxFlankingSize;
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

			using namespace std::placeholders;
			std::random_shuffle(shuffle.begin(), shuffle.end());
#ifdef _DEBUG_OUT_
			//OutputMissing("test/test7/segment.txt", "missing");
			OutputSubgraph("alignment.maf", "subgraph.txt", "pair.txt");
			exit(0);
#endif
			srand(time(0));
			time_t mark = time(0);
			count_ = 0;
			std::cout << '[' << std::flush;
			progressPortion_ = shuffle.size() / progressCount_;
			tbb::task_scheduler_init init(threads);
			tbb::parallel_for(tbb::blocked_range<size_t>(0, shuffle.size()), ProcessVertex(*this, shuffle));
			std::cout << ']' << std::endl;
			//std::cout << "Time: " << time(0) - mark << std::endl;
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


		void GenerateOutput(const std::string & outDir, bool genSeq) const
		{
			BlockList instance;
			std::vector<std::vector<bool> > covered(storage_.GetChrNumber());
			for (size_t i = 0; i < covered.size(); i++)
			{
				covered[i].assign(storage_.GetChrSequence(i).size(), false);
			}

			for (size_t chr = 0; chr < blockId_.size(); chr++)
			{
				for (size_t i = 0; i < blockId_[chr].size();)
				{
					if (storage_.GetIterator(chr, i).IsUsed())
					{
						int64_t bid = blockId_[chr][i].block;
						size_t j = i;
						for (; j < blockId_[chr].size() && blockId_[chr][i] == blockId_[chr][j]; j++);
						j--;
						int64_t cstart = storage_.GetIterator(chr, i, bid > 0).GetPosition();
						int64_t cend = storage_.GetIterator(chr, j, bid > 0).GetPosition() + (bid > 0 ? k_ : -k_);
						int64_t start = min(cstart, cend);
						int64_t end = max(cstart, cend);
						instance.push_back(BlockInstance(bid, chr, start, end));
						i = j + 1;
					}
					else
					{
						++i;
					}
				}
			}

			std::cout.setf(std::cout.fixed);
			std::cout.precision(2);
			std::cout << "Blocks found: " << blocksFound_ << std::endl;
			std::cout << "Total coverage: " << CalculateCoverage(instance) << std::endl;


			CreateOutDirectory(outDir);
			std::string blocksDir = outDir + "/blocks";
			ListBlocksIndicesGFF(instance, outDir + "/" + "blocks_coords.gff");
			if (genSeq)
			{
				CreateOutDirectory(blocksDir);
				ListBlocksSequences(instance, blocksDir);
			}

	//		GenerateReport(instance, outDir + "/" + "coverage_report.txt");*/
	//		ListBlocksIndices(instance, outDir + "/" + "blocks_coords.txt");
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
		void GenerateReport(const BlockList & block, const std::string & fileName) const;
		void CalculateCoverage(GroupedBlockList::const_iterator start, GroupedBlockList::const_iterator end, std::vector<double> & ret) const;
		std::string OutputIndex(const BlockInstance & block) const;
		void OutputBlocks(const std::vector<BlockInstance>& block, std::ofstream& out) const;
		void ListBlocksIndices(const BlockList & block, const std::string & fileName) const;
		void ListBlocksIndicesGFF(const BlockList & blockList, const std::string & fileName) const;
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

		bool TryFinalizeBlock(const Path & currentPath, Path & finalizer, size_t bestRightSize, size_t bestLeftSize)
		{
			bool ret = false;
			std::vector<Path::InstanceSet::const_iterator> lockInstance;
			for (auto it : currentPath.GoodInstancesList())
			{
				lockInstance.push_back(it);
			}
			
			std::sort(lockInstance.begin(), lockInstance.end(), Path::CmpInstance);
			{
				std::pair<size_t, size_t> idx(SIZE_MAX, SIZE_MAX);
				for (auto & instance : lockInstance)
				{
					if (instance->Front().IsPositiveStrand())
					{
						storage_.LockRange(instance->Front(), instance->Back(), idx);
					}
					else
					{
						storage_.LockRange(instance->Back().Reverse(), instance->Front().Reverse(), idx);
					}
				}
			}
		
			finalizer.Init(currentPath.Origin());
			for (size_t i = 0; i < bestRightSize - 1 && finalizer.PointPushBack(currentPath.RightPoint(i).GetEdge()); i++);
			for (size_t i = 0; i < bestLeftSize - 1 && finalizer.PointPushFront(currentPath.LeftPoint(i).GetEdge()); i++);
			if (finalizer.Score() > 0 && finalizer.GoodInstances() > 1)
			{
				ret = true;
				int64_t instanceCount = 0;
				int64_t currentBlock = ++blocksFound_;	
				for (auto jt : finalizer.AllInstances())
				{
					if (finalizer.IsGoodInstance(*jt))
					{
						auto it = jt->Front();
						do
						{
							it.MarkUsed();							
							int64_t idx = it.GetIndex();
							int64_t maxidx = storage_.GetChrVerticesCount(it.GetChrId());
							blockId_[it.GetChrId()][it.GetIndex()].block = it.IsPositiveStrand() ? +currentBlock : -currentBlock;
							blockId_[it.GetChrId()][it.GetIndex()].instance = instanceCount;

						} while (it++ != jt->Back());

						instanceCount++;
					}
				}
			}
				
			finalizer.Clear();
			std::pair<size_t, size_t> idx(SIZE_MAX, SIZE_MAX);
			for (auto & instance : lockInstance)
			{
				if (instance->Front().IsPositiveStrand())
				{
					storage_.UnlockRange(instance->Front(), instance->Back(), idx);
				}
				else
				{
					storage_.UnlockRange(instance->Back().Reverse(), instance->Front().Reverse(), idx);
				}
			}

			return ret;
		}

		struct NextVertex
		{
			int32_t diff;
			int32_t count;
			JunctionStorage::JunctionSequentialIterator origin;
			NextVertex() : count(0)
			{

			}

			NextVertex(int64_t diff, JunctionStorage::JunctionSequentialIterator origin) : origin(origin), diff(diff), count(1)
			{

			}
		};

		std::pair<int32_t, NextVertex> MostPopularVertex(const Path & currentPath, bool forward, std::vector<uint32_t> & count, std::vector<uint32_t> & data)
		{
			NextVertex ret;
			int32_t bestVid = 0;
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
					for (size_t d = 1; it.Valid() && (d < lookingDepth_  || abs(it.GetPosition() - origin.GetPosition()) <= maxBranchSize_); d++)
					{
						int32_t vid = it.GetVertexId();
						if (!currentPath.IsInPath(vid) && !it.IsUsed())
						{
							auto adjVid = vid + storage_.GetVerticesNumber();
							if (count[adjVid] == 0)
							{
								data.push_back(adjVid);
							}

							count[adjVid] += weight;
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
			std::vector<uint32_t> & data,
			size_t & bestRightSize,
			int64_t & bestScore,
			int64_t & nowScore)
		{
			bool success = false;
			int64_t origin = currentPath.Origin();
			std::pair<int32_t, NextVertex> nextForwardVid;
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
						}
					}
				}
			}

			return success;
		}

		bool ExtendPathBackward(Path & currentPath,
			std::vector<uint32_t> & count,
			std::vector<uint32_t> & data,
			size_t & bestLeftSize,
			int64_t & bestScore,
			int64_t & nowScore)
		{
			bool success = false;
			std::pair<int32_t, NextVertex> nextBackwardVid;
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
						}
					}
				}
			}

			return success;
		}

		int64_t k_;
		size_t progressCount_;
		size_t progressPortion_;
		std::atomic<int64_t> count_;
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
		std::ofstream debugOut_;
		std::vector<std::vector<Edge> > syntenyPath_;
		std::vector<std::vector<Assignment> > blockId_;
#ifdef _DEBUG_OUT_
		bool debug_;
		std::set<int64_t> missingVertex_;
#endif
	};
}

#endif
