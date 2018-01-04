#include "blocksfinder.h"


namespace Sibelia
{
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef _WIN32
#include <direct.h>
#endif

	void CreateOutDirectory(const std::string & path)
	{
		int result = 0;
#ifdef _WIN32
		result = _mkdir(path.c_str());
#else
		result = mkdir(path.c_str(), 0755);
#endif
		if (result != 0 && errno != EEXIST)
		{
			throw std::runtime_error(("Cannot create dir " + path).c_str());
		}
	}

	JunctionStorage * JunctionStorage::this_;

	bool compareById(const BlockInstance & a, const BlockInstance & b)
	{
		return CompareBlocks(a, b, &BlockInstance::GetBlockId);
	}

	bool compareByChrId(const BlockInstance & a, const BlockInstance & b)
	{
		return CompareBlocks(a, b, &BlockInstance::GetChrId);
	}

	bool compareByStart(const BlockInstance & a, const BlockInstance & b)
	{
		return CompareBlocks(a, b, &BlockInstance::GetChrId);
	}

	const std::string DELIMITER(80, '-');

	int BlockInstance::GetSignedBlockId() const
	{
		return id_;
	}

	bool BlockInstance::GetDirection() const
	{
		return id_ > 0;
	}

	int BlockInstance::GetSign() const
	{
		return GetSignedBlockId() > 0 ? +1 : -1;
	}

	int BlockInstance::GetBlockId() const
	{
		return abs(id_);
	}

	size_t BlockInstance::GetChrId() const
	{
		return chr_;
	}

	size_t BlockInstance::GetStart() const
	{
		return start_;
	}

	size_t BlockInstance::GetEnd() const
	{
		return end_;
	}

	size_t BlockInstance::GetConventionalStart() const
	{
		if (GetDirection())
		{
			return start_ + 1;
		}

		return end_;
	}

	size_t BlockInstance::GetConventionalEnd() const
	{
		if (GetDirection())
		{
			return end_;
		}

		return start_ + 1;
	}

	std::pair<size_t, size_t> BlockInstance::CalculateOverlap(const BlockInstance & instance) const
	{
		if (GetChrId() == instance.GetChrId())
		{
			size_t overlap = 0;
			if (GetStart() >= instance.GetStart() && GetStart() <= instance.GetEnd())
			{
				return std::pair<size_t, size_t>(GetStart(), min(GetEnd(), instance.GetEnd()));
			}

			if (instance.GetStart() >= GetStart() && instance.GetStart() <= GetEnd())
			{
				return std::pair<size_t, size_t>(instance.GetStart(), min(GetEnd(), instance.GetEnd()));
			}
		}

		return std::pair<size_t, size_t>(0, 0);
	}

	bool BlockInstance::operator == (const BlockInstance & toCompare) const
	{
		return start_ == toCompare.start_ && end_ == toCompare.end_ && GetChrId() == toCompare.GetChrId() && id_ == toCompare.id_;
	}

	bool BlockInstance::operator != (const BlockInstance & toCompare) const
	{
		return !(*this == toCompare);
	}

	void BlockInstance::Reverse()
	{
		id_ = -id_;
	}

	size_t BlockInstance::GetLength() const
	{
		return end_ - start_;
	}

	bool BlockInstance::operator < (const BlockInstance & toCompare) const
	{
		return std::make_pair(GetBlockId(), std::make_pair(GetChrId(), GetStart())) < std::make_pair(toCompare.GetBlockId(), std::make_pair(toCompare.GetChrId(), toCompare.GetStart()));
	}

	void BlocksFinder::GenerateReport(const BlockList & block, const std::string & fileName) const
	{
		std::ofstream out;
		TryOpenFile(fileName, out);
		GroupedBlockList sepBlock;
		std::vector<IndexPair> group;
		BlockList blockList = block;
		GroupBy(blockList, compareById, std::back_inserter(group));
		for (std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
		{
			sepBlock.push_back(std::make_pair(it->second - it->first, std::vector<BlockInstance>(blockList.begin() + it->first, blockList.begin() + it->second)));
		}

		ListChrs(out);
		out << "Degree\tCount\tTotal";
		for (size_t i = 0; i < storage_.GetChrNumber(); i++)
		{
			out << "\tSeq " << i + 1;
		}

		out << std::endl;
		group.clear();
		std::vector<bool> cover;
		GroupBy(sepBlock, ByFirstElement, std::back_inserter(group));
		group.push_back(IndexPair(0, sepBlock.size()));
		for (std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
		{
			if (it != group.end() - 1)
			{
				out << sepBlock[it->first].first << '\t' << it->second - it->first << '\t';
			}
			else
			{
				out << "All\t" << it->second - it->first << "\t";
			}

			out.precision(2);
			out.setf(std::ostream::fixed);
			std::vector<double> coverage = CalculateCoverage(sepBlock.begin() + it->first, sepBlock.begin() + it->second, cover);
			std::copy(coverage.begin(), coverage.end(), std::ostream_iterator<double>(out, "%\t"));
			out << std::endl;
		}

		out << DELIMITER << std::endl;
	}

	std::vector<double> BlocksFinder::CalculateCoverage(GroupedBlockList::const_iterator start, GroupedBlockList::const_iterator end, std::vector<bool> & cover) const
	{
		std::vector<double> ret;
		double totalBp = 0;
		double totalCoveredBp = 0;
		for (size_t chr = 0; chr < storage_.GetChrNumber(); chr++)
		{
			totalBp += storage_.GetChrSequence(chr).size();
			cover.assign(storage_.GetChrSequence(chr).size(), 0);
			for (GroupedBlockList::const_iterator it = start; it != end; ++it)
			{
				for (size_t i = 0; i < it->second.size(); i++)
				{
					if (it->second[i].GetChrId() == chr)
					{
						std::fill(cover.begin() + it->second[i].GetStart(), cover.begin() + it->second[i].GetEnd(), COVERED);
					}
				}
			}

			double nowCoveredBp = static_cast<double>(std::count(cover.begin(), cover.end(), COVERED));
			ret.push_back(nowCoveredBp / cover.size() * 100);
			totalCoveredBp += nowCoveredBp;
		}

		ret.insert(ret.begin(), totalCoveredBp / totalBp * 100);
		return ret;
	}


	std::string BlocksFinder::OutputIndex(const BlockInstance & block) const
	{
		std::stringstream out;
		out << block.GetChrId() + 1 << '\t' << (block.GetSignedBlockId() < 0 ? '-' : '+') << '\t';
		out << block.GetConventionalStart() << '\t' << block.GetConventionalEnd() << '\t' << block.GetEnd() - block.GetStart();
		return out.str();
	}


	void BlocksFinder::OutputBlocks(const std::vector<BlockInstance>& block, std::ofstream& out) const
	{
		std::vector<IndexPair> group;
		std::vector<BlockInstance> blockList = block;
		GroupBy(blockList, compareById, std::back_inserter(group));
		for (std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
		{
			size_t length = it->second - it->first;
			std::sort(blockList.begin() + it->first, blockList.begin() + it->second, compareByChrId);
			out << "Block #" << blockList[it->first].GetBlockId() << std::endl;
			out << "Seq_id\tStrand\tStart\tEnd\tLength" << std::endl;
			for (auto jt = blockList.begin() + it->first; jt < blockList.begin() + it->first + length; ++jt)
			{
				out << OutputIndex(*jt) << std::endl;
			}

			out << DELIMITER << std::endl;
		}
	}


	void BlocksFinder::ListBlocksIndices(const BlockList & block, const std::string & fileName) const
	{
		std::ofstream out;
		TryOpenFile(fileName, out);
		ListChrs(out);
		OutputBlocks(block, out);
	}

	void BlocksFinder::ListChromosomesAsPermutations(const BlockList & block, const std::string & fileName) const
	{
		std::ofstream out;
		TryOpenFile(fileName, out);
		std::vector<IndexPair> group;
		BlockList blockList = block;
		GroupBy(blockList, compareByChrId, std::back_inserter(group));
		for (std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
		{
			out.setf(std::ios_base::showpos);
			size_t length = it->second - it->first;
			size_t chr = blockList[it->first].GetChrId();
			out << '>' << storage_.GetChrDescription(chr) << std::endl;
			std::sort(blockList.begin() + it->first, blockList.begin() + it->second);
			for (auto jt = blockList.begin() + it->first; jt < blockList.begin() + it->first + length; ++jt)
			{
				out << jt->GetSignedBlockId() << " ";
			}

			out << "$" << std::endl;
		}
	}


	void BlocksFinder::TryOpenFile(const std::string & fileName, std::ofstream & stream) const
	{
		stream.open(fileName.c_str());
		if (!stream)
		{
			throw std::runtime_error(("Cannot open file " + fileName).c_str());
		}
	}

	void BlocksFinder::ListChrs(std::ostream & out) const
	{
		out << "Seq_id\tSize\tDescription" << std::endl;
		for (size_t i = 0; i < storage_.GetChrNumber(); i++)
		{
			out << i + 1 << '\t' << storage_.GetChrSequence(i).size() << '\t' << storage_.GetChrDescription(i) << std::endl;
		}

		out << DELIMITER << std::endl;
	}
}