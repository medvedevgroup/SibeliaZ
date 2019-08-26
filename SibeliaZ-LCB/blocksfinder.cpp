#include "blocksfinder.h"


namespace Sibelia
{
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef _WIN32
#include <direct.h>
#endif

	bool operator < (const Template & a, const Template & b)
	{
		int64_t dist1 = a.second.GetPosition() < a.first.GetPosition();
		int64_t dist2 = b.second.GetPosition() < b.first.GetPosition();
		return dist1 < dist2;
	}

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

	const std::string DELIMITER(80, '-');

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

	JunctionStorage * JunctionStorage::this_;
	extern const std::string VERSION = "1.0.0";

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
	template<class It>
	std::string Join(It start, It end, const std::string & delimiter)
	{
		It last = --end;
		std::stringstream ss;
		for (; start != end; ++start)
		{
			ss << *start << delimiter;
		}

		ss << *last;
		return ss.str();
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


	void BlocksFinder::ListBlocksIndices(const BlockList & block, const std::string & fileName) const
	{
		std::ofstream out;
		TryOpenFile(fileName, out);
		ListChrs(out);
		OutputBlocks(block, out);
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

	void BlocksFinder::ListBlocksIndicesGFF(const BlockList & blockList, const std::string & fileName) const
	{
		std::ofstream out;
		TryOpenFile(fileName, out);
		BlockList block(blockList);
		std::sort(block.begin(), block.end(), compareById);
		const std::string header[] =
		{
			"##gff-version 2",
			std::string("##source-version SibeliaZ ") + VERSION,
			"##Type DNA"
		};

		out << Join(header, header + 3, "\n") << std::endl;
		for (BlockList::const_iterator it = block.begin(); it != block.end(); ++it)
		{
			size_t start = it->GetStart() + 1;
			size_t end = it->GetEnd();
			const std::string record[] =
			{
				storage_.GetChrDescription(it->GetChrId()), 
				"SibeliaZ",
				"LCB_copy",
				IntToStr(start),
				IntToStr(end),
				".",
				(it->GetDirection() ? "+" : "-"),
				".",
				"id=" + IntToStr(static_cast<size_t>(it->GetBlockId()))
			};

			out << Join(record, record + sizeof(record) / sizeof(record[0]), "\t") << std::endl;
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
}
