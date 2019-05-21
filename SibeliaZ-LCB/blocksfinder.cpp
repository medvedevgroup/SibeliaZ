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
	extern const std::string VERSION = "1.0.0";

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

	int64_t BlockInstance::GetStart() const
	{
		return start_;
	}

	int64_t BlockInstance::GetEnd() const
	{
		return end_;
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

	int64_t BlockInstance::GetLength() const
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
