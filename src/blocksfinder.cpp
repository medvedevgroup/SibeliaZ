#include "blocksfinder.h"

namespace Sibelia
{
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
				return std::pair<size_t, size_t>(GetStart(), std::min(GetEnd(), instance.GetEnd()));
			}

			if (instance.GetStart() >= GetStart() && instance.GetStart() <= GetEnd())
			{
				return std::pair<size_t, size_t>(instance.GetStart(), std::min(GetEnd(), instance.GetEnd()));
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
}