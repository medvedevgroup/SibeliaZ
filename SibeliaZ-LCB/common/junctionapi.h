#ifndef _JUNCTION_POSITION_API_H_
#define _JUNCTION_POSITION_API_H_

#include <fstream>
#include <cstdint>
#include <exception>

namespace TwoPaCo
{
	struct JunctionPosition
	{
	public:
		JunctionPosition() : chr_(UINT32_MAX), pos_(UINT32_MAX), bifId_(INT64_MAX) {}
		JunctionPosition(uint32_t chr, uint32_t pos, int64_t bifId) :
			chr_(chr), pos_(pos), bifId_(bifId) {}
		uint32_t GetPos() const
		{
			return pos_;
		}

		uint32_t GetChr() const
		{
			return chr_;
		}

		int64_t GetId() const
		{
			return bifId_;
		}

	private:
		uint32_t chr_;
		uint32_t pos_;
		int64_t bifId_;
		static const int64_t SEPARATOR_BIF = INT64_MAX;
		static const uint32_t SEPARATOR_POS = -1;		
		friend class JunctionPositionReader;
		friend class JunctionPositionWriter;
	};

	class JunctionPositionReader
	{
	public:
		JunctionPositionReader(const std::string & inFileName) : nowChr_(0), in_(inFileName.c_str(), std::ios::binary)
		{
			if (!in_)
			{
				throw std::runtime_error("Can't read the input file");
			}
		}

		void RestoreVector(std::vector<bool> & mark, size_t chr)
		{
			JunctionPosition pos;
			mark.assign(mark.size(), false);
			while (NextJunctionPosition(pos))
			{
				if (pos.GetChr() == chr)
				{
					mark[pos.GetPos()] = true;
				}
				else
				{
					size_t cur = in_.tellg();
					in_.seekg(cur - sizeof(pos.pos_) - sizeof(pos.bifId_), in_.beg);
					break;					
				}
			}
		}

		void RestoreAllVectors(std::vector<std::vector<bool> > & mark)
		{
			JunctionPosition pos;
			while (NextJunctionPosition(pos))
			{
				mark[pos.GetChr()][pos.GetPos()] = true;
			}
		}

		bool NextJunctionPosition(JunctionPosition & pos)
		{
			for (;; nowChr_++)
			{
				pos = JunctionPosition(nowChr_, 0, 0);
				in_.read(reinterpret_cast<char*>(&pos.pos_), sizeof(pos.pos_));
				in_.read(reinterpret_cast<char*>(&pos.bifId_), sizeof(pos.bifId_));

				if (!in_)
				{
					return false;
				}

				if (pos.pos_ != JunctionPosition::SEPARATOR_POS && pos.bifId_ != JunctionPosition::SEPARATOR_BIF)
				{
					return true;
				}
			}
		}

	private:
		uint32_t nowChr_;
		std::ifstream in_;
	};
	

	class JunctionPositionWriter
	{
	public:
		JunctionPositionWriter(const std::string & outFileName) : nowChr_(0), out_(outFileName.c_str(), std::ios::binary)
		{			
			if (!out_)
			{
				throw std::runtime_error("Can't create the output file");
			}
		}

		void WriteJunction(JunctionPosition pos)
		{
			for (; pos.chr_ > nowChr_; ++nowChr_)
			{
				WriteJunction(JunctionPosition(nowChr_, JunctionPosition::SEPARATOR_POS, JunctionPosition::SEPARATOR_BIF));
			}

			out_.write(reinterpret_cast<const char*>(&pos.pos_), sizeof(pos.pos_));
			out_.write(reinterpret_cast<const char*>(&pos.bifId_), sizeof(pos.bifId_));

			if (!out_)
			{
				throw std::runtime_error("Can't write to the output file");
			}
		}

	private:
		uint32_t nowChr_;
		std::ofstream out_;
	};
}

#endif