#include <cassert>
#include "dnachar.h"

namespace TwoPaCo
{
	bool DnaChar::isValid_[CHAR_SIZE];
	bool DnaChar::isDefinite_[CHAR_SIZE];
	char DnaChar::reverseTable_[CHAR_SIZE];
	const std::string DnaChar::LITERAL = "ACGT";
	const std::string DnaChar::EXT_LITERAL = "ACGTN";
	const std::string DnaChar::VALID_CHARS = "ACGTURYKMSWBDHWNXV";
	
	namespace
	{
		DnaChar helper;
	}
	
	size_t DnaChar::MakeUpChar(char ch)
	{
		switch (ch)
		{
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		}

		return -1;
	}

	char DnaChar::UnMakeUpChar(size_t ch)
	{
		switch (ch)
		{
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
		}

		return 'N';
	}

	DnaChar::DnaChar()
	{
		std::fill(reverseTable_, reverseTable_ + CHAR_SIZE, 'N');
		reverseTable_['A'] = 'T';
		reverseTable_['T'] = 'A';
		reverseTable_['C'] = 'G';
		reverseTable_['G'] = 'C';
		std::fill(isValid_, isValid_ + CHAR_SIZE, false);
		std::fill(isDefinite_, isDefinite_ + CHAR_SIZE, false);
		for (char ch : LITERAL)
		{
			isDefinite_[ch] = true;
		}

		for (char ch : VALID_CHARS)
		{
			isValid_[ch] = true;
		}
	}

	bool DnaChar::IsValid(char ch)
	{
		return isValid_[ch];
	}

	bool DnaChar::IsDefinite(char ch)
	{
		return isDefinite_[ch];
	}

	char DnaChar::ReverseChar(char ch)
	{
		return reverseTable_[ch];
	}

	std::string DnaChar::ReverseCompliment(const std::string & str)
	{
		std::string ret;
		for (auto it = str.rbegin(); it != str.rend(); ++it)
		{
			ret.push_back(DnaChar::ReverseChar(*it));
		}

		return ret;
	}

	bool DnaChar::LessSelfReverseComplement(std::string::const_iterator pit, size_t size)
	{
		std::string::const_reverse_iterator nit(pit + size);
		for (size_t i = 0; i < size; i++)
		{
			char reverse = DnaChar::ReverseChar(*nit);
			if (*pit != reverse)
			{
				return *pit < reverse;
			}

			++nit;
			++pit;
		}

		return false;
	}
}