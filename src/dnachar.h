#ifndef _DNA_CHAR_
#define _DNA_CHAR_

#include <string>
#include <cstdint>

namespace TwoPaCo
{
	class DnaChar
	{
	public:
		DnaChar();
		static const std::string LITERAL;
		static const std::string EXT_LITERAL;
		static const std::string VALID_CHARS;
		static bool IsValid(char ch);
		static bool IsDefinite(char ch);
		static char ReverseChar(char ch);		
		static std::string ReverseCompliment(const std::string & str);
		
		static size_t MakeUpChar(char ch);
		static char UnMakeUpChar(size_t ch);
		static bool LessSelfReverseComplement(std::string::const_iterator pit, size_t vertexSize);
	private:
		static const size_t CHAR_SIZE = 1 << 8;
		static bool isValid_[CHAR_SIZE];
		static bool isDefinite_[CHAR_SIZE];	
		static char reverseTable_[CHAR_SIZE];
	};
}

#endif
