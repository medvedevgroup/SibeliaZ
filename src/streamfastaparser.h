#ifndef _STREAM_FASTA_PARSER_H_
#define _STREAM_FASTA_PARSER_H_

#include <vector>
#include <fstream>
#include <stdexcept>
#include <algorithm>

#include "dnachar.h"

namespace TwoPaCo
{
	class StreamFastaParser
	{
	public:
		class Exception: public std::runtime_error
		{
		public:
			Exception(const std::string & msg);
		};
		
		bool ReadRecord();
		~StreamFastaParser();
		bool GetChar(char & ch);		
		std::string GetErrorMessage() const;
		std::string GetCurrentHeader() const;
		StreamFastaParser(const std::string & fileName);
	private:				
		static const size_t BUF_SIZE = 1 << 20;

		bool Peek(char & ch);
		bool GetCh(char & ch);		

		std::ifstream stream_;
		std::string errorMessage_;
		std::string currentHeader_;
		char * buffer_;
		size_t bufferSize_;
		size_t bufferPos_;
	};
}

#endif