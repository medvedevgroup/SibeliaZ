#ifndef _STREAM_FASTA_PARSER_H_
#define _STREAM_FASTA_PARSER_H_

#include <vector>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <tbb/mutex.h>
#include <iostream>
#include <memory>

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

	struct NewTask
	{
#ifdef _DEBUG
		static const size_t TASK_SIZE = 32;
#else
		static const size_t TASK_SIZE = 1 << 20;
#endif	
		bool isFinal;
		size_t seqId;
		size_t read;
		size_t start;		
		size_t piece;
		std::string str;
		std::string overlap;		
		char buffer[TASK_SIZE + 1];

		void Commence()
		{
			str = overlap;
			str.insert(str.end(), buffer, buffer + read);
		}
	};

	class GenomeReader
	{
	public:
		GenomeReader(size_t overlapSize) : pieceId_(0), start_(0), seqId_(0), overlapSize_(overlapSize), in_("parsed.txt")
		{

		}

		bool Read(NewTask & task)
		{			
			mutex_.lock();
			task.isFinal = true;
			bool ret = true;
			task.overlap = overlapBuffer_;
			overlapBuffer_.clear();			
						
			char ch;
			for (task.read = 0; task.read < NewTask::TASK_SIZE && in_.get(ch); task.read++)
			{
				if (isspace(ch))
				{
					break;
				}

				task.buffer[task.read] = ch;
			}

			if (task.read == 0)
			{				
				ret = false;
			}
			else
			{
				task.seqId = seqId_;
				task.start = start_;
				
				task.piece = pieceId_++;
				std::copy(task.buffer + task.read - overlapSize_, task.buffer + task.read, std::back_inserter(overlapBuffer_));
				if (task.read == NewTask::TASK_SIZE)
				{
					task.isFinal = false;
					start_ += task.read;
					overlapBuffer_.resize(overlapSize_);
					
				}
				else
				{
					task.isFinal = true;
					start_ = 0;
					++seqId_;
				}
			}

			mutex_.unlock();
			if (ret)
			{
				task.Commence();
			}
			
			return ret;
		}

	private:
		size_t seqId_;
		size_t start_;
		size_t pieceId_;
		tbb::mutex mutex_;
		size_t overlapSize_;
		std::ifstream in_;
		std::string overlapBuffer_;
	};

	class ChrReader
	{
	public:
		ChrReader(const std::vector<std::string> & fileName) : currentFile_(0), fileName_(fileName)
		{
			if (fileName.size() > 0)
			{
				parser_.reset(new TwoPaCo::StreamFastaParser(fileName[0]));
			}
		}

		bool NextChr(std::string & buf)
		{
			buf.clear();
			while (currentFile_ < fileName_.size())
			{
				if (parser_->ReadRecord())
				{
					char ch;
					while (parser_->GetChar(ch))
					{
						buf.push_back(ch);
					}

					return true;
				}
				else
				{
					if (++currentFile_ < fileName_.size())
					{
						parser_.reset(new TwoPaCo::StreamFastaParser(fileName_[currentFile_]));
					}
				}
			}

			return false;
		}

	private:
		size_t currentFile_;
		std::vector<std::string> fileName_;
		std::unique_ptr<TwoPaCo::StreamFastaParser> parser_;
	};

}

#endif
