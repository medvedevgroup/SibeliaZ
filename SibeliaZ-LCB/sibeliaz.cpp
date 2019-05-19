#include <tclap/CmdLine.h>

#include "blocksfinder.h"

size_t Atoi(const char * str)
{
	size_t ret;
	std::stringstream ss(str);
	ss >> ret;
	return ret;
}

class OddConstraint : public TCLAP::Constraint < unsigned int >
{
public:
	~OddConstraint()
	{

	}

	std::string description() const
	{
		return "value of K must be odd";
	}

	std::string shortID() const
	{
		return "oddc";
	}

	bool check(const unsigned & value) const
	{
		return (value % 2) == 1;
	}
};

int main(int argc, char * argv[])
{
	OddConstraint constraint;

	try
	{
		TCLAP::CmdLine cmd("SibeliaZ-LCB, a program for construction of locally-collinear blocks from complete genomes", ' ', Sibelia::VERSION);

		TCLAP::ValueArg<unsigned int> kvalue("k",
			"kvalue",
			"Value of k",
			false,
			25,
			&constraint,
			cmd);

		TCLAP::ValueArg<unsigned int> maxBranchSize("b",
			"branchsize",
			"Maximum branch size",
			false,
			200,
			"integer",
			cmd);

		TCLAP::ValueArg<unsigned int> minBlockSize("m",
			"blocksize",
			"Minimum block size",
			false,
			50,
			"integer",
			cmd);

		TCLAP::ValueArg<unsigned int> threads("t",
			"threads",
			"Number of worker threads",
			false,
			1,
			"integer",
			cmd);

		TCLAP::ValueArg<unsigned int> abundanceThreshold("a",
			"abundance",
			"Max abundance of a junction",
			false,
			150,
			"integer",
			cmd);

		TCLAP::ValueArg<std::string> inFileName("",
			"graph",
			"Binary file containing the graph",
			true,
			"de_bruijn.bin",
			"file name",
			cmd);

		TCLAP::ValueArg<std::string> outDirName("o",
			"outdir",
			"Output dir for blocks sequences",
			false,
			"",
			"directory name",
			cmd);

		TCLAP::SwitchArg noSeq("",
			"noseq",
			"Do not output blocks sequences",
			cmd,
			false);

		TCLAP::UnlabeledMultiArg<std::string> genomesFileName("filenames",
			"FASTA file(s) with nucleotide sequences.",
			true,
			"fasta files with genomes",
			cmd);

		cmd.parse(argc, argv);

		std::cout << "Loading the graph..." << std::endl;
		Sibelia::JunctionStorage storage(inFileName.getValue(),
			genomesFileName.getValue(),
			kvalue.getValue(),
			threads.getValue(),
			abundanceThreshold.getValue(),
			0);

		std::cout << "Analyzing the graph..." << std::endl;
		Sibelia::BlocksFinder finder(storage, kvalue.getValue());
		finder.FindBlocks(minBlockSize.getValue(),
			maxBranchSize.getValue(),
			maxBranchSize.getValue(),
			8,
			0,
			threads.getValue(),
			outDirName.getValue() + "/paths.txt");

		std::cout << "Generating the output..." << std::endl;
		finder.GenerateOutput(outDirName.getValue(), !noSeq.getValue());
	}
	catch (TCLAP::ArgException & e)
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		return 1;
	}
	catch (std::runtime_error & e)
	{
		std::cerr << "error: " << e.what() << std::endl;
		return 1;
	}

	return 0;
}
