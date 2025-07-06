[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/sibeliaz/README.html)

SibeliaZ 1.2.7
===============

Release date: 6th July 2025
================================

Authors
=======
* Ilia Minkin (Pennsylvania State University)
* Paul Medvedev (Pennsylvania State University)

Introduction
============
SibeliaZ is a whole-genome alignment and locally-coliinear blocks construction
pipeline. The blocks coordinates are output in GFF format and the alignment is
in MAF.

SibeliaZ was designed for the inputs consisting of multiple similar genomes,
like different strains of the same species. The tool works best for the datasets
with the distance from a leaf to the most recent common ancestor not exceeding
0.09 substitutions per site, or 9 PAM units.

Currently SibeliaZ does not support chromosomes in the input longer than
4294967296 bp, this will be fixed in the future releases.

Compilation and installation
============================
The easiset way to install SibeliaZ, is to use [bioconda](https://bioconda.github.io/).
Once you have bioconda environment installed, install package sibeliaz:

	conda install sibeliaz

To compile the code yourself, you will need recent installations of the following
software (Linux only):

* Git
* CMake 
* A GCC compiler supporting C++11

SibeliaZ was tested with CMake 3.5.1 and gcc 5.4.0.
The easiest way to install the dependencies is to use a package management
system. For APT on Debian systems they can be installed by the following
command, which should take only a couple of minutes:

	sudo apt-get install git cmake g++

Once you installed the things above, do the following:

Clone the repository by running:

	git clone https://github.com/medvedevgroup/SibeliaZ 

Go to the root directory of the project and create the "build" folder by
executing:

	cd SibeliaZ
	mkdir build

Initialize dependencies by executing:

	git submodule update --init --recursive

Go to the "build" directory and compile and install the project by running:
	
	cd build
	cmake .. -DCMAKE_INSTALL_PREFIX=<path to install the binaries>
	make install

The make run will produce and install the executables of twopaco, sibeliaz-lcb,
spoa and a wrapper script sibeliaz which implements the pipeline.

SibeliaZ usage
==============
SibeliaZ takes a collection of FASTA file as an input. The simplest way to run
SibeliaZ is to run the following command:

	sibeliaz <input FASTA files>

For example:

	sibeliaz genome1.fa genome2.fa 

The alignment will be reported relative to the sequence ids, so all the input
sequences should have a unique id in the fasta header. By default, the output
will be written in the directory "sibeliaz_out" in the current working
directory. It will contain a GFF file "blocks_coords.gff" containing
coordinates of the found blocks, and file "alignment.maf" with the actual
alignment. The subdirectory "examples" contains an example of running
SibeliaZ and the output it produces. Running SibeliaZ on this example should
require less than 5 minutes on a typical machine. SibeliaZ has several
parameters that affect the accuracy and performance, they are described below.

Building synteny blocks
=======================
You can also construct longer synteny blocks required for certain analyses. By
default, sibeliaz also installs program [maf2synteny](https://github.com/fenderglass/maf2synteny)
written by [Mikhail Kolmogorov](https://github.com/fenderglass). If you already
have output of SibeliaZ, you can run maf2synteny on it:

	maf2synteny <GFF or MAF file created by SibeliaZ>

Otherwise, it is best to run sibeliaz without the alignment step saving time
and memory:
	
	sibeliaz -n <output_directory>
	maf2synteny <output_directory>/blocks_coords.gff

For further information on synteny blocks construction, please refer to the 
documentation of maf2synteny.

Output description
==================
The output directory will contain:

1) A GFF file with coordinates of the locally-collinear blocks. Lines that
have identical id fields correspond to different copies of the same block.
The file name is "blocks_coords.gff"
2) A MAF file with the whole-genome alignment of the input. The file name
is "alignment.maf".

Note: the actual alignment is produced by globally aliging the locally-colliner
blocks, which is memory-hungry. It could be impossible to align certain blocks,
especially if they have a lot of copies and/or long due to the aligner running
out of memory, even on machines with large RAM. The output directory will have
a subdirectory "blocks" that will contain FASTA files with blocks that were
impossible to align. Each file correspond to a block and contains its copies.
FASTA headers contain the coordinates of all copies of the block in the same
format as MAF records, except that fields are separated by a semicolon.

It is possible to skip the alignment (use the -n switch) step and produce only
coordinates of the blocks if the alignment is not needed for downstream analysis.
In this case SibeliaZ will not produce the "alignment.maf" file and "blocks"
subdirectory.

Parameters affecting accuracy
=============================

The value of k
--------------
This parameter defines the order of the de Bruijn graph being used and controls
the tradeoff between the sensitivity on one hand, and speed and memory usage
on the other. The parameter is set by the key

	-k <an odd integer>

In general the lower the k, the slower and more sensitive the alignment is. For
small datasets, like bacteria, we recommend k=15, and for mammalian-sized
genomes k=25. The default is 25.

Vertices frequency threshold
----------------------------
Input genomes often contain repeated elements that make the graph large and
convoluted. To deal with this issue, SibeliaZ removes all k-mers with frequency
more than a threshold, which is controlled by the option:

	-a <integer>

We recommend setting it to twice the maximum number of copies of a homologous
block the underlying input genome collection has. For example, assume a single 
input genome in the collection of the input genomes has a genomic element that
is duplicated **D** times and the user is interested in this element being properly
aligned across the input. This element could be a gene that was duplicated multiple
times or a transposon; we assume that **D** is the maximum frequency of occurrence
across all such genomic elements of interest in a single genome.

If there are **N** input genomes containing the genomic element,
then the parameter -a should be set to at least **D * N * 2**, otherwise,
the genomic element will likely be missing in the resulting alignment. 
If the resources allow, this parameter could be set to a higher value, but setting
it too high could significantly delay the computation.
For example, suppose we align 10 genomes, and each genome has a gene that is duplicated 4 times.
Hence, the parameter -a should be set to at least 10 * 4 * 2 = 80 or higher.

The default parameter value is 150 and is optimized for a collection of mammalian genomes.
This parameter must be adjusted accordingly in case the input contains hundreds
of short genomes like bacteria or viruses. 


Bubble size threshold
---------------------
SibeliaZ analyzes the graph by looking for long chains of bubbles in it. A bubble
is a pair of paths having the same endpoints. A long chain of bubbles is likely
to be generated by a pair of homologous sequences, which SibeliaZ looks for.
However, if the paths between endpoints of a bubble is too long, it may arise
through the spurious similarity. To avoid this, SibeliaZ discards bubbles 
with paths longer than the threshold -b, which can be set by:

	-b <integer>

The default value of -b is 200. Increasing value may increase recall of divergent
sequences, but if -b is too high, it will decrease accuracy as well.

Locally-collinear block size
----------------------------
SibelaZ only output blocks longer than a specified threshold, which is set by

	-m <integer>

The default value is 50. Warning: increasing this parameter may significantly
slow down the computation.

Technical parameters
====================

Skipping the alignment
----------------------
To skip the alignment and only output coordinates of the blocks, use the
switch

	-n

Threads number
--------------
The maximum number of thread for SibeliaZ to use. This parameter is set by 

	-t <integer>

By default SibeliaZ tries to use as much threads as possible. You can limit
this number by using the above switch. Note that different stages of the
pipeline have different scalabilities. TwoPaCo will not use more than
16 threads, while graph analyzer sibeliaz-lcb and the global aligner will use
as much as possible.

Memory allocation
-----------------
The graph constructor TwoPaCo preallocates memory for Bloom filter. By default,
the Bloom filter size is thrice of the size of the input files. The Bloom
filter size can be set manually with the option:

	-f <memory amount in GB>

Output directory
----------------
The directory for the output files can be set by the argument

	-o <directory>

The default is "sibeliaz_out" in the current working directory.

A note about the repeat masking
==============================
SibeliaZ and TwoPaCo currently do not recognize soft-masked characters (i.e. using
lowercase characters), so please convert soft-masked repeats to hard-maksed ones
(with Ns) if you would like to mask the repeats explicitly. However, it is not
necessary as SibeliaZ uses the abundance parameter -a to filter out high-copy
repeats.

Difference between Sibelia and SibeliaZ
=======================================
SibeliaZ is the future developement of synteny-finder Sibelia. The key difference
is that old Sibelia was designed to produce long synteny blocks, while SibeliaZ
produces shorter locally-collinear blocks or LCBs. Output of SibeliaZ is very
similar to Sibelia's when it is run in a single stage. At the same time, SibeliaZ
is much faster and can handle longer genomes.

Export to GFA1 (experimental)
=============================
The script located at Sibeliaz-LCB/maf_to_gfa1.py lets you convert a MAF file
produced by SibeliaZ to a GFA1 file representing a graph induced by the alignment.
The GFA1 file then can be imported into [vg](https://github.com/vgteam/vg) or visualized. Usage:

	python maf_to_gfa1.py <MAF alignment file> <input FASTA files>

Conversion to XMFA
==================
The script located at Sibeliaz-LCB/maf_to_gfa1.py lets you convert a MAF file
to XMFA format. Requires BioPython of version >= 1.6.9. Usage:

        python maf_to_xmfa.py < <MAF alignment file>


Troubleshooting
===============
It could be that SibeliaZ runs out of memory on large inputs. Possible reasons
include:

* TwoPaCo having the Bloom filter too small. To increase its size, use the -f switch

* SibeliaZ-LCB running out of memory. You can try to reduce the abundance parameter
-a to prune the internal data structure and reduce its size

Citation
========
If you use SibeliaZ in your research, please cite:

	Minkin, I., Medvedev, P. Scalable multiple whole-genome alignment and locally collinear block construction with SibeliaZ.
	Nat Commun 11, 6327 (2020). https://doi.org/10.1038/s41467-020-19777-8

If you also used maf2synteny, please cite the Ragout [paper](https://github.com/fenderglass/maf2synteny#citation).

License
=======
See LICENSE.txt

Contacts
========
E-mail your feedback at ivminkin@gmail.com.

Datasets used for analyses in the paper
======================================
See: https://github.com/medvedevgroup/SibeliaZ/blob/master/DATA.txt
