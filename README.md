SibeliaZ 1.2.0
===============

Release date: 30th August 2019
=================

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
0.085 substitutions per site, or 9 PAM units.

Currently SibeliaZ does not support chromosomes in the input longer than
4294967296 bp, this will be fixed in the future releases.

Compilation and installation
============================
To compile the code, you need recent installations of the following software
(Linux only):

* Git
* CMake 
* A GCC compiler supporting C++11
* Intel TBB library properly installed on your system. In other words, G++
  should be able to find TBB libs (future releases will not depend on TBB)

The easiest way to install the dependencies is to use a package management
system, for APT on Debian systems they can be installed by the following:

	sudo apt-get install git cmake g++ libtbb-dev

Once you installed the things above, do the following:

Clone the repository https://github.com/medvedevgroup/SibeliaZ by
running:

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
SibeliaZ and the output it produces. SibeliaZ has several parameters that
affect the accuracy and performance, they are described below.

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
Mammalian genomes contain many repeated elements that make the graph large and
convoluted. To deal with this issue, SibeliaZ removes all k-mers with frequency
more than a threshold, which is controlled by the option:

	-a <integer>

We recommend to set it to the twice of the maximum number of copies a homologous
block in the input genomes has. For example, if the largest gene family of the input
genomes has N members, set -a to at least N * 2. However, increasing this value may
significantly slow down the computation. The default value is 150.

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
The GFA1 file then can be imported into vg (https://github.com/vgteam/vg) or
visualized. Usage:

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

	Scalable multiple whole-genome alignment and locally collinear block construction with SibeliaZ
	Ilia Minkin, Paul Medvedev
	bioRxiv 548123; doi: https://doi.org/10.1101/548123

License
=======
See LICENSE.txt

Contacts
========
E-mail your feedback at ivminkin@gmail.com.

Datasets used of analyses in the paper
======================================
See: https://github.com/medvedevgroup/SibeliaZ/blob/master/DATA.txt
