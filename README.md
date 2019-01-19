SibeliaZ 1.0.0
===============

Authors
=======
* Ilia Minkin (Pennsylvania State University)
* Paul Medvedev (Pennsylvania State University)

Introduction
============
SibeliaZ is a whole-genome alignment pipeline based on analysis of the 
compacted de Bruijn graph. SibeliaZ first constructs the graph, finds
homologous blocks free of rearrangements and then aligns them globally.
The graph is constructed by TwoPaCo and the the blocks are aligned by spoa.
The blocks coordinates are output in GFF format and the alignment is in MAF.

SibeliaZ works best for the inputs consisting of multiple similar genomes,
like different strains of the same species. The tool is designed for the
datasets with the distance from a leaf to the most recent common ancestor
not exceeding 0.085 substitutions per site, or 9 PAM units. It can be 
used for more divergent datasets, but it will have worse recall of divergent
blocks.

Difference between Sibelia and SibeliaZ
=======================================
SibeliaZ is the future developement of synteny-finder Sibelia. The key difference
is that old Sibelia was designed to produce long synteny blocks, while SibeliaZ 
produces shorter locally-collinear blocks or LCBs. Output of SibeliaZ is very
similar to Sibelia's when it is run in a single stage. Locally collinear blocks can
be globally aligned to produce a whole genome alignment, or they can be chained to
get longer synteny blocks. At the same time, SibeliaZ is much faster and can
handle longer genomes.

Compilation
===========
To compile the code, you need the following (Linux only):

* CMake 
* A GCC compiler supporting C++11
* Intel TBB library properly installed on your system. In other words, G++
  should be able to find TBB libs (future releases will not depend on TBB)

Once you installed the things above, do the following:

* Clone the repository https://github.com/medvedevgroup/SibeliaZ
* Go to the root directory of the project and create the "build" folder
* Type git submodule update --init --recursive
* Go to the "build" directory
* Run cmake .. -DCMAKE_INSTALL_PREFIX=<path to install the binaries>
* Run make install

The make run will produce the executables of twopaco, sibeliaz-lcb, spoa and
a wrapper script sibeliaz which implements the pipeline.

SibeliaZ usage
===============
SibeliaZ takes a FASTA file as an input. The simplest way to run SibeliaZ
is to enter the following command:

	sibeliaz -t <number of threads> -f <memory amount in GB> <input FASTA file>

which will run the pipeline using at most -t threads and allocating -f GBs of 
memory. The limit specified by -f is "soft", sibeliaz may actually use more
memory (see detailed description below). For small datasets, it at least should
be several times of the input file size, for large genomes it is best to use
all memory available. 

The output will be written in the directory "sibeliaz_out" in the current
working directory. It will contain a GFF file "blocks_coords.gff" containing
coordinates of the found blocks, a directory "block" with FASTA file containing
the sequences of the blocks, and file "alignment.maf" with the actual alignment.
It is possible to skip the alignment (use the -n switch) step and produce only
coordinates of the blocks if the alignment is not needed for downstream analysis.

SibeliaZ has other parameters affecting the running time, sensitivity and the
output, which are described below.

Output description
==================
The output directory will contain:

1) A GFF file with coordinates of the locally-collinear blocks. Lines that
have identical id fields correspond to different copies of the same block.
The file name is "blocks_coords.gff"
2) A MAF file with the whole-genome alignment of the input. The file name
is "alignment.maf"
3) A directory with sequences of blocks in FASTA format. Each file correspond
to a block and contains its copies. FASTA headers contain the coordinates
of all copies of the block in the same format as MAF records, except that
fields are separated by a semicolon.

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
convoluted. To deal with this issue, SibeliaZ removes all vertices corresponding
to k-mers with frequence more than a threshold, which is controlled by the option:

	-a <integer>

The default value is 150. Increasing this value may significantly slow down the
alignment.

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
sequences, but if -b is too high, it will decrease accuracy as well. The value
of -b aslo should not exceed the value of -m due to the properties of our 
algorithm.

Locally-collinear block size
----------------------------
SibelaZ only output blocks longer than a specified threshold, which is set by

	-m <integer>

The default value is 250. This value can be increased to filter out short
homologous blocks if it is necessary for the downstream analysis. The value
of -m should be larger than -b due to algorithmic reasons.


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

For the best performance allocate as many threads as possible (but no more 
than the number of cores in the machine). Note that different stages of 
the pipeline have different scalabilities. TwoPaCo will not use more than
16 threads, graph analyzer sibeliaz-lcb will not use more than 4, while
the global aligner will use as much as possible.

Memory allocation
-----------------
Some prgorams in the pipeline require an amount of memory specified beforehand.
For example, the graph constructor TwoPaCo uses it for Bloom filter, and the
global aligner requires a memory limit to be set to work correctly. This can
be set usin the option:

	-f <memory amount in GB>

Roughly speaking the memory should be at least 2-3 times of the size of the 
input genomes, and for large dataset it is adviced to allocate as much as
possible.


Output directory
----------------
The directory for the output files can be set by the argument

	-o <directory>

The default is "sibeliaz_out" in the current working directory.

A note about the repeat masking
==============================
SibeliaZ TwoPaCo currently does not recognize soft-masked characters (i.e. using
lowercase characters), so please convert soft-masked repeats to hard-maksed ones
(with Ns). 

Citation
========
If you use SibeliaZ in your research, please cite:


License
=======

Contacts
========
E-mail your feedback at ivminkin@gmail.com.

Datasets used of analyses in the paper
======================================
1) Simulated bacterial genomes:
2) Mice gene alignment:
