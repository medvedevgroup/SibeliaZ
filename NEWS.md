SibeliaZ 1.2.0
==============
* The results are now reproducible regardless to the number of threads
being used
* Changed behavior of "-f" flag: now it is optional and only controls size of
the Bloom filter allocated by TwoPaCo. By default, it is now thrice of the
total input file size
* Fixed a bug in handling corrupt FASTA files
* Slightly improved accuracy
* Added a script to convert MAF alignment to XMFA (maf_to_xmfa.py)
                                                                  
SibeliaZ 1.1.0
==============
* Added a script for converting MAF alignment to GFA1 (maf_to_gf1.py)
* Input genomes can be provided in multiple FASTA files now
* Fixed a bug causing some blocks to overlap by k basepairs
* Fixed a bug resulting some exact matching blocks to be missing

SibeliaZ 1.0.0
==============
* Hello, World!
