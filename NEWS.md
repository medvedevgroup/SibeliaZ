SibeliaZ 1.2.4
==============
* Made GFF output compliant with the [GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).
* Added a pragma to GFF output showing the lengths of input sequences
* Fixed the incorrect path to spoa in bioconda version

SibeliaZ 1.2.3
==============
* Reduced the number of temporary files created during the global alignment phase
* Removed dependency on GNU time from the main SibeliaZ script

SibeliaZ 1.2.2
==============
* Updated the version of spoa
* Added submodule maf2synteny, a program for constructing synteny blocks from
the output of Sibeliaz
* Fixed a bug in SibeliaZ script resulting in failed global alignment calls in
certain cases

SibeliaZ 1.2.1
==============
* Updated the version of TwoPaCo
* Fixed a bug resulting in missing of some collinear blocks
* Added script "glue_gfa1.py" that glues back GFA1 produced by "maf_to_gfa1.py"
to the original genomes
* Fixed inconsistent use of spaces and tabs in "maf_to_gfa1.py"
* MAF file now contains the version number and the arguments used

SibeliaZ 1.2.0
==============
* The results are now deterministic and do not depend on the number of threads
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
