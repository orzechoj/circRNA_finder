circRNA_finder
==============

Scripts required for running the pipeline to find circular RNAs from RNA-seq data, as used in (Westholm et al, Cell Reports, 2014).

Contailn the following files:
- filterCirc.awk
- filterSpliceSiteCircles.pl
- postProcessStarAlignment.pl
- runStar.pl
- starCirclesToBed.pl


These scripts have been tested on various Linux distributions. Before they can be run, make sure that the following prerequisites are installed:
 - perl
 - awk
 - STAR (version 2.3.1o and 2.3.1s.t have worked, not version 2.3.0)
 - samtools


To run the scripts to identify circular RNAs, first run STAR, once for each data set:

./runStar.pl <R1 fastq> <R2 fastq> <path to STAR genome> <output directory and prefix>


Next, run the post processing scripts. If there are STAR outputs for many data sets in the same folder, this command will process each of these in turn:

./postProcessStarAlignment.pl <directory with STAR results> <output directory>


