circRNA_finder
==============

Scripts required for running the pipeline to find circular RNAs from RNA-seq data, as used in

Jakub O. Westholm, Pedro Miura, Sara Olson, Sol Shenker, Brian Joseph, Piero Sanfilippo, Susan E. Celniker, Brenton R. Graveley, and Eric C. Lai. Genome-wide Analysis of Drosophila Circular RNAs Reveals Their Structural and Sequence Properties and Age-Dependent Neural Accumulation. Westholm et al. Cell Reports, 2014.

Contains the following files:
- filterCirc.awk
- filterSpliceSiteCircles.pl
- nrForwardSplicedReads.pl
- postProcessStarAlignment.pl
- runStar.pl
- starCirclesToBed.pl


These scripts have been tested on various Linux distributions. Before they can be run, make sure that the following prerequisites are installed:
 - perl
 - awk
 - STAR (versions 2.4.1c, 2.3.1o and 2.3.1s.t have worked, not version 2.3.0)
 - samtools


To run the scripts to identify circular RNAs, first run STAR, once for each data set:

./runStar.pl [R1 fastq] [R2 fastq] [path to STAR genome] [output directory and prefix]


Next, run the post processing scripts. If there are STAR outputs for many data sets in the same folder, this command will process each of these in turn:

./postProcessStarAlignment.pl [directory with STAR results] [output directory]


For each library the following output files are produced:

a) <lib name>_filteredJunctions.bed: A bed file with all circular junctions found by the pipeline. The score column indicates the  number reads spanning each junction.

b) <lib name>_s_filteredJunctions.bed: A bed file with those juction in (a) that are flanked by GT-AG splice sites. The score column indicates the  number reads spanning each junction.

c) <lib name>_s_filteredJunctions_fw.bed: A bed file with the same circular junctions as in file (b), but here the score column gives the average number of forward spliced reads at both splice sites around each circular junction.
