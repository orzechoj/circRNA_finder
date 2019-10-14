circRNA_finder
==============

Scripts required for running the pipeline to find circular RNAs from RNA-seq data, as used in

Jakub O. Westholm, Pedro Miura, Sara Olson, Sol Shenker, Brian Joseph, Piero Sanfilippo, Susan E. Celniker, Brenton R. Graveley, and Eric C. Lai. [Genome-wide Analysis of Drosophila Circular RNAs Reveals Their Structural and Sequence Properties and Age-Dependent Neural Accumulation](https://www.cell.com/cell-reports/abstract/S2211-1247(14)00931-0) Westholm et al. Cell Reports, 2014.


These scripts have been tested on various Linux distributions. Before they can be run, make sure that the following prerequisites are installed:
 - perl
 - awk
 - STAR (tesed on version 2.7.2d)
 - samtools


To run the scripts to identify circular RNAs, first run STAR, once for each data set. For paired end data, the command is

```bash
./runStar.pl --inFile1 [R1 fastq] --inFile2 [R2 fastq] --genomeDir [path to STAR genome] --maxMismatch [max mismatches realtive to read length, default 0.02] --outPrefix [output directory and prefix]
```

Here `--maxMismatch` sets the `outFilterMismatchNoverReadLmax` parameter in STAR.



Next, run the post processing scripts. If there are STAR outputs for many data sets in the same folder, this command will process each of these in turn:

```bash
./postProcessStarAlignment.pl --starDir [directory with STAR results] --minLen [minimum length of circular RNAs] --outDir [output directory]
```

For each library the following output files are produced:

a) <lib name>_filteredJunctions.bed: A bed file with all circular junctions found by the pipeline. The score column indicates the  number reads spanning each junction.

b) <lib name>_s_filteredJunctions.bed: A bed file with those juction in (a) that are flanked by GT-AG splice sites. The score column indicates the  number reads spanning each junction.

c) <lib name>_s_filteredJunctions_fw.bed: A bed file with the same circular junctions as in file (b), but here the score column gives the average number of forward spliced reads at both splice sites around each circular junction.

d) (Sorted and indexed) bam file with all chimeric reads identified by STAR. The circRNA junction spanning reads are a subset of these.
