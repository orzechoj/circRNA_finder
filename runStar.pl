#!/usr/bin/perl

use Getopt::Long;
use strict;

###################
## Run STAR to finding circular RNAs
##
## This command runs STAR (version 2.6.0 or later), to find candidates for chimeric transcripts.
## The results of running this script are then post-processed with postProcessStarAlignment.pl

## in:
## - fastq file with R1 reads
## - fastq file with R2 reads
## - directory with genome formatted for STAR
## - directory and prefix of output files (e.g. "~/star_out/lib1_")

## out:
## - results of running STAR are printed to ~/star_out/lib1_<suffix>


#############
## Arguments
my $inFile1 = "";
my $inFile2 = "";
my $genomeDir = "";
my $outPrefix = "";
my $maxMismatchFraction = 0.02;

GetOptions("inFile1=s" => \$inFile1,
           "inFile2=s" => \$inFile2,
           "genomeDir=s" => \$genomeDir,
           "outPrefix=s" => \$outPrefix,
	   "maxMismatch=f" => \$maxMismatchFraction);


#############
## Hard coded parameters
my $starCmd = "/proj/uppstore2017171/private/jacke/circrna_benchmark/STAR/source/STAR";
my $nThreads = 4;
my $chimSegMin = 20;
my $alignIntronMax = 1000000;
my $alignTxPerReadMax = 100000;
my $outFilterMmMax = 4;

my $fullCmd = "$starCmd --genomeDir $genomeDir --readFilesCommand gunzip -c --readFilesIn $inFile1 $inFile2 --runThreadN $nThreads --chimSegmentMin $chimSegMin --chimScoreMin 1 --alignIntronMax $alignIntronMax --outFilterMismatchNoverReadLmax $maxMismatchFraction --alignTranscriptsPerReadNmax $alignTxPerReadMax --twopassMode Basic --outSAMtype BAM SortedByCoordinate --chimOutType Junctions SeparateSAMold --outFilterMultimapNmax 2 --outFileNamePrefix $outPrefix";

system($fullCmd);
