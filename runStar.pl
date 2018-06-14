#!/usr/bin/perl

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
my $inFile1 = $ARGV[0];
my $inFile2 = $ARGV[1];
my $genomeDir = $ARGV[2];
my $outPrefix = $ARGV[3];

#############
## Hard coded parameters
my $starCmd = "STAR";
my $nThreads = 2;
my $chimSegMin = 20;
my $alignIntronMax = 500000;
my $alignTxPerReadMax = 100000;
my $outFilterMmMax = 4;

my $fullCmd = "$starCmd --genomeDir $genomeDir --readFilesCommand gunzip -c --readFilesIn $inFile1 $inFile2 --runThreadN $nThreads --chimSegmentMin $chimSegMin --chimScoreMin 1 --alignIntronMax $alignIntronMax --outFilterMismatchNmax $outFilterMmMax --alignTranscriptsPerReadNmax $alignTxPerReadMax --twopassMode Basic --outSAMtype BAM SortedByCoordinate --chimOutType Junctions SeparateSAMold --outFilterMultimapNmax 2 --outFileNamePrefix $outPrefix";

system($fullCmd);
