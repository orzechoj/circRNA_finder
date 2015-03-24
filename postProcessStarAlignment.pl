#!/usr/bin/perl

use strict;
use File::Basename;
use FindBin;

#######################################################################################################
## Script for post processing the results from the STAR alignment, to filter and format circular RNAs.
## Takes the following arguments:
## - a directory with STAR output
## - a directory for the bed files with circular RNAs
##
## For each library, creates the following files
## - bed file with all circular junctions from STAR
## - bed file with circular junction supported by splice sites
## - indexed bam files with all chimeric reads


#############
## Arguments
my $inDir = $ARGV[0];
my $outDir = $ARGV[1];

## The scripts used to process STAR output
my $scriptDir = $FindBin::Bin;
my $filterStarProg = $scriptDir."/filterCirc.awk";
my $makeBedFileProg = $scriptDir."/starCirclesToBed.pl";
my $filterBedFileProg = $scriptDir."/filterSpliceSiteCircles.pl";


################
## 'Main'

## Get all names of files and libraries to process
my @samFiles = `ls $inDir*Chimeric.out.sam`;
my @junctionFiles = `ls $inDir*Chimeric.out.junction`;
my @libnames = ();
foreach my $currentLib (@junctionFiles){
  chomp $currentLib;
  my $newName = basename($currentLib);
  $newName =~ s/Chimeric.out.junction//g;
  push @libnames, $newName;
}
my $nrLibs = @libnames;

foreach my $i (1..$nrLibs){
  print STDERR "****** Processing ".$libnames[$i-1]." (".$i."/".$nrLibs.") ******\n";

  #####
  ## Filter to get the circular junctions from the chimeric STAR output
  my $filteredStarFile = $outDir.$libnames[$i-1]."filteredJunctions.txt";
  my $filterJunctionsCmd = "cat $junctionFiles[$i-1] | awk -f $filterStarProg | sort | uniq -c | sort -k1,1rn > $filteredStarFile";
  print STDERR "Filtering chimeric transcripts...\n";
  system($filterJunctionsCmd);



  ####
  ## Comvert circular junctions to bed format
  my $allCirclesBedFile = $outDir.$libnames[$i-1]."filteredJunctions.bed";
  my $makeBedFileCmd = "$makeBedFileProg $filteredStarFile > $allCirclesBedFile";
  print STDERR "Converting to .bed file...\n";
  system($makeBedFileCmd);


  #####
  ## Make bed file with only circular junctions supported by splice sites
  my $spliceCirclesBedFile = $outDir.$libnames[$i-1]."s_filteredJunctions.bed";
  my $filterBedFileCmd = "$filterBedFileProg $allCirclesBedFile > $spliceCirclesBedFile";
  print STDERR "Filtering circle junctions on splice sites...\n";
  system($filterBedFileCmd);


  ######
  ## Make indexed bed files
  chomp $samFiles[$i-1];
  my $junctionBamFile = $samFiles[$i-1];
  $junctionBamFile =~ s/\.sam/\.bam/g;
  my $junctionSortedBamFile = $junctionBamFile;
  $junctionSortedBamFile =~ s/\.bam/\.sorted/g;
  my $makeBamCmd = "samtools view -bS -o  $junctionBamFile $samFiles[$i-1]";
  my $sortBamCmd = "samtools sort $junctionBamFile $junctionSortedBamFile";
  my $indexBamCmd = "samtools index $junctionSortedBamFile.bam";

  print STDERR "Creating .bam file...\n";
  print $makeBamCmd."\n";
  system($makeBamCmd);
  print STDERR "Sorting .bam file...\n";
  system($sortBamCmd);
  print STDERR "Indexing .bam file...\n";
  system($indexBamCmd);

  ######
  ## Clean up
  system("rm ".$filteredStarFile);
  system("rm ".$junctionBamFile);
}
