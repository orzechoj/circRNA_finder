#!/usr/bin/perl

use strict;
use File::Basename;

##############################################################################################################
## Script for parsing STAR output to get number of reads supporting normal splice events around circular RNAs.
## Takes the following arguments
## - A bed file with circular RNAs
## - A "SJ.out.tab" file created by STAR
## 
## Prints the input bed table to STDOUT, but with a new score which is the total nr of forward spliced reads
## at both splice sites of the circular RNA / 2.
##
## ex: ./getNrSpliceReads.pl ~/projects/feuk_2015/data/circle_pipeline_out/cyto/up_063_2_s_filteredJunctions.bed ~/tmp/SJ/up_063_2_SJ.out.tab > ~/tmp/fwSplice.txt


#############
## Arguments
my $circleBedFile = $ARGV[0];
my $starSjFile = $ARGV[1];

my %spliceBefore = (); ## Map chr_pos_strand -> nr reads
my %spliceAfter = ();  ## Map chr_pos_strand -> nr reads


###############################
## Read forward spliced reads from STAR output
open SJ_FILE, $starSjFile or die "Cannot open $starSjFile\n";
while(<SJ_FILE>){
  chomp;
  my @fields = split;
  my $chr = $fields[0];
  my $start = $fields[1];
  my $end = $fields[2];
  my $strand = $fields[3];
  my $reads= $fields[6];
  
  ## Convert strand from +/- to 1/2
  if($strand eq "1"){ $strand = "+"; } 
  if($strand eq "2"){ $strand = "-"; } 
  
  ## Adjust coordinates
  my $start = $start -1;
  my $end = $end;

  ## Get spliced reads ending at this coordinate
  my $beforeStr = $chr."_".$end."_".$strand;
  $spliceBefore{$beforeStr} = $spliceBefore{$beforeStr} + $reads;

  ## Get spliced reads starting at this coordinate
  my $afterStr = $chr."_".$start."_".$strand;
  $spliceAfter{$afterStr} = $spliceBefore{$afterStr} + $reads;
}

#################################
## Read bed file with circular RNAs, and for each get the nr of forward spliced reads
open BED_FILE, $circleBedFile or die "Cannot open $circleBedFile\n";
while(<BED_FILE>){
  chomp;
  my @fields = split;
  my $chr = $fields[0];
  my $start = $fields[1];
  my $end = $fields[2];
  my $id = $fields[3];
  my $score = $fields[4];
  my $strand = $fields[5];

  my $beforeStr = $chr."_".$start."_".$strand;
  my $afterStr = $chr."_".$end."_".$strand;
  my $newScore = $spliceAfter{$afterStr} + $spliceBefore{$beforeStr} / 2;

  print $chr."\t".$start."\t".$end."\t".$id."\t".$newScore."\t".$strand."\n";
}



