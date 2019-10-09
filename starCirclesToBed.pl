#!/usr/bin/perl

## use strict;

#####################################################################################################
## Converts the post processed circular RNA data from STAR + filterCirc.awk to bed format:
##
## - chromosome
## - start (first nucleotide of circle, 0 based)
## - end (last nucleotide of circle, 0 based)
## - name (unique nr + "s" if the circle matches a splice junction)
## - score (nr of supporting reads)
## - strand
##
## Also takes an argument for minimum genomic lengh. Cadidate circular RNAs shorter than this are removed.
##
## ex: prog/starCirclesToBed.pl MHcircCollapsed.txt 200


#############
## Arguments
my $starPostFile = $ARGV[0];
my $minLen = $ARGV[1];

################
## 'Main'
open IN, $starPostFile or die "Cannot open $starPostFile\n";
my $i = 0;

while(<IN>){
  chomp;
  my @fields = split /\s+/;

  my $score = $fields[1];
  my $chr = $fields[2];
  my $start = $fields[3];
  my $end = $fields[4] - 1;
  my $strand = $fields[5];
  my $hasSplice = $fields[6];

  $i++;
  my $useId = $i;
  if($hasSplice){ $useId = $useId."s"; }

  if(abs($end-$start) >= $minLen){
    print $chr."\t".$start."\t".$end."\t".$useId."\t".$score."\t".$strand."\n";
  }
}
