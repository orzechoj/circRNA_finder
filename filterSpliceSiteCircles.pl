#!/usr/bin/perl

use strict;

#####################################################################################################
## Perl script for filtering bed files with circular RNAs.
## Circular RNAs with "s" in the name (e.g. "24s") are kept and printed to result file. All other are removed.

## ex: prog/filterSpliceSiteCircles.pl in.bed > out.bed


#############
## Arguments
my $inBedFile = $ARGV[0];


################
## 'Main'
open IN, $inBedFile or die "Cannot open $inBedFile\n";
my $i = 0;

while(<IN>){
  chomp;
  my @fields = split /\s+/;

  my $chr = $fields[0];
  my $start = $fields[1];
  my $end = $fields[2];
  my $name = $fields[3];
  my $score = $fields[4];
  my $strand = $fields[5];

  if($name =~ /s/){
    print $chr."\t".$start."\t".$end."\t".$name."\t".$score."\t".$strand."\n";
  }
}
