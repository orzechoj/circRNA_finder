#!/usr/bin/perl

use strict;
use File::Basename;
use Getopt::Long;
use FindBin;

#######################################################################################################
## Script for post processing the results from the STAR alignment, to filter and format circular RNAs.
## Takes the following arguments:
## - a directory with STAR output
## - the minimum genomic length of the circular RNAs (optional, defualt is 200nt).
## - a directory for the bed files with circular RNAs
##
## For each library, creates the following files
## a) a bed file with all circular junctions from STAR, where score is the nr
##    of junction spanning reads
## b) a bed file with circular junction supported by splice sites, where score
##    in the nr of junction spanning reads
## c) a bed file with circular junction supported by splice sites, where score
##    in the nr of forward spliced reads
## d) an indexed bam file with all chimeric reads


#############
## Arguments
my $inDir = "";
my $minLen = 200;
my $outDir = "";

GetOptions("starDir=s" => \$inDir,
           "outDir=s" => \$outDir,
           "minLen=f" => \$minLen);

if (not defined $inDir or not defined $outDir) {
    die "Run the post processung script like this: postProcessStarAlignment.pl <directory with STAR data> <output directory>\n";
}

## The scripts used to process STAR output
my $scriptDir = $FindBin::Bin;
my $filterStarProg = $scriptDir."/filterCirc.awk";
my $makeBedFileProg = $scriptDir."/starCirclesToBed.pl";
my $filterBedFileProg = $scriptDir."/filterSpliceSiteCircles.pl";
my $nrFwSplicedReadsProg = $scriptDir."/nrForwardSplicedReads.pl";



################
## 'Main'

## Get all names of files and libraries to process
my @samFiles = `ls $inDir*Chimeric.out.sam`;
my @chimericFiles = `ls $inDir*Chimeric.out.junction`;
my @junctionFiles = `ls $inDir*SJ.out.tab`;

my @libnames = ();
foreach my $currentLib (@chimericFiles){
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
    my $filterJunctionsCmd = "cat $chimericFiles[$i-1] | awk -f $filterStarProg | sort | uniq -c | sort -k1,1rn > $filteredStarFile";
    print STDERR "Filtering chimeric transcripts...\n";
    system($filterJunctionsCmd);


    ####
    ## Convert circular junctions to bed format and filter on length
    my $allCirclesBedFile = $outDir.$libnames[$i-1]."filteredJunctions.bed";
    my $makeBedFileCmd = "$makeBedFileProg $filteredStarFile $minLen > $allCirclesBedFile";
    print STDERR "Converting to .bed file...\n";
    system($makeBedFileCmd);


    #####
    ## Make bed file with only circular junctions supported by splice sites
    my $spliceCirclesBedFile = $outDir.$libnames[$i-1]."s_filteredJunctions.bed";
    my $filterBedFileCmd = "$filterBedFileProg $allCirclesBedFile > $spliceCirclesBedFile";
    system($filterBedFileCmd);


    #####
    ## Make bed file for circular junctions supported by splice sites, with nr of forward spliced reads
    chomp $junctionFiles[$i-1];
    my $tmpJunctionFile = $junctionFiles[$i-1];
    my $fwSplicedBedFile = $outDir.$libnames[$i-1]."s_filteredJunctions_fw.bed";
    my $getFwSpliceCmd = "$nrFwSplicedReadsProg $spliceCirclesBedFile $tmpJunctionFile > $fwSplicedBedFile";
    print STDERR "Getting forward spliced reads for circle junctions with splice sites...\n";
    system($getFwSpliceCmd);


    ######
    ## Make indexed bam files
    chomp $samFiles[$i-1];
    my $junctionBamFile = $outDir."/".basename($samFiles[$i-1]);
    $junctionBamFile =~ s/\.sam/\.bam/g;
    my $junctionSortedBamFile = $junctionBamFile;
    $junctionSortedBamFile =~ s/\.bam/\.sorted.bam/g;
    my $makeBamCmd = "samtools view -bS -o $junctionBamFile $samFiles[$i-1]";
    my $sortBamCmd = "samtools sort -o $junctionSortedBamFile $junctionBamFile";
    my $indexBamCmd = "samtools index $junctionSortedBamFile";

    print STDERR "Creating .bam file...\n";
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
