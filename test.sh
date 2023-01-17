#!/bin/bash
set -euo pipefail


# Generate STAR index
mkdir mbl_locus_star_genome
../star_2.7.2d/STAR-2.7.2d/source/STAR --runThreadN 4 \
     --runMode genomeGenerate \
     --genomeDir mbl_locus_star_genome \
     --genomeSAindexNbases 8 \
     --genomeFastaFiles tests/data/mbl_locus.fa \
     --sjdbGTFfile tests/data/mbl_locus.gtf \
     --sjdbOverhang 100

# Test STAR
mkdir STAR_SRR1197328_mbl
./runStar.pl --inFile1 tests/data/SRR1197328_mbl_1.fastq.gz --inFile2 tests/data/SRR1197328_mbl_2.fastq.gz --genomeDir mbl_locus_star_genome --outPrefix STAR_SRR1197328_mbl/

# Test post-processing of paired end data
./postProcessStarAlignment.pl --starDir STAR_SRR1197328_mbl/ --minLen 200 --outDir ./
diff -u filteredJunctions.bed tests/expected/filteredJunctions.bed 
