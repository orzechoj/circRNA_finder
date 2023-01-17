#!/bin/bash
set -euo pipefail

# Test post-processing of paired end data
./postProcessStarAlignment.pl --starDir tests/data/STAR_SRR1197328/ --minLen 200 --outDir ./
diff -u filteredJunctions.bed tests/expected/filteredJunctions.bed 
