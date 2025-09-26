#!/bin/bash
set -e

echo "Bioinformatics container started"
echo "Working directory: /work"
echo "EggNOG DB directory: /data/eggnog"
echo "Running Pangenomic Automated Workflow..."
echo ""
cd /work/Pangenomics-Pipeline/codes
exec bash workflow1.sh "$@"
