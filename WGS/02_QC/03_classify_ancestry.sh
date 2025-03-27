#!/bin/bash
#$ -cwd
#$ -pe smp 4 -binding linear:4
#$ -l h_rt=12:00:00,h_vmem=25G
source /broad/software/scripts/useuse
reuse -q R-4.3

pcadir=$1
pcplot0=$2
npc=$3
pred_prob=$4
batch=$5
scrdir="/stanley/huang_lab/shared/data/sc-asia/script"

Rscript $scrdir/03_classify_ancestry.R $pcadir $pcplot0 $npc $pred_prob $batch
