#!/bin/bash
#$ -cwd
#$ -pe smp 4 -binding linear:4
#$ -l h_rt=12:00:00,h_vmem=5G

source /broad/software/scripts/useuse
reuse -q R-4.3

datadir=$1
preqcdir=${datadir}/preimp_qc
pcadir=${datadir}/pca
pop=$2
pop_preqcdir=${datadir}/pop_qc
scrdir="/stanley/huang_lab/shared/data/sc-asia/script"
batch=$3

PLINK="/home/unix/yu/software/00install_softwares/bin/plink"

cd $pop_preqcdir

# Perform LD pruning: chrX
$PLINK \
--bfile ${pop}_${batch}_btqc_mgqc \
--chr 23 \
--geno 0.02 \
--maf 0.05 \
--snps-only just-acgt \
--indep-pairwise 200 100 0.1 \
--out ${pop}_${batch}_btqc_mgqc-ldpr-chrx

# Check sex
$PLINK \
--bfile ${pop}_${batch}_btqc_mgqc \
--extract ${pop}_${batch}_btqc_mgqc-ldpr-chrx.prune.in \
--check-sex \
--out ${pop}_${batch}_btqc_mgqc-chrx


# Plot F-statistic
Rscript $scrdir/05_pop_check_sex.R $pop $pop_preqcdir $batch
