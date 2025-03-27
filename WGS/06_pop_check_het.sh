#!/bin/bash
#$ -cwd
#$ -pe smp 4 -binding linear:4
#$ -l h_rt=12:00:00,h_vmem=25G

source /broad/software/scripts/useuse
reuse -q R-4.3

datadir=$1
preqcdir=${datadir}/preimp_qc
pcadir=${datadir}/pca
pop=$2
pop_preqcdir=${datadir}/pop_qc
fhighld_region="/stanley/huang_lab/home/ychen/resources/long_range_LD_intervals.txt"
scrdir="/stanley/huang_lab/shared/data/sc-asia/script"
batch=$3

PLINK="/home/unix/yu/software/00install_softwares/bin/plink"
cd $pop_preqcdir

# Perform LD pruning: autosomal SNPs
$PLINK \
--bfile ${pop}_${batch}_btqc_mgqc \
--autosome \
--geno 0.02 \
--maf 0.05 \
--snps-only just-acgt \
--extract $pcadir/${batch}_btqc-rsid.non-atgc.snplist \
--exclude range $fhighld_region \
--indep-pairwise 200 100 0.1 \
--out ${pop}_${batch}_btqc_mgqc-ldpr-aut


# Estimate autosomal heterozygosity rate/inbreeding coeff
# also use pruned SNPs! see: http://zzz.bwh.harvard.edu/plink/ibdibs.shtml
$PLINK \
--bfile ${pop}_${batch}_btqc_mgqc \
--extract ${pop}_${batch}_btqc_mgqc-ldpr-aut.prune.in \
--het \
--out ${pop}_${batch}_btqc_mgqc-inbr

# Calculate hetorozygosity rate (i.e., the proportion of heterozygous genotypes for a given individual.)
awk 'NR>1{print ($5-$3)/$5}' ${pop}_${batch}_btqc_mgqc-inbr.het | sed '1i HetRate' | paste ${pop}_${batch}_btqc_mgqc-inbr.het - > ${pop}_${batch}_btqc_mgqc-inbr.hetrate

# Plot het rate distribution
Rscript $scrdir/06_pop_check_het.R $pop $pop_preqcdir $batch


