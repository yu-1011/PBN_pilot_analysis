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
scrdir="/stanley/huang_lab/shared/data/sc-asia/script"
batch=$3
fhighld_region="/stanley/huang_lab/home/ychen/resources/long_range_LD_intervals.txt"

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
--out ${pop}_${batch}_btqc_mgqc-ldpr-aut-tmp


# if [ !-f ${pop}_pbk_btqc_mgqc-ibd-min0.125.genome.gz ]; then
# Estimate IBD
$PLINK \
--bfile ${pop}_${batch}_btqc_mgqc \
--extract ${pop}_${batch}_btqc_mgqc-ldpr-aut-tmp.prune.in \
--genome \
--min 0.125 \
--out ${pop}_${batch}_btqc_mgqc-ibd-min0.125

gzip -f ${pop}_${batch}_btqc_mgqc-ibd-min0.125.genome

# Plot IBD to infer relatedness and extract one from each pair of related's (pihat>0.2) for removal
Rscript $scrdir/07_pop_check_ibd.R $pop $pop_preqcdir $batch


rm ${pop}_${batch}_btqc_mgqc-ldpr-aut-tmp*
