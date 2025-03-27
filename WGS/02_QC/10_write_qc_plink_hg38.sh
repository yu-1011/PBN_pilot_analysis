#!/bin/bash
#$ -cwd
#$ -pe smp 1 -binding linear:1
#$ -l h_rt=2:00:00,h_vmem=25G

datadir=$1
pop_preqcdir=${datadir}/pop_qc
pop_pcadir=${datadir}/pop_pca
preqcdir=${datadir}/preimp_qc
pcadir=${datadir}/pca
sexqcdir=${datadir}/chrx_qc
finalqc=${datadir}/final_qc
pop=$2
batch=$3

PLINK="/home/unix/yu/software/00install_softwares/bin/plink"

cd $pop_preqcdir

if [ ! -d $finalqc ]; then
	mkdir $finalqc
fi

# Rmove samples not passing sex, het, and ibd filter
$PLINK \
--bfile $pop_preqcdir/${pop}_${batch}_btqc_mgqc \
--remove ${pop}_${batch}_sexcheck_het_ibd.remove.indlist \
--make-bed \
--out ${pop}_${batch}_unrel-tmp

# Perform final SNP-level QC on predicted pop-specific samples 
$PLINK \
--bfile ${pop}_${batch}_unrel-tmp \
--geno 0.02 \
--write-snplist \
--out ${pop}_${batch}_unrel_geno02

$PLINK \
--bfile ${pop}_${batch}_unrel-tmp \
--hardy \
--out ${pop}_${batch}_unrel-hardy
# --hwe

awk '$9<1e-10{print $2}' ${pop}_${batch}_unrel-hardy.hwe > ${pop}_${batch}_unrel_hwe_pe-10.snplist

# Remove SNPs with call rate < 0.98 and p-hwe < 1e-10
$PLINK \
--bfile ${pop}_${batch}_unrel-tmp \
--extract ${pop}_${batch}_unrel_geno02.snplist \
--exclude ${pop}_${batch}_unrel_hwe_pe-10.snplist \
--make-bed \
--out ${pop}_${batch}_unrel_qc

rm ${pop}_${batch}_unrel-tmp*

# Save non-pseudo-autosomal regions, non-indel, and polymorphic SNPs
$PLINK \
--bfile ${pop}_${batch}_unrel_qc \
--exclude $sexqcdir/pars.txt \
--range \
--make-bed \
--out ${pop}_${batch}_unrel_qc_nopar

# Remove SNPs with sex-specific MAF & HWE & assoc
$PLINK \
--bfile  ${pop}_${batch}_unrel_qc_nopar \
--exclude ${pop}_${batch}_sex_diff_snps.diff \
--make-bed \
--out ${pop}_${batch}_unrel_qc_final

# Calc freq
$PLINK \
--bfile ${pop}_${batch}_unrel_qc_final \
--freq \
--out ${pop}_${batch}_unrel_qc_final-af
# # when not using --mac 1 (not excluding monomorphic SNPs):
# awk '$5==0' ${pop}_${batch}_unrel_qc_aut-af.frq | wc -l #-> still ~80K monomorphic sites; to remove

mv ${pop}_${batch}_unrel_qc_final* $finalqc/
