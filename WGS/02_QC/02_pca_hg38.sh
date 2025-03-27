#!/bin/bash
#$ -cwd
#$ -pe smp 4 -binding linear:4
#$ -l h_rt=12:00:00,h_vmem=25G,

datadir=$1
preqcdir=${datadir}/preimp_qc
batch=$2
pcadir=${datadir}/pca
ref="/stanley/huang_lab/home/ychen/resources/1KG_hg38/kgpAll.maf5"
fhighld_region="/stanley/huang_lab/home/ychen/resources/long_range_LD_intervals.txt"
scrdir="/stanley/huang_lab/shared/data/sc-asia/script"

# ref=/data/js95/yfeng/utility/1kG/ALL.1KG_phase3.20130502.genotypes.maf005

PLINK="/home/unix/yu/software/00install_softwares/bin/plink"

if [ ! -d $pcadir ]; then
	mkdir $pcadir
fi

cd $pcadir

###############################
# [ PCA with reference data ] #
###############################

# Find SNPs in common between study sample and ref sample
awk 'NR==FNR{a[$2];next} ($2 in a)' $preqcdir/${batch}_geno05_mind02_geno02_rsid.bim $ref.bim > $pcadir/${batch}_ref.common.snplist

#awk 'NR==FNR {a[$2] = $0; next} {key = $2; if (key in a) print a[key]; else print $0}' $ref.bim $preqcdir/${batch}_geno05_mind02_geno02.bim > $preqcdir/${batch}_geno05_mind02_geno02_update.bim
#mv $preqcdir/${batch}_geno05_mind02_geno02_update.bim $preqcdir/${batch}_geno05_mind02_geno02.bim

# Retain only overlapping SNPs (do this separately for each file; necessary step before merging b/c bmerge does not retain only overlapping SNPs)
$PLINK \
--bfile $preqcdir/${batch}_geno05_mind02_geno02_rsid \
--extract ${batch}_ref.common.snplist \
--allow-no-sex \
--make-bed \
--out ${batch}-tmp

$PLINK \
--bfile $ref \
--extract ${batch}_ref.common.snplist \
--make-bed \
--out ref-tmp

# Merge study sample with ref panel
$PLINK --bfile ${batch}-tmp \
--keep-allele-order \
--bmerge ref-tmp \
--allow-no-sex \
--make-bed \
--out ${batch}_ref

rm *tmp*

# If there are strand-flipping or multi-allelic SNPs...exclude them and repeat the steps (shouldn't be too many)
if [ -f ${batch}_ref-merge.missnp ]; then
    
    mv ${batch}_ref.log ${batch}_ref-merge.log

    $PLINK \
    --bfile $preqcdir/${batch}_geno05_mind02_geno02_rsid \
    --extract ${batch}_ref.common.snplist \
    --exclude ${batch}_ref-merge.missnp \
    --allow-no-sex \
    --make-bed \
    --out ${batch}-tmp

    $PLINK \
    --bfile $ref \
    --extract ${batch}_ref.common.snplist \
    --exclude ${batch}_ref-merge.missnp \
    --make-bed \
    --out ref-tmp

    $PLINK \
    --bfile ${batch}-tmp \
    --keep-allele-order \
    --bmerge ref-tmp \
    --allow-no-sex \
    --make-bed \
    --out ${batch}_ref

    rm *tmp*
fi


# Find strand ambiguous SNPs
python $scrdir/find_atgc_snps.py ${batch}_ref.bim > ${batch}_ref.atgc.snplist

# Write a list of non-strand ambiguous SNPs to keep
awk 'NR==FNR{a[$1];next} !($2 in a) {print $2}' ${batch}_ref.atgc.snplist ${batch}_ref.bim > ${batch}_ref.non-atgc.snplist


# Perform LD pruning
$PLINK \
--bfile ${batch}_ref \
--autosome \
--geno 0.02 \
--maf 0.05 \
--snps-only just-acgt \
--extract ${batch}_ref.non-atgc.snplist \
--exclude range $fhighld_region \
--indep-pairwise 200 100 0.1 \
--out ${batch}_ref-ldpr


# Run pca
$PLINK \
--bfile ${batch}_ref \
--extract ${batch}_ref-ldpr.prune.in \
--pca 20 header tabs \
--out ${batch}_ref.pca



#############################
#  [ PCA on study sample ]  #
#############################

# Find strand ambiguous SNPs
python $scrdir/find_atgc_snps.py $preqcdir/${batch}_geno05_mind02_geno02_rsid.bim > ${batch}_btqc-rsid.atgc.snplist

# Write a list of non-strand ambiguous SNPs to keep
awk 'NR==FNR{a[$1];next} !($2 in a) {print $2}' ${batch}_btqc-rsid.atgc.snplist $preqcdir/${batch}_geno05_mind02_geno02_rsid.bim > ${batch}_btqc-rsid.non-atgc.snplist

# Perform LD pruning
$PLINK \
--bfile $preqcdir/${batch}_geno05_mind02_geno02_rsid \
--autosome \
--geno 0.02 \
--maf 0.05 \
--snps-only just-acgt \
--extract ${batch}_btqc-rsid.non-atgc.snplist \
--exclude range $fhighld_region \
--indep-pairwise 200 100 0.1 \
--make-bed \
--out ${batch}_btqc-rsid

# Run pca
$PLINK \
--bfile ${batch}_btqc-rsid \
--pca 20 header tabs \
--out ${batch}_btqc-rsid.pca
