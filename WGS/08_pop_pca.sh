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
pop_pcadir=${datadir}/pop_pca
scrdir="/stanley/huang_lab/shared/data/sc-asia/script"
ref="/stanley/huang_lab/home/ychen/resources/1kG/ALL.1KG_phase3.20130502.genotypes.maf005"
fhighld_region="/stanley/huang_lab/home/ychen/resources/long_range_LD_intervals.txt"
batch=$3

PLINK="/home/unix/yu/software/00install_softwares/bin/plink"

if [ ! -d $pop_pcadir ]; then
	mkdir $pop_pcadir
fi

cd $pop_pcadir

if [ ! -f ${pop}_pbk_sexcheck_het_ibd.remove.indlist ]; then
    cat ${pop_preqcdir}/${pop}_${batch}_sex_mismatch_F025_M075.indlist ${pop_preqcdir}/${pop}_${batch}_het_outlier_5sd.indlist ${pop_preqcdir}/${pop}_${batch}_ibd_pihat02.indlist | sort | uniq > ${pop}_${batch}_sexcheck_het_ibd.remove.indlist
	cat ${pop}_${batch}_sexcheck_het_ibd.remove.indlist > $pop_preqcdir/${pop}_${batch}_sexcheck_het_ibd.remove.indlist
fi

#############################
#  [ PCA on study sample ]  #
#############################

# Find strand ambiguous SNPs
python $scrdir/find_atgc_snps.py $pop_preqcdir/${pop}_${batch}_btqc_mgqc-rsid.bim > ${pop}_${batch}_btqc_mgqc-rsid.atgc.snplist

cat ${pop}_${batch}_btqc_mgqc-rsid.atgc.snplist| sort | uniq > ${pop}_${batch}_btqc_mgqc-rsid.remove.snplist

awk 'NR==FNR{a[$1];next} !($2 in a) {print $2}' ${pop}_${batch}_btqc_mgqc-rsid.remove.snplist $pop_preqcdir/${pop}_${batch}_btqc_mgqc-rsid.bim > ${pop}_${batch}_btqc_mgqc-rsid.keep.snplist

# Perform LD pruning
$PLINK \
--bfile $pop_preqcdir/${pop}_${batch}_btqc_mgqc-rsid \
--autosome \
--geno 0.02 \
--maf 0.05 \
--snps-only just-acgt \
--remove ${pop}_${batch}_sexcheck_het_ibd.remove.indlist \
--extract ${pop}_${batch}_btqc_mgqc-rsid.keep.snplist \
--exclude range $fhighld_region \
--indep-pairwise 200 100 0.1 \
--out ${pop}_${batch}_unrel-ldpr

# Run pca
$PLINK \
--bfile $pop_preqcdir/${pop}_${batch}_btqc_mgqc-rsid \
--remove ${pop}_${batch}_sexcheck_het_ibd.remove.indlist \
--extract ${pop}_${batch}_unrel-ldpr.prune.in \
--pca 20 header tabs \
--out ${pop}_${batch}_unrel.pca


rm ${pop}_${batch}_btqc_mgqc-rsid.remove.snplist
rm ${pop}_${batch}_btqc_mgqc-rsid.keep.snplist




###############################
# [ PCA with reference data ] #
###############################

# Subset to pop-specific samples in the ref panel
if [ $pop == "EUR" ]; then 
    awk '$1=="CEU" || $1=="TSI" || $1=="FIN" || $1=="GBR" || $1=="IBS"' $ref.fam | awk '{print $1"\t"$2}' > $ref.${pop}.fam
elif [ $pop == "AFR" ]; then 
    awk '$1=="YRI" || $1=="LWK" || $1=="GWD" || $1=="MSL" || $1=="ESN" || $1=="ASW" || $1=="ACB"' $ref.fam | awk '{print $1"\t"$2}' > $ref.${pop}.fam 
elif [ $pop == "EAS" ]; then 
    awk '$1=="CHB" || $1=="JPT" || $1=="CHS" || $1=="CDX" || $1=="KHV"' $ref.fam | awk '{print $1"\t"$2}' > $ref.${pop}.fam
elif [ $pop == "AMR" ]; then
    awk '$1=="MXL" || $1=="PUR" || $1=="CLM" || $1=="PEL"' $ref.fam | awk '{print $1"\t"$2}' > $ref.${pop}.fam
elif [ $pop == "SAS" ]; then
    awk '$1=="GIH" || $1=="PJL" || $1=="BEB" || $1=="STU" || $1=="ITU"' $ref.fam | awk '{print $1"\t"$2}' > $ref.${pop}.fam
fi


# Find SNPs in common between study sample and ref sample
awk 'NR==FNR{a[$0];next} ($0 in a)' $pop_preqcdir/${pop}_${batch}_btqc_mgqc-rsid.bim $ref.bim > ${pop}_${batch}_ref.common.snplist

# Merge with specific subpopulation Retain only overlapping SNPs (do this separately for each file; necessary step before merging b/c bmerge does not retain only overlapping SNPs)
$PLINK \
--bfile $pop_preqcdir/${pop}_${batch}_btqc_mgqc-rsid \
--extract ${pop}_${batch}_ref.common.snplist \
--remove ${pop}_${batch}_sexcheck_het_ibd.remove.indlist \
--make-bed \
--out ${batch}-tmp

$PLINK \
--bfile $ref \
--keep $ref.${pop}.fam \
--extract ${pop}_${batch}_ref.common.snplist \
--make-bed \
--out ref-tmp

$PLINK \
--bfile $ref \
--extract ${pop}_${batch}_ref.common.snplist \
--make-bed \
--out all-1kg-ref-tmp

# Merge study sample with ref panel
$PLINK --bfile ${batch}-tmp \
--keep-allele-order \
--allow-no-sex \
--bmerge ref-tmp \
--make-bed \
--out ${pop}_${batch}_ref

$PLINK --bfile ${batch}-tmp \
--keep-allele-order \
--bmerge all-1kg-ref-tmp \
--allow-no-sex \
--make-bed \
--out ${pop}_${batch}_all_1kg_ref

rm *tmp*

# If there are strand-flipping or multi-allelic SNPs...exclude them and repeat the steps (shouldn't be too many)
if [ -f ${pop}_${batch}_ref-merge.missnp ]; then
    
    mv ${pop}_${batch}_ref.log ${pop}_${batch}_ref-merge.log

    $PLINK \
    --bfile $pop_preqcdir/${pop}_${batch}_btqc_mgqc-rsid \
    --remove ${pop}_${batch}_sexcheck_het_ibd.remove.indlist \
    --extract ${pop}_${batch}_ref.common.snplist \
    --exclude ${pop}_${batch}_ref-merge.missnp \
    --make-bed \
    --out ${batch}-tmp

    $PLINK \
    --bfile $ref \
    --keep $ref.${pop}.fam \
    --extract ${pop}_${batch}_ref.common.snplist \
    --exclude ${pop}_${batch}_ref-merge.missnp \
    --make-bed \
    --out ref-tmp

    $PLINK \
    --bfile ${batch}-tmp \
    --keep-allele-order \
    --bmerge ref-tmp \
    --make-bed \
    --allow-no-sex \
    --out ${pop}_${batch}_ref

    rm *tmp*
fi

if [ -f ${pop}_${batch}_all_1kg_ref-merge.missnp ]; then
    
    mv ${pop}_${batch}_all_1kg_ref.log ${pop}_${batch}_all_1kg_ref-merge.log

    $PLINK \
    --bfile $pop_preqcdir/${pop}_${batch}_btqc_mgqc-rsid \
    --remove ${pop}_${batch}_sexcheck_het_ibd.remove.indlist \
    --extract ${pop}_${batch}_ref.common.snplist \
    --exclude ${pop}_${batch}_all_1kg_ref-merge.missnp \
    --make-bed \
    --out ${batch}_all_1kg_ref-tmp

    $PLINK \
    --bfile $ref \
    --extract ${pop}_${batch}_ref.common.snplist \
    --exclude ${pop}_${batch}_ref-merge.missnp \
    --make-bed \
    --out all-1kg-ref-tmp

	$PLINK \
	--bfile ${batch}_all_1kg_ref-tmp \
	--keep-allele-order \
	--bmerge all-1kg-ref-tmp \
	--allow-no-sex \
	--make-bed \
	--out ${pop}_${batch}_all_1kg_ref

    rm *tmp*
fi

# Find strand ambiguous SNPs
python $scrdir/find_atgc_snps.py ${pop}_${batch}_ref.bim > ${pop}_${batch}_ref.atgc.snplist

# Write a list of non-strand ambiguous SNPs to keep
awk 'NR==FNR{a[$1];next} !($2 in a) {print $2}' ${pop}_${batch}_ref.atgc.snplist ${pop}_${batch}_ref.bim > ${pop}_${batch}_ref.non-atgc.snplist


# Perform LD pruning
$PLINK \
--bfile ${pop}_${batch}_ref \
--autosome \
--geno 0.02 \
--maf 0.05 \
--snps-only just-acgt \
--extract ${pop}_${batch}_ref.non-atgc.snplist \
--exclude range $fhighld_region \
--indep-pairwise 200 100 0.1 \
--out ${pop}_${batch}_ref-ldpr

$PLINK \
--bfile ${pop}_${batch}_all_1kg_ref \
--autosome \
--geno 0.02 \
--maf 0.05 \
--snps-only just-acgt \
--extract ${pop}_${batch}_ref.non-atgc.snplist \
--exclude range $fhighld_region \
--indep-pairwise 200 100 0.1 \
--out ${pop}_${batch}_all_1kg_ref-ldpr

# Run pca
$PLINK \
--bfile ${pop}_${batch}_ref \
--extract ${pop}_${batch}_ref-ldpr.prune.in \
--pca 20 header tabs \
--out ${pop}_${batch}_ref.pca

$PLINK \
--bfile ${pop}_${batch}_all_1kg_ref \
--extract ${pop}_${batch}_all_1kg_ref-ldpr.prune.in \
--pca 20 header tabs \
--out ${pop}_${batch}_all_1kg_ref.pca

# Plot PCs
Rscript $scrdir/08_pop_pca.R $pop $pop_pcadir $batch