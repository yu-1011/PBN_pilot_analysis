#!/bin/bash
#$ -cwd
#$ -pe smp 4 -binding linear:4
#$ -l h_rt=12:00:00,h_vmem=25G

source /broad/software/scripts/useuse

reuse Python-3.9

datadir=$1
pop=$2
batch=$3
pop_preqcdir=${datadir}/pop_qc
pop_pcadir=${datadir}/pop_pca
preqcdir=${datadir}/preimp_qc
pcadir=${datadir}/pca
sexqcdir=${datadir}/chrx_qc

alpha=0.05
mindthresh=0.02
minmaf=0.05
genothresh=0.02

scrdir="/stanley/huang_lab/shared/data/sc-asia/script"

PLINK="/home/unix/yu/software/00install_softwares/bin/plink"
xwasloc="/stanley/huang_lab/home/ychen/software/xwas-3.0/bin"
xwasQ="xwas"
PYTHON="/broad/software/free/Linux/redhat_7_x86_64/pkgs/anaconda3_2022.10/bin/python3.9"


if [ ! -d $sexqcdir ]; then
	mkdir $sexqcdir
fi

cd $sexqcdir

# Set PAR locations based on hg 19 build
echo  '23 60001 2699520 par1\n' > pars.txt
echo  '23 154931044 155260560 par2' >> pars.txt

# remove pseudo-autosomal regions
echo -e "Removing pseudo-autosomal regions"
$PLINK \
--bfile $pop_preqcdir/${pop}_${batch}_btqc_mgqc-rsid \
--remove $pop_preqcdir/${pop}_${batch}_sexcheck_het_ibd.remove.indlist \
--make-bed \
--out ${pop}_${batch}_unrel-ldpr_no_par \
--exclude pars.txt \
--range

totsnps=$(wc -l ${pop}_${batch}_unrel-ldpr_no_par.bim | awk '{print $1}')
bonf=$(echo "scale=20;$alpha/$totsnps" | bc)

# QC separately for male and female
malenum=`awk '$5 == 1' ${pop}_${batch}_unrel-ldpr_no_par.fam | wc -l`;
femalenum=`awk '$5 == 2' ${pop}_${batch}_unrel-ldpr_no_par.fam | wc -l`;

if [ $malenum -eq 0 ] || [ $femalenum -eq 0 ]; then
    echo "Only one gender exist. Male # is $malenum; Female # is $femalenum"
    echo "QC is only performed in one gender"
    echo -e "Quality Control for Binary Traits"

	echo "HWE"
	$PLINK \
	--bfile ${pop}_${batch}_unrel-ldpr_no_par \
	--hardy \
	--out ${pop}_${batch}_unrel-ldpr_no_par_1gender_hwe \
	--filter-controls
	cat ${pop}_${batch}_unrel-ldpr_no_par_1gender_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${fname}_1gender_snp.exclude

	echo "Correlation between missingness and phenotype"
	$PLINK \
	--bfile ${pop}_${batch}_unrel-ldpr_no_par \
	--test-missing \
	--out ${pop}_${batch}_unrel-ldpr_no_par_1gender_mcc
	awk -v bf=$bonf '$5<bf {print $2}' ${pop}_${batch}_unrel-ldpr_no_par_1gender_mcc.missing >> ${pop}_${batch}_unrel-ldpr_no_par_1gender_snp.exclude
	
	echo "MAF, missingness per SNP, missingness per individual"
	$PLINK \
	--bfile ${pop}_${batch}_unrel-ldpr_no_par \
	--make-bed \
	--out ${pop}_${batch}_unrel-ldpr_no_par_1gender_ft \
	--mind ${mindthresh} \
	--maf ${minmaf} \
	--geno ${genothresh} \
	--exclude ${pop}_${batch}_unrel-ldpr_no_par_1gender_snp.exclude

else

	echo "QC is performed separately for each gender"
	echo -e "Quality Control for Binary Traits"

	echo "Quality Control for Male"
	$PLINK \
	--bfile ${pop}_${batch}_unrel-ldpr_no_par \
	--filter-males \
	--make-bed \
	--out ${pop}_${batch}_unrel-ldpr_no_par_male

	echo "HWE"
	$PLINK \
	--bfile ${pop}_${batch}_unrel-ldpr_no_par_male \
	--hardy \
	--out ${pop}_${batch}_unrel-ldpr_no_par_male_hwe \
	--filter-controls
	cat ${pop}_${batch}_unrel-ldpr_no_par_male_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${pop}_${batch}_unrel-ldpr_no_par_male_snp.exclude

	echo "Correlation between missingness and phenotype"
	$PLINK \
	--bfile ${pop}_${batch}_unrel-ldpr_no_par_male \
	--test-missing \
	--out ${pop}_${batch}_unrel-ldpr_no_par_male_mcc
	awk -v bf=$bonf '$5<bf {print $2}' ${pop}_${batch}_unrel-ldpr_no_par_male_mcc.missing >> ${pop}_${batch}_unrel-ldpr_no_par_male_snp.exclude

	echo "MAF, missingness per SNP, missingness per individual"
	$PLINK \
	--bfile ${pop}_${batch}_unrel-ldpr_no_par_male \
	--make-bed \
	--out ${pop}_${batch}_unrel-ldpr_no_par_male_ft \
	--mind ${mindthresh} \
	--maf ${minmaf} \
	--geno ${genothresh} \
	--exclude ${pop}_${batch}_unrel-ldpr_no_par_male_snp.exclude

	echo "Quality Control for Female"
	$PLINK \
	--bfile ${pop}_${batch}_unrel-ldpr_no_par \
	--filter-females \
	--make-bed \
	--out ${pop}_${batch}_unrel-ldpr_no_par_female

	echo "HWE"
	$PLINK \
	--bfile ${pop}_${batch}_unrel-ldpr_no_par_female \
	--hardy \
	--out ${pop}_${batch}_unrel-ldpr_no_par_female_hwe \
	--filter-controls
	cat ${pop}_${batch}_unrel-ldpr_no_par_female_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${pop}_${batch}_unrel-ldpr_no_par_female_snp.exclude

	echo "Correlation between missingness and phenotype"
	$PLINK \
	--bfile ${pop}_${batch}_unrel-ldpr_no_par_female \
	--test-missing \
	--out ${pop}_${batch}_unrel-ldpr_no_par_female_mcc
	awk -v bf=$bonf '$5<bf {print $2}' ${pop}_${batch}_unrel-ldpr_no_par_female_mcc.missing >> ${pop}_${batch}_unrel-ldpr_no_par_female_snp.exclude

	echo "MAF, missingness per SNP, missingness per individual"
	$PLINK \
	--bfile ${pop}_${batch}_unrel-ldpr_no_par_female \
	--make-bed \
	--out ${pop}_${batch}_unrel-ldpr_no_par_female_ft \
	--mind ${mindthresh} --maf ${minmaf} --geno ${genothresh} \
	--exclude ${pop}_${batch}_unrel-ldpr_no_par_female_snp.exclude

	echo -e "Merging male and female QC"
	$PYTHON $scrdir/sex_snp_compare.py ${pop}_${batch}_unrel-ldpr_no_par_female_ft.bim  ${pop}_${batch}_unrel-ldpr_no_par_male_ft.bim  ${pop}_${batch}_unrel-ldpr_no_par_sex_diff_snps.diff
	#cat ${pop}_${batch}_unrel-ldpr_no_par_sex_diff_snps.diff > $pop_preqcdir/${pop}_${batch}_unrel-ldpr_no_par_sex_diff_snps.diff
	
	$PLINK \
	--bfile ${pop}_${batch}_unrel-ldpr_no_par_female_ft \
	--exclude ${pop}_${batch}_unrel-ldpr_no_par_sex_diff_snps.diff \
	--make-bed \
	--out ${pop}_${batch}_unrel-ldpr_no_par_female_ft_tmp
	
	$PLINK \
	--bfile ${pop}_${batch}_unrel-ldpr_no_par_male_ft \
	--exclude ${pop}_${batch}_unrel-ldpr_no_par_sex_diff_snps.diff \
	--make-bed \
	--out ${pop}_${batch}_unrel-ldpr_no_par_male_ft_tmp
	
	$PLINK \
	--bfile ${pop}_${batch}_unrel-ldpr_no_par_female_ft_tmp \
	--bmerge ${pop}_${batch}_unrel-ldpr_no_par_male_ft_tmp.bed ${pop}_${batch}_unrel-ldpr_no_par_male_ft_tmp.bim ${pop}_${batch}_unrel-ldpr_no_par_male_ft_tmp.fam \
	--make-bed \
	--out ${pop}_${batch}_unrel-ldpr_no_par_sexqc_merge
fi

# Significant diff in MAF between males and females only for binary traits and only in controls
echo -e "Significant difference in MAF between males and females (only for binary traits and only in controls)"
echo "Extracting X chromosome"
totx=$(awk '$1=="X" || $1==23 {print $0}' ${pop}_${batch}_unrel-ldpr_no_par_sexqc_merge.bim | wc -l)
bonfx=$(echo "scale=20;$alpha/$totx" | bc)
${xwasloc}/$xwasQ \
--bfile ${pop}_${batch}_unrel-ldpr_no_par_sexqc_merge \
--xwas \
--noweb \
--make-bed \
--out ${pop}_${batch}_unrel-ldpr_no_par_sexqc_final_qcx  \
--chr X \
--freqdiff-x ${bonfx}

awk '$NF<0.05{print $2}' ${pop}_${batch}_unrel-ldpr_no_par_sexqc_final_qcx.xtest| cat - ${pop}_${batch}_unrel-ldpr_no_par_sex_diff_snps.diff |awk '!a[$1]++{print}' > $pop_preqcdir/${pop}_${batch}_sex_diff_snps.diff

rm *tmp*
