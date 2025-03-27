# /bin/bash
#$ -cwd
#$ -pe smp 4 -binding linear:4
#$ -l h_rt=12:00:00,h_vmem=25G,


datadir=$1
batch=$2
ori_bfile_name=$3 		# e.g., MEGA
ref="/stanley/huang_lab/home/ychen/resources/1KG_hg38/kgpAll.maf5"
preqcdir=${datadir}/preimp_qc


######################
# test:
# datadir=/stanley/huang_lab/shared/data/sc-asia/Release/Hong/kor1/Genotype_QCed/qc/imputation
# batch=scz_kor1_asn_yc-qc1

PLINK="/home/unix/yu/software/00install_softwares/bin/plink"

if [ ! -d $datadir/preimp_qc ]; then
	mkdir $datadir/preimp_qc
fi

cd $datadir/preimp_qc
######################


# 1. Output a list of SNPs with call rate > 0.95
$PLINK \
--bfile $datadir/$ori_bfile_name \
--geno 0.05 \
--write-snplist \
--out ${batch}_geno05


# 2. Filter to SNPs in the previous step, output a list of samples with call rate > 0.98
$PLINK \
--bfile $datadir/$ori_bfile_name \
--extract ${batch}_geno05.snplist \
--missing \
--out ${batch}_geno05
awk 'NR>1{if($6<0.02){print $1,$2}}' ${batch}_geno05.imiss > ${batch}_geno05_mind02.indlist


# 3. Filter to samples in the previous step, output a list of SNPs with call rate > 0.98
$PLINK \
--bfile $datadir/$ori_bfile_name \
--extract ${batch}_geno05.snplist \
--keep ${batch}_geno05_mind02.indlist \
--geno 0.02 \
--write-snplist \
--out ${batch}_geno05_mind02_geno02


# 4. Keep samples and SNPs passing previous filters, calculate SNP-level missing rate
$PLINK \
--bfile $datadir/$ori_bfile_name \
--extract ${batch}_geno05_mind02_geno02.snplist \
--keep ${batch}_geno05_mind02.indlist \
--missing \
--out ${batch}_geno05_mind02_geno02

# 5. make bed by batch
$PLINK \
--bfile $datadir/$ori_bfile_name \
--extract ${batch}_geno05_mind02_geno02.snplist \
--keep ${batch}_geno05_mind02.indlist \
--make-bed \
--out ${batch}_geno05_mind02_geno02

$PLINK \
--bfile $datadir/$ori_bfile_name \
--extract ${batch}_geno05_mind02_geno02.snplist \
--keep ${batch}_geno05_mind02.indlist \
--make-bed \
--out ${batch}_geno05_mind02_geno02_tmp

awk 'NR==FNR {a[$1""$4""$5""$6] = $0; next} {key = $1""$4""$5""$6; if (key in a) print a[key]; else print $0}' $ref.bim $preqcdir/${batch}_geno05_mind02_geno02_tmp.bim > $preqcdir/${batch}_geno05_mind02_geno02_update.bim
mv $preqcdir/${batch}_geno05_mind02_geno02_update.bim $preqcdir/${batch}_geno05_mind02_geno02_tmp.bim

$PLINK \
--bfile ${batch}_geno05_mind02_geno02_tmp \
--make-bed \
--out ${batch}_geno05_mind02_geno02_rsid

rm *tmp*

######
# Summarize number of samples and SNPs removed at each step
# touch qc_summary.tsv
# ninds=$(cat $datadir/$batch.fam | wc -l)
# nsnps=$(cat $datadir/$batch.bim | wc -l)
# echo "Ninds "$ninds > qc_summary.tsv
# echo "Nsnps "$nsnps >> qc_summary.tsv
# 
# echo "Nsnps with call rate >0.95 "$(cat ${batch}_geno05.snplist | wc -l) >> qc_summary.tsv
# echo "Ninds with call rate >0.98 "$(cat ${batch}_geno05_mind02.indlist | wc -l) >> qc_summary.tsv
# echo "Nsnps with call rate >0.98 "$(cat ${batch}_geno05_mind02_geno02.snplist | wc -l) >> qc_summary.tsv
# 
