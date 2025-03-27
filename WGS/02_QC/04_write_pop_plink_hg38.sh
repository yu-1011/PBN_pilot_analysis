#!/bin/bash
#$ -cwd
#$ -pe smp 4 -binding linear:4
#$ -l h_rt=12:00:00,h_vmem=25G,

datadir=$1
preqcdir=${datadir}/preimp_qc
pcadir=${datadir}/pca
pop=$2
pop_preqcdir=${datadir}/pop_qc
pred_prob=$3
batch=$4

PLINK="/home/unix/yu/software/00install_softwares/bin/plink"

if [ ! -d $pop_preqcdir ]; then
	mkdir $pop_preqcdir
fi

cd $pop_preqcdir

POP=$(echo $pop | tr '[a-z]' '[A-Z]')
awk -v x=$POP -v y=$pred_prob '$25==x && $26>y' $pcadir/${batch}_ref.PC.predPop.tsv | awk '{print $2"\t"$1}' > $pcadir/${batch}_ref.PC.predPop${pred_prob}.${pop}.indlist
awk '{print $2}' $preqcdir/${batch}_geno05_mind02_geno02.bim > ${batch}_geno05_mind02_geno02_rsid.snp

$PLINK \
--bfile $preqcdir/${batch}_geno05_mind02_geno02_rsid \
--make-bed \
--extract ${batch}_geno05_mind02_geno02_rsid.snp \
--out ${pop}_${batch}_btqc_mgqc-rsid

$PLINK \
--bfile ${pop}_${batch}_btqc_mgqc-rsid \
--freq \
--out ${pop}_${batch}_btqc_mgqc-rsid-af

#######
$PLINK \
--bfile $preqcdir/${batch}_geno05_mind02_geno02 \
--make-bed \
--out ${pop}_${batch}_btqc_mgqc
