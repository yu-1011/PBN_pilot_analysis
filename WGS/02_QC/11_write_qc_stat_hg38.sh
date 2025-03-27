#!/bin/bash
#$ -cwd
#$ -pe smp 1 -binding linear:1
#$ -l h_rt=2:00:00,h_vmem=5G,
# Define variable for job, to be stored in ${SGE_TASK_ID}

datadir=$1
pop_preqcdir=${datadir}/pop_qc
pop_pcadir=${datadir}/pop_pca
preqcdir=${datadir}/preimp_qc
pcadir=${datadir}/pca
sexqcdir=${datadir}/chrx_qc
finalqc=${datadir}/final_qc
pop=$2
batch=$3
ori_bfile_name=$4

if [ ! -d $finalqc/qc_stat ]; then
	mkdir $finalqc/qc_stat
fi

cd $finalqc/qc_stat

# Summarize number of samples and SNPs removed at each step
touch qc_summary.tsv
ninds=$(cat $datadir/$ori_bfile_name.fam | wc -l)

ninds_m=$(awk '$5==1{print}' $datadir/$ori_bfile_name.fam | wc -l)
ninds_f=$(awk '$5==2{print}' $datadir/$ori_bfile_name.fam | wc -l)

ninds_con=$(awk '$6==1{print}' $datadir/$ori_bfile_name.fam | wc -l)
ninds_case=$(awk '$6==2{print}' $datadir/$ori_bfile_name.fam | wc -l)

nsnps=$(cat $datadir/$ori_bfile_name.bim | wc -l)

echo "cohorts_id "$batch > qc_summary.tsv
echo "population "$pop >> qc_summary.tsv

echo "initial_Ninds "$ninds >> qc_summary.tsv
echo "initial_Ninds_Male "$ninds_m >> qc_summary.tsv
echo "initial_Ninds_Female "$ninds_f >> qc_summary.tsv

echo "initial_Ninds_Case "$ninds_case >> qc_summary.tsv
echo "initial_Ninds_Controls "$ninds_con >> qc_summary.tsv

echo "initial_Nsnps "$nsnps >> qc_summary.tsv

echo "Nsnps_with_call_rate_>0.95 "$(cat $preqcdir/${batch}_geno05.snplist | wc -l) >> qc_summary.tsv
echo "Nsnps_with_call_rate_>0.98 "$(cat $preqcdir/${batch}_geno05_mind02_geno02.snplist | wc -l) >> qc_summary.tsv
echo "Nsnps_filtered_by_HWE_check "$(cat $pop_preqcdir/${pop}_${batch}_unrel_hwe_pe-10.snplist | wc -l) >> qc_summary.tsv
echo "Nsnps_filtered_by_sex_check "$(cat $pop_preqcdir/${pop}_${batch}_sex_diff_snps.diff | wc -l) >> qc_summary.tsv
echo "Nsnps_non-pseudo-autosomal_regions_non-indel_polymorphic_SNPs "$(cat $pop_preqcdir/${pop}_${batch}_unrel_qc_nopar.bim | wc -l) >> qc_summary.tsv
echo "Nsnps_in_HRC "$(cat $finalqc/${pop}_${batch}_unrel_qc_final.bim | wc -l) >> qc_summary.tsv

echo "Ninds_with_call_rate_>0.98 "$(cat $preqcdir/${batch}_geno05_mind02.indlist | wc -l) >> qc_summary.tsv
echo "Ninds_in_pop_pred "$(cat $pcadir/${batch}_ref.PC.predPop0.9.${pop}.indlist | wc -l) >> qc_summary.tsv
echo "Ninds_filtered_by_sex_check "$(cat $pop_preqcdir/${pop}_${batch}_sex_mismatch_F025_M075.indlist | wc -l) >> qc_summary.tsv
echo "Ninds_filtered_by_het_check "$(cat $pop_preqcdir/${pop}_${batch}_het_outlier_5sd.indlist | wc -l) >> qc_summary.tsv
echo "Ninds_filtered_by_IBD_check "$(cat $pop_preqcdir/${pop}_${batch}_ibd_pihat02.indlist  | wc -l) >> qc_summary.tsv

ninds_final=$(cat $finalqc/${pop}_${batch}_unrel_qc_final.fam | wc -l)
nsnps_final=$(cat $finalqc/${pop}_${batch}_unrel_qc_final.bim | wc -l)
ninds_final_m=$(awk '$5==1{print}' $finalqc/${pop}_${batch}_unrel_qc_final.fam | wc -l)
ninds_final_f=$(awk '$5==2{print}' $finalqc/${pop}_${batch}_unrel_qc_final.fam | wc -l)
ninds_final_sex_un=$ninds - $ninds_m - $ninds_f

ninds_final_con=$(awk '$6==1{print}' $finalqc/${pop}_${batch}_unrel_qc_final.fam | wc -l)
ninds_final_case=$(awk '$6==2{print}' $finalqc/${pop}_${batch}_unrel_qc_final.fam | wc -l)
ninds_final_phe_un=$ninds - $ninds_con - $ninds_case

echo "final_Ninds "$ninds_final >> qc_summary.tsv

echo "final_Ninds_Male "$ninds_final_m >> qc_summary.tsv
echo "final_Ninds_Female "$ninds_final_f >> qc_summary.tsv

echo "final_Ninds_Case "$ninds_final_case >> qc_summary.tsv
echo "final_Ninds_Controls "$ninds_final_con >> qc_summary.tsv

echo "final_Nsnps "$nsnps_final >> qc_summary.tsv

mv qc_summary.tsv ${pop}_${batch}_unrel_qc_summary.tsv

cp $pop_pcadir/*.pdf .
cp $pop_preqcdir/*.pdf .
