scr="/stanley/huang_lab/shared/data/sc-asia/script"
datadir="/stanley/huang_lab/home/ychen/proj-PBN/00_raw_data/02QCed_file"
batch="scz_pbn_pilot_yc"  #dis_cohortid_pop_author
ori_bfile_name="PBN_pilot_final" 
pop="EUR"
pcplot0=True
npc=5
pred_prob=0.9

qsub -N $batch"_01_callrate_qc" -e $datadir/${batch}_01_e.file -o $datadir/${batch}_01_o.file $scr/01_callrate_qc_by_batch_hg38.sh $datadir $batch $ori_bfile_name

qsub -hold_jid $batch"_01_callrate_qc" -N $batch"_02_pca" -e $datadir/${batch}_02_e.file -o $datadir/${batch}_02_o.file $scr/02_pca_hg38.sh $datadir $batch

qsub -hold_jid $batch"_02_pca" -N $batch"_03_classify_ancestry" -e $datadir/${batch}_03_e.file -o $datadir/${batch}_03_o.file $scr/03_classify_ancestry.sh ${datadir}"/pca" $pcplot0 $npc $pred_prob $batch 

qsub -hold_jid $batch"_03_classify_ancestry" -N $batch"_04_write_pop_plink" -e $datadir/${batch}_04_e.file -o $datadir/${batch}_04_o.file $scr/04_write_pop_plink_hg38.sh ${datadir} $pop $pred_prob $batch 

qsub -hold_jid $batch"_04_write_pop_plink" -N $batch"_05_pop_check_sex" -e $datadir/${batch}_05_e.file -o $datadir/${batch}_05_o.file $scr/05_pop_check_sex.sh ${datadir} $pop $batch 

qsub -hold_jid $batch"_05_pop_check_sex" -N $batch"_06_pop_check_het" -e $datadir/${batch}_06_e.file -o $datadir/${batch}_06_o.file $scr/06_pop_check_het.sh ${datadir} $pop $batch 

qsub -hold_jid $batch"_06_pop_check_het" -N $batch"_07_pop_filter_ibd" -e $datadir/${batch}_07_e.file -o $datadir/${batch}_07_o.file $scr/07_pop_check_ibd.sh ${datadir} $pop $batch 

qsub -hold_jid $batch"_07_pop_filter_ibd" -N $batch"_08_pop_pca" -e $datadir/${batch}_08_e.file -o $datadir/${batch}_08_o.file $scr/08_pop_pca_hg38.sh ${datadir} $pop $batch 

qsub -hold_jid $batch"_08_pop_pca" -N $batch"_09_chrx_qc" -e $datadir/${batch}_09_e.file -o $datadir/${batch}_09_o.file $scr/09_chr_x_qc.sh ${datadir} $pop $batch 

qsub -hold_jid $batch"_09_chrx_qc" -N $batch"_10_final_qc" -e $datadir/${batch}_10_e.file -o $datadir/${batch}_10_o.file $scr/10_write_qc_plink_hg38.sh ${datadir} $pop $batch 

qsub -hold_jid $batch"_10_final_qc" -N $batch"_11_final_qc_sum_stat" -e $datadir/${batch}_11_e.file -o $scr/log/${batch}_11_o.file $scr/11_write_qc_stat_hg38.sh ${datadir} $pop $batch $ori_bfile_name
