cores=4
memory_type="standard"
threads=4
sample_size=200
case_or_control="case"
terra_path="fc-3eddf4b4-1ae3-4d30-9235-4a20971bac43"
mount_folder="/home/kimsoyeo/mount_folder/"

python3 phase_impute.glimpse2_hb.py $cores $memory_type $threads $sample_size \
	$case_or_control $terra_path $mount_folder
