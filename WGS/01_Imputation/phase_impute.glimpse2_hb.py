import hailtop.fs as hfs
import hailtop.batch as hb
import pandas as pd
import re
import sys

cores = sys.argv[1]
memory_type = sys.argv[2]
threads = sys.argv[3]
sample_size = sys.argv[4]
case_or_control = sys.argv[5]
terra_path = sys.argv[6]
mount_folder = sys.argv[7]

backend = hb.ServiceBackend(billing_project='huang-SC-Asia-schizophrenia', 
        remote_tmpdir='gs://soyeon-sga/Hail_Batch_huangSCAsiaSCZ', 
        regions=['us-central1']) 

b = hb.Batch(backend=backend, name=f'phase_index: c{cores}-{memory_type}-t{threads}: sample{sample_size}', default_image='skim212/glimpse2_hb:1.0')
fasta_file = b.read_input_group(**{'fasta': 'gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta', 
	'fasta.fai': 'gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai'})

ref_file = hfs.open('gs://glimpse2-datasharing/2023-09-4cM-splitRef/sg5k.split_ref.bin_list.full.rmSingle.gcpPath.txt') 

for line in ref_file.readlines():
    ref=line.strip("\n")
    reference_bin = ref
    reference_bin = b.read_input(reference_bin) 
    ref_name=ref.split("chr")[-1].split(".bin")[0]
    print(ref_name)
    sample_batch_list = hfs.open(f"gs://soyeon-sga/ReferencePanels/biox18_imputed/Cram_list_{sample_size}split/{case_or_control}/sample_batch_{case_or_control}.list") 
    for sample_batch in sample_batch_list.readlines():
        sample_batch=sample_batch.strip("\n") 
        sample_batch_name=sample_batch.split("sample_")[-1].split(".")[0]
        # print(sample_batch_name)
        # submit batch job
        j = b.new_job(name='biox18k-sg5k')
        j._preemptible = True
        #j.memory('16Gi')
        j.cpu(f'{cores}').memory(f'{memory_type}')
        # mount a bucket
        j.cloudfuse(f'{terra_path}', f'{mount_folder}')
        sample_batch = b.read_input(sample_batch)
        refsam_name = ref_name + "." + sample_batch_name
        print(refsam_name)
        j.declare_resource_group(**{refsam_name: {'bcf': '{root}.bcf', 'csi':'{root}.csi'}})
        j.command(f'GLIMPSE2_phase_static --bam-list {sample_batch} --fasta {fasta_file.fasta} --reference {reference_bin} --output {j[refsam_name].bcf} --threads {threads}')
        j.command(f'bcftools index -f {j[refsam_name].bcf} -o {j[refsam_name].csi} --threads {threads}')
        # save output
        b.write_output(j[refsam_name].bcf, f"gs://soyeon-sga/ReferencePanels/biox18_imputed/20231203_{sample_size}sample/c{cores}-{memory_type}-t{threads}-in2loop/{sample_batch_name}/{refsam_name}.bcf")
        b.write_output(j[refsam_name].csi, f"gs://soyeon-sga/ReferencePanels/biox18_imputed/20231203_{sample_size}sample/c{cores}-{memory_type}-t{threads}-in2loop/{sample_batch_name}/{refsam_name}.csi")
b.run(disable_progress_bar=True, wait=False)
