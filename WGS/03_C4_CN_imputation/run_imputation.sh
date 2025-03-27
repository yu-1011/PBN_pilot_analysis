#!/bin/bash
#$ -cwd
#$ -pe smp 4 -binding linear:4
#$ -l h_rt=12:00:00,h_vmem=25G

# Example of running imputation using Osprey version 0.1-r10 (May, 2022) with the "universal" 1000 Genomes reference panel.
# There are two input arguments:
#  inputPanel: A file containing snp genotypes into which C4 CNVs/PSVs should be imputed
#    The inputPanel can be in either .vcf, .vcf.gz or .bcf format.
#    The input panel must be a valid vcf/bcf on hg38 (see below).
#  outDir: An output directory (will be created if it doesn't exist) for the output and intermediate files

inputPanel="$1"
outDir="$2"
if [ -z "${outDir}" ]; then
    echo "Usage: inputPanel outDir" && exit 1
fi

source /broad/software/scripts/useuse

reuse .bcftools-1.10.2
reuse .shapeit4-4.2.1
reuse R-4.0
reuse .tabix-0.2.6

# The following variables may need to be set based on the input panel and local installation diretories.
# Note that the input panel is assumed to be in bgzip'd vcf format (with tabix index) or a bcf file (with csi index).
# Currently only hg38 is supported

genome=hg38
refPanelDir=osprey_C4_1000G_Universal_Nov2022

# Shapeit 4.2 is used in the new pipeline
shapeit=/humgen/cnp04/bobh/shapeit/shapeit4-4.2.0/bin/shapeit4.2
shapeit_geneticMap=/humgen/cnp04/bobh/shapeit/shapeit4-4.2.0/maps/b38/chr6.b38.gmap.gz
shapeit_threads=4

ospDir=${refPanelDir}/osprey
ospIters=100
ospThreads=4
# Example of doing scattering for psvs to speed things up, see below.
ospScatter=false

snpReferencePanel=${refPanelDir}/C4_panel_1000G_ext_Universal_${genome}.bcf
cnvReferencePanel=${refPanelDir}/C4_panel_1000G_ext_Universal_${genome}.C4AB_cnvs.vcf.gz
psvReferencePanel=${refPanelDir}/C4_panel_1000G_ext_Universal_${genome}.C4_psvs.vcf.gz
refPanelSampleList=${refPanelDir}/C4_panel_1000G_ext_Universal.samples.list

declare -A targetIntervals=([hg19]=6:31948000-32014000 [hg38]=chr6:31980500-32046500)
targetInterval=${targetIntervals[${genome}]}

mkdir -p ${outDir} || exit 1

inputBaseName="$(basename ${inputPanel} | sed 's/.bcf$//' | sed 's/.vcf.gz$//' | sed 's/.vcf$//')"

mergedPanel="${outDir}/${inputBaseName}.merged.vcf.gz"
mergedPanelPhased="${outDir}/${inputBaseName}.merged.phased.vcf.gz"
ibsMatrix="${outDir}/${inputBaseName}.merged.ibs.gz"
cnvOutputFile="${outDir}/${inputBaseName}.cnvs.imputed.vcf.gz"
psvOutputFile="${outDir}/${inputBaseName}.psvs.imputed.vcf.gz"

# Create merged panel
# This can be quite slow/expensive.
if [ ! -e "${mergedPanel}" ]; then
    echo $(date) "Merging input panel and reference panel ..."
    isecDir=${outDir}/isec
    mkdir -p ${isecDir} || exit 1
    bcftools view -Ob ${inputPanel} > ${isecDir}/input_panel.tmp.bcf || exit 1
    bcftools index ${isecDir}/input_panel.tmp.bcf || exit 1
    bcftools isec -Ob -p ${isecDir} -w 1,2 -n=2 ${isecDir}/input_panel.tmp.bcf ${snpReferencePanel} || exit 1
    bcftools merge -Ou ${isecDir}/0000.bcf ${isecDir}/0001.bcf \
        | bcftools annotate -Oz -x 'INFO,^FORMAT/GT' > ${mergedPanel} || exit 1
    tabix -p vcf ${mergedPanel} || exit 1
    # Leave this for debugging for now
    #rm -rf ${isecDir}
fi

# Phase the merged panel with shapeit4.2
if [ ! -e "${mergedPanelPhased}" ]; then
    echo $(date) "Phasing merged panel ${mergedPanel} ..."
    # Note: consider adding --sequencing if using WGS data
    logFile="$(echo ${mergedPanelPhased} | sed 's/.vcf.gz$//' | sed 's/$/.log/')"
    chr="$(echo ${targetInterval} | awk -F : '{ print $1 }')"
    ${shapeit} \
        --input ${mergedPanel} \
        --map ${shapeit_geneticMap} \
        --region ${chr} \
        --sequencing \
        --output ${mergedPanelPhased} \
        --thread ${shapeit_threads} \
        --log ${logFile} \
        || exit 1
    tabix -p vcf ${mergedPanelPhased} || exit 1
fi

# Compute IBS matrix from merged panel
if [ ! -e "${ibsMatrix}" ]; then
    echo $(date) "Computing IBS matrix for merged panel ${mergedPanelPhased} ..."
    ${ospDir}/ospreyIBS \
        -vcf ${mergedPanelPhased} \
        -ibs ${ibsMatrix} \
        -gmap ${shapeit_geneticMap} \
        -L ${targetInterval} \
        -t ${ospThreads} \
        -rs ${refPanelSampleList} \
        || exit 1
fi

# Cross-validation of CNV imputation from merged panel
if [ ! -e ${outDir}/crossval/1000G.C4AB_cnvs.CNFR2.summary.txt ]; then
    # The results from cross validation in ${outDir}/crossval should be inspected for QC
    ${refPanelDir}/scripts/osprey_cross_validate.sh \
        ${mergedPanelPhased} ${refPanelDir} || exit 1
fi

# CNV imputation
if [ ! -e ${cnvOutputFile} ]; then
    echo $(date) "Imputing CNVs into ${inputPanel} ..."
    ${ospDir}/osprey --version > ${outDir}/osprey_version.txt || exit 1
    ${ospDir}/osprey \
        -vcf ${cnvReferencePanel} \
        -o ${cnvOutputFile} \
        -ibs ${ibsMatrix} \
        -t ${ospThreads} \
        -iter ${ospIters} \
        || exit 1
fi

# PSV imputation
psvOutputFile=${outDir}/${inputBaseName}.psvs.imputed.vcf.gz
if [ ! -e ${psvOutputFile} ]; then
    if [ "${ospScatter}" == "true" ]; then
        # Osprey supports scattering by site using the --siteIndex argument.
        # For example --siteIndex 1@10 means the first batch of 10 sites, 2@10 means sites 11-20, etc.
        # The code below uses tools that run on our local cluster, you will have to adapt to your own environment.
        # In particular, run_queue_scatter.sh scatters a script (first argument) based on the second argument and waits for all jobs to complete.
        # bcftools concat can be used to merge the results as shown.
        psvBatchSize=10
        nSites="$(zcat ${psvReferencePanel} | grep -v ^# | wc -l)"
        nSiteBatches="$(echo ${nSites} ${psvBatchSize} | awk '{ print int(($1 + ($2-1))/$2) }')"

        echo $(date) "Imputing PSVs into ${inputPanel} (${nSites} sites in ${nSiteBatches} batches of ${psvBatchSize} sites each) ..."
        scatterArgs="$(seq 1 ${nSiteBatches} | awk -v B=${psvBatchSize} '{ print $1 "@" B }')"
        scatterArgs="$(echo ${scatterArgs} | sed 's/ /,/g')"
        scatterDir=${outDir}/scatter/psvs
        scatterLogDir=${scatterDir}/logs
        scatterOutputFilePrefix=${scatterDir}/"$(basename ${psvOutputFile} | sed 's/.vcf.gz$//')"

        # Note: Deleting the log files forces redo
        rm -rf ${scatterLogDir}
        mkdir -p ${scatterLogDir} || exit 1

        # This script (not supplied) scatters the jobs by siteIndex and waits for completion.
        # Each output file name is of the form ${scatterOutputFilePrefix}.k@N.vcf.gz
        /stanley/genome_strip/tools/scripts/run_queue_scatter.sh \
            -jobLogDir ${scatterLogDir} \
            scripts/osprey_scatter.sh \
                ${scatterArgs} \
                ${scatterOutputFilePrefix} \
                -vcf ${psvReferencePanel} \
                -ibs ${ibsMatrix} \
                -t ${ospThreads} \
                -iter ${ospIters} \
                || exit 1

        echo $(date) "Merging PSV output files ..."
        mergeFiles="$(ls ${scatterDir}/$(basename ${scatterOutputFilePrefix}).*.vcf.gz | sort -V)"
        bcftools concat -Oz ${mergeFiles} > ${psvOutputFile} || exit 1
        tabix -p vcf ${psvOutputFile} || exit 1
    else
        echo $(date) "Imputing PSVs into ${inputPanel} ..."
        ${ospDir}/osprey \
            -vcf ${psvReferencePanel} \
            -o ${psvOutputFile} \
            -ibs ${ibsMatrix} \
            -t ${ospThreads} \
            -iter ${ospIters} \
            || exit 1
    fi
fi

if [ ! -e ${outDir}/refpanel_phasing_summary.txt ]; then
    (bcftools query -f 'ID[\t%SAMPLE]\n' ${cnvOutputFile} | head -n 1;
     bcftools query -f 'PST[\t%PST]\n' ${cnvOutputFile} | head -n 1;
     bcftools query -f '%ID[\t%PCN]\n' ${cnvOutputFile}) \
        > ${outDir}/${inputBaseName}.cnvs.pcns.txt || exit 1

    Rscript ${refPanelDir}/scripts/format_phased_haps.R ${outDir}/${inputBaseName}.cnvs.pcns.txt P \
        > ${outDir}/refpanel_sample_haps_phased.txt || exit 1
    Rscript ${refPanelDir}/scripts/format_phased_haps.R ${outDir}/${inputBaseName}.cnvs.pcns.txt I \
        > ${outDir}/refpanel_sample_haps_imputed.txt || exit 1

    Rscript ${refPanelDir}/scripts/compute_phasing_accuracy.R \
        ${outDir}/refpanel_sample_haps_phased.txt \
        ${refPanelDir}/1000G_sample_diplotypes.txt \
        "C4,C4A,C4B,HERV" \
        ${refPanelDir}/1000G_sample_populations.txt \
        > ${outDir}/refpanel_phasing_summary.txt \
        || exit 1
fi

echo $(date) "Script completed successfully."
