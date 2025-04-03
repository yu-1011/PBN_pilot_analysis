# PBN WGS analysis

## Imputation
Imputation is implemented in the Terra workspace. It hosts workflows and example data for imputation of variants in low coverage sequencing data. The workflows utilize GLIMPSE2 for imputation. GLIMPSE2 is designed for reference panels containing hundreds of thousands of reference samples, with a special focus on rare variants.

The input to this workflow can be either single-sample or multi-sample VCFs with existing GT and PL annotations or CRAMs, in which case GLIMPSE2 will calculate the PLs using pileup calling.

### Output
- File imputed_vcf: Single, multi-sample imputed VCF that covers all regions defined in the contig_regions input in GlimpseSplitReference. The name of the file is the basename of input_vcf with .imputed.vcf.gz added.
- File imputed_vcf_index: Index to imputed_vcf.

## WGS QC
This repository details the quality control (QC) pipeline of the Psychiatric Biomarker Network (PBN), which largely follows the recommendations described in:

For autosome: Peterson, R. E., Kuchenbaecker, K., Walters, R. K., Chen, C.-Y., Popejoy, A. B., Periyasamy, S., et al. (2019). Genome-wide Association Studies in Ancestrally Diverse Populations: Opportunities, Methods, Pitfalls, and Recommendations. Cell, 179(3), 589–603. http://doi.org/10.1016/j.cell.2019.08.051

For X chromosome: Khramtsova, Ekaterina A., et al. "Quality control and analytic best practices for testing genetic models of sex differences in large populations." Cell 186.10 (2023): 2044-2061. https://www.sciencedirect.com/science/article/pii/S0092867423004105


### Sample QC
- Sample-level call rate: > 0.98
- Remove samples that fail sex check: `--check-sex`
- Absolute value of autosomal heterozygosity rate deviating from the mean: e.g., 3SD; `--het`
- Identify unrelated individuals (Pi_hat < 0.2) within Asian samples

### Variants QC
- SNP-level call rate: > 0.98
- HWE: > 1e-10
- Retain only SNPs: Excluding indels and monomorphic SNPs (for imputation; HRC: SNP only)
- For chromosome X:
  - Remove PAR region: The PARs are frequently removed because, although on the sex chromosomes, they do not behave as strictly X or Y regions.
  - Remove sex-specificity MAF SNPs: Variants are filtered if they have significantly different MAF between male and female controls.
  - Remove sex-specificity missingness SNPs: Variants are filtered if they have significantly different missingness between male and female controls.
  - Remove sex-assoc SNPs: Variants are filtered if they have significantly different minor allele frequency between male and female controls.

### Population Assignment
- Select common, high-quality SNPs for population inference:
  - SNP-level call rate > 0.98
  - Remove strand ambiguous SNPs and long-range LD regions (chr6:25-35Mb; chr8:7-13Mb inversion)
  - Prune to < 100K independent SNPs
  - Merge with 1KG reference

### Prepare Data for TopMed/Singapore 50k Imputation Server
- Harmonize study data with HRC data
- Convert plink to vcf by chromosome

### Send Unrelated East Asian Samples to TOPMed Server for Imputation
- Chromosome X:
  - Ploidy Check: Verifies if all variants in the nonPAR region are either haploid or diploid.
  - Mixed Genotypes Check: Verifies if the number of mixed genotypes (e.g., 1/.) is < 10 %.

### Post-imputation QC (Converting vcf dosages to plink hard-call genotypes)
- INFO score/Imputation R2: > 0.8
- MAF (based on imputed dosages): > 1%
- HWE: > 1e-10
- SNP-level call rate: > 0.98

### Population Assignment After imputation
- Select common, high-quality SNPs for population inference:
  - SNP-level call rate > 0.98
  - Remove strand ambiguous SNPs and long-range LD regions (chr6:25-35Mb; chr8:7-13Mb inversion)

### Output for QC
- Sex check figure
- IBD check figure
- Het check figure
-	Population assignment figures before imputation 
-	Population assignment figures after imputation 
-	QC statistics table 

## C4 Copy number imputation

C4 CN imputation is performed by Osprey (r0.1-r10)

### Output for C4 imputation
There are two main outputs of the imputation algorithm:

- CNV Imputation: cnvOutputFile="${outDir}/${inputBaseName}.cnvs.imputed.vcf.gz"
  - This VCF file contains phased/imputed genotypes for the "classic" C4 CNV features, including the copy number of
total C4, separate copy numbers for C4A and C4B, and HERV copy number.

  The phased genotypes are represented using the following application-defined FORMAT tags in the VCF:
  - PCN: This indicates the phased copy numbers in the format "n|n".
  - PCNF: Phased copy number dosage estimates in the format "n.nn|n.nn".
  - PCNQ: A Phred-scaled estimate of the confidence in the phasing/imputation quality, one value for each haploptype, in the format "q1|q2".
  - PCNL: Two vectors of copy number likelihoods, one for each haplotype, separated by vertical bar.

  These fields are analogous to the Genome STRiP CN,CNF,CNQ,CNL/CNP fields for diploid copy number estimates,
except that there are two values in each tag, one for each phased haplotype.

  In addition, there is an INFO tag OSPR2 on each variant which gives an estimate of the imputation quality
for each variant. This may be useful for filtering out poorly-imputed variants.

  - The information in the VCF are further summarized in the text file ${outDir}/refpanel_sample_haps_imputed.txt
  in a format consistent with earlier versios of the pipeline:
    This is a tab-delimited text file with three columns.
    - SAMPLE: The sample identifier from the VCF for each target sample.
    - H1,H2: Phased C4 alleles encoded as (H_n_n_n_n).
  
      The four numbers correspond in the H_n_n_n_n notation correspond to the number of (haploid) copies
    of C4, C4A, C4B and the HERV, respectively. For example, the haplotype H_2_1_1_1 represents a
    haplotype with 2 copies of C4, one of which are C4A and one of which is C4B, along with one HERV.
      We typically write this "AL-BS", although the relationship of the HERV to C4A or C4B is not necessarily accurate.
    The H1 and H2 haplotypes correspond to the two haplotypes in the phased VCF created by eagle.
  
  - There is also an additional text file produced ${outDir}/refpanel_phasing_summary.txt.
    This is a tab-delimited file containing statistics on how well the phasing/imputation preserved the diploid copy number measurements of the reference panel (the 1000 Genomes samples).
    This may be useful for debugging. The R2 numbers should generally be above 0.9.
    This is not a cross-validation result.

    - Cross-Validation of CNV Imputation
  
    The script osprey_cross_validate.sh is called by run_imputation_example to perform cross-validation using the
  1000 Genomes samples to provide information on the accuracy of the imputation based on the SNPs in your input
  panel. The results of the cross-validation are written to ${outDir}/crossval.
  
    The cross-validation iteratively holds out subsets of the input samples (default batch size is 51 samples)
  and imputes the held-out samples using the remainder of the 1000 Genomes cohort in the reference panel.
  This is done separately for each continental population. Then the imputed values for all samples are compared
  to the diploid copy number dosages from the reference panel.
  
    - The main output of the cross-validation consists of two summary tables:
      - ${outDir}/crossval/1000G.C4AB_cnvs.CNR2.summary.txt
      - ${outDir}/crossval/1000G.C4AB_cnvs.CNFR2.summary.txt
  
    These tables give the squared Pearson correlation of either integer copy number (CNR2) or fractional copy
  number (CNFR2) between the imputed diploid dosage and the estimated dosage for each population.
  The correlation across all samples in the reference panels is shown as population "ALL".
  
    These summary tables also list the number of SNPs in common between the reference panel and your input panel.
  In general, we would expect the imputation to work better when there are more overlapping SNPs.


  - Haplotype Refinement

    An optional but recommended step after performing CNV imputation is to run scripts/refine_C4_haplotypes.sh.
    This script takes the VCF file produced by CNV imputation as input.
    It does a C4-specific refinement of the phased haplotypes by recomputing the likelihoods over a set of possible
    haplotypes that take into account features of the C4 haplotype encoding, including the dependency between total
    C4 as the sum of C4A+C4B and the constraint that the HERV copy number cannot exceed the C4 copy number.
    
    This script produces two output files:
    - .phased_haps_refined.txt:
      There are four columns for each sample.
      - H1 and H2 have the refined haplotypes, encoded as H_n_n_n_n as above.
      - Q1 and Q2 have phred-scaled quality values for the haplothpe as a whole.
  
    - .phased_hap_refinement.txt:
      There are two rows for each sample, corresponding to haplotype 1 and haplotype 2, indicated by the HAP column.
      For each haplotype, there are four columns:
      ORIGINAL and REFINED contain the original haplotype and the refined haplotype, in H_n_n_n_n encoding as above.
      OQUAL and RQUAL contain the original quality and the refined quality scores, phred-scaled.

- PSV Imputation

  - psvOutputFile="${outDir}/${inputBaseName}.psvs.imputed.vcf.gz"
  This VCF file contains phased/imputed values for a set of paralogous sequence variants (PSVs) ascertained from
  the 1000 Genomes sequencing data. The PSVs include rare forms of C4 (C4A1 and C4A2), some rare putative C4 LoF
  variants, and other coding and non-coding PSVs.
  
  The fields in the PSV output file are the same as for the CNV output file.
  Each PSV is modeled as a copy number variant giving the number of copies of the PSV estimated to be carried on
  each haplotype.


