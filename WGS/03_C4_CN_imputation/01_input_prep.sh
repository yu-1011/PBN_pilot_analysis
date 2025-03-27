tabix -p vcf PBN_pilot.imputed.vcf.gz

plink --bfile PBN_pilot_final \
      --chr 6 \
      --from-bp 24000001 \
      --to-bp 34000000 \
      --make-bed \
      --out PBN_pilot_final_chr6_24M_34M
      
plink --bfile PBN_pilot_final_chr6_24M_34M \
      --recode vcf bgz \
      --out PBN_pilot_final_chr6_24M_34M

zcat PBN_pilot_final_chr6_24M_34M.vcf.gz | \
awk 'BEGIN{OFS="\t"} /^#/ {print} !/^#/ {$1="chr"$1; print}' | \
bgzip > PBN_pilot_final_chr6_24M_34M_chr_updated.vcf.gz
bcftools index PBN_pilot_final_chr6_24M_34M_chr_updated.vcf.gz

