
##Creating uniq ids based on Hg38
for i in {1..22} X XY; do 

awk -F"\t" 'BEGIN{OFS="\t"}{print $1,$1"_"$4"_"$6"_"$5,$3,$4,$5,$6}' /edgehpc/dept/compbio/human_genetics/users/bsun/UKB_PPP_pQTLs/analysis/interim/ukb_genetic_variants/ukb_imp_chr${i}_mac50_info07_b0_7_patched_bfiles.bim > ukb_imp_chr${i}_mac50_info07_b0_7_patched_bfiles.bim

cp /edgehpc/dept/compbio/human_genetics/users/bsun/UKB_PPP_pQTLs/analysis/interim/ukb_genetic_variants/ukb_imp_chr${i}_mac50_info07_b0_7_patched_bfiles.bed  ukb_imp_chr${i}_mac50_info07_b0_7_patched_bfiles.bed
cp /edgehpc/dept/compbio/human_genetics/users/bsun/UKB_PPP_pQTLs/analysis/interim/ukb_genetic_variants/ukb_imp_chr${i}_mac50_info07_b0_7_patched_bfiles.fam  ukb_imp_chr${i}_mac50_info07_b0_7_patched_bfiles.fam

done 


##Mergee the plink files
for i in {2..22} X XY; do
    echo ukb_imp_chr1_mac50_info07_b0_7_patched_bfiles | sed 's/\.bed$//g' >> all_plinkfiles.txt
done

plink --bfile ukb_imp_chr1_mac50_info07_b0_7_patched_bfiles --merge-list all_plinkfiles.txt --make-bed --out ukb_imp_chr1_22XXY_mac50_info07_b0_7_patched_bfiles


####Few duplicate variants presents that needs to be removeed
plink --bfile ukb_imp_chr1_22XXY_mac50_info07_b0_7_patched_bfiles --list-duplicate-vars
plink --bfile ukb_imp_chr1_22XXY_mac50_info07_b0_7_patched_bfiles --exclude suppress-first plink.dupvar --out ukb_imp_chr1_22XXY_mac50_info07_b0_7_patched_bfiles_uniq

