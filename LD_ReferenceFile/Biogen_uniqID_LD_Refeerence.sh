
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


for i in {1..22} X XY; do 

rm ukb_imp_chr${i}_mac50_info07_b0_7_patched_bfiles.bim
rm ukb_imp_chr${i}_mac50_info07_b0_7_patched_bfiles.bed
rm ukb_imp_chr${i}_mac50_info07_b0_7_patched_bfiles.fam

done 



####Few duplicate variants presents that needs to be removeed
#plink --bfile ukb_imp_chr1_22XXY_mac50_info07_b0_7_patched_bfiles --list-duplicate-vars
#plink --bfile ukb_imp_chr1_22XXY_mac50_info07_b0_7_patched_bfiles --exclude suppress-first plink.dupvar --out ukb_imp_chr1_22XXY_mac50_info07_b0_7_patched_bfiles_uniq











##########Not completed

for chr in {1..22} X XY; do 

plink --bfile /edgehpc/dept/human_genetics/users/jjohn1/Referencefile_h19/ukb_imp_chr1_22_mac50_info07_b0_7_patched_bfiles_Hg38/ukb_imp_chr1_22XXY_mac50_info07_b0_7_patched_bfiles \
    --make-bed  --snps-only --chr ${chr} --out ukb_imp_chr1_22XXY_mac50_info07_b0_7_patched_bfiles_SNP_Only_Chr${chr} 

done 



for chr in {1..22} X XY; do 

/edgehpc/dept/human_genetics/users/jjohn1/Software/gotcloud/bin/vcfCooker \
    --write-vcf --bgzf \
    --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --in-bfile ukb_imp_chr1_22XXY_mac50_info07_b0_7_patched_bfiles_SNP_Only_Chr${chr} \
    --out ukb_imp_chr1_22XXY_mac50_info07_b0_7_patched_bfiles_SNP_Only_Chr${chr}_vcfCooker

done


for chr in {1..22} X XY; do 

plink  --vcf ukb_imp_chr1_22XXY_mac50_info07_b0_7_patched_bfiles_SNP_Only_Chr${chr}_vcfCooker.vcf.gz --keep-allele-order \
    --make-bed --out ukb_imp_chr1_22XXY_mac50_info07_b0_7_patched_bfiles_SNP_Only_Chr${chr}_vcfCooker 

done 


source /home/jjohn1/modulefiles/anaconda3/bin/activate
module load SAMtools/1.15-GCC-11.2.0
module load PLINK/1.9b_6.21-x86_64
module load R/4.1.3-foss-2021b


