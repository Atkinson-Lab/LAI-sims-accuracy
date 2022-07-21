

cd /storage/atkinson/home/u242335/jointcall_phasing

#maf, triallelic, indel filter

for i in {22..1}; do \
/storage/atkinson/shared_resources/software/plink2/plink2 --vcf \
/storage/atkinson/shared_resources/reference/1kG_HGDP_jointcall/raw/genotype/gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.bgz \
--max-alleles 2 --snps-only just-acgt --maf 0.005 --make-bed --out filter_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}; done 

# 3 - remove dup ids

#per id

#check first
for i in {22..1}; do \
/storage/atkinson/shared_resources/software/plink2/plink2 --bfile filter_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i} \
--set-missing-var-ids @:#[b38] --rm-dup --make-bed --out jointcall_dups_chr${i}; done

#exclude

for i in {22..1}; do \
/storage/atkinson/shared_resources/software/plink2/plink2 --bfile filter_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i} \
--set-missing-var-ids @:#[b38] --rm-dup exclude-all --geno 0.1 --make-bed --out QCed_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}; done

## for shapeit 4

for i in {22..1}; do \
/storage/atkinson/shared_resources/software/plink2/plink2 --bfile filter_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i} \
--set-missing-var-ids @:#[b38] --rm-dup exclude-all --geno 0.1 --export vcf --out QCed_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}; done

for i in {22..1}; do bgzip -c QCed_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf > QCed_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.gz; done

for i in {22..5}; do bcftools index QCed_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.gz; done

#running shapeit4
for i in {22..19}; do \
/mnt/genetica_1/shapeit4.2 \
--input QCed_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.gz \
--map /mnt/genetica_1/Jessica/admixsimu/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt  \
--output phased_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.gz \
--region ${i} \
--sequencing \
--log phased_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.log \
--thread 12 ; done

## old: shapeit2

#../phasing/shapeit \
#-B QCed_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i} \
#--input-map ../phasing/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt \ #check
#--output-max phased_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.haps.gz \
#phased_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.sample \
#--thread 16

