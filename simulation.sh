
# 1kg .sample file from https://mathgen.stats.ox.ac.uk/impute/1000GP%20Phase%203%20haplotypes%206%20October%202014.html

##### Select which indivs go to simulation
### ADMIX-SIMU
# 31 NAT (half of Surui, Karitiana, Maya, Pima) = SIM
egrep "Surui|Karitiana|Maya|Pima|Colombian" BHRC_KIDS > NAT_HGDP_SIMU
nano NAT_HGDP_SIMU

awk '{print $1}' NAT_HGDP_SIMU > NAT_HGDP_SIMU_id

# 30 IBS = SIM
grep IBS ./1000GP_Phase3/1000GP_Phase3.sample | awk '{print $1}' | head -n 30 > IBS_1KG_SIMU

# 30 YRI = SIM
grep IBS ./1000GP_Phase3/1000GP_Phase3.sample | awk '{print $1}' | head -n 30 > YRI_1KG_SIMU

##### make hapsample files out of phased dataset with the individuals for SIMU.
for i in {1..22}; do \
/storage/atkinson/home/u242335/phasing/shapeit -convert --input-haps \
../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--include-ind NAT_HGDP_SIMU_id \
--output-haps NAT_HGDP_SIMU_chr${i}.haps NAT_HGDP_SIMU_chr${i}.sample ; done

for i in {1..22}; do \
/storage/atkinson/home/u242335/phasing/shapeit -convert --input-haps \
../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--include-ind IBS_1KG_SIMU \
--output-haps IBS_1KG_SIMU_chr${i}.haps IBS_1KG_SIMU_chr${i}.sample ; done

for i in {1..22}; do \
/storage/atkinson/home/u242335/phasing/shapeit -convert --input-haps \
../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--include-ind YRI_1KG_SIMU \
--output-haps YRI_1KG_SIMU_chr${i}.haps YRI_1KG_SIMU_chr${i}.sample ; done

#### Make .phgeno files out of hapsample for SIMU
for i in {1..22};do \
cut -d' ' -f 6- NAT_HGDP_SIMU_chr${i}.haps | sed 's/\s//g' > NAT${i}.phgeno; \
cut -d' ' -f 6- YRI_1KG_SIMU_chr${i}.haps | sed 's/\s//g' > AFR${i}.phgeno; \
cut -d' ' -f 6- IBS_1KG_SIMU_chr${i}.haps | sed 's/\s//g' > EUR${i}.phgeno;done


##### make a VCF file out of phased dataset with the individuals for SIMU (for .snp file).
for i in {1..22}; do \
/storage/atkinson/home/u242335/phasing/shapeit -convert --input-haps \
../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--include-ind NAT_HGDP_SIMU_id \
--output-vcf 1KG_HGDP_ADMIXSIMU_NAT_chr${i}.vcf; done


## create the .snp files
for i in {1..22};do \
/storage/atkinson/shared_resources/software/plink2/plink2 --vcf 1KG_HGDP_ADMIXSIMU_NAT_chr${i}.vcf \
--const-fid 0 --max-alleles 2 --make-bed --out 1kg_hgdp_pBim2_${i};done

for i in {1..22};do perl ./admix-simu-master/insert-map.pl 1kg_hgdp_pBim2_${i}.bim ../../phasing/genetic_map/recomb-hg38/genetic_map_chr${i}_hg38.txt > chr${i}.pos; \
awk -F' ' '{ print $1":"$4"_"$5"_"$6, $1, $3, $4, $5, $6 }' chr${i}.pos > 1kg_hgdp.chr${i}.snp; done


## Edit here the ancestry orders and which population you will test, remember to manually write the .dat file.
export  anc1=NAT
export  anc2=EUR
export  anc3=AFR
export  AdmixPOP=Brasa


# edit .dat file with number of generations and ancestry proportions
# 3 way admixture
for i in {1..22}; do ../admix-simu-master/simu-mix.pl ./$AdmixPOP.dat ../1kg_hgdp.chr${i}.snp $AdmixPOP${i} -$anc1 ../$anc1${i}.phgeno -$anc2 ../$anc2${i}.phgeno -$anc3 ../$anc3${i}.phgeno ; done


# convert .bp from the simulation to .hanc
for i in {1..22}; do ../admix-simu-master/bp2anc.pl $AdmixPOP${i}.bp > $AdmixPOP${i}.hanc; done

# hanc1
for i in {1..22}; do sed 's/./& /g' $AdmixPOP${i}.hanc > $AdmixPOP${i}.hanc1 ;done

# transpose hanc1
for i in {1..22};do python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < $AdmixPOP${i}.hanc1 > $AdmixPOP${i}.hanc2; done

head $AdmixPOP${i}.hanc2
rm *.hanc1

# From Elizabeth's scripts: Now I need to separate the .phgeno file to have a space between each character
for i in {1..22}; do sed 's/./& /g' ./$AdmixPOP${i}.phgeno > temp_Genos_chr${i}; done
for i in {1..22}; do awk -F' ' '{ print $2, $1, $4, $5, $6 }' ../1kg_hgdp.chr${i}.snp > temp_SNP_chr${i} ; done
for i in {1..22}; do paste temp_SNP_chr${i} temp_Genos_chr${i} > temp_admix_chr${i}.haps ; done
for i in {1..22}; do sed 's/\t/ /g' temp_admix_chr${i}.haps > $AdmixPOP.chr${i}.haps ; done
head $AdmixPOP.chr1.haps

rm temp_*

# Manually created the .sample file with 30 ID names (e.g. SIM1, SIM2, SIM3)

for i in {1..22}; do cp $AdmixPOP.sample.txt $AdmixPOP.chr${i}.sample; done

