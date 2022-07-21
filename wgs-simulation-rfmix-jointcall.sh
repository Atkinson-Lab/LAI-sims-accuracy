### RFMIX
# 31 NAT (other half of Surui, Karitiana, Maya, Pima) = RFMIX
egrep "Surui|Karitiana|Maya|Pima|Colombian" BHRC_KIDS > NAT_HGDP_RFMIX
nano NAT_HGDP_RFMIX # edited manually to be the counterpart of NAT_HGDP_SIMU

mv NAT_HGDP_RFMIX temp_NAT_HGDP_RFMIX
awk '{print $2}' tepmp_NAT_HGDP_RFMIX > NAT_HGDP_RFMIX

rm temp_NAT_HGDP_RFMIX

# 77 IBS = RFMIX . using the same indivs as used for the prev simulations = unrelated samples.
grep IBS ./1000GP_Phase3/1000GP_Phase3.sample | awk '{print $1}' | tail -n 77 > IBS_1KG_RFMIX

# 77 YRI = RFMIX. using the same indivs as used for the prev simulations = unrelated samples.
grep YRI ./1000GP_Phase3/1000GP_Phase3.sample | awk '{print $1}' | tail -n 77 > YRI_1KG_RFMIX

##################
# RFmix - preparo dos arquivos
#################

################# NAT-HGDP IBS YRI #######################

#from /storage/atkinson/home/u242335/lai/simu-jointcall/wgs
##### make hapsample files out of phased dataset with the individuals for RFMix REF

for i in {1..22}; do \
/storage/atkinson/shared_resources/software/plink2/plink2 --vcf \
/storage/atkinson/home/u242335/1kg_hg38_phased/QC_phased_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.gz \
--keep /storage/atkinson/home/u242335/lai/simu-jointcall/NAT_HGDP_RFMIX \
--export haps --out REF_NAT_HGDP_chr${i}; done

for i in {1..22}; do \
/storage/atkinson/shared_resources/software/plink2/plink2 --vcf \
/storage/atkinson/home/u242335/1kg_hg38_phased/QC_phased_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.gz \
--keep /storage/atkinson/home/u242335/lai/simu-jointcall/IBS_1KG_RFMIX \
--export haps --out REF_IBS_1KG_chr${i}; done

for i in {1..22}; do \
/storage/atkinson/shared_resources/software/plink2/plink2 --vcf \
/storage/atkinson/home/u242335/1kg_hg38_phased/QC_phased_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.gz \
--keep /storage/atkinson/home/u242335/lai/simu-jointcall/YRI_1KG_RFMIX \
--export haps --out REF_YRI_1KG_chr${i}; done


# printar as amostras ref em um mesmo arquivo.

sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_NAT_HGDP_chr1.sample | awk '{print $2}' > /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_RFMIX-SIMU_NAT_EUR_AFR.ref
sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_IBS_1KG_chr1.sample | awk '{print $2}' >> /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_RFMIX-SIMU_NAT_EUR_AFR.ref
sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_YRI_1KG_chr1.sample | awk '{print $2}' >> /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_RFMIX-SIMU_NAT_EUR_AFR.ref

#### 2 way admixture

sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_NAT_HGDP_chr1.sample | awk '{print $2}' > /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_RFMIX-SIMU_NAT_EUR.ref
sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_IBS_1KG_chr1.sample | awk '{print $2}' >> /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_RFMIX-SIMU_NAT_EUR.ref

# arrumar o .sample do plink (O PLINK MUDA A ORDEM DOS INDIVS DA REF então não dá pra copiar da ultima run)

head -n 2 /storage/atkinson/home/u242335/lai/simu-jointcall/REF_NAT_HGDP_chr1.sample > header_sample.txt

for i in {1..22}; do \
sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_NAT_HGDP_chr${i}.sample | awk '{print $1, $2, $3, 0, 0, 0, -9}' > temp_REF_NAT_HGDP_chr${i}; \
rm /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_NAT_HGDP_chr${i}.sample ; \
cat header_sample.txt temp_REF_NAT_HGDP_chr${i} > REF_NAT_HGDP_chr${i}.sample; done

for i in {1..22}; do \
sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_IBS_1KG_chr${i}.sample | awk '{print $1, $2, $3, 0, 0, 0, -9}' > temp_REF_IBS_1KG_chr${i}; \
rm /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_IBS_1KG_chr${i}.sample ; \
cat header_sample.txt temp_REF_IBS_1KG_chr${i} > REF_IBS_1KG_chr${i}.sample; done

for i in {1..22}; do \
sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_YRI_1KG_chr${i}.sample | awk '{print $1, $2, $3, 0, 0, 0, -9}' > temp_REF_YRI_1KG_chr${i}; \
rm /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_YRI_1KG_chr${i}.sample ; \
cat header_sample.txt temp_REF_YRI_1KG_chr${i} > REF_YRI_1KG_chr${i}.sample; done


################# PEL IBS YRI #######################
##### make hapsample files out of phased dataset with the individuals for RFMix REF

for i in {1..1}; do \
/storage/atkinson/shared_resources/software/plink2/plink2 --vcf \
/storage/atkinson/home/u242335/1kg_hg38_phased/QC_phased_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.gz \
--keep /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/PEL_1KG \
--export haps --out REF_NAT_PEL_1KG_chr${i}; done

# printar as amostras ref em um mesmo arquivo.

sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_NAT_PEL_1KG_chr1.sample | awk '{print $2}' > /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_RFMIX-SIMU_PEL_EUR_AFR.ref
sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_IBS_1KG_chr1.sample | awk '{print $2}' >> /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_RFMIX-SIMU_PEL_EUR_AFR.ref
sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_YRI_1KG_chr1.sample | awk '{print $2}' >> /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_RFMIX-SIMU_PEL_EUR_AFR.ref

#### 2 way admixture

sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_NAT_PEL_1KG_chr1.sample | awk '{print $2}' > /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_RFMIX-SIMU_PEL_EUR.ref
sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_IBS_1KG_chr1.sample | awk '{print $2}' >> /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_RFMIX-SIMU_PEL_EUR.ref

# arrumar o .sample do plink (O PLINK MUDA A ORDEM DOS INDIVS DA REF então não dá pra copiar da ultima run)

for i in {1..1}; do \
sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_NAT_PEL_1KG_chr${i}.sample | awk '{print $1, $2, $3, 0, 0, 0, -9}' > temp_REF_PEL_1KG_chr${i}; \
rm /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_NAT_PEL_1KG_chr${i}.sample ; \
cat header_sample.txt temp_REF_PEL_1KG_chr${i} > REF_NAT_PEL_1KG_chr${i}.sample; done


################# PEL-EAS IBS YRI #######################

for i in {1..1}; do \
/storage/atkinson/shared_resources/software/plink2/plink2 --vcf \
/storage/atkinson/home/u242335/1kg_hg38_phased/QC_phased_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.gz \
--keep /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/PEL_EAS_1KG_ids \
--export haps --out REF_NAT_PEL_EAS_1KG_chr${i}; done

# printar as amostras ref em um mesmo arquivo.

sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_NAT_PEL_EAS_1KG_chr1.sample | awk '{print $2}' > /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_RFMIX-SIMU_REF_PEL_EAS_EUR_AFR.ref
sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_IBS_1KG_chr1.sample | awk '{print $2}' >> /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_RFMIX-SIMU_REF_PEL_EAS_EUR_AFR.ref
sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_YRI_1KG_chr1.sample | awk '{print $2}' >> /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_RFMIX-SIMU_REF_PEL_EAS_EUR_AFR.ref

#### 2 way admixture

sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_NAT_PEL_EAS_1KG_chr1.sample | awk '{print $2}' > /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_RFMIX-SIMU_REF_PEL_EAS_EUR.ref
sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_IBS_1KG_chr1.sample | awk '{print $2}' >> /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_RFMIX-SIMU_REF_PEL_EAS_EUR.ref

# arrumar o .sample do plink (O PLINK MUDA A ORDEM DOS INDIVS DA REF então não dá pra copiar da ultima run)

for i in {1..1}; do \
sed '1,2d' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_NAT_PEL_EAS_1KG_chr${i}.sample | awk '{print $1, $2, $3, 0, 0, 0, -9}' > temp_REF_PELEAS_1KG_chr${i}; \
rm /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_NAT_PEL_EAS_1KG_chr${i}.sample ; \
cat header_sample.txt temp_REF_PELEAS_1KG_chr${i} > REF_NAT_PEL_EAS_1KG_chr${i}.sample; done


##################
# SIMU - preparo dos arquivos
#################

awk '{print 0, $1}' ../NAT_HGDP_SIMU_id > NAT_HGDP_SIMU_plink
awk '{print 0, $1}' ../IBS_1KG_SIMU > IBS_1KG_SIMU_plink
awk '{print 0, $1}' ../YRI_1KG_SIMU > YRI_1KG_SIMU_plink

for i in {1..22}; do \
/storage/atkinson/shared_resources/software/plink2/plink2 --vcf \
/storage/atkinson/home/u242335/1kg_hg38_phased/phased_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.gz \
--exclude paraexcluir_plink_${i}.txt \
--export vcf --out /storage/atkinson/home/u242335/1kg_hg38_phased/QC_phased_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i} ; done

for i in {1..22}; do \
bgzip -c /storage/atkinson/home/u242335/1kg_hg38_phased/QC_phased_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf > /storage/atkinson/home/u242335/1kg_hg38_phased/QC_phased_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.gz ; done


##### make hapsample files out of phased dataset with the individuals for SIMU. ran from /storage/atkinson/home/u242335/lai/simu-jointcall/wgs
for i in {1..22}; do \
/storage/atkinson/shared_resources/software/plink2/plink2 --vcf \
/storage/atkinson/home/u242335/1kg_hg38_phased/QC_phased_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.gz \
--keep /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/NAT_HGDP_SIMU_plink \
--export haps --out NAT_HGDP_SIMU_chr${i}; done

for i in {1..22}; do \
/storage/atkinson/shared_resources/software/plink2/plink2 --vcf \
/storage/atkinson/home/u242335/1kg_hg38_phased/QC_phased_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.gz \
--keep /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/IBS_1KG_SIMU_plink \
--export haps --out IBS_1KG_SIMU_chr${i}; done

for i in {1..22}; do \
/storage/atkinson/shared_resources/software/plink2/plink2 --vcf \
/storage/atkinson/home/u242335/1kg_hg38_phased/QC_phased_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.gz \
--keep /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/YRI_1KG_SIMU_plink \
--export haps --out YRI_1KG_SIMU_chr${i}; done

#### Make .phgeno files out of hapsample for SIMU
for i in {1..22}; do \
cut -d' ' -f 6- NAT_HGDP_SIMU_chr${i}.haps | sed 's/\s//g' > NAT${i}.phgeno; \
cut -d' ' -f 6- YRI_1KG_SIMU_chr${i}.haps | sed 's/\s//g' > AFR${i}.phgeno; \
cut -d' ' -f 6- IBS_1KG_SIMU_chr${i}.haps | sed 's/\s//g' > EUR${i}.phgeno ; done


##### make a VCF file out of phased dataset with the individuals for SIMU (for .snp file).
for i in {1..22}; do \
/storage/atkinson/shared_resources/software/plink2/plink2 --vcf \
/storage/atkinson/home/u242335/1kg_hg38_phased/QC_phased_gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.gz \
--keep /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/NAT_HGDP_SIMU_plink \
--const-fid 0 --max-alleles 2 \
--make-bed --out 1kg_hgdp_pBim2_${i}; done


for i in {1..22};do perl /storage/atkinson/home/u242335/lai/simu-jointcall/admix-simu-master/insert-map.pl \
1kg_hgdp_pBim2_${i}.bim /storage/atkinson/home/u242335/phasing/genetic_map/wgs/genetic_map_chr${i}_hg38_firstLine.txt > chr${i}.pos; \
awk -F' ' '{ print $1":"$4"_"$5"_"$6, $1, $3, $4, $5, $6 }' chr${i}.pos > 1kg_hgdp.chr${i}.snp; done

#2 way NAT/EUR
#70/30 PEL
#50/50 MAM
#05/95 EXEUR
#95/05 EXNAT

# 3 way NAT/EUR/AFR
# 33/33/34 = Mix
# 15/60/25 = Brasa

################################################
# SIMU 9 Gen, chr 1, MODELS wgs/gen9/
################################################

## Edit here the ancestry orders and which population will you test, remember to manually write the .dat file.
export  anc1=NAT
export  anc2=EUR
export  anc3=AFR
export  AdmixPOP=EXEUR

# editar o arquivo Baiano.dat

# running from /storage/atkinson/home/u242335/lai/simu-jointcall/gen12

# 2 way admixture
for i in {1..1}; do /storage/atkinson/home/u242335/lai/simu-jointcall/admix-simu-master/simu-mix.pl \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.dat \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/1kg_hgdp.chr${i}.snp $AdmixPOP${i} \
-$anc1 /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/$anc1${i}.phgeno \
-$anc2 /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/$anc2${i}.phgeno ; done

# 3 way admixture
for i in {1..1}; do /storage/atkinson/home/u242335/lai/simu-jointcall/admix-simu-master/simu-mix.pl \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.dat \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/1kg_hgdp.chr${i}.snp $AdmixPOP${i} \
-$anc1 /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/$anc1${i}.phgeno \
-$anc2 /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/$anc2${i}.phgeno \
-$anc3 /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/$anc3${i}.phgeno ; done


# convert .bp from the simulation to .hanc
for i in {1..1}; do /storage/atkinson/home/u242335/lai/simu-jointcall/admix-simu-master/bp2anc.pl $AdmixPOP${i}.bp > $AdmixPOP${i}.hanc; done

# hanc1
for i in {1..1}; do sed 's/./& /g' $AdmixPOP${i}.hanc > $AdmixPOP${i}.hanc1 ;done

# transpose hanc1
for i in {1..1};do python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < $AdmixPOP${i}.hanc1 > $AdmixPOP${i}.hanc2; done

head $AdmixPOP${i}.hanc2
rm *.hanc1

# From Elizabeth's scripts: Now I need to separate the .phgeno file to have a space between each character
for i in {1..1}; do sed 's/./& /g' ./$AdmixPOP${i}.phgeno > temp_Genos_chr${i}; done
for i in {1..1}; do awk -F' ' '{ print $2, $1, $4, $5, $6 }' ../1kg_hgdp.chr${i}.snp > temp_SNP_chr${i} ; done
for i in {1..1}; do paste temp_SNP_chr${i} temp_Genos_chr${i} > temp_admix_chr${i}.haps ; done
for i in {1..1}; do sed 's/\t/ /g' temp_admix_chr${i}.haps > $AdmixPOP.chr${i}.haps ; done
head $AdmixPOP.chr1.haps

rm temp_*

# Manually created the .sample file with fake names
cp ../gen12/Brasa.chr1.sample ./
cp ../gen12/Brasa.notref ./

# Mix sample
for i in {1..1}; do sed 's/Brasa/EXEUR/g' Brasa.chr${i}.sample > EXEUR.chr${i}.sample ; done

sed 's/Brasa/EXEUR/g' Brasa.notref > EXEUR.notref


##################
# RFmix
#################

######## MUDAR PRO ALVO From /storage/atkinson/home/u242335/lai/simu-jointcall/gen9

export  AdmixPOP=Mix

#### 2 way admixture

# shapeit output to rfmix. From /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9
for i in {1..1}; do \
python /storage/atkinson/home/u242335/lai/ancestry_pipeline-master/shapeit2rfmix.py \
--shapeit_hap_ref ../REF_NAT_HGDP_chr${i}.haps,../REF_IBS_1KG_chr${i}.haps \
--shapeit_hap_admixed ./$AdmixPOP.chr${i}.haps \
--shapeit_sample_ref ../REF_NAT_HGDP_chr${i}.sample,../REF_IBS_1KG_chr${i}.sample \
--shapeit_sample_admixed ./$AdmixPOP.chr${i}.sample \
--chr ${i} \
--genetic_map /storage/atkinson/home/u242335/phasing/genetic_map/wgs/genetic_map_chr${i}_hg38_firstLine.txt \
--ref_keep /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/REF_RFMIX-SIMU_NAT_EUR.ref \
--admixed_keep ./$AdmixPOP.notref \
--out $AdmixPOP.gen9_RefHGDPNAT_1kgIBS; done

# running RFMIx
cd /storage/atkinson/home/u242335/lai/RFMix_v1.5.4/

for i in {1..1}; do \
python RunRFMix.py \
-e 2 \
-w 0.2 \
-n 5 \
-G 9 \
--num-threads 12 \
--use-reference-panels-in-EM \
--forward-backward \
TrioPhased \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefHGDPNAT_1kgIBS_chr${i}.alleles \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefHGDPNAT_1kgIBS.classes \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefHGDPNAT_1kgIBS_chr${i}.snp_locations \
-o /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefHGDPNAT_1kgIBS_chr${i}.rfmix; done



# 3 way admixture
# shapeit output to rfmix. From /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9
for i in {1..1}; do \
python /storage/atkinson/home/u242335/lai/ancestry_pipeline-master/shapeit2rfmix.py \
--shapeit_hap_ref ../REF_NAT_HGDP_chr${i}.haps,../REF_IBS_1KG_chr${i}.haps,../REF_YRI_1KG_chr${i}.haps \
--shapeit_hap_admixed ./$AdmixPOP.chr${i}.haps \
--shapeit_sample_ref ../REF_NAT_HGDP_chr${i}.sample,../REF_IBS_1KG_chr${i}.sample,../REF_YRI_1KG_chr${i}.sample \
--shapeit_sample_admixed ./$AdmixPOP.chr${i}.sample \
--chr ${i} \
--genetic_map /storage/atkinson/home/u242335/phasing/genetic_map/wgs/genetic_map_chr${i}_hg38_firstLine.txt \
--ref_keep ../REF_RFMIX-SIMU_NAT_EUR_AFR.ref \
--admixed_keep ./$AdmixPOP.notref \
--out $AdmixPOP.gen9_RefHGDPNAT_1kgIBS_1kgYRI; done


# running RFMIx
cd /storage/atkinson/home/u242335/lai/RFMix_v1.5.4/

for i in {1..1}; do \
python RunRFMix.py \
-e 2 \
-w 0.2 \
-n 5 \
-G 9 \
--num-threads 20 \
--use-reference-panels-in-EM \
--forward-backward \
TrioPhased \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.alleles \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefHGDPNAT_1kgIBS_1kgYRI.classes \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.snp_locations \
-o /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix; done

### ACCURACY PREP ###

export  AdmixPOP=Brasa

# preparar Simu Real para verificar a acurácia. De /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9

for i in {1..1}; do awk '{print $1, $4}' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/chr${i}.pos > chr${i}.posmini; \
paste chr${i}.posmini ./Brasa${i}.hanc2 | sed 's/\t/ /g' > SimuReal_Brasa_chr${i}.AncCalls; done


#####
# Brasa
for i in {1..1}; do sed 's/1/0/g' $AdmixPOP.gen9_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt > Brasa_1_Lat2_${i}; done
for i in {1..1}; do sed 's/2/1/g' Brasa_1_Lat2_${i} > Brasa_1_Lat2.1_${i}; done
for i in {1..1}; do sed 's/3/2/g' Brasa_1_Lat2.1_${i} > Brasa_1_Lat2.2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' $AdmixPOP.gen9_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.map | sed 's/:.*//g' > Brasa_1_mappin_${i}_Lat; done
for i in {1..1}; do paste Brasa_1_mappin_${i}_Lat Brasa_1_Lat2.2_${i} | sed 's/\t/ /g' > Brasa_HGDP_Lat3_${i}; done


################# PEL IBS YRI #######################

######## MUDAR PRO ALVO From /storage/atkinson/home/u242335/lai/simu-jointcall/gen9

export  AdmixPOP=EXNAT

# 2 way admixture
# shapeit output to rfmix. From /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9
for i in {1..1}; do \
python /storage/atkinson/home/u242335/lai/ancestry_pipeline-master/shapeit2rfmix.py \
--shapeit_hap_ref ../REF_NAT_PEL_1KG_chr${i}.haps,../REF_IBS_1KG_chr${i}.haps \
--shapeit_hap_admixed ./$AdmixPOP.chr${i}.haps \
--shapeit_sample_ref ../REF_NAT_PEL_1KG_chr${i}.sample,../REF_IBS_1KG_chr${i}.sample \
--shapeit_sample_admixed ./$AdmixPOP.chr${i}.sample \
--chr ${i} \
--genetic_map /storage/atkinson/home/u242335/phasing/genetic_map/wgs/genetic_map_chr${i}_hg38_firstLine.txt \
--ref_keep ../REF_RFMIX-SIMU_PEL_EUR.ref \
--admixed_keep ./$AdmixPOP.notref \
--out $AdmixPOP.gen9_RefPEL_1kgIBS; done


# running RFMIx
cd /storage/atkinson/home/u242335/lai/RFMix_v1.5.4/

for i in {1..1}; do \
python RunRFMix.py \
-e 2 \
-w 0.2 \
-n 5 \
-G 9 \
--num-threads 20 \
--use-reference-panels-in-EM \
--forward-backward \
TrioPhased \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefPEL_1kgIBS_chr${i}.alleles \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefPEL_1kgIBS.classes \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefPEL_1kgIBS_chr${i}.snp_locations \
-o /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefPEL_1kgIBS_chr${i}.rfmix; done





# 3 way admixture
# shapeit output to rfmix. From /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9
for i in {1..1}; do \
python /storage/atkinson/home/u242335/lai/ancestry_pipeline-master/shapeit2rfmix.py \
--shapeit_hap_ref ../REF_NAT_PEL_1KG_chr${i}.haps,../REF_IBS_1KG_chr${i}.haps,../REF_YRI_1KG_chr${i}.haps \
--shapeit_hap_admixed ./$AdmixPOP.chr${i}.haps \
--shapeit_sample_ref ../REF_NAT_PEL_1KG_chr${i}.sample,../REF_IBS_1KG_chr${i}.sample,../REF_YRI_1KG_chr${i}.sample \
--shapeit_sample_admixed ./$AdmixPOP.chr${i}.sample \
--chr ${i} \
--genetic_map /storage/atkinson/home/u242335/phasing/genetic_map/wgs/genetic_map_chr${i}_hg38_firstLine.txt \
--ref_keep ../REF_RFMIX-SIMU_PEL_EUR_AFR.ref \
--admixed_keep ./$AdmixPOP.notref \
--out $AdmixPOP.gen9_RefPEL_1kgIBS_1kgYRI; done


# running RFMIx
cd /storage/atkinson/home/u242335/lai/RFMix_v1.5.4/

for i in {1..1}; do \
python RunRFMix.py \
-e 2 \
-w 0.2 \
-n 5 \
-G 9 \
--num-threads 20 \
--use-reference-panels-in-EM \
--forward-backward \
TrioPhased \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefPEL_1kgIBS_1kgYRI_chr${i}.alleles \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefPEL_1kgIBS_1kgYRI.classes \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefPEL_1kgIBS_1kgYRI_chr${i}.snp_locations \
-o /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefPEL_1kgIBS_1kgYRI_chr${i}.rfmix; done

### ACCURACY PREP ###

export  AdmixPOP=Brasa

# Brasa
for i in {1..1}; do sed 's/1/0/g' $AdmixPOP.gen9_RefPEL_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt > Brasa_1_Lat2_${i}; done
for i in {1..1}; do sed 's/2/1/g' Brasa_1_Lat2_${i} > Brasa_1_Lat2.1_${i}; done
for i in {1..1}; do sed 's/3/2/g' Brasa_1_Lat2.1_${i} > Brasa_1_Lat2.2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' $AdmixPOP.gen9_RefPEL_1kgIBS_1kgYRI_chr${i}.map | sed 's/:.*//g' > Brasa_1_mappin_${i}_Lat; done
for i in {1..1}; do paste Brasa_1_mappin_${i}_Lat Brasa_1_Lat2.2_${i} | sed 's/\t/ /g' > Brasa_PEL_Lat3_${i}; done


################# PEL-EAS IBS YRI #######################

######## MUDAR PRO ALVO From /storage/atkinson/home/u242335/lai/simu-jointcall/gen9

export  AdmixPOP=EXNAT

# 2 way admixture
# shapeit output to rfmix. From /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9
for i in {1..1}; do \
python /storage/atkinson/home/u242335/lai/ancestry_pipeline-master/shapeit2rfmix.py \
--shapeit_hap_ref ../REF_NAT_PEL_EAS_1KG_chr${i}.haps,../REF_IBS_1KG_chr${i}.haps \
--shapeit_hap_admixed ./$AdmixPOP.chr${i}.haps \
--shapeit_sample_ref ../REF_NAT_PEL_EAS_1KG_chr${i}.sample,../REF_IBS_1KG_chr${i}.sample \
--shapeit_sample_admixed ./$AdmixPOP.chr${i}.sample \
--chr ${i} \
--genetic_map /storage/atkinson/home/u242335/phasing/genetic_map/wgs/genetic_map_chr${i}_hg38_firstLine.txt \
--ref_keep ../REF_RFMIX-SIMU_REF_PEL_EAS_EUR.ref \
--admixed_keep ./$AdmixPOP.notref \
--out $AdmixPOP.gen9_RefPEL_EAS_1kgIBS; done


# running RFMIx
cd /storage/atkinson/home/u242335/lai/RFMix_v1.5.4/

for i in {1..1}; do \
python RunRFMix.py \
-e 2 \
-w 0.2 \
-n 5 \
-G 9 \
--num-threads 20 \
--use-reference-panels-in-EM \
--forward-backward \
TrioPhased \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefPEL_EAS_1kgIBS_chr${i}.alleles \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefPEL_EAS_1kgIBS.classes \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefPEL_EAS_1kgIBS_chr${i}.snp_locations \
-o /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefPEL_EAS_1kgIBS_chr${i}.rfmix; done




# 3 way admixture
# shapeit output to rfmix. From /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9
for i in {1..1}; do \
python /storage/atkinson/home/u242335/lai/ancestry_pipeline-master/shapeit2rfmix.py \
--shapeit_hap_ref ../REF_NAT_PEL_EAS_1KG_chr${i}.haps,../REF_IBS_1KG_chr${i}.haps,../REF_YRI_1KG_chr${i}.haps \
--shapeit_hap_admixed ./$AdmixPOP.chr${i}.haps \
--shapeit_sample_ref ../REF_NAT_PEL_EAS_1KG_chr${i}.sample,../REF_IBS_1KG_chr${i}.sample,../REF_YRI_1KG_chr${i}.sample \
--shapeit_sample_admixed ./$AdmixPOP.chr${i}.sample \
--chr ${i} \
--genetic_map /storage/atkinson/home/u242335/phasing/genetic_map/wgs/genetic_map_chr${i}_hg38_firstLine.txt \
--ref_keep ../REF_RFMIX-SIMU_REF_PEL_EAS_EUR_AFR.ref \
--admixed_keep ./$AdmixPOP.notref \
--out $AdmixPOP.gen9_RefPEL_EAS_1kgIBS_1kgYRI; done


# running RFMIx
cd /storage/atkinson/home/u242335/lai/RFMix_v1.5.4/

for i in {1..1}; do \
python RunRFMix.py \
-e 2 \
-w 0.2 \
-n 5 \
-G 9 \
--num-threads 20 \
--use-reference-panels-in-EM \
--forward-backward \
TrioPhased \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefPEL_EAS_1kgIBS_1kgYRI_chr${i}.alleles \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefPEL_EAS_1kgIBS_1kgYRI.classes \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefPEL_EAS_1kgIBS_1kgYRI_chr${i}.snp_locations \
-o /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen9/$AdmixPOP.gen9_RefPEL_EAS_1kgIBS_1kgYRI_chr${i}.rfmix; done

### ACCURACY PREP ###

export  AdmixPOP=Brasa

# Brasa
for i in {1..1}; do sed 's/1/0/g' $AdmixPOP.gen9_RefPEL_EAS_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt > Brasa_1_Lat2_${i}; done
for i in {1..1}; do sed 's/2/1/g' Brasa_1_Lat2_${i} > Brasa_1_Lat2.1_${i}; done
for i in {1..1}; do sed 's/3/2/g' Brasa_1_Lat2.1_${i} > Brasa_1_Lat2.2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' $AdmixPOP.gen9_RefPEL_EAS_1kgIBS_1kgYRI_chr${i}.map | sed 's/:.*//g' > Brasa_1_mappin_${i}_Lat; done
for i in {1..1}; do paste Brasa_1_mappin_${i}_Lat Brasa_1_Lat2.2_${i} | sed 's/\t/ /g' > Brasa_PELEAS_Lat3_${i}; done



################################################
# SIMU 12 Gen, all chromosomes wgs/gen12/
################################################

## Edit here the ancestry orders and which population will you test, remember to manually write the .dat file.
export  anc1=NAT
export  anc2=EUR
export  anc3=AFR
export  AdmixPOP=Brasa

# editar o arquivo Baiano.dat

# running from /storage/atkinson/home/u242335/lai/simu-jointcall/gen12

# 3 way admixture
for i in {7..22}; do /storage/atkinson/home/u242335/lai/simu-jointcall/admix-simu-master/simu-mix.pl \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen12/$AdmixPOP.dat \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/1kg_hgdp.chr${i}.snp $AdmixPOP${i} \
-$anc1 /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/$anc1${i}.phgeno \
-$anc2 /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/$anc2${i}.phgeno \
-$anc3 /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/$anc3${i}.phgeno ; done

# convert .bp from the simulation to .hanc
for i in {1..22}; do  /storage/atkinson/home/u242335/lai/simu-jointcall/admix-simu-master/bp2anc.pl $AdmixPOP${i}.bp > $AdmixPOP${i}.hanc; done

# hanc1
for i in {1..22}; do sed 's/./& /g' $AdmixPOP${i}.hanc > $AdmixPOP${i}.hanc1 ;done

# transpose hanc1
for i in {1..22};do python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < $AdmixPOP${i}.hanc1 > $AdmixPOP${i}.hanc2; done

head $AdmixPOP${i}.hanc2
rm *.hanc1

# From Elizabeth's scripts: Now I need to separate the .phgeno file to have a space between each character
for i in {1..22}; do sed 's/./& /g' ./$AdmixPOP${i}.phgeno > temp_Genos_chr${i}; done
for i in {1..22}; do awk -F' ' '{ print $2, $1, $4, $5, $6 }' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/1kg_hgdp.chr${i}.snp > temp_SNP_chr${i} ; done
for i in {1..22}; do paste temp_SNP_chr${i} temp_Genos_chr${i} > temp_admix_chr${i}.haps ; done
for i in {1..22}; do sed 's/\t/ /g' temp_admix_chr${i}.haps > $AdmixPOP.chr${i}.haps ; done
head $AdmixPOP.chr1.haps

rm temp_*

# Manually created the .sample file with fake names
# Jessica cuidado aqui
sed 's/PEL/Brasa/g' ../../Jessica/admixsimu/PEL.sample.txt > $AdmixPOP.sample.txt
for i in {1..22}; do cp $AdmixPOP.sample $AdmixPOP.chr${i}.sample; done

##################
# RFmix
#################

######## MUDAR PRO ALVO From /storage/atkinson/home/u242335/lai/simu-jointcall/gen12

# alvos, arquivo criado na mão apenas com nomes ficticios da nossa amostra alvo. editado p tirar o cabecalho
awk '{print $2}' ./Brasa.sample > Brasa.notref
nano Brasa.notref

# OPTIONAL for other runs
#cp Brasa.notref $AdmixPOP.no
#sed 's/Brasa/Baiano/g' $AdmixPOP.no > $AdmixPOP.notref
export  pasta=/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen12
export  AdmixPOP=Brasa

# shapeit output to rfmix. From /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen12
for i in {1..22}; do \
python /storage/atkinson/home/u242335/lai/ancestry_pipeline-master/shapeit2rfmix.py \
--shapeit_hap_ref ../REF_NAT_HGDP_chr${i}.haps,../REF_IBS_1KG_chr${i}.haps,../REF_YRI_1KG_chr${i}.haps \
--shapeit_hap_admixed ./$AdmixPOP.chr${i}.haps \
--shapeit_sample_ref ../REF_NAT_HGDP_chr${i}.sample,../REF_IBS_1KG_chr${i}.sample,../REF_YRI_1KG_chr${i}.sample \
--shapeit_sample_admixed ./$AdmixPOP.chr${i}.sample \
--chr ${i} \
--genetic_map /storage/atkinson/home/u242335/phasing/genetic_map/wgs/genetic_map_chr${i}_hg38_firstLine.txt \
--ref_keep ../REF_RFMIX-SIMU_NAT_EUR_AFR.ref \
--admixed_keep ./$AdmixPOP.notref \
--out Brasa_gen12_RefHGDPNAT_1kgIBS_1kgYRI; done

# running RFMIx
cd /storage/atkinson/home/u242335/lai/RFMix_v1.5.4/

for i in {1..10}; do \
python RunRFMix.py \
-e 2 \
-w 0.2 \
-n 5 \
-G 12 \
--num-threads 12 \
--use-reference-panels-in-EM \
--forward-backward \
TrioPhased \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen12/Brasa_gen12_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.alleles \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen12/Brasa_gen12_RefHGDPNAT_1kgIBS_1kgYRI.classes \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen12/Brasa_gen12_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.snp_locations \
-o /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen12/Brasa_gen12_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix; done


### ACCURACY PREP ###

# preparar Simu Real para verificar a acurácia. De /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen12


#####
# Brasa
for i in {1..22}; do sed 's/1/0/g' Brasa_gen12_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt > Brasa_1_Lat2_${i}; done
for i in {1..22}; do sed 's/2/1/g' Brasa_1_Lat2_${i} > Brasa_1_Lat2.1_${i}; done
for i in {1..22}; do sed 's/3/2/g' Brasa_1_Lat2.1_${i} > Brasa_1_Lat2.2_${i}; done
for i in {1..22}; do awk '{print $1, $3}' Brasa_gen12_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.map | sed 's/:.*//g' > Brasa_1_mappin_${i}_Lat; done
for i in {1..22}; do paste Brasa_1_mappin_${i}_Lat Brasa_1_Lat2.2_${i} | sed 's/\t/ /g' > Brasa_1_Lat3_${i}; done


## whole genome combining

for i in {1..22}; do sed 's/1/0/g' Brasa_GSAIMPUT_gen12_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt > Brasa_1_Lat2_${i}; done
for i in {1..22}; do sed 's/2/1/g' Brasa_1_Lat2_${i} > Brasa_1_Lat2.1_${i}; done
for i in {1..22}; do sed 's/3/2/g' Brasa_1_Lat2.1_${i} > Brasa_1_Lat2.2_${i}; done
for i in {1..22}; do awk '{print $1, $3}' Brasa_GSAIMPUT_gen12_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.map | sed 's/:.*//g' > Brasa_1_mappin_${i}_Lat; done
for i in {1..22}; do paste Brasa_1_mappin_${i}_Lat Brasa_1_Lat2.2_${i} | sed 's/\t/ /g' > Brasa_1_Lat3_${i}; done



#simu
cat SimuReal_Brasa_chr1.AncCalls SimuReal_Brasa_chr2.AncCalls SimuReal_Brasa_chr3.AncCalls SimuReal_Brasa_chr4.AncCalls \
SimuReal_Brasa_chr5.AncCalls SimuReal_Brasa_chr6.AncCalls SimuReal_Brasa_chr7.AncCalls SimuReal_Brasa_chr8.AncCalls \
SimuReal_Brasa_chr9.AncCalls SimuReal_Brasa_chr10.AncCalls SimuReal_Brasa_chr11.AncCalls SimuReal_Brasa_chr12.AncCalls \
SimuReal_Brasa_chr13.AncCalls SimuReal_Brasa_chr14.AncCalls SimuReal_Brasa_chr15.AncCalls SimuReal_Brasa_chr16.AncCalls \
SimuReal_Brasa_chr17.AncCalls SimuReal_Brasa_chr18.AncCalls SimuReal_Brasa_chr19.AncCalls SimuReal_Brasa_chr20.AncCalls \
SimuReal_Brasa_chr21.AncCalls SimuReal_Brasa_chr22.AncCalls > SimuReal_Brasa_genome.AncCalls 

# lai
# cuidado aqui com o double carrot
for i in {1..22}; do \
paste -d' ' <(cut -d ' ' -f -2 Brasa_1_Lat3_${i}) <(cut -d ' ' -f 373- Brasa_1_Lat3_${i}) >> Brasa_HGDP_Lat3_genome; done


################################################
# SIMU 17 Gen, chr 1, BRASA wgs/gen17/
################################################

## Edit here the ancestry orders and which population will you test, remember to manually write the .dat file.
export  anc1=NAT
export  anc2=EUR
export  anc3=AFR
export  AdmixPOP=Mix

# editar o arquivo Baiano.dat

# running from /storage/atkinson/home/u242335/lai/simu-jointcall/gen12

# 3 way admixture
for i in {1..1}; do /storage/atkinson/home/u242335/lai/simu-jointcall/admix-simu-master/simu-mix.pl \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen17/$AdmixPOP.dat \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/1kg_hgdp.chr${i}.snp $AdmixPOP${i} \
-$anc1 /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/$anc1${i}.phgeno \
-$anc2 /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/$anc2${i}.phgeno \
-$anc3 /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/$anc3${i}.phgeno ; done

# convert .bp from the simulation to .hanc
for i in {1..1}; do /storage/atkinson/home/u242335/lai/simu-jointcall/admix-simu-master/bp2anc.pl $AdmixPOP${i}.bp > $AdmixPOP${i}.hanc; done

# hanc1
for i in {1..1}; do sed 's/./& /g' $AdmixPOP${i}.hanc > $AdmixPOP${i}.hanc1 ;done

# transpose hanc1
for i in {1..1};do python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < $AdmixPOP${i}.hanc1 > $AdmixPOP${i}.hanc2; done

head $AdmixPOP${i}.hanc2
rm *.hanc1

# From Elizabeth's scripts: Now I need to separate the .phgeno file to have a space between each character
for i in {1..1}; do sed 's/./& /g' ./$AdmixPOP${i}.phgeno > temp_Genos_chr${i}; done
for i in {1..1}; do awk -F' ' '{ print $2, $1, $4, $5, $6 }' /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/1kg_hgdp.chr${i}.snp > temp_SNP_chr${i} ; done
for i in {1..1}; do paste temp_SNP_chr${i} temp_Genos_chr${i} > temp_admix_chr${i}.haps ; done
for i in {1..1}; do sed 's/\t/ /g' temp_admix_chr${i}.haps > $AdmixPOP.chr${i}.haps ; done
head $AdmixPOP.chr1.haps

rm temp_*

# Manually created the .sample file with fake names
cp ../gen12/Brasa.chr1.sample ./
cp ../gen12/Brasa.notref ./

# Mix sample
for i in {1..1}; do sed 's/Brasa/Mix/g' Brasa.chr${i}.sample > Mix.chr${i}.sample ; done

sed 's/Brasa/Mix/g' Brasa.notref > Mix.notref

##################
# RFmix
#################

######## MUDAR PRO ALVO From /storage/atkinson/home/u242335/lai/simu-jointcall/gen17

export  AdmixPOP=Brasa

# shapeit output to rfmix. From /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen17
for i in {1..1}; do \
python /storage/atkinson/home/u242335/lai/ancestry_pipeline-master/shapeit2rfmix.py \
--shapeit_hap_ref ../REF_NAT_HGDP_chr${i}.haps,../REF_IBS_1KG_chr${i}.haps,../REF_YRI_1KG_chr${i}.haps \
--shapeit_hap_admixed ./$AdmixPOP.chr${i}.haps \
--shapeit_sample_ref ../REF_NAT_HGDP_chr${i}.sample,../REF_IBS_1KG_chr${i}.sample,../REF_YRI_1KG_chr${i}.sample \
--shapeit_sample_admixed ./$AdmixPOP.chr${i}.sample \
--chr ${i} \
--genetic_map /storage/atkinson/home/u242335/phasing/genetic_map/wgs/genetic_map_chr${i}_hg38_firstLine.txt \
--ref_keep ../REF_RFMIX-SIMU_NAT_EUR_AFR.ref \
--admixed_keep ./$AdmixPOP.notref \
--out $AdmixPOP.gen17_RefHGDPNAT_1kgIBS_1kgYRI; done

# running RFMIx
cd /storage/atkinson/home/u242335/lai/RFMix_v1.5.4/

for i in {1..1}; do \
python RunRFMix.py \
-e 2 \
-w 0.2 \
-n 5 \
-G 17 \
--num-threads 20 \
--use-reference-panels-in-EM \
--forward-backward \
TrioPhased \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen17/$AdmixPOP.gen17_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.alleles \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen17/$AdmixPOP.gen17_RefHGDPNAT_1kgIBS_1kgYRI.classes \
/storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen17/$AdmixPOP.gen17_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.snp_locations \
-o /storage/atkinson/home/u242335/lai/simu-jointcall/wgs/gen17/$AdmixPOP.gen17_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix; done

