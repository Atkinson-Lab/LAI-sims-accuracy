
# 1kg .sample file from https://mathgen.stats.ox.ac.uk/impute/1000GP%20Phase%203%20haplotypes%206%20October%202014.html

##### Select which indivs go to simulation and which go to rfmix (same as used before)
### ADMIX-SIMU
# 31 NAT (half of Surui, Karitiana, Maya, Pima) = SIM
egrep "Surui|Karitiana|Maya|Pima|Colombian" BHRC_KIDS > NAT_HGDP_SIMU
nano NAT_HGDP_SIMU

awk '{print $1}' NAT_HGDP_SIMU > NAT_HGDP_SIMU_id

# 30 IBS = SIM
grep IBS ./1000GP_Phase3/1000GP_Phase3.sample | awk '{print $1}' | head -n 30 > IBS_1KG_SIMU

# 30 YRI = SIM
grep IBS ./1000GP_Phase3/1000GP_Phase3.sample | awk '{print $1}' | head -n 30 > YRI_1KG_SIMU


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


##### make hapsample files out of phased dataset with the individuals for SIMU. ran from /storage/atkinson/home/u242335/lai/simu-jointcall
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


## Edit here the ancestry orders and which population will you test, remember to manually write the .dat file.
export  anc1=NAT
export  anc2=EUR
export  anc3=AFR
export  AdmixPOP=Brasa

# 33/33/34 = Mix
# 15/60/25 = Brasa

# editar o arquivo Baiano.dat

# running from /storage/atkinson/home/u242335/lai/simu-jointcall/gen12

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

# Manually created the .sample file with fake names
# Jessica cuidado aqui
sed 's/PEL/Brasa/g' ../../Jessica/admixsimu/PEL.sample.txt > $AdmixPOP.sample.txt
for i in {1..22}; do cp $AdmixPOP.sample.txt $AdmixPOP.chr${i}.sample; done




##################
# RFmix
#################

#from /storage/atkinson/home/u242335/lai/simu-jointcall/
##### make hapsample files out of phased dataset with the individuals for RFMix REF

for i in {1..22}; do \
/storage/atkinson/home/u242335/phasing/shapeit -convert --input-haps ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--include-ind NAT_HGDP_RFMIX \
--output-haps REF_NAT_HGDP_chr${i}.haps REF_NAT_HGDP_chr${i}.sample ; done

for i in {1..22}; do \
/storage/atkinson/home/u242335/phasing/shapeit -convert --input-haps ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--include-ind IBS_1KG_RFMIX \
--output-haps REF_IBS_1KG_chr${i}.haps REF_IBS_1KG_chr${i}.sample ; done

for i in {1..22}; do \
/storage/atkinson/home/u242335/phasing/shapeit -convert --input-haps ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--include-ind YRI_1KG_RFMIX \
--output-haps REF_YRI_1KG_chr${i}.haps REF_YRI_1KG_chr${i}.sample ; done


# printar as amostras ref em um mesmo arquivo.
cat NAT_HGDP_RFMIX IBS_1KG_RFMIX YRI_1KG_RFMIX > ./REF_RFMIX-SIMU_NAT_EUR_AFR.ref


# alvos, arquivo criado na mão apenas com nomes ficticios da nossa amostra alvo. editado p tirar o cabecalho
awk '{print $2}' ./gen12/Brasa.sample > Brasa.notref
nano Brasa.notref

# OPTIONAL for other runs
#cp Brasa.notref $AdmixPOP.no
#sed 's/Brasa/Baiano/g' $AdmixPOP.no > $AdmixPOP.notref


# shapeit output to rfmix. From /storage/atkinson/home/u242335/lai/simu-jointcall/gen12
for i in {1..22}; do \
python ../../ancestry_pipeline-master/shapeit2rfmix.py \
--shapeit_hap_ref ../REF_NAT_HGDP_chr${i}.haps,../REF_IBS_1KG_chr${i}.haps,../REF_YRI_1KG_chr${i}.haps \
--shapeit_hap_admixed ./Brasa.chr${i}.haps \
--shapeit_sample_ref ../REF_NAT_HGDP_chr${i}.sample,../REF_IBS_1KG_chr${i}.sample,../REF_YRI_1KG_chr${i}.sample \
--shapeit_sample_admixed ./Brasa.chr${i}.sample \
--chr ${i} \
--genetic_map /storage/atkinson/home/u242335/phasing/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt \
--ref_keep ../REF_RFMIX-SIMU_NAT_EUR_AFR.ref \
--admixed_keep ./Brasa.notref \
--out Brasa_gen12_RefHGDPNAT_1kgIBS_1kgYRI; done

# running RFMIX
cd /storage/atkinson/home/u242335/lai/RFMix_v1.5.4/

for i in {19..22}; do \
python RunRFMix.py \
-e 2 \
-w 0.2 \
-n 5 \
-G 12 \
--num-threads 12 \
--use-reference-panels-in-EM \
--forward-backward \
TrioPhased \
/storage/atkinson/home/u242335/lai/simu-jointcall/gen12/Brasa_gen12_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.alleles \
/storage/atkinson/home/u242335/lai/simu-jointcall/gen12/Brasa_gen12_RefHGDPNAT_1kgIBS_1kgYRI.classes \
/storage/atkinson/home/u242335/lai/simu-jointcall/gen12/Brasa_gen12_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.snp_locations \
-o /storage/atkinson/home/u242335/lai/simu-jointcall/gen12/Brasa_gen12_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix; done

### ACCURACY PREP ###

# preparar Simu Real para verificar a acurácia. De /storage/atkinson/home/u242335/lai/simu-jointcall/gen12

for i in {1..22}; do awk '{print $1, $4}' /storage/atkinson/home/u242335/lai/simu-jointcall/chr${i}.pos > chr${i}.posmini; \
paste chr${i}.posmini ./Brasa${i}.hanc2 | sed 's/\t/ /g' > SimuReal_Brasa_chr${i}.AncCalls; done


#####
# Brasa
for i in {1..22}; do sed 's/1/0/g' Brasa_gen12_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt > Brasa_1_Lat2_${i}; done
for i in {1..22}; do sed 's/2/1/g' Brasa_1_Lat2_${i} > Brasa_1_Lat2.1_${i}; done
for i in {1..22}; do sed 's/3/2/g' Brasa_1_Lat2.1_${i} > Brasa_1_Lat2.2_${i}; done
for i in {1..22}; do awk '{print $1, $3}' Brasa_gen12_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.map | sed 's/:.*//g' > Brasa_1_mappin_${i}_Lat; done
for i in {1..22}; do paste Brasa_1_mappin_${i}_Lat Brasa_1_Lat2.2_${i} | sed 's/\t/ /g' > Brasa_1_Lat3_${i}; done

