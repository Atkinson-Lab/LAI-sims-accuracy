#lendo os simureal
listBrasasimureal  <- c() # creates a list
simureal <- list() # creates a list
simureal[1] <- "vazio"
for (i in 1:22){
  listBrasasimureal[i] <- paste("SimuReal_Brasa_chr", i, ".AncCalls", sep='')
}

for (k in 1:length(listBrasasimureal)){
  simureal[[k]] <- read.table(listBrasasimureal[k], h=F)
}

for (k in 1:length(simureal)){
  for (i in c(3:62)){
    simureal[[k]][,i] <- as.numeric(gsub(1,3,simureal[[k]][,i]))
  }}

# lendo os rfmixed
listBrasarfmix  <- c() # creates a list
ldf <- list() # creates a list
ldf[1] <- "vazio"

for (i in 1:22){
  listBrasarfmix[i] <- paste("Brasa_1_Lat3_", i, sep="")
}
for (k in 1:length(listBrasarfmix)){
  ldf[[k]] <- read.table(listBrasarfmix[k], h=F)
}

for (k in 1:length(ldf)){
  ldf[[k]] <- ldf[[k]][,c(2,1,(ncol(ldf[[k]])-59):ncol(ldf[[k]]))]
  for (i in c(3:62)){
    ldf[[k]][,i] <- as.numeric(gsub(1,3,ldf[[k]][,i]))
  }}

simureal2 <- list()
simureal2[1] <- "vazio"
for (k in 1:length(ldf)){
  simureal2[[k]] <- simureal[[k]][simureal[[k]]$V2 %in% ldf[[k]]$V1,]
}

#preservar espaco
rm(simureal)

install.packages("ggplot2")
library(ggplot2)
a_Brasa <- list()
b_Brasa <- list()
c_Brasa <- list()
d_Brasa <- list()
e_Brasa <- list()
f_Brasa <- list()
g_Brasa <- list()

mean_AFR_Brasa <- list()
sd_AFR_Brasa <- list()
mean_EUR_Brasa <- list()
sd_EUR_Brasa <- list()
mean_NAT_Brasa <- list()
sd_NAT_Brasa <- list()

a_Brasa[1] <- "vazio"
b_Brasa[1] <- "vazio"
c_Brasa[1] <- "vazio"
d_Brasa[1] <- "vazio"
e_Brasa[1] <- "vazio"
f_Brasa[1] <- "vazio"
g_Brasa[1] <- "vazio"

mean_AFR_Brasa[1] <- "vazio"
sd_AFR_Brasa[1] <- "vazio"
mean_EUR_Brasa[1] <- "vazio"
sd_EUR_Brasa[1] <- "vazio"
mean_NAT_Brasa[1] <- "vazio"
sd_NAT_Brasa[1] <- "vazio"


for (k in 1:22){
  a_Brasa[[k]] <- simureal2[[k]] == ldf[[k]]
  b_Brasa[[k]] <- as.data.frame(1*a_Brasa[[k]])
  c_Brasa[[k]] <-  colSums(b_Brasa[[k]])/dim(b_Brasa[[k]])[1]
  d_Brasa[[k]] <- simureal2[[k]] + ldf[[k]]
  
  # Accu AFR
  e_Brasa[[k]] <- colSums(d_Brasa[[k]]==4)/colSums(simureal2[[k]]==2)
  mean_AFR_Brasa[[k]] <- mean(e_Brasa[[k]][c(3:62)])
  sd_AFR_Brasa[[k]] <- sd(e_Brasa[[k]][c(3:62)])
  
  # Accu EUR
  f_Brasa[[k]] <- colSums(d_Brasa[[k]]==6)/colSums(simureal2[[k]]==3)
  mean_EUR_Brasa[[k]] <- mean(f_Brasa[[k]][3:62])
  sd_EUR_Brasa[[k]] <- sd(f_Brasa[[k]][3:62])
  
  # Accu NAT
  g_Brasa[[k]] <- colSums(d_Brasa[[k]]==0)/colSums(simureal2[[k]]==0)
  mean_NAT_Brasa[[k]] <- mean(g_Brasa[[k]][3:62])
  sd_NAT_Brasa[[k]] <- sd(g_Brasa[[k]][3:62])
}

rm(ldf)

#GRAFICOS
TabBrasa <- list()
TabBrasa[1] <- "vazio"

for (k in 1:22){
  TabBrasa[[k]] <- as.data.frame(c(e_Brasa[[k]][3:62],f_Brasa[[k]][3:62],g_Brasa[[k]][3:62]))
  TabBrasa[[k]]$Variable[1:60] <- c("AFR")
  TabBrasa[[k]]$Variable[61:120] <-c("EUR")
  TabBrasa[[k]]$Variable[121:180] <-c("NAT")
  TabBrasa[[k]]$Ref <- "12 Gen - 15% NAT / 60% EUR / 25% AFR"
  names(TabBrasa[[k]])[1] <- "Accuracy" 
  rownames(TabBrasa[[k]])[1:60] <- paste("AFR", c(1:60), sep="")
  rownames(TabBrasa[[k]])[61:120] <- paste("EUR", c(1:60), sep="")
  rownames(TabBrasa[[k]])[121:180] <- paste("NAT", c(1:60), sep="")
}

anc = c("AFR", "EUR", "NAT")
vivid = c("#5D69B1", "#FCA315", "#1BB6AF")

library(ggthemr)
ggthemr('fresh', layout = "clean")

i=1
ggplot(TabBrasa[[i]], aes(x=Ref, fill=Variable, y=Accuracy)) +
  ggtitle(paste("LAI sensitivity by ancestry for GSA -","chromosome",i,sep=" ")) +
  geom_boxplot() +
  scale_fill_manual(breaks = anc,
                    values = vivid, 
                    name = "Ancestry", 
                    labels = c("African", "European", "Native American")) +
  xlab("Admixture model") + ylab("True Positive Rate") + 
  labs(fill = "Ancestry", size = 14) +
  coord_cartesian(ylim=c(0.30,1)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), 
        title = element_text(size=12), legend.text=element_text(size=12)) +
  stat_summary(aes(label = , 1))

### GETTING OUTLIER NAMES
library(car)

Boxplot(TabBrasa[[1]]$Accuracy ~ TabBrasa[[1]]$Variable, labels = rownames(TabBrasa[[1]]))

library(dplyr)

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

library(tibble)

dat <- c()
dat[[1]] <- "vazio"
tab_narm <- c()
tab_narm[[1]] <- "vazio"

for (k in 1:22){
  tab_narm[[k]] <- na.omit(TabBrasa[[k]], Accuracy)
  dat[[k]] <- tab_narm[[k]] %>% tibble::rownames_to_column(var="outlier") %>% group_by(Variable) %>% mutate(is_outlier=ifelse(is_outlier(Accuracy), Accuracy, as.numeric(NA)))
  dat[[k]]$outlier[which(is.na(dat[[k]]$is_outlier))] <- as.numeric(NA)
}

i=6

ggplot(dat[[i]], aes(x=Variable, fill=Variable, y=Accuracy)) +
  ggtitle(paste("LAI sensitivity by ancestry portion for GSA -","chromosome",i,sep=" ")) +
  geom_boxplot() +
  scale_fill_manual(breaks = anc,
                    values = vivid, 
                    name = "Ancestry", 
                    labels = c("African", "European", "Native American")) +
  xlab("Admixture model") + ylab("True Positive Rate") + 
  labs(fill = "Ancestry", size = 14) +
  coord_cartesian(ylim=c(0.30,1)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), 
        title = element_text(size=12), legend.text=element_text(size=12)) + 
  geom_text(aes(label=outlier),na.rm=TRUE, hjust=1.25) 


#optional:
#library(ggbeeswarm)
#add + geom_quasirandom() to ggplot

####### GLOBAL ACC ######

SimuReal_Brasa <- read.table("SimuReal_Brasa_genome.AncCalls", header = F)
for (i in c(3:62)){
  SimuReal_Brasa[,i] <- as.numeric(gsub(1,3,SimuReal_Brasa[,i]))
}

Rfmixed_refHGDP <-read.table("Brasa_HGDP_Lat3_genome")
Rfmixed_refHGDP <- Rfmixed_refHGDP[,c(2,1,(ncol(Rfmixed_refHGDP)-59):ncol(Rfmixed_refHGDP))]
for (i in c(3:62)){
  Rfmixed_refHGDP[,i] <- as.numeric(gsub(1,3,Rfmixed_refHGDP[,i]))
}

rfmix_global <- Rfmixed_refHGDP
rm(Rfmixed_refHGDP)
simureal2_global <- SimuReal_Brasa[SimuReal_Brasa$V2 %in% rfmix_global$V1,]

a_Brasa_global <- simureal2_global == rfmix_global
b_Brasa_global <- as.data.frame(1*a_Brasa_global)
c_Brasa_global <-  colSums(b_Brasa_global)/dim(b_Brasa_global)[1]
d_Brasa_global <- simureal2_global + rfmix_global

rm(a_Brasa_global)
rm(b_Brasa_global)
rm(c_Brasa_global)

# Accu AFR
e_Brasa_global <- colSums(d_Brasa_global==4)/colSums(simureal2_global==2)
mean_AFR_Brasa_global <- mean(e_Brasa_global[c(3:62)])
sd_AFR_Brasa_global <- sd(e_Brasa_global[c(3:62)])

# Accu EUR
f_Brasa_global <- colSums(d_Brasa_global==6)/colSums(simureal2_global==3)
mean_EUR_Brasa_global <- mean(f_Brasa_global[3:62])
sd_EUR_Brasa_global <- sd(f_Brasa_global[3:62])

# Accu NAT
g_Brasa_global <- colSums(d_Brasa_global==0)/colSums(simureal2_global==0)
mean_NAT_Brasa_global <- mean(g_Brasa_global[3:62])
sd_NAT_Brasa_global <- sd(g_Brasa_global[3:62])

TabBrasa_global <- as.data.frame(c(e_Brasa_global[3:62],f_Brasa_global[3:62],g_Brasa_global[3:62]))
TabBrasa_global$Variable[1:60] <- c("AFR")
TabBrasa_global$Variable[61:120] <-c("EUR")
TabBrasa_global$Variable[121:180] <-c("NAT")
TabBrasa_global$Ref <- "12 Gen - 15% NAT / 60% EUR / 25% AFR"
names(TabBrasa_global)[1] <- "Accuracy" 
rownames(TabBrasa_global)[1:60] <- paste("AFR", c(1:60), sep="")
rownames(TabBrasa_global)[61:120] <- paste("EUR", c(1:60), sep="")
rownames(TabBrasa_global)[121:180] <- paste("NAT", c(1:60), sep="")

write.table(TabBrasa_global, "Brasa_12gen_wgs_genome.txt", col.names = T, row.names = F, quote = F, sep="\t")

i=1
ggplot(TabBrasa_global, aes(x=Ref, fill=Variable, y=Accuracy)) +
  ggtitle("LAI sensitivity by ancestry portion for GSA - All chromosomes") +
  geom_boxplot() +
  scale_fill_manual(breaks = anc,
                    values = neovivid, 
                    name = "Ancestry", 
                    labels = c("African", "European", "Native American")) +
  xlab("Admixture model") + ylab("True Positive Rate") + 
  labs(fill = "Ancestry", size = 14) +
  coord_cartesian(ylim=c(0.30,1)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), 
        title = element_text(size=12), legend.text=element_text(size=12)) 


## outlier names -global
tab_narm_global <- na.omit(TabBrasa_global, Accuracy)
dat_global <- tab_narm_global %>% tibble::rownames_to_column(var="outlier") %>% group_by(Variable) %>% mutate(is_outlier=ifelse(is_outlier(Accuracy), Accuracy, as.numeric(NA)))
dat_global$outlier[which(is.na(dat_global$is_outlier))] <- as.numeric(NA)


ggplot(dat_global, aes(x=Variable, fill=Variable, y=Accuracy)) +
  ggtitle(paste("LAI sensitivity by ancestry portion for GSA - all chromosomes",sep=" ")) +
  geom_boxplot() +
  scale_fill_manual(breaks = anc,
                    values = vivid, 
                    name = "Ancestry", 
                    labels = c("African", "European", "Native American")) +
  xlab("Admixture model") + ylab("True Positive Rate") + 
  labs(fill = "Ancestry", size = 14) +
  coord_cartesian(ylim=c(0.80,1)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), 
        title = element_text(size=12), legend.text=element_text(size=12)) + 
  geom_text(aes(label=outlier),na.rm=TRUE, hjust=1.25) + geom_quasirandom()



###### direcao do erro - PER CHR ########
tabela_concat <- list()
tabela_concat[1] <- "vazio"
an <- list()
an[1] <- "vazio"
na <- list()
na[1] <- "vazio"
en <- list()
en[1] <- "vazio"
ne <- list()
ne[1] <- "vazio"
ae <- list()
ae[1] <- "vazio"
ea <- list()
ea[1] <- "vazio"

d <- list()
d[1] <- "vazio"
e <- list()
e[1] <- "vazio"

for (k in 2:22){
    tabela_concat[[k]] <- as.data.frame(cbind(simureal2[[k]]$V1, simureal2[[k]]$V2))
    for (i in 3:dim(simureal2[[k]])[2]){
        tabela_concat[[k]][,i] <- paste(simureal2[[k]][,i],ldf[[k]][,i], sep=" ")
    }
tabela_concat[[k]]$WC_erroAN <- rowSums(tabela_concat[[k]][,3:62] == "2 0")
tabela_concat[[k]]$WC_erroNA <- rowSums(tabela_concat[[k]][,3:62] == "0 2")
tabela_concat[[k]]$WC_erroEN <- rowSums(tabela_concat[[k]][,3:62] == "3 0")
tabela_concat[[k]]$WC_erroNE <- rowSums(tabela_concat[[k]][,3:62] == "0 3")
tabela_concat[[k]]$WC_erroAE <- rowSums(tabela_concat[[k]][,3:62] == "2 3")
tabela_concat[[k]]$WC_erroEA <- rowSums(tabela_concat[[k]][,3:62] == "3 2")

an[[k]] <- tabela_concat[,c("V2", "WC_erroAN")]
na[[k]] <- tabela_concat[,c("V2", "WC_erroNA")]
en[[k]] <- tabela_concat[,c("V2", "WC_erroEN")]
ne[[k]] <- tabela_concat[,c("V2", "WC_erroNE")]
ae[[k]] <- tabela_concat[,c("V2", "WC_erroAE")]
ea[[k]] <- tabela_concat[,c("V2", "WC_erroEA")]

an[[k]]$grupo <- "AN"
colnames(an[[k]])[2] <- "erro"
na[[k]]$grupo <- "NA"
colnames(na[[k]])[2] <- "erro"

en[[k]]$grupo <- "EN"
colnames(en[[k]])[2] <- "erro"
ne[[k]]$grupo <- "NE"
colnames(ne[[k]])[2] <- "erro"

ae[[k]]$grupo <- "AE"
colnames(ae[[k]])[2] <- "erro"
ea[[k]]$grupo <- "EA"
colnames(ea[[k]])[2] <- "erro"

d[[k]] <- rbind(an[[k]], na[[k]], en[[k]], ne[[k]], ae[[k]], ea[[k]])
e[[k]] <- d[[k]][order(d[[k]]$V2),]

## all chromosomes

  tabela_concat <- as.data.frame(cbind(simureal2_global$V1, simureal2_global$V2))
  for (i in 3:dim(simureal2_global)[2]){
    tabela_concat[,i] <- paste(simureal2_global[,i],rfmix_global[,i], sep=" ")
  }
  tabela_concat$WC_erroAN <- rowSums(tabela_concat[,3:62] == "2 0")
  tabela_concat$WC_erroNA <- rowSums(tabela_concat[,3:62] == "0 2")
  tabela_concat$WC_erroEN <- rowSums(tabela_concat[,3:62] == "3 0")
  tabela_concat$WC_erroNE <- rowSums(tabela_concat[,3:62] == "0 3")
  tabela_concat$WC_erroAE <- rowSums(tabela_concat[,3:62] == "2 3")
  tabela_concat$WC_erroEA <- rowSums(tabela_concat[,3:62] == "3 2")
  
  an <- tabela_concat[,c("V2", "WC_erroAN")]
  na <- tabela_concat[,c("V2", "WC_erroNA")]
  en <- tabela_concat[,c("V2", "WC_erroEN")]
  ne <- tabela_concat[,c("V2", "WC_erroNE")]
  ae <- tabela_concat[,c("V2", "WC_erroAE")]
  ea <- tabela_concat[,c("V2", "WC_erroEA")]
  
  an$grupo <- "AN"
  colnames(an)[2] <- "erro"
  na$grupo <- "NA"
  colnames(na)[2] <- "erro"
  
  en$grupo <- "EN"
  colnames(en)[2] <- "erro"
  ne$grupo <- "NE"
  colnames(ne)[2] <- "erro"
  
  ae$grupo <- "AE"
  colnames(ae)[2] <- "erro"
  ea$grupo <- "EA"
  colnames(ea)[2] <- "erro"
  
  d <- rbind(an, na, en, ne, ae, ea)
  e <- d[order(d$V2),]
  
}

# an[[k]], na[[k]], en[[k]], ne[[k]], ae[[k]], ea[[k]] são as tabelas que vão pro
#plotgardener.

erros_global <- list()
erros_global[1] <- "vazio"

for (k in 1:22){
  d_Brasa[[k]]$WC <- rowSums(d_Brasa[[k]][,2:62] == 2 | d_Brasa[[k]][,2:62] == 3| d_Brasa[[k]][,2:62] == 5)
  as.data.frame(cbind((d_Brasa[[k]]$V2/2), d_Brasa[[k]]$WC, (d_Brasa[[k]]$V1/2))) -> erros_global[[k]]
  colnames(erros_global[[k]]) <- c("bp", "WC", "chr")
}

##### HEATMAP/CARIOGRAMA
# 3 é eur. 6 é correto soma eur-eur
# 0 é nat. 0 é correto soma nat-nat
# 2 é afr. 4 é correto soma afr-afr
# 2 erro afr-nat
# 3 erro eur-nat
# 5 erro eur-afr


## REGIOES MAIOR ERRO
hotspot <- list()
hotspot[1] <- "vazio"
for (k in 2:22){
  hotspot[[k]] <- d_Brasa[[k]][d_Brasa[[k]]$WC > 6, c(1,2)]
  write.table((hotspot[[k]]/2), paste("Wrongcallvariants_chr",k,".txt", sep=""), col.names = F, row.names = F, quote = F)
}
