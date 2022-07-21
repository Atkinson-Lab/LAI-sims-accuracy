#lendo os simureal
listBrasasimureal  <- c() # creates a list
simureal <- list() # creates a list
simureal[1] <- "vazio"
for (i in 1:22){
  listBrasasimureal[i] <- paste("/media/i123/genetica_1/Marq/allchr/SimuReal_Brasa_chr", i, ".AncCalls", sep='')
}

for (k in 2:length(listBrasasimureal)){
  simureal[[k]] <- read.table(listBrasasimureal[k], h=F)
}

for (k in 2:length(simureal)){
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
for (k in 2:length(listBrasarfmix)){
  ldf[[k]] <- read.table(listBrasarfmix[k], h=F)
}

for (k in 2:length(ldf)){
  ldf[[k]] <- ldf[[k]][,c(2,1,(ncol(ldf[[k]])-59):ncol(ldf[[k]]))]
  for (i in c(3:62)){
    ldf[[k]][,i] <- as.numeric(gsub(1,3,ldf[[k]][,i]))
  }}

simureal2 <- list()
simureal2[1] <- "vazio"
for (k in 2:length(ldf)){
  simureal2[[k]] <- simureal[[k]][simureal[[k]]$V2 %in% ldf[[k]]$V1,]
}

#preservar espaco
rm(simureal)

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


for (k in 2:22){
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
  
#GRAFICOS
TabBrasa <- list()
TabBrasa[1] <- "vazio"

for (k in 2:22){
  TabBrasa[[k]] <- as.data.frame(c(c_Brasa[[k]][3:62],e_Brasa[[k]][3:62],f_Brasa[[k]][3:62],g_Brasa[[k]][3:62]))
  TabBrasa[[k]]$Variable[1:60] <- c("TotalAccu")
  TabBrasa[[k]]$Variable[61:120] <- c("AFR")
  TabBrasa[[k]]$Variable[121:180] <-c("EUR")
  TabBrasa[[k]]$Variable[181:240] <-c("NAT")
  TabBrasa[[k]]$Ref <- "HGDP_NAT+IBS+YRI"
  names(TabBrasa[[k]])[1] <- "Accuracy" 
}

pdf("acc_allchr.pdf")
for (k in 2:22){
# mudar o titulo
  ggplot(TabBrasa[[k]], aes(x=Ref, fill=Variable, y=Accuracy)) +
    ggtitle(paste("True positive rate for WGS, 15% NAT / 60% EUR / 25% AFR proportions","chr",k,sep=" ")) +
    geom_boxplot() +
    xlab("Reference Panel") + ylab("True Positive Rate") + labs(fill = "Ancestry", size = 12) + theme_bw() + coord_cartesian(ylim=c(0.30,1)) + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), title = element_text(size=12))
}
dev.off()

#library(shades)
#library(RColorBrewer)
#LaCroixColoR Berry and Orange scale_fill_manual(values=gradient("viridis",3)
#RColorBrewer::display.brewer.pal(3,"Accent") 
# lightness(scale_fill_brewer(type="qual"), delta(-20))
#rcartocolor::display_carto_pal(4,"Vivid")

#Inspect if colorblind friendly
library(colorblindr)
cvd_grid(fig)

# if not, use okabe ito's orange #E69F00 bluish green #009E73 and reddish purple
# #CC79A7

#maybe like
ggplot(all_fst_chr) + geom_smooth(aes(x=POS,y=FST_eurafr, col="royalblue")) 
  + geom_smooth(aes(x=POS,y=FST_eurnat,col="forestgreen")) + 
  geom_smooth(aes(x=POS,y=FST_natafr, col = "coral1")) + 
  xlab(paste("Bp Chromosome ",i,sep="") ) + ylab("True Positive Rate") + 
  ggtitle(paste("LAI Sensitivity", i, sep="")) + 
  scale_colour_discrete(name = "Population", labels = c("NAT", "EUR", "AFR")) + 
  theme(plot.title=element_text(size=12))

###### direcao do erro ########
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

}

pdf("stacked_brasa_12g_allchr.pdf")
for (k in 2:22){
    ggplot(e[[k]], aes(x=V2, y=erro, fill=grupo)) + geom_area() + ggtitle(paste("Direction of error - 15/60/25 12 gen simulation - chr", k, sep=" "))
}
dev.off()

pdf("histooverlap_brasa_12g_allchr.pdf")
for (k in 2:22){
    ggplot(e[[k]], aes(x=V2, y=erro)) + geom_col(alpha=0.5, position="dodge", aes(color=grupo)) + ggtitle(paste("Direction of error - 15/60/25 12 gen simulation - chr", k, sep=" "))
}
dev.off()


##### HEATMAP/CARIOGRAMA
# 3 é eur. 6 é correto soma eur-eur
# 0 é nat. 0 é correto soma nat-nat
# 2 é afr. 4 é correto soma afr-afr
# 2 erro afr-nat
# 3 erro eur-nat
# 5 erro eur-afr

tabela_paraheat <- list()
tabela_paraheat[1] <- "vazio"

for (k in 2:22){
  d_Brasa[[k]]$WC <- rowSums(d_Brasa[[k]][,2:62] == 2 | d_Brasa[[k]][,2:62] == 3| d_Brasa[[k]][,2:62] == 5)
  as.data.frame(cbind((d_Brasa[[k]]$V2/2), d_Brasa[[k]]$WC, (d_Brasa[[k]]$V1/2))) -> tabela_paraheat[[k]]
  colnames(tabela_paraheat[[k]]) <- c("bp", "WC", "chr")
  write.csv(tabela_paraheat[[k]], paste("Brasa_errosglobal",k,".csv", sep=""), col.names = T, row.names = F, quote = F)
}

## REGIOES MAIOR ERRO
hotspot <- list()
hotspot[1] <- "vazio"
for (k in 2:22){
  hotspot[[k]] <- d_Brasa[[k]][d_Brasa[[k]]$WC > 6, c(1,2)]
  write.table((hotspot[[k]]/2), paste("Wrongcallvariants_chr",k,".txt", sep=""), col.names = F, row.names = F, quote = F)
}