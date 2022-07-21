SimuReal_Brasa2 <- read.table("SimuReal_Brasa_chr1.AncCalls", header = F)
SimuReal_Brasa <- SimuReal_Brasa2 
for (i in c(3:62)){
  SimuReal_Brasa[,i] <- as.numeric(gsub(1,3,SimuReal_Brasa2[,i]))
}

SimuReal_Baiano <- read.table("SimuReal_Baiano_chr1.AncCalls", header = F)

Rfmixed_refHGDP <-read.table("Brasa_HGDP_Lat3_1")
Rfmixed_refHGDP <- Rfmixed_refHGDP[,c(2,1,(ncol(Rfmixed_refHGDP)-59):ncol(Rfmixed_refHGDP))]
for (i in c(3:62)){
  Rfmixed_refHGDP[,i] <- as.numeric(gsub(1,3,Rfmixed_refHGDP[,i]))
}

Rfmixed_refPelEurs <-read.table("Brasa_PEL_Lat3_1")
Rfmixed_refPelEurs <- Rfmixed_refPelEurs[,c(2,1,(ncol(Rfmixed_refPelEurs)-59):ncol(Rfmixed_refPelEurs))]
for (i in c(3:62)){
  Rfmixed_refPelEurs[,i] <- as.numeric(gsub(1,3,Rfmixed_refPelEurs[,i]))
}

Rfmixed_refPelEasEurs <-read.table("Brasa_PELEAS_Lat3_1")
Rfmixed_refPelEasEurs <- Rfmixed_refPelEasEurs[,c(2,1,(ncol(Rfmixed_refPelEasEurs)-59):ncol(Rfmixed_refPelEasEurs))]
for (i in c(3:62)){
  Rfmixed_refPelEasEurs[,i] <- as.numeric(gsub(1,3,Rfmixed_refPelEasEurs[,i]))
}

Rfmixed_Baiano <-read.table("Baiano_Lat3_1")
Rfmixed_Baiano <- Rfmixed_Baiano[,c(2,1,(ncol(Rfmixed_Baiano)-59):ncol(Rfmixed_Baiano))]

############
# MIX
#HGDP only
SimuReal2_HGDP <- SimuReal_Brasa[SimuReal_Brasa$V2 %in% Rfmixed_refHGDP$V1,]

a_refHGDP <- SimuReal2_HGDP == Rfmixed_refHGDP
b_refHGDP <- as.data.frame(1*a_refHGDP)
c_refHGDP <-  colSums(b_refHGDP)/dim(b_refHGDP)[1]
d_refHGDP <- SimuReal2_HGDP + Rfmixed_refHGDP
# Accu afr
e_refHGDP <- colSums(d_refHGDP==4)/colSums(SimuReal2_HGDP==2)
mean_AFR_refHGDP <- mean(e_refHGDP[c(3:62)])
sd_AFR_refHGDP <- sd(e_refHGDP[c(3:62)])
# Accu eur
f_refHGDP <- colSums(d_refHGDP==6)/colSums(SimuReal2_HGDP==3)
mean_EUR_refHGDP <- mean(f_refHGDP[3:62])
sd_EUR_refHGDP <- sd(f_refHGDP[3:62])
# Accu nat
g_refHGDP <- colSums(d_refHGDP==0)/colSums(SimuReal2_HGDP==0)
mean_NAT_refHGDP <- mean(g_refHGDP[3:62])
sd_NAT_refHGDP <- sd(g_refHGDP[3:62])

TabHGDP <- as.data.frame(c(e_refHGDP[3:62],f_refHGDP[3:62],g_refHGDP[3:62]))
TabHGDP$Variable[1:60] <- c("AFR")
TabHGDP$Variable[61:120] <-c("EUR")
TabHGDP$Variable[121:180] <-c("NAT")
TabHGDP$Ref <- "HGDP NAT | IBS | YRI"
names(TabHGDP)[1] <- "Accuracy" 
rownames(TabHGDP)[1:60] <- paste("AFR", c(1:60), sep="")
rownames(TabHGDP)[61:120] <- paste("EUR", c(1:60), sep="")
rownames(TabHGDP)[121:180] <- paste("NAT", c(1:60), sep="")

#################
#Pel and Eurs 1KG
SimuReal2_peleurs <- SimuReal_Brasa[SimuReal_Brasa$V2 %in% Rfmixed_refPelEurs$V1,]

a_refPelEurs <- SimuReal2_peleurs == Rfmixed_refPelEurs
b_refPelEurs <- as.data.frame(1*a_refPelEurs)
c_refPelEurs <-  colSums(b_refPelEurs)/dim(b_refPelEurs)[1]
d_refPelEurs <- SimuReal2_peleurs + Rfmixed_refPelEurs

# Accu afr
e_refPelEurs <- colSums(d_refPelEurs==4)/colSums(SimuReal2_peleurs==2)
mean_AFR_refPelEurs <- mean(e_refPelEurs[c(3:62)])
sd_AFR_refPelEurs <- sd(e_refPelEurs[c(3:62)])
# Accu EUR
f_refPelEurs <- colSums(d_refPelEurs==6)/colSums(SimuReal2_peleurs==3)
mean_EUR_refPelEurs <- mean(f_refPelEurs[3:62])
sd_EUR_refPelEurs <- sd(f_refPelEurs[3:62])
# Accu Nat
g_refPelEurs <- colSums(d_refPelEurs==0)/colSums(SimuReal2_peleurs==0)
mean_NAT_refPelEurs <- mean(g_refPelEurs[3:62])
sd_NAT_refPelEurs <- sd(g_refPelEurs[3:62])

TabrefPelEurs <- as.data.frame(c(e_refPelEurs[3:62],f_refPelEurs[3:62],g_refPelEurs[3:62]))
TabrefPelEurs$Variable[1:60] <- c("AFR")
TabrefPelEurs$Variable[61:120] <-c("EUR")
TabrefPelEurs$Variable[121:180] <-c("NAT")
TabrefPelEurs$Ref <- "PEL | IBS | YRI"
names(TabrefPelEurs)[1] <- "Accuracy" 
rownames(TabrefPelEurs)[1:60] <- paste("AFR", c(1:60), sep="")
rownames(TabrefPelEurs)[61:120] <- paste("EUR", c(1:60), sep="")
rownames(TabrefPelEurs)[121:180] <- paste("NAT", c(1:60), sep="")

#################
#Pel + Eas and Eurs 1KG
SimuReal2_PelEasEurs <- SimuReal_Brasa[SimuReal_Brasa$V2 %in% Rfmixed_refPelEasEurs$V1,]
a_refPelEasEurs <- SimuReal2_PelEasEurs == Rfmixed_refPelEasEurs
b_refPelEasEurs <- as.data.frame(1*a_refPelEasEurs)
c_refPelEasEurs <-  colSums(b_refPelEasEurs)/dim(b_refPelEasEurs)[1]
d_refPelEasEurs <- SimuReal2_PelEasEurs + Rfmixed_refPelEasEurs
# Accu afr
e_refPelEasEurs <- colSums(d_refPelEasEurs==4)/colSums(SimuReal2_PelEasEurs==2)
mean_AFR_refPelEasEurs <- mean(e_refPelEasEurs[c(3:62)])
sd_AFR_refPelEasEurs <- sd(e_refPelEasEurs[c(3:62)])
# Accu EUR
f_refPelEasEurs <- colSums(d_refPelEasEurs==6)/colSums(SimuReal2_PelEasEurs==3)
mean_EUR_refPelEasEurs <- mean(f_refPelEasEurs[3:62])
sd_EUR_refPelEasEurs <- sd(f_refPelEasEurs[3:62])
# Accu NAT
g_refPelEasEurs <- colSums(d_refPelEasEurs==0)/colSums(SimuReal2_PelEasEurs==0)
mean_NAT_refPelEasEurs <- mean(g_refPelEasEurs[3:62])
sd_NAT_refPelEasEurs <- sd(g_refPelEasEurs[3:62])

TabrefPelEasEurs <- as.data.frame(c(e_refPelEasEurs[3:62],f_refPelEasEurs[3:62],g_refPelEasEurs[3:62]))
TabrefPelEasEurs$Variable[1:60] <- c("AFR")
TabrefPelEasEurs$Variable[61:120] <-c("EUR")
TabrefPelEasEurs$Variable[121:180] <-c("NAT")
TabrefPelEasEurs$Ref <- "PEL + EAS | IBS | YRI"
names(TabrefPelEasEurs)[1] <- "Accuracy" 
rownames(TabrefPelEasEurs)[1:60] <- paste("AFR", c(1:60), sep="")
rownames(TabrefPelEasEurs)[61:120] <- paste("EUR", c(1:60), sep="")
rownames(TabrefPelEasEurs)[121:180] <- paste("NAT", c(1:60), sep="")

######
# Ploting the simulations
#####




library(reshape)
# pode nao ser necessário se for só uma tab
a <- as.data.frame(rbind(TabHGDP,TabrefPelEurs,TabrefPelEasEurs))
mdata <- melt(a)

write.table(mdata, "refpaneltab_brasa9.txt", col.names = T, row.names = F, quote = F)

level_order = c("HGDP NAT | IBS | YRI", "PEL | IBS | YRI", "PEL + EAS | IBS | YRI")

# mudar o titulo
p2 <- ggplot(mdata, aes(x=factor(Ref, level = level_order), fill=Variable, y=value)) +
  ggtitle("Accuracy for WGS, 15% NAT/60% EUR/25% AFR proportions") +
  geom_boxplot() +
  xlab("Reference Panel") + ylab("Accuracy") + labs(fill = "AccuracyType") + theme_bw() + coord_cartesian(ylim=c(0,1))
p2


anc = c("AFR", "EUR", "NAT")
vivid = c("#5D69B1", "#FCA315", "#1BB6AF")
juice = c("#B25D91", "#FCA315", "#1BB6AF")
okabe = c("#CC79A7", "#E69F00", "#009E73")

library(ggthemr)
ggthemr('fresh', layout = "clean")

library(ggplot2)
i=1
ggplot(mdata, aes(x=Ref, fill=Variable, y=Accuracy)) +
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






#################
## BAIANO

#Afr (YRI) and Eurs (IBR) 1KG 50/50
SimuReal2_Baiano <- SimuReal_Baiano[SimuReal_Baiano$V2 %in% Rfmixed_Baiano$V1,]
a_refBaiano <- SimuReal2_Baiano == Rfmixed_Baiano
b_refBaiano <- as.data.frame(1*a_refBaiano)
c_refBaiano <-  colSums(b_refBaiano)/dim(b_refBaiano)[1]
d_refBaiano <- SimuReal2_Baiano + Rfmixed_Baiano

# Accu afr
e_refBaiano <- colSums(d_refBaiano==2)/colSums(SimuReal2_Baiano==1)
mean_AFR_refBaiano <- mean(e_refBaiano[3:62])
sd_AFR_refBaiano <- sd(e_refBaiano[3:62])

# Accu Eur
g_refBaiano <- colSums(d_refBaiano==0)/colSums(SimuReal2_Baiano==0)
mean_EUR_refBaiano <- mean(g_refBaiano[3:62])
sd_EUR_refBaiano <- sd(g_refBaiano[3:62])