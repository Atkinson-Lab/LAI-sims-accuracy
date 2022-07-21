######
# Ploting the true positive rates - reference panel comparison
#####

library(reshape)
library(ggplot2)
i=1
a <- as.data.frame(rbind(TabHGDP,TabrefPelEurs,TabrefPelEasEurs))
mdata <- melt(a)
level_order = c("HGDP AMR | IBS | YRI", "PEL | IBS | YRI", "PEL + EAS | IBS | YRI")

anc = c("AFR", "EUR", "AMR")
vivid = c("#5D69B1", "#FCA315", "#1BB6AF")

library(ggthemr)
ggthemr('fresh', layout = "clean")

ggplot(mdata, aes(x=Ref, fill=Variable, y=Accuracy)) +
  ggtitle(paste("LAI sensitivity by ancestry for GSA -","chromosome",i,sep=" ")) +
  geom_boxplot() +
  scale_fill_manual(breaks = anc,
                    values = vivid, 
                    name = "Ancestry", 
                    labels = c("African", "European", "Amerindigenous")) +
  xlab("Admixture model") + ylab("True Positive Rate") + 
  labs(fill = "Ancestry", size = 14) +
  coord_cartesian(ylim=c(0.30,1)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), 
        title = element_text(size=12), legend.text=element_text(size=12)) +
  stat_summary(aes(label = , 1))
  
#################
# Miscalls by error modes plots
#################
  
# Load libraries and datasets
library("plotgardener")
library("ggplot2")
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("plotgardenerData")
library("AnnotationHub")

# Formatting miscall tables
#alter chromosome number and run following commands

i = 1

# overall (whole genome) error counts

read.table("errosglobal_chr1.txt", h=F) -> erro
erro_reform <- erro[,c(3,1,2)]
colnames(erro_reform) <- c("chrom","bp", "score")
paste("chr",erro_reform$chrom,sep='') -> erro_reform$chrom
scale(erro_reform$score) -> erro_reform$wc_scaled
as.numeric(erro_reform$wc_scaled) -> erro_reform$wc_scaled
erro_reform$end <- c(erro_reform$bp[2:length(erro_reform$bp)] - 1, erro_reform$bp[length(erro_reform$bp)] + 1)
erro_reform <- erro_reform[,c(1,2,5,3)]
colnames(erro_reform) <- c("chrom", "start", "end", "score")

# error mode counts
# en[[k]], ne[[k]], an[[k]], na[[k]], ae[[k]], ea[[k]]
### EN 
read.table("/Users/mauerjh/Documents/atkinsonlab/figures/errormodes_brasa12gen/error_EN_chr1.txt", h=F) -> erroEN
erroEN_reform <- erroEN[,c(3,1,2)]
colnames(erroEN_reform) <- c("chrom","bp", "score")
paste("chr",erroEN_reform$chrom,sep='') -> erroEN_reform$chrom
scale(erroEN_reform$score) -> erroEN_reform$wc_scaled
as.numeric(erroEN_reform$wc_scaled) -> erroEN_reform$wc_scaled
erroEN_reform$end <- c(erroEN_reform$bp[2:length(erroEN_reform$bp)] - 1, erroEN_reform$bp[length(erroEN_reform$bp)] + 1)
erroEN_reform <- erroEN_reform[,c(1,2,5,3)]
colnames(erroEN_reform) <- c("chrom", "start", "end", "score")
### NE
read.table("/Users/mauerjh/Documents/atkinsonlab/figures/errormodes_brasa12gen/error_NE_chr1.txt", h=F) -> erroNE
erroNE_reform <- erroNE[,c(3,1,2)]
colnames(erroNE_reform) <- c("chrom","bp", "score")
paste("chr",erroNE_reform$chrom,sep='') -> erroNE_reform$chrom
scale(erroNE_reform$score) -> erroNE_reform$wc_scaled
as.numeric(erroNE_reform$wc_scaled) -> erroNE_reform$wc_scaled
erroNE_reform$end <- c(erroNE_reform$bp[2:length(erroNE_reform$bp)] - 1, erroNE_reform$bp[length(erroNE_reform$bp)] + 1)
erroNE_reform <- erroNE_reform[,c(1,2,5,3)]
colnames(erroNE_reform) <- c("chrom", "start", "end", "score")
### AN
read.table("/Users/mauerjh/Documents/atkinsonlab/figures/errormodes_brasa12gen/error_AN_chr1.txt", h=F) -> erroAN
erroAN_reform <- erroAN[,c(3,1,2)]
colnames(erroAN_reform) <- c("chrom","bp", "score")
paste("chr",erroAN_reform$chrom,sep='') -> erroAN_reform$chrom
scale(erroAN_reform$score) -> erroAN_reform$wc_scaled
as.numeric(erroAN_reform$wc_scaled) -> erroAN_reform$wc_scaled
erroAN_reform$end <- c(erroAN_reform$bp[2:length(erroAN_reform$bp)] - 1, erroAN_reform$bp[length(erroAN_reform$bp)] + 1)
erroAN_reform <- erroAN_reform[,c(1,2,5,3)]
colnames(erroAN_reform) <- c("chrom", "start", "end", "score")
### NA
read.table("/Users/mauerjh/Documents/atkinsonlab/figures/errormodes_brasa12gen/error_NA_chr1.txt", h=F) -> erroNA
erroNA_reform <- erroNA[,c(3,1,2)]
colnames(erroNA_reform) <- c("chrom","bp", "score")
paste("chr",erroNA_reform$chrom,sep='') -> erroNA_reform$chrom
scale(erroNA_reform$score) -> erroNA_reform$wc_scaled
as.numeric(erroNA_reform$wc_scaled) -> erroNA_reform$wc_scaled
erroNA_reform$end <- c(erroNA_reform$bp[2:length(erroNA_reform$bp)] - 1, erroNA_reform$bp[length(erroNA_reform$bp)] + 1)
erroNA_reform <- erroNA_reform[,c(1,2,5,3)]
colnames(erroNA_reform) <- c("chrom", "start", "end", "score")
### EA
read.table("/Users/mauerjh/Documents/atkinsonlab/figures/errormodes_brasa12gen/error_EA_chr1.txt", h=F) -> erroEA
erroEA_reform <- erroEA[,c(3,1,2)]
colnames(erroEA_reform) <- c("chrom","bp", "score")
paste("chr",erroEA_reform$chrom,sep='') -> erroEA_reform$chrom
scale(erroEA_reform$score) -> erroEA_reform$wc_scaled
as.numeric(erroEA_reform$wc_scaled) -> erroEA_reform$wc_scaled
erroEA_reform$end <- c(erroEA_reform$bp[2:length(erroEA_reform$bp)] - 1, erroEA_reform$bp[length(erroEA_reform$bp)] + 1)
erroEA_reform <- erroEA_reform[,c(1,2,5,3)]
colnames(erroEA_reform) <- c("chrom", "start", "end", "score")
### AE
read.table("/Users/mauerjh/Documents/atkinsonlab/figures/errormodes_brasa12gen/error_AE_chr1.txt", h=F) -> erroAE
erroAE_reform <- erroAE[,c(3,1,2)]
colnames(erroAE_reform) <- c("chrom","bp", "score")
paste("chr",erroAE_reform$chrom,sep='') -> erroAE_reform$chrom
scale(erroAE_reform$score) -> erroAE_reform$wc_scaled
as.numeric(erroAE_reform$wc_scaled) -> erroAE_reform$wc_scaled
erroAE_reform$end <- c(erroAE_reform$bp[2:length(erroAE_reform$bp)] - 1, erroAE_reform$bp[length(erroAE_reform$bp)] + 1)
erroAE_reform <- erroAE_reform[,c(1,2,5,3)]
colnames(erroAE_reform) <- c("chrom", "start", "end", "score")

# plotting parameters

maxScore <- max(c(erro_reform$score, erroNE_reform$score, erroNA_reform$score, erroEN_reform$score, 
                  erroEA_reform$score, erroAN_reform$score, erroAE_reform$score))
print(maxScore)
params <- pgParams(
    chrom = "chr1", chromstart = 13504, chromend = 248945522,
    assembly = "hg38",
    x = 0.5, just = c("left", "top"),
    width = 4, length = 4, default.units = "inches",
    range = c(0, maxScore)
)

## plot

pageCreate(
    width = 5, height = 6, default.units = "inches", xgrid = 0, ygrid = 0
)
ideogramPlot <- plotIdeogram(params= params,  
    y = 0.5,height = 0.2,
)
plotText(
    label = "Global", fonsize = 4, fontcolor = "darkgray",
    x = 0.5, y = 0.85, just = c("left", "top"),
    default.units = "inches"
)
signalplot <- plotSignal(
    data = erro_reform,
    linecolor = "darkgray",
    params = params,
    y = "0.4b", height = 0.5,
)
## Annotate y-axis
annoYaxis(
    plot = signalplot, at = c(0, 3, 6, 9),label=c("0", "5", "10", "15"),
    axisLine = TRUE, fontsize = 9
)
## Plot y-axis label
plotText(
    label = "% wrong calls", x = 0.10, y = 1.37, rot = 90,
    fontsize = 9, fontface = "plain", just = "center",
    default.units = "inches"
)
## Plot and place signal plots
signalplotNA <- plotSignal(
    data = erroNA_reform,
    linecolor = "#1BB6AF",
    params = params,
    y = "0.03b", height = 0.5,
)
signalplotNE <- plotSignal(
    data = erroNE_reform,
    linecolor = "#1BB6AF",
    params = params,
    y = "0.03b", height = 0.5,
)
signalplotEA <- plotSignal(
    data = erroEA_reform,
    linecolor = "#FCA315",
    params = params,
    y = "0.03b", height = 0.5,
)
signalplotEN <- plotSignal(
    data = erroEN_reform,
    linecolor = "#FCA315",
    params = params,
    y = "0.03b", height = 0.5,
)
signalplotAN <- plotSignal(
    data = erroAN_reform,
    linecolor = "#7481E7",
    params = params,
    y = "0.03b", height = 0.5,
)
signalplotAE <- plotSignal(
    data = erroAE_reform,
    linecolor = "#7481E7",
    params = params,
    y = "0.03b", height = 0.5,
)
plotGenomeLabel(
    params = params,
    y = "0.01b", height = 0.2, scale = "Mb",
)
#NA NE EA EN AN AE
legendPlot <- plotLegend(
    legend = c("NAT[truth] to AFR[infer]", "NAT[truth] to EUR[infer]",
               "EUR[truth] to AFR[infer]", "EUR[truth] to NAT[infer]",
               "AFR[truth] to NAT[infer]","AFR[truth] to EUR[infer]"),
    fill = c("#1BB6AF","#1BB6AF", "#FCA315", "#FCA315", "#7481E7", "#7481E7"),
    border = FALSE,
    x = 5, y = 2.5, width = 1.5, height = 2,
    just = c("left", "top"),
    default.units = "inches"
)
pageGuideHide()

##
### Plot example of miscall hotspots and genes in this region
```{r}
#creating a page and plotting ideogram:
pageCreate(
    width = 6.25, height = 2.25, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
ideogramPlot <- plotIdeogram(
    chrom = "chr20", assembly = "hg38",
    orientation = "h",
    x = 0.25, y = 0.5, width = 5.75, height = 0.3, just = "left"
)
# We can then use annoHighlight() to highlight our genomic region of interest with a box of our desired height:
# 21: 13,187,535-13,337,741
# 20: 58,415-112,164
region <- pgParams(chrom = "chr20", chromstart = 58415, chromend = 112164)
annoHighlight(
    plot = ideogramPlot, params = region,
    fill = "red",
    y = 0.25, height = 0.5, just = c("left", "top"), default.units = "inches"
)
#To make it clearer that we are zooming in on a genomic region, we can then use annoZoomLines() to add zoom lines from the genomic region we highlighted:
annoZoomLines(
    plot = ideogramPlot, params = region,
    y0 = 0.75, x1 = c(0.25, 6), y1 = 1.25, default.units = "inches"
)
#Finally, we can add our zoomed-in gene data within the zoom lines (indicate same reg as highlight):
plotGenes(params = region, assembly = "hg38",
   x = 0.25, y = 1.5, width = 5.75, height = 1, just = "left"
)
## Plot genome label
plotGenomeLabel(params = region,
    assembly = "hg38",
    x = 0.25, y = 1.68, length = 5.75,
    default.units = "inches"
)
