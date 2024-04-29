library('pheatmap')
library(edgeR)
library(limma)
library(ggplot2)
library(tidyverse) 
library(RColorBrewer)
library(ggrepel)

setwd('/Users/ian/Desktop/GItR_Capstone/CMDA_Capstone')
countData <-read.table(header = TRUE, sep = "", file = "counts.txt")

colnames(countData) <- c("Geneid","Chr", "Start" ,"End", "Strand","Length","0hr1","0hr2","0hr3","Avr1","Avr2","Avr3", "EV1", "EV2","EV3")
sampleInfo <- read.table(header = TRUE, sep = "", file = "counts.txt.summary")
#countData$condition <- substr(rownames(sampleInfo), 1,nchar(rownames(sampleInfo))-1)
groupNames <- c("0hr", "0hr","0hr","AvrRxo1", "AvrRxo1","AvrRxo1","EV","EV","EV")
geneSymbols <- read.csv("geneSymbols (1).csv")
geneTranslator <- read.table(header = TRUE, sep ="",file="Nbe.v1.1_vs_Nbe.v1.txt")

dge <- DGEList(counts = countData[7:15],group = factor(groupNames))

keep<-rowSums(cpm(dge$counts)>10) >= 2
dge$counts <- dge$counts[keep, ]
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge) 
##############################################################################

dgeDispersion <- estimateCommonDisp(dge, verbose=T)

dgeDispersion <- estimateTagwiseDisp(dgeDispersion)
names(dgeDispersion$span)

et0vsAvr <- exactTest(dgeDispersion, pair=c("0hr","AvrRxo1")) # compare groups 1 and 2
et0hrvsEV <- exactTest(dgeDispersion, pair=c("0hr","EV")) # compare groups 1 and 3
etEVvsAvr <- exactTest(dgeDispersion, pair=c("EV","AvrRxo1")) # compare groups 2 and 3
topTags(etEVvsAvr, n=10)
topTags(et0hrvsEV, n=10)
topTags(et0vsAvr, n=10)

sigCount <- 10
log2Fchange <- 1
pval <- .01




######################################################################
indexs = as.double(rownames(dge$counts))

zscoreMat <- data.frame(Zerohr1=double(),Zerohr2=double(),Zerohr3=double(),Avr1=double(),Avr2=double(),Avr3=double(), EV1=double(),
                     EV2=double(),EV3=double(), index=double())
for(i in seq(1,3098)){
  zscoreMat[nrow(zscoreMat) + 1,] <- c((dge$counts[i,]-mean(dge$counts[i,]))/sd(dge$counts[i,]), indexs[i])
}





just_scores <- zscoreMat %>% select(-index)


p1 <- pheatmap(just_scores)

