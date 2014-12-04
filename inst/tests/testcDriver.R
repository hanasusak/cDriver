
rm(list=ls()) 
library(devtools)
library(testthat)
setwd('/Users/hanasusak/Desktop/projects/cDriver')
document()
check()
install()
#install.packages('data.table')
#build()
library(data.table)
#library(cDriver)
#data()
?CCF
# calculate CCF 
head(sample.genes.mutect)
sample.genes.mutect <- CCF(sample.genes.mutect)

#df <- read.table('../PR_TCGA_HNSC_PAIR_Capture_All_Pairs_QCPASS_v4.aggregated.capture.tcga.uuid.automated.somatic_509_hana_cleaned_luis.maf.cadd.mod.added.correct', sep='\t', quote="", header=T)
#df <- CCF(df)

# without header
b <- CCF(sample.genes.mutect, VAF = 7, ploidy = 18, CCF_CNV = 19)
head(b)
rm(b)

# load genes  of interest as vector
#bcgr.prob <- bcgr.combine(sample.genes.mutect, length.genes$Hugo_Symbol, length.genes$Coverd_len)
bcgr.prob <- bcgr.combine(sample.genes.mutect)
#bcgr.prob <- bcgr.combine(df)
head(bcgr.prob)

#df1 <- bayes.risk(sample.genes.mutect, bcgr.prob, genes=length.genes$Hugo_Symbol, prior.sick = 0.00007) 
df1 <- bayes.risk(sample.genes.mutect, bcgr.prob, prior.sick = 0.00007) 
#df1 <- bayes.risk(df, bcgr.prob, prior.sick = 0.00007) 
head(df1)

df2 <- bayes.driver(sample.genes.mutect, bcgr.prob, prior.driver = 0.001) 
#df2 <- bayes.driver(df, bcgr.prob, prior.driver = 0.001) 
head(df2)

df.final <- combine.ranking(list(df1, df2),  min.mut = 2 )
head(df.final)

plotSamplesMut(sample.genes.mutect)
plotMutChange(sample.genes.mutect)
boxplotCCF.mutations(sample.genes.mutect, df.final, color='Variant_Type')
boxplotCCF.patients(sample.genes.mutect, df.final)
boxplotCCF.patients(sample.genes.mutect, df.final, color='gender')

plotStaircase(sample.genes.mutect, df.final)


####
df22 <- read.table('../test.final.BLCA.maf', header=T)
df22 <- df22[df22$Hugo_Symbol != 'Unknown',]

source("/Users/hanasusak/Desktop/projects/genLov/genLov/r/Common_script.R")
setwd('/Users/hanasusak/Desktop/projects/cDriver')

# problem 
df.mapping[2200,] 
df.mapping[2220,] # problem here 

cor.nm <- correct.names(df22$Hugo_Symbol)

df22$Hugo_Symbol <- cor.nm$new.names

#prev double
df22[df22$Hugo_Symbol == 'AK3L1' & df22$Chromosome == 1, 'Hugo_Symbol'] <- 'AK4'
df22[df22$Hugo_Symbol == 'ODZ3' & df22$Chromosome == 4, 'Hugo_Symbol'] <- 'TENM3'
df22[df22$Hugo_Symbol == 'PRSSL1' & df22$Chromosome == 19, 'Hugo_Symbol'] <- 'PRSS57'
#syno double
df22[df22$Hugo_Symbol == 'MLL4' & df22$Chromosome == 19, 'Hugo_Symbol'] <- 'KMT2B'
# prev and syno
df22[df22$Hugo_Symbol == 'DBC1' & df22$Chromosome == 9, 'Hugo_Symbol'] <- 'BRINP1'
df22[df22$Hugo_Symbol == 'MLL2' & df22$Chromosome == 12, 'Hugo_Symbol'] <- 'KMT2D'
df22[df22$Hugo_Symbol == 'CCRL1' & df22$Chromosome == 3, 'Hugo_Symbol'] <- 'ACKR4'
df22[df22$Hugo_Symbol == 'CCBP2' & df22$Chromosome == 3, 'Hugo_Symbol'] <- 'ACKR2'
df22[df22$Hugo_Symbol == 'RAB39' & df22$Chromosome == 11, 'Hugo_Symbol'] <- 'RAB39A'
df22[df22$Hugo_Symbol == 'C11orf2' & df22$Chromosome == 11, 'Hugo_Symbol'] <- 'VPS51'
    
g2 <- unique(df22$Hugo_Symbol)
sum(!g2 %in% length.genes$Hugo_Symbol)


patDat <- read.table('../padData.txt', header=T, sep='\t')
head(patDat)
patDat <- patDat[,1:12]
df23 <- merge(x=df22,y=patDat, by.x='Tumor_Sample_Barcode', by.y='Sample.ID', all.x=TRUE, all.y=FALSE)

df23 <- CCF(df23)
bcgr.prob <- bcgr(df23)
head(bcgr.prob)

bcgr.lw <- bcgr.lawrence(df23)
head(bcgr.lw)

plot(bcgr.prob, bcgr.lw)
text(bcgr.prob, bcgr.lw, names(bcgr.prob), cex=0.75)
cor(bcgr.prob, bcgr.lw)

bcgr.2 <- bcgr.combine(df23)
head(bcgr.2)


df1.2 <- bayes.risk(df23,  bcgr.2, prior.sick =  250/(10^6)) 
head(df1.2)

df2.2 <- bayes.driver(df23, bcgr.2, prior.driver = 0.001) 
head(df2.2)

df.final <- combine.ranking(list(df1.2, df2.2))
head(df.final)

df23$Tumor.grade <- as.factor(df23$Tumor.grade)
df23$Tumor.stage <- as.factor(df23$Tumor.stage)
df23$Date.of.initial.pathologic.diagnosis <- as.factor(df23$Date.of.initial.pathologic.diagnosis)

boxplotCCF(df23, df.final)
boxplotCCF(df23, df.final, topGenes=30, color='Tumor.stage', shape=('Gender'))
plotStaircase(df23, df.final, allSamples=T)
plotStaircase(df23, df.final, allSamples=F)

boxplotCCF(df23, df.final, topGenes=30, color='Tumor.stage', shape='Gender')
boxplotCCF(df23, df.final, topGenes=30, color='Tumor.stage')
