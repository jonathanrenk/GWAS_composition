## Script written by Jonathan Renk
## 6 June 2020

## This script is setting up GAPIT to run GWAS

## Clearing the global environment
rm(list=ls(all=TRUE))

## Setting up the working directory
getwd()
setwd("/Users/jonathanrenk/Desktop/EMS /GWAS_Christine_SNPs/")

#source("http://www.bioconductor.org/biocLite.R")
#biocLite("multtest")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("multtest")

### These packages are now installed
#install.packages("gplots")
#install.packages("LDheatmap")
#install.packages("genetics")
#install.packages("ape")
#install.packages("EMMREML")
#install.packages("scatterplot3d")

## Loading in packages for GAPIT
library("multtest")
library("gplots")
library("LDheatmap")
library("genetics")
library("ape")
library("EMMREML")
library("compiler")
library("scatterplot3d")

source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

## Loading in packages for FarmCPU
library("bigmemory")
library("biganalytics")

source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")

#Step 1: Set working directory and import data
myY <- read.table("gwas_blup_traits_christine.txt", head = TRUE)
#myG <- read.table("widiv_447g_mona_SNPs.hmp.txt", head=FALSE)
myGD <- read.big.matrix("GAPIT.Genotype.Numerical.christine.txt", type="char", sep="\t", head = TRUE)
myGM <- read.table("GAPIT.Genotype.map.christine.txt", head = TRUE)

# Converting HapMap format to numerical for FarmCPU
#myGAPIT <- GAPIT(G=myG, output.numerical=TRUE)

#Step 2: Run GAPIT 
#myGAPIT <- GAPIT( 
#  Y=myY[,c(1,18)], #Trait Starch_As_is
#  G=myG,
#  PCA.total=5,
#  method.bin="optimum",
#  model="FarmCPU"
#)

PCA <- read.csv("/Users/jonathanrenk/Desktop/EMS /GWAS_Christine_SNPs/GAPIT_christine/GAPIT.PCA.christine.csv", stringsAsFactors=FALSE) # might need to have header TRUE
PCA <- PCA[,-1]

pvals <- read.table("FarmCPU.p.threshold.optimize.Fat_As_is.txt", head = FALSE) # change this for each trait
pval <- quantile(pvals$V1, 0.01)
print(paste0("Pvalue threshold is: ", pval))

#Step 3: Run FarmCPU
myFarmCPU <- FarmCPU(
  Y=myY[,c(1,15)], #phenotype ***change this for every trait***
  GD=myGD, #Genotype matrix
  GM=myGM, #Genotypic map
  CV=PCA, #Covariate variables (First 5 PCAs from GAPIT)
  threshold.output=0.01, #P value smaller than threshold.output will be output in GWAS table
  p.threshold=pval, 
  MAF.calculate=TRUE, #Calculate minor allele frequency (MAF) or not, if set to TRUE, the SNPs with a lower MAF (<maf.threshold) will be deleted
  method.bin="optimum",
  maf.threshold=0.05, #When MAF.calculate=TRUE, the SNPs with a lower MAF (<maf.threshold) will be deleted
  maxLoop=50 #Maximum number of iterations allowed
)

#-log10(0.01/2412791) # *This is the default signifiance threshold for FarmCPU Bonferroni (8.38252) *
#-log10(0.05/2412791) # (7.68355)
#-log10(0.01/2386666) # (8.377792)
#-log10(0.05/2386666) # (7.678822)
site_sum <- read.table("widiv_447g_christine_SNPs_SiteSummary.txt", header=T, sep = "\t")
site_sum <- site_sum[,c(-5:-14,-16:-37)]
hist(site_sum$Minor.Allele.Frequency, main="Per Site", xlab="Minor Allele Frequency", col="cadetblue", breaks = 50)
threshold_maf <- sum(site_sum$Minor.Allele.Frequency < 0.05) #25928 SNPs less than 0.05 maf
min(site_sum$Minor.Allele.Frequency)
