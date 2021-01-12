## Script written by Jonathan Renk
## 6 June 2020

## This script is setting up GAPIT to generate PCA covariate file and then perform GWAS through FarmCPU

## Clearing the global environment
rm(list=ls(all=TRUE))

## Setting up the working directory
getwd()
setwd("/Users/jonathanrenk/Desktop/EMS /GWAS_Christine_SNPs/christine_redo/")

#source("http://www.bioconductor.org/biocLite.R")
#biocLite("multtest")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("multtest")
BiocManager::install("snpStats")

install.packages("gplots")
install.packages("LDheatmap")
install.packages("genetics")
install.packages("ape")
install.packages("EMMREML")
install.packages("scatterplot3d")

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

############################################################################################################
#### Generating the numerical format genotype files used in FarmCPU #####
############################################################################################################

myY <- read.table("gwas_blup_traits_christine_v2.txt", head = TRUE)
myG <- read.table("widiv_446g_christine_SNPs.hmp.txt", head=FALSE)

# Converting HapMap format to numerical for FarmCPU
myGAPIT <- GAPIT(G=myG, output.numerical=TRUE)

############################################################################################################
#### Generating the PCA file used in FarmCPU #####
############################################################################################################

myY <- read.table("gwas_blup_traits_christine_v2.txt", head = TRUE)
myGD <- read.big.matrix("GAPIT.Genotype.Numerical.txt", type="char", sep="\t", head = TRUE)
myGM <- read.table("GAPIT.Genotype.map.txt", head = TRUE)

myGAPIT <- GAPIT( 
  Y=myY[,c(1,2)], 
  GD=myGD,
  GM=myGM,
  PCA.total=5,
  method.bin="optimum",
  model="FarmCPU"
)

############################################################################################################
#### Generating the pvalue threshold used in FarmCPU #####
############################################################################################################

myY <- read.table("gwas_blup_traits_christine_v2.txt", head = TRUE)
myGD <- read.big.matrix("GAPIT.Genotype.Numerical.txt", type="char", sep="\t", head = TRUE)
myGM <- read.table("GAPIT.Genotype.map.txt", head = TRUE)

### NOTE! this section is for computing permuted model entry thresholds. Only needs to be run once per trait! ###
pvals <- FarmCPU.P.Threshold(
  Y=myY[,c(1,2)], #only two columns allowed, the first column is taxa name and the second is phenotype
  GD=myGD,
  GM=myGM,
  trait="Ankom.Crude.Fiber", #name of the trait, only used for the output file name
  theRep=100 #number of permutation times 
)

############################################################################################################
#### Conducting GWAS in FarmCPU #####
############################################################################################################

myY <- read.table("gwas_blup_traits_christine_v2.txt", head = TRUE)
myGD <- read.big.matrix("/Users/jonathanrenk/Desktop/EMS /GWAS_Christine_SNPs/christine_redo/GWAS_materials/GAPIT.Genotype.Numerical.txt", type="char", sep="\t", head = TRUE)
myGM <- read.table("/Users/jonathanrenk/Desktop/EMS /GWAS_Christine_SNPs/christine_redo/GWAS_materials/GAPIT.Genotype.map.txt", head = TRUE)

PCA <- read.csv("/Users/jonathanrenk/Desktop/EMS /GWAS_Christine_SNPs/christine_redo/GWAS_materials/GAPIT.PCA.csv", stringsAsFactors=FALSE) # might need to have header TRUE
PCA <- PCA[,c(-1,-7)]

pvals <- read.table("/Users/jonathanrenk/Desktop/EMS /GWAS_Christine_SNPs/christine_redo/GWAS_materials/p.threshold/FarmCPU.p.threshold.optimize.Starch_As_is.txt", head = FALSE) # change this for each trait
pval <- quantile(pvals$V1, 0.01)
print(paste0("Pvalue threshold is: ", pval))

myFarmCPU <- FarmCPU(
  Y=myY[,c(1,17)], #phenotype ***change this for every trait***
  GD=myGD, #Genotype matrix
  GM=myGM, #Genotypic map
  CV=PCA, #Covariate variables (First 5 PCAs from GAPIT)
  threshold.output=1, #P value smaller than threshold.output will be output in GWAS table
  p.threshold=pval, 
  MAF.calculate=TRUE, #Calculate minor allele frequency (MAF) or not, if set to TRUE, the SNPs with a lower MAF (<maf.threshold) will be deleted
  method.bin="optimum",
  maf.threshold=0.05, #When MAF.calculate=TRUE, the SNPs with a lower MAF (<maf.threshold) will be deleted
  maxLoop=50 #Maximum number of iterations allowed
)
