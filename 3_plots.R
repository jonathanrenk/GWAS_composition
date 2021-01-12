## Script written by Jonathan Renk
## 13 Apr 2020

## This script is contains code for all of the figures and supplemental figures used in the paper

## Clearing the global environment
rm(list=ls(all=TRUE))

## Loading in packages
library(scatterplot3d)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(cowplot)
library(reshape2)
library(plyr)
library(PerformanceAnalytics)
library(here)
library(corrplot)
library(ggcorrplot)
library(car)
library(ggpubr)
library(ggfortify)
library(cluster)
library(RColorBrewer)

############################################################################################################
#### Figure 1 #####
############################################################################################################
## Loading in the data
data <- read.csv("YC_final.csv")

## Removing SampleID, Analysis.Time, and M.Distance columns
data <- data[,c(-1, -6, -24:-29)]
## Removing Amylopectin 
data <- data[,-5]
## Removal of entry 22 and 653 being outliers from residuals
data <- data[c(-22,-653),]
## Removal of Genotype, Block, and Rep
data <- data[,c(-1:-3)]

## Setting variables as factors for Env
data[,1] <- as.factor(data[,1])

str(data)
data$Env <- revalue(data$Env, c("1"="MN16", "2"="IA16", "3"="MN17", "4"="IA17", "5"="MO17" ))
data_long <- melt(data)
head(data_long)

ggplot(data = data_long, aes(x=value, fill=Env)) +
  geom_density(alpha=0.4) + 
  labs(y = "Density", x = "Value") + 
  facet_wrap( ~ variable, scales="free", 
              labeller = labeller(variable = c("Amylopectin_average" = "Amylopectin","Ankom_Crude_Fiber_average" = " Ankom Fiber",
                                               "Ash" = "Ash", "Crude.fat" = "Crude Fat", "Crude.fiber" = "Crude Fiber",
                                               "Fructose" = "Fructose", "Glucose" = "Glucose", "N_combustion_average" = "N Combustion",
                                               "N_Kjeltec_average" = "N Kjeltec", "Sucrose_average" = "Sucrose", "Total_Sugars_average" = "Total Sugars",
                                               "Moisture" = "Moisture", "Protein_As_is" = "Protein As is", "Fat_As_is" = "Fat As is", "Fiber_As_is" = "Fiber As is",
                                               "Ash_As_is" = "Ash As is", "Starch_As_is" = "Starch As is"))) +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Blues")

############################################################################################################
#### Figure 2 #####
############################################################################################################
## Loading in the data
data <- read.csv("variance_components_v2.csv", header=T, stringsAsFactors=F)
colnames(data) <- c("Variance Component", "Ankom Crude Fiber", "Ash As Is", "Ash", "Crude Fat", "Crude Fiber", 
                    "Fat As Is", "Fiber As Is", "Fructose", "Glucose", "Moisture", "N Combustion", "N Kjeltec", 
                    "Protein As Is", "Starch As Is", "Sucrose", "Total Sugars")

meltedVarComp <- melt(data, id = "Variance Component")
meltedVarComp$value <- as.numeric(meltedVarComp$value)

ggplot(meltedVarComp, aes(variable, value)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  facet_grid(~ meltedVarComp[,"Variance Component"]) +
  scale_x_discrete(limits = rev(levels(meltedVarComp$variable))) +
  labs(x = "Trait", y = "Percent phenotypic variance explained") +
  theme_bw() + 
  theme(strip.text = element_text(size = 8), strip.background = element_rect(color = "black", fill = "white"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        legend.position = 'none')

############################################################################################################
#### Figure 3 #####
############################################################################################################
## Loading in the data
data <- read.csv("christine_blups_het.csv", header=T, stringsAsFactors=F)
str(data)
summary(data$heterotic_group)
data[,18] <- as.factor(data[,18])
data[,19] <- as.factor(data[,19])

dent <- data.frame(data [ which(data$heterotic_group_2 == 'dent'),])
other <- data.frame(data [ which(data$heterotic_group_2 == 'other'),])

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

#Dent
colnames(dent) <- c("Taxa", "Ankom Crude Fiber", "Ash", "Crude Fat", "Crude Fiber", 
                    "Fructose", "Glucose", "N Combustion", "N Kjeltec", "Sucrose", "Total Sugars", "Moisture",
                    "Protein As Is", "Fat As Is", "Fiber As Is", "Ash As Is", "Starch As Is", 
                    "Heterotic.Group", "heterotic_group_2")

dent <- dent[,c(-1,-18:-19)]
str(dent)
cormat <- cor(dent, use = "p", method = "pearson")

# matrix of the p-value of the correlation
p.mat <- cor.mtest(dent)
head(p.mat[, 1:5])

ggcorrplot(cor(dent), p.mat = p.mat, type = "upper",
                        colors = c("darkred", "white", "steelblue"),
                        ggtheme = ggplot2::theme_classic, legend.title = "Correlation\nCoefficient") + 
  labs(title = "Dent (n = 338)", x = NULL, y = NULL) +
  theme(text = element_text(),
        title =element_text(size=8, face='bold'),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 6, angle = 90, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))

#Other
colnames(other) <- c("Taxa", "Ankom Crude Fiber", "Ash", "Crude Fat", "Crude Fiber", 
                     "Fructose", "Glucose", "N Combustion", "N Kjeltec", "Sucrose", "Total Sugars", "Moisture",
                     "Protein As Is", "Fat As Is", "Fiber As Is", "Ash As Is", "Starch As Is", 
                     "Heterotic.Group", "heterotic_group_2")

other <- other[,c(-1,-18:-20)]
str(other)
cormat <- cor(other, use = "p", method = "pearson")

# matrix of the p-value of the correlation
p.mat <- cor.mtest(other)
head(p.mat[, 1:5])

ggcorrplot(cor(other), p.mat = p.mat, type = "upper",
                         colors = c("darkred", "white", "steelblue"),
                         ggtheme = ggplot2::theme_classic, legend.title = "Correlation\nCoefficient") + 
  labs(title = "Popcorn, Sweet Corn, Flint (n = 17)", x = NULL, y = NULL) +
  theme(text = element_text(),
        title =element_text(size=8, face='bold'),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 6, angle = 90, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))

############################################################################################################
#### Figure 4 #####
############################################################################################################
#### A ####
## Loading in the data
data <- read.csv("/Users/jonathanrenk/Desktop/EMS /GWAS_Christine_SNPs/christine_redo/GWAS_materials/GAPIT.PCA.csv", header=T, stringsAsFactors=F)
colnames(data) <- c("taxa", "PC1", "PC2", "PC3", "PC4", "PC5", "Maize Type")
str(data)
summary(data)
unique(data$`Maize Type`)
print(levels(data$`Maize Type`))

p <- ggplot(data) +
  geom_point(aes(x = PC1, y = PC2, color = `Maize Type`), alpha = 0.5, size = 1) +
  geom_point(data = subset(data, `Maize Type` == 'Flint'),
             aes(x = PC1, y = PC2, color = `Maize Type`),  alpha = 0.5, size = 1) +
  geom_point(data = subset(data, `Maize Type` == 'Popcorn'),
             aes(x = PC1, y = PC2, color = `Maize Type`),  alpha = 0.5, size = 1) +
  geom_point(data = subset(data, `Maize Type` == 'Sweet Corn'),
             aes(x = PC1, y = PC2, color = `Maize Type`),  alpha = 0.5, size = 1) +
  theme_classic() +
  scale_color_manual(values = c("#F8766D","#DE9200","#93AA00","#00BA38","#000000","#00B9E3","#619CFF","#DB72FB","#FF61C3"))

p2 <- ggplot(data) +
  geom_point(aes(x = PC1, y = PC3, color = `Maize Type`),  alpha = 0.5, size = 1) +
  geom_point(data = subset(data, `Maize Type` == 'Flint'),
             aes(x = PC1, y = PC3, color = `Maize Type`),  alpha = 0.5, size = 1) +
  geom_point(data = subset(data, `Maize Type` == 'Popcorn'),
             aes(x = PC1, y = PC3, color = `Maize Type`),  alpha = 0.5, size = 1) +
  geom_point(data = subset(data, `Maize Type` == 'Sweet Corn'),
             aes(x = PC1, y = PC3, color = `Maize Type`),  alpha = 0.5, size = 1) +
  theme_classic() +
  scale_color_manual(values = c("#F8766D","#DE9200","#93AA00","#00BA38","#000000","#00B9E3","#619CFF","#DB72FB","#FF61C3")) +
  theme(legend.position="none")

p3 <- ggplot(data) +
  geom_point(aes(x = PC2, y = PC3, color = `Maize Type`),  alpha = 0.5, size = 1) +
  geom_point(data = subset(data, `Maize Type` == 'Flint'),
             aes(x = PC2, y = PC3, color = `Maize Type`),  alpha = 0.5, size = 1) +
  geom_point(data = subset(data, `Maize Type` == 'Popcorn'),
             aes(x = PC2, y = PC3, color = `Maize Type`),  alpha = 0.5, size = 1) +
  geom_point(data = subset(data, `Maize Type` == 'Sweet Corn'),
             aes(x = PC2, y = PC3, color = `Maize Type`),  alpha = 0.5, size = 1) +
  theme_classic() +
  scale_color_manual(values = c("#F8766D","#DE9200","#93AA00","#00BA38","#000000","#00B9E3","#619CFF","#DB72FB","#FF61C3")) +
  theme(legend.position="none")

#Legend
legend <- get_legend(p)
#Remove the legend
p <- p + theme(legend.position="none")
#Combined figure
p4 <- grid.arrange(p, p2, p3, legend, ncol=4, widths=c(1.5, 1.5, 1.5, 1.0)) 

#### B ####
## Loading in the data
SNP <- read.csv("significant_snp_locations_v2.csv", header=T, stringsAsFactors=F)
str(SNP)
#SNP <- SNP[c(-1:-2),]

# make chromosome and end1 as numeric, trait as factor
SNP$chr <- as.factor(SNP$chr)
SNP$pos <- as.numeric(SNP$pos)
SNP$Trait <- as.factor(SNP$Trait)
SNP$Trait <- factor(SNP$Trait, levels = c("Protein As Is", "Ankom Crude Fiber", "Ash As Is", "Fat As Is", 
                                          "Fiber As Is", "Fructose", "Sucrose", "Starch As Is", "Crude Fiber", 
                                          "N Combustion", "Ash", "N Kjeltec"))

#create a new table of maize chromosome length and make it a data frame
maize_chromosomes <- cbind(chromosome = c(1:10), start = c(rep(0,10)), end = c(307041717, 243907191, 235667834, 246994605, 223902240, 174033170, 182381542, 181122637, 159769782, 150982314))
maize_chromosomes <- data.frame(maize_chromosomes)
str(maize_chromosomes)

ggplot(SNP, aes(as.integer(chr), pos)) +
  geom_segment(data = maize_chromosomes, aes(x = chromosome, xend = chromosome, y = start, yend = end), lineend = "round", color = "black", size = 16, inherit.aes = FALSE) +
  geom_segment(data = maize_chromosomes, aes(x = chromosome, xend = chromosome, y = start, yend = end), lineend = "round", color = "white", size = 15, inherit.aes = FALSE) +
  scale_y_reverse(breaks = seq(3.5e8, 0, -50e6), labels = c(350, seq(300, 0, -50)), limits = c(3.5e8, 0)) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10)) +
  ylab("Genomic positions (Mb)") + xlab ("Chromosome") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12), axis.title.x = element_text(size= 12), axis.text.y = element_text(size=11), axis.title.y = element_text(size=12)) +
  theme(legend.text = element_text(size=12), legend.title = element_text(size=12)) +
  geom_point(aes(color= Trait), position = position_dodge(width = 0.4), size = 3, alpha = 1) +
  scale_color_brewer(palette = "Paired") +
  theme(legend.title = element_blank(), legend.position = c(0.55,0.15), legend.direction = "horizontal", legend.text=element_text(size=10)) +
  scale_color_manual(values = c("#6A3D9A", "#A6CEE3", "#B2DF8A", "#FB9A99", "#E31A1C", "#FDBF6F", "#B15928", "#FFCC33", "#33A02C", "#FF7F00", "#1F78B4", "#CAB2D6"))

#### C ####
## Loading in the data
data_all <- read.csv("pve_snps.csv", header=T, stringsAsFactors=F)
str(data_all)
data_all[,1] <- as.factor(data_all[,1])
data_all$Trait <- factor(data_all$Trait, levels = c("Sucrose", "Glucose", "Total Sugars", "Fructose", "Crude Fat", "Moisture", "Starch As Is", "N Kjeltec", "Crude Fiber", "Protein As Is", "Ash", "Fat As Is", "N Combustion", "Fiber As Is", "Ankom Crude Fiber", "Ash As Is"))

ggplot(NULL, aes(Trait, PVE)) + 
  geom_bar(aes(fill = "All SNPs"), data = data, alpha = 1, stat="identity", color = "black") +
  geom_bar(aes(fill = "Significant GWAS SNPs"), data = data2, alpha = 0.8, stat="identity", color = "black") +
  labs(y = "PVE (%)", fill = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 90), axis.ticks = element_blank(), axis.text.y = element_text(size = 12)) +
  scale_fill_manual(values=c("steelblue", "lightgray"), labels = c("All SNPs", "Significant GWAS SNPs")) +
  geom_errorbar(data=data_all, aes(ymin=PVE-SE, ymax=PVE+SE)) +
  geom_errorbar(data=data_all, aes(ymin=PVE_GWAS-SE_GWAS, ymax=PVE_GWAS+SE_GWAS)) 

############################################################################################################
#### Supplemental Figure 1 #####
############################################################################################################
## Loading in the data
data <- read.table("christine_common_snp_dist.txt", header=F, sep = "\t")
str(data)
## Removing taxa name column
data <- data[,c(-1)]
colnames(data) <- NULL # removing column headers
data <- as.matrix(data) # converting to matrix 
data[upper.tri(data)] <- NA # make the upper half of the matrix as NAs
diag(data)=NA # make the diagonal NAs since they are zero
data <- cbind(which(!is.na(data), arr.ind = TRUE), na.omit(as.vector(data))) # getting into a vector of values for plotting
data <- as.data.frame(data)

## Loading in the data
data2 <- read.table("mona_common_snp_dist.txt", header=F, sep = "\t")
str(data2)
## Removing taxa name column
data2 <- data2[,c(-1)]
colnames(data2) <- NULL # removing column headers
data2 <- as.matrix(data2) # converting to matrix 
data2[upper.tri(data2)] <- NA # make the upper half of the matrix as NAs
diag(data2)=NA # make the diagonal NAs since they are zero
data2 <- cbind(which(!is.na(data2), arr.ind = TRUE), na.omit(as.vector(data2))) # getting into a vector of values for plotting
data2 <- as.data.frame(data2)

### Making a new combined data frame ###
combined_data <- data[,c(1:3)]
combined_data$Mona_GD <- data2[,c(3)]
names(combined_data)[3] <- "Christine_GD"
str(combined_data)

plot(combined_data$Christine_GD, combined_data$Mona_GD, xlab = "Genetic Distance (O'Conner et al., 2020)", ylab = "Genetic Distance (Mazaheri et al., 2019)", col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3)

############################################################################################################
#### Supplemental Figure 3 #####
############################################################################################################
lm_eqn <- function(data){
  m <- lm(Prediction ~ Observed, data);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(sqrt(summary(m)$r.squared), digits = 3)))
  as.character(as.expression(eq));
}

data <- read.csv("ankom_crude_fiber.csv", header=T, stringsAsFactors=F) #change for every trait
str(data)

p <- ggplot(data = data, aes(x = Observed, y = Prediction)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() +
  ggtitle("Ankom Crude Fiber") +
  labs(x = "Wet Lab Observed Value") +
  theme_classic() 
p
p1 <- p + annotate("text", x = 1.6, y = 2.40, label = lm_eqn(data), parse = TRUE, size = 3)

############################################################################################################
#### Supplemental Figure 4 #####
############################################################################################################
## Loading in the data
#data <- read.csv("christine_YC_final_common_447.csv", header=T, stringsAsFactors=F)
data <- read.csv("YC_final.csv")

## Removing SampleID, Analysis.Time, and M.Distance columns
data <- data[,c(-1, -6, -24:-29)]
## Removing Amylopectin 
data <- data[,-5]
## Removal of entry 22 and 653 being outliers from residuals
data <- data[c(-22,-653),]

## Setting variables as factors for genotype, block, rep, and env
data[,1] <- as.factor(data[,1])
data[,2] <- as.factor(data[,2])
data[,3] <- as.factor(data[,3])
data[,4] <- as.factor(data[,4])
str(data)

## Change for each trait ##
# Run a random effects model
model.1 <- lmer(Ankom_Crude_Fiber_average ~ (1|Genotype) + (1|Env) + (1|Env/Rep) + (1|Rep/Block) + (1|Genotype:Env), data = data, REML = TRUE)
# Decreasing stopping tolerances
strict_tol <- lmerControl(optCtrl=list(xtol_abs=1e-8, ftol_abs=1e-8))
if (all(model.1@optinfo$optimizer=="nloptwrap")) {
  model <- update(model.1, control=strict_tol)
}
summary(model, correlation=FALSE)
# qq plot
qqPlot(residuals(model), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Ankom Crude Fiber", id = FALSE) 

############################################################################################################
#### Supplemental Figure 6 #####
############################################################################################################
## Loading in the data
data <- read.csv("loc_year_variance.csv", header=T, stringsAsFactors=F)
colnames(data) <- c("Variance Component", "Ankom Crude Fiber", "Ash As Is", "Ash", "Crude Fat", "Crude Fiber", 
                    "Fat As Is", "Fiber As Is", "Fructose", "Glucose", "Moisture", "N Combustion", "N Kjeltec", 
                    "Protein As Is", "Starch As Is", "Sucrose", "Total Sugars")
meltedVarComp <- melt(data, id = "Variance Component")
meltedVarComp$value <- as.numeric(meltedVarComp$value)

ggplot(meltedVarComp, aes(variable, value)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  facet_grid(~ meltedVarComp[,"Variance Component"]) +
  scale_x_discrete(limits = rev(levels(meltedVarComp$variable))) +
  labs(x = "Trait", y = "Percent phenotypic variance explained") +
  theme_bw() + 
  theme(strip.text = element_text(size = 4), strip.background = element_rect(color = "black", fill = "white"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 4),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        legend.position = 'none')

############################################################################################################
#### Supplemental Figure 7 #####
############################################################################################################
## Loading in packages
library("qqman")

## Loading in the data
data <- read.csv("FarmCPU.Total.Sugars.GWAS.Results.csv")
#removing maf and effect
data <- data[,c(-5:-6)]
#renaming columns
names(data)[1] <- "SNP"
names(data)[2] <- "CHR"
names(data)[3] <- "BP"
names(data)[4] <- "P"
str(data)
unique(data$CHR)
NAs <- which(is.na(data$P)) #may need to run if the Error in plot.window(...) : need finite 'ylim' values appears when plotting
data <- data[c(-NAs),]
new.names <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

manhattan(data, chr="CHR", bp="BP", snp="SNP", p="P", main="Total Sugars", chrlabs = new.names, suggestiveline = F, genomewideline = F, ylim = c(0,7))
