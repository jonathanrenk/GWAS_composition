## Script written by Jonathan Renk
## 13 Apr 2020

## This script is setting up GAPIT to run GWAS

## Clearing the global environment
rm(list=ls(all=TRUE))

## Setting up the working directory
getwd()
setwd("/Users/jonathanrenk/Desktop/EMS /GWAS_Christine_SNPs/GAPIT_christine/")

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

## Loading in the data
data <- read.csv("GAPIT.PCA.christine.csv", header=T, stringsAsFactors=F)

str(data)
summary(data)
unique(data$heterotic_group)
print(levels(data$heterotic_group))

## Setting variable of heterotic group as factor
data[,7] <- as.factor(data[,7])
head(data)
## 2D PCA plots
## Formatting the plot
#theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
## plotting
#PCA1 vs PCA2
p<-ggplot(data,aes(x=PC1,y=PC2,color=heterotic_group ))
p<-p+geom_point()+ theme_classic()
p
#PCA1 vs PCA3
p2<-ggplot(data,aes(x=PC1,y=PC3,color=heterotic_group ))
p2<-p2+geom_point()+ theme_classic() + theme(legend.position="none")
p2
#PCA2 vs PCA3
p3<-ggplot(data,aes(x=PC2,y=PC3,color=heterotic_group ))
p3<-p3+geom_point()+ theme_classic() + theme(legend.position="none")
p3
#Legend
legend <- get_legend(p)
#Remove the legend
p <- p + theme(legend.position="none")
#Combined figure
p4 <- grid.arrange(p, p2, p3, legend, ncol=4, widths=c(1.5, 1.5, 1.5, 1.0)) 
ggsave("pc_column.2.pdf", p4, width = 7.5, height = 3.5, units = c("in"))

## 3D Scatterplot
colors <- c("#F8766D","#DE9200","#93AA00","#00BA38","#00C19F","#00B9E3","#619CFF","#DB72FB","#FF61C3")
colors <- colors[as.numeric(data$heterotic_group)]
s3d <- scatterplot3d(data$PC1,data$PC2,data$PC3,xlab="PC1",ylab="PC2",zlab="PC3",pch=16,color=colors,cex=1,cex.lab=1.4, cex.axis=1.0,lwd=3,angle=55,scale.y=0.7) 
legend("bottomleft", inset=0.1, legend = levels(data$heterotic_group),
       col =  c("#F8766D","#DE9200","#93AA00","#00BA38","#00C19F","#00B9E3","#619CFF","#DB72FB","#FF61C3"), pch = 16, cex = 0.5)

### Manuscript figures - distribution of raw phenotypic values
## Clearing the global environment
rm(list=ls(all=TRUE))

## Setting up the working directory
getwd()
setwd("/Users/jonathanrenk/Desktop/EMS /GWAS_Christine_SNPs/")

## Loading in the data
data <- read.csv("christine_YC_final_common_447.csv", header=T, stringsAsFactors=F)

## Removing SampleID, Analysis.Time, and M.Distance columns
data <- data[,c(-1:-4,-6, -24:-29)]
data <- data[,c(-2)] # removing amylopectin
data[,1] <- as.factor(data[,1])
str(data)
data$Env <- revalue(data$Env, c("1"="MN16", "2"="IA16", "3"="MN17", "4"="IA17", "5"="MO17" ))
data_long <- melt(data)
head(data_long)

p <- ggplot(data = data_long, aes(x=value, fill=Env)) +
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
p

### Manuscript figures - Spearman correlations of BLUPs
## Clearing the global environment
rm(list=ls(all=TRUE))

## Setting up the working directory
getwd()
setwd("/Users/jonathanrenk/Desktop/EMS /GWAS_Christine_SNPs/")

## Loading in the data
data <- read.csv("gwas_blup_traits_christine.csv", header=T, stringsAsFactors=F)
data <- data[,c(-1)]

chart.Correlation(data, histogram = FALSE, method = c("spearman"), pch=19)
chart.Correlation(data, histogram = FALSE, method = c("pearson"), pch=19)

### Making barplot of variance components ###
## Clearing the global environment
rm(list=ls(all=TRUE))

## Setting up the working directory
getwd()
setwd("/Users/jonathanrenk/Desktop/EMS /GWAS_Christine_SNPs/phenos/")

## Loading in the data
data <- read.csv("variance_components_2.csv", header=T, stringsAsFactors=F)

meltedVarComp <- melt(data, id = "Variance.Component")
meltedVarComp$value <- as.numeric(meltedVarComp$value)

ggplot(meltedVarComp, aes(variable, value)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  facet_grid(~ meltedVarComp[,"Variance.Component"]) +
  scale_x_discrete(limits = rev(levels(meltedVarComp$variable))) +
  labs(x = "Trait", y = "Percent phenotypic variance explained") +
  theme_bw() + 
  theme(strip.text = element_text(size = 8), strip.background = element_rect(color = "black", fill = "white"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        legend.position = 'none')

### Making the chromosome view plot of significant SNPs
## Clearing the global environment
rm(list=ls(all=TRUE))

## Setting up the working directory
getwd()
setwd("/Users/jonathanrenk/Desktop/EMS /GWAS_Christine_SNPs/")

## Loading in the data
SNP <- read.csv("significant_snp_locations.csv", header=T, stringsAsFactors=F)
str(SNP)
SNP <- SNP[c(-1:-2),]

# make chromosome and end1 as numeric, trait as factor
SNP$chr <- as.factor(SNP$chr)
SNP$pos <- as.numeric(SNP$pos)
SNP$Trait <- as.factor(SNP$Trait)
SNP$Trait_1 <- as.factor(SNP$Trait_1)

#create a new table of maize chromosome length and make it a data frame
maize_chromosomes <- cbind(chromosome = c(1:10), start = c(rep(0,10)), end = c(307041717, 243907191, 235667834, 246994605, 223902240, 174033170, 182381542, 181122637, 159769782, 150982314))
maize_chromosomes <- data.frame(maize_chromosomes)
str(maize_chromosomes)

#making the figure
plot2 <- ggplot() +
  geom_segment(data = maize_chromosomes, aes(x = chromosome, xend = chromosome, y = start, yend = end), lineend = "round", color = "lightgrey", size = 5) +
  geom_segment(data = SNP, aes(x = as.integer(chr) - 0.05, xend = as.integer(chr) + 0.05, y = pos, yend = pos, color = Trait), size = 1) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10)) + scale_y_continuous(breaks = seq(0, 3.5e8, 50e6), labels = c(0, seq(50, 350, 50)), limits = c(0, 3.5e8)) + ylab("Genomic positions (Mb)") + xlab ("Chromosome") + 
  theme_classic() +
  theme(axis.text.x = element_text(size=11), axis.title.x = element_text(size= 12), axis.text.y = element_text(size=11), axis.title.y = element_text(size=12)) +
  theme(legend.text = element_text(size=12), legend.title = element_text(size=12))

plot3 <- ggplot() +
  geom_segment(data = maize_chromosomes, aes(x = chromosome, xend = chromosome, y = start, yend = end), lineend = "round", color = "black", size = 6) +
  geom_segment(data = maize_chromosomes, aes(x = chromosome, xend = chromosome, y = start, yend = end), lineend = "round", color = "white", size = 5) +
  geom_segment(data = SNP, aes(x = as.integer(chr) - 0.05, xend = as.integer(chr) + 0.05, y = pos, yend = pos, color = Trait), size = 1) +
  scale_y_reverse(breaks = seq(3.5e8, 0, -50e6), labels = c(350, seq(300, 0, -50)), limits = c(3.5e8, 0)) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10)) +
  ylab("Genomic positions (Mb)") + xlab ("Chromosome") +
  theme_classic() +
  theme(axis.text.x = element_text(size=11), axis.title.x = element_text(size= 12), axis.text.y = element_text(size=11), axis.title.y = element_text(size=12)) +
  theme(legend.text = element_text(size=12), legend.title = element_text(size=12))

plot <- ggplot(SNP, aes(as.integer(chr), pos)) +
  geom_segment(data = maize_chromosomes, aes(x = chromosome, xend = chromosome, y = start, yend = end), lineend = "round", color = "black", size = 16, inherit.aes = FALSE) +
  geom_segment(data = maize_chromosomes, aes(x = chromosome, xend = chromosome, y = start, yend = end), lineend = "round", color = "white", size = 15, inherit.aes = FALSE) +
  scale_y_reverse(breaks = seq(3.5e8, 0, -50e6), labels = c(350, seq(300, 0, -50)), limits = c(3.5e8, 0)) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10)) +
  ylab("Genomic positions (Mb)") + xlab ("Chromosome") +
  theme_classic() +
  theme(axis.text.x = element_text(size=10), axis.title.x = element_text(size= 10), axis.text.y = element_text(size=11), axis.title.y = element_text(size=12)) +
  theme(legend.text = element_text(size=10), legend.title = element_text(size=10)) +
  geom_point(aes(color= Trait), position = position_dodge(width = 0.4), size = 3) +
  theme(legend.position = "bottom")

ggsave(plot, filename = "karotype.pdf", width = 7.5, height = 7, units = c("in"))

### Making plots for Pearson correlations of BLUPs per heterotic group ###
## Clearing the global environment
rm(list=ls(all=TRUE))

## Setting up the working directory
getwd()
setwd("/Users/jonathanrenk/Desktop/EMS /GWAS_Christine_SNPs/")

## Loading in the data
data <- read.csv("gwas_blup_traits_christine.csv", header=T, stringsAsFactors=F)
str(data)
data <- data[,c(-2)]
summary(data$heterotic_group)
data[,18] <- as.factor(data[,18])

## Subsetting the data by heterotic group
stiff_stalk <- data.frame(data [ which(data$heterotic_group == 'stiff_stalk'),])
non_stiff_stalk <- data.frame(data [ which(data$heterotic_group == 'non_stiff_stalk'),])
unknown <- data.frame(data [ which(data$heterotic_group == 'unknown'),])
iodent <- data.frame(data [ which(data$heterotic_group == 'iodent'),])
mixed <- data.frame(data [ which(data$heterotic_group == 'mixed'),])
popcorn <- data.frame(data [ which(data$heterotic_group == 'popcorn'),])
tropical <- data.frame(data [ which(data$heterotic_group == 'tropical'),])
sweet_corn <- data.frame(data [ which(data$heterotic_group == 'sweet_corn'),])
flint <- data.frame(data [ which(data$heterotic_group == 'flint'),])

#Stiff stalk
stiff_stalk <- stiff_stalk[,c(-1,-18)]
str(stiff_stalk)
cormat <- cor(stiff_stalk, use = "p", method = "pearson")

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

# matrix of the p-value of the correlation
p.mat <- cor.mtest(stiff_stalk)
head(p.mat[, 1:5])

stiff_stalk_plot <- ggcorrplot(cor(stiff_stalk), p.mat = p.mat, type = "upper",
                               colors = c("darkred", "white", "steelblue"),
                               ggtheme = ggplot2::theme_classic, legend.title = "Correlation\nCoefficient") + 
  labs(title = "Stiff Stalk (n = 135)", x = NULL, y = NULL) +
  theme(text = element_text(),
        title =element_text(size=8, face='bold'),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 6, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))

#Non stiff stalk
non_stiff_stalk <- non_stiff_stalk[,c(-1,-18)]
str(non_stiff_stalk)
cormat <- cor(non_stiff_stalk, use = "p", method = "pearson")

# matrix of the p-value of the correlation
p.mat <- cor.mtest(non_stiff_stalk)
head(p.mat[, 1:5])

non_stiff_stalk_plot <- ggcorrplot(cor(non_stiff_stalk), p.mat = p.mat, type = "upper",
                                   colors = c("darkred", "white", "steelblue"),
                                   ggtheme = ggplot2::theme_classic, legend.title = "Correlation\nCoefficient") + 
  labs(title = " Non Stiff Stalk (n = 123)", x = NULL, y = NULL) +
  theme(text = element_text(),
        title =element_text(size=8, face='bold'),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 6, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))

#Unknown
unknown <- unknown[,c(-1,-18)]
str(unknown)
cormat <- cor(unknown, use = "p", method = "pearson")

# matrix of the p-value of the correlation
p.mat <- cor.mtest(unknown)
head(p.mat[, 1:5])

unknown_plot <- ggcorrplot(cor(unknown), p.mat = p.mat, type = "upper",
                           colors = c("darkred", "white", "steelblue"),
                           ggtheme = ggplot2::theme_classic, legend.title = "Correlation\nCoefficient") + 
  labs(title = "Unknown (n = 71)", x = NULL, y = NULL) +
  theme(text = element_text(),
        title =element_text(size=8, face='bold'),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 6, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))

#Iodent
iodent <- iodent[,c(-1,-18)]
str(iodent)
cormat <- cor(iodent, use = "p", method = "pearson")

# matrix of the p-value of the correlation
p.mat <- cor.mtest(iodent)
head(p.mat[, 1:5])

iodent_plot <- ggcorrplot(cor(iodent), p.mat = p.mat, type = "upper",
                          colors = c("darkred", "white", "steelblue"),
                          ggtheme = ggplot2::theme_classic, legend.title = "Correlation\nCoefficient") + 
  labs(title = "Iodent (n = 46)", x = NULL, y = NULL) +
  theme(text = element_text(),
        title =element_text(size=8, face='bold'),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 6, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))

#Mixed
mixed <- mixed[,c(-1,-18)]
str(mixed)
cormat <- cor(mixed, use = "p", method = "pearson")

# matrix of the p-value of the correlation
p.mat <- cor.mtest(mixed)
head(p.mat[, 1:5])

mixed_plot <- ggcorrplot(cor(mixed), p.mat = p.mat, type = "upper", 
                         colors = c("darkred", "white", "steelblue"),
                         ggtheme = ggplot2::theme_classic, legend.title = "Correlation\nCoefficient") + 
  labs(title = "Mixed (n = 46)", x = NULL, y = NULL) +
  theme(text = element_text(),
        title =element_text(size=8, face='bold'),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 6, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))

#Popcorn
popcorn <- popcorn[,c(-1,-18)]
str(popcorn)
cormat <- cor(popcorn, use = "p", method = "pearson")

# matrix of the p-value of the correlation
p.mat <- cor.mtest(popcorn)
head(p.mat[, 1:5])

popcorn_plot <- ggcorrplot(cor(popcorn), p.mat = p.mat, type = "upper", 
                           colors = c("darkred", "white", "steelblue"),
                           ggtheme = ggplot2::theme_classic, legend.title = "Correlation\nCoefficient") + 
  labs(title = "Popcorn (n = 11)", x = NULL, y = NULL) +
  theme(text = element_text(),
        title =element_text(size=8, face='bold'),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 6, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))

#Tropical
tropical <- tropical[,c(-1,-18)]
str(tropical)
cormat <- cor(tropical, use = "p", method = "pearson")

# matrix of the p-value of the correlation
p.mat <- cor.mtest(tropical)
head(p.mat[, 1:5])

tropical_plot <- ggcorrplot(cor(tropical), p.mat = p.mat, type = "upper", insig = "blank", 
                            colors = c("darkred", "white", "steelblue"),
                            ggtheme = ggplot2::theme_classic, legend.title = "Correlation\nCoefficient") + 
  labs(title = "Tropical (n = 9)", x = NULL, y = NULL) +
  theme(text = element_text(),
        title =element_text(size=8, face='bold'),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 6, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))

#Sweet corn
sweet_corn <- sweet_corn[,c(-1,-18)]
str(sweet_corn)
cormat <- cor(sweet_corn, use = "p", method = "pearson")

# matrix of the p-value of the correlation
p.mat <- cor.mtest(sweet_corn)
head(p.mat[, 1:5])

sweet_corn_plot <- ggcorrplot(cor(sweet_corn), p.mat = p.mat, type = "upper", insig = "blank", 
                              colors = c("darkred", "white", "steelblue"),
                              ggtheme = ggplot2::theme_classic, legend.title = "Correlation\nCoefficient") + 
  labs(title = "Sweet Corn (n = 5)", x = NULL, y = NULL) +
  theme(text = element_text(),
        title =element_text(size=8, face='bold'),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 6, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))

#Flint cannot due because of one genotype
flint <- flint[,c(-1,-19)]
str(flint)
cormat <- cor(flint, use = "p", method = "pearson")

# matrix of the p-value of the correlation
p.mat <- cor.mtest(flint)
head(p.mat[, 1:5])

ggcorrplot(cor(flint), p.mat = p.mat, type = "upper",
           colors = c("darkred", "white", "steelblue"),
           ggtheme = ggplot2::theme_classic, legend.title = "Correlation\nCoefficient") + 
  labs(title = "flint\nn = 1", x = NULL, y = NULL) +
  theme(text = element_text(),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))

combined <- ggarrange(stiff_stalk_plot, non_stiff_stalk_plot, unknown_plot, iodent_plot, mixed_plot, popcorn_plot,
                      ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")

ggsave(combined, filename = "corrs.pdf", width = 11, height = 7.5, units = c("in"))

### Making the plot for Christine's SNPs of percent heterozygosity per genotype ###
## Clearing the global environment
rm(list=ls(all=TRUE))

## Setting up the working directory
getwd()
setwd("/Users/jonathanrenk/Desktop/EMS /GWAS_Christine_SNPs/")

## Loading in the data
data <- read.table("taxa_summary.txt", header=T, sep = "\t")
str(data)

# Number of Heterozygous SNPs per genotype
hist(data$Number.Heterozygous, main="Histogram of number of heterozygous SNPs", xlab="Heterozygous SNPs", col="cadetblue")
boxplot(data$Number.Heterozygous)

# Proportion of Heterozygous SNPs per genotype
hist(data$Proportion.Heterozygous, main="Histogram of proportion of heterozygous SNPs", xlab="Heterozygous SNPs", col="cadetblue", breaks = 50)
boxplot(data$Proportion.Heterozygous)

# Number of missing data per genotype
hist(data$Gametes.Missing, main="Histogram of number of missing SNPs", xlab=" Missing SNPs", col="cadetblue")
boxplot(data$Gametes.Missing)

# Proportion of missing data per genotype
hist(data$Proportion.Missing, main="Histogram of proportion of missing SNPs", xlab=" Missing SNPs", col="cadetblue")
boxplot(data$Proportion.Missing)

newdata <- data[order(data$Proportion.Heterozygous),]
