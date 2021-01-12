## Script written by Jonathan Renk
## 7 June 2020

## Clearing the global environment
rm(list=ls(all=TRUE))

## Setting up the working directory
getwd()
setwd("/Users/jonathanrenk/Desktop/EMS /GWAS_Christine_SNPs/")

## Loading packages
library("car")
library("lme4")
library("dfoptim")
library("lmerTest")

## Loading in the data
data <- read.csv("christine_YC_final_common_447.csv", header=T, stringsAsFactors=F)

## Removing SampleID, Analysis.Time, and M.Distance columns
data <- data[,c(-1, -6, -24:-29)]

## Setting variables as factors for genotype, block, rep, and env
data[,1] <- as.factor(data[,1])
data[,2] <- as.factor(data[,2])
data[,3] <- as.factor(data[,3])
data[,4] <- as.factor(data[,4])

str(data)

traits <- colnames(data)[5:ncol(data)]

for(trait in traits){
  pdf(paste0("plots_", trait, ".pdf"), width = 15, height = 5)
  par(mfrow=c(1,2))
  
  # Histogram of raw data
  hist(data[,trait],
       main=paste0("Histogram of raw\n", trait, " values"),
       xlab=trait,
       col="cadetblue")
  
  # Stripchart of raw data across environments
  stripchart(data[,trait] ~data$Env,
             main=paste0("Stripchart of raw\n", trait, " values"),
             xlab="Environment",
             ylab=trait,
             vertical=TRUE, 
             method="jitter",
             bg="cadetblue",
             pch=21
  )
  
  dev.off()
  
  if(paste0("stats_", trait, ".txt") %in% list.files()){
    system(paste0("rm stats_", trait, ".txt"))}
  # Summary statistics of the trait
  summary <- summary(data[,trait], )
  out <- capture.output(summary)
  cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  
  # Test for normality
  normality <- shapiro.test(data[,trait])
  out <- capture.output(normality)
  cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  
  # Run an ANOVA (switched : in Env:Rep and Rep:Block for nesting)
  model <- lm(get(trait) ~ Genotype + Env + Env/Rep + Rep/Block + Genotype:Env, data=data)
  anova <- anova(model)
  out <- capture.output(anova)
  cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  
  # Run a random effects model
  model.1 <- lmer(get(trait) ~ (1|Genotype) + (1|Env) + (1|Env/Rep) + (1|Rep/Block) + (1|Genotype:Env), data = data, REML = TRUE)
  # Decreasing stopping tolerances
  strict_tol <- lmerControl(optCtrl=list(xtol_abs=1e-8, ftol_abs=1e-8))
  if (all(model.1@optinfo$optimizer=="nloptwrap")) {
    model <- update(model.1, control=strict_tol)
  }
  summary(model, correlation=FALSE)
  random_effects <- ranef(model)
  # Write out BLUPs for Genotypes
  write.table(random_effects$Genotype, paste0("blups_", trait, ".csv"), col.names=F, row.names=F, sep=",")
  # Summary of random effects
  summary <- summary(model, correlation=FALSE)
  out <- capture.output(summary)
  cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  # Write out residuals from ANOVA
  write.table(resid(model), paste0("resids_", trait, ".csv"), col.names=F, row.names=F, sep=",")
  # Calculate hertiability 
  model_variances <- as.data.frame(VarCorr(model))
  h2 <- model_variances$vcov[2]/(model_variances$vcov[2]+(model_variances$vcov[1]/5)+(model_variances$vcov[8]/10))
  out <- capture.output(h2)
  cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  
  pdf(paste0("assumptions_", trait, ".pdf"), width = 15, height = 5)
  par(mfrow=c(1,3))
  
  # Model Fit with REML
  plot(fitted(model), residuals(model), pch=19, col="dark blue", ylab="Residuals", xlab="Predicted")
  abline(h=0,col="red", lwd=1, lty=1)
  # histogram of residuals
  hist(residuals(model),main="Histogram of residuals",freq=F, xlab="Residuals", ylab= "Freq", col="palegreen", col.main="darkblue")
  x=seq(-5e-15,9e-15,5e-15)
  curve(dnorm(x,mean(residuals(model)),sd(residuals(model))),add=T,lwd=2, col="red", lty=1)
  # qq plot
  qqPlot(residuals(model), pch=19, col="dark blue", col.lines="red", xlab="Pred quantiles", ylab="Obs quantiles") 
  
  dev.off()
  
}
