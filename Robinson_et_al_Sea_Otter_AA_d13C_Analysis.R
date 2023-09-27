#####################################################################################################################
## R script for plots, linear mixed effects models, linear regressions, linear discriminant analysis  ###############
## Robinson et al. 2023 Oecologia                                                                     ###############
## Title: "Tissue-Specific Carbon Isotope Patterns of Amino Acids in Southern Sea Otters"             ###############
## This R code runs all analyses presented in the manuscript                                          ###############
## Authors: A Robinson, alanalea@unm.edu; EA Elliott Smith; AC Besser 09/10/2023                      ###############
#####################################################################################################################


# load the libraries into R
library(lme4)
library(lmerTest)
library(reshape2)
library (emmeans)
library(multcomp)
library(ggplot2)
library(MASS)
library(plyr)
library(reshape2)
library(stats)

# load data
AA.Tissue.Ott <- read.csv("Robinson_et_al_Sea_Otter_AA_d13C_Data.csv", header=T)

str(AA.Tissue.Ott)
# removing missing data and data with high standard deviations
AA.Tissue.Ott$Tyr13C <- NULL
AA.Tissue.Ott$Lys13C <- NULL

# removing metadata columns that will not be used for statistical analyses
AA.Tissue.Ott$Date <- NULL
AA.Tissue.Ott$Condition <- NULL
AA.Tissue.Ott$Age <- NULL

# structuring the data frame
AA.Tissue.Ott2 <- melt(AA.Tissue.Ott, id=c("Individual","Tissue", "Sex"))
colnames(AA.Tissue.Ott2) <- c("Indiv","Tissue", "Sex", "Amino", "d13C")

# set sex and tissue as factors
AA.Tissue.Ott2$Sex <- as.factor(AA.Tissue.Ott2$Sex)
AA.Tissue.Ott2$Tissue <- as.factor(AA.Tissue.Ott2$Tissue)
str(AA.Tissue.Ott2) # checking structure


###########################################
##   Run Linear Mixed Effects Models     ##
###########################################


# Do we have a multi-level model? 
Ott.Null <- lmer(d13C ~ 1 + (1|Amino:Tissue), 
                 data=AA.Tissue.Ott2, REML=FALSE)
summary(Ott.Null)

# checking the intraclass correlation (ICC) to determine if multi-level modeling is the correct choice for our analysis
# the ICC measures the degree of clustering in our data 
# An ICC greater than 0 indicates a multi-level study
ICC.Model<-function(Ott.Null) {
  tau.Null<-as.numeric(lapply(summary(Ott.Null)$varcor, diag))
  sigma.Null <- as.numeric(attr(summary(Ott.Null)$varcor, "sc")^2)
  ICC.Null <- tau.Null/(tau.Null+sigma.Null)
  return(ICC.Null)
}
round(ICC.Model(Ott.Null), 2)

# Model 1 - including sex, AA, and tissue as fixed effects
# (1|ID)... how individual means vary from the global mean
# (1 + Amino|Indiv)... looking at the effect of AA for each individual and how this varies from the global effect of AA
# (1 + Tissue|Indiv)... looking at the tissue effect for each individual and how this deviates from the global effect of tissue
Model.Ott <- lmer(d13C ~ Sex*Amino*Tissue + (1|Indiv) + (1|Amino:Indiv) + (1|Tissue:Indiv), 
                   data=AA.Tissue.Ott2, REML=FALSE)
summary(Model.Ott, ddf = "Satterthwaite") # this adds the p-values (if you have installed and loaded the "lmerTest" package)

hist(resid(Model.Ott))
shapiro.test(resid(Model.Ott))

# choose a model by AIC using a stepwise algorithm
step(Model.Ott)

# step function removes sex from the model
# model found:
  # d13C ~ Amino + Tissue + (1 | Indiv) + (1 | Amino:Indiv) + (1 | Tissue:Indiv) + Amino:Tissue


# Model 2 - setting AA and tissue as fixed effects
# (1|ID)... how individual means vary from the global mean
# (1 + Amino|Indiv)... looking at the effect of AA for each individual and how this varies from the global effect of AA
# (1 + Tissue|Indiv)... looking at the tissue effect for each individual and how this deviates from the global effect of tissue

Model.Ott2 <- lmer(d13C ~ Amino*Tissue + (1|Indiv) + (1|Indiv:Amino) + (1|Indiv:Tissue), 
                  data=AA.Tissue.Ott2, REML=F)
summary(Model.Ott2, ddf = "Satterthwaite") # this adds the p-values (if you have installed and loaded the "lmerTest" package)

# choose a model by AIC using a stepwise algorithm
step(Model.Ott2)

# model found:
 # d13C ~ Amino * Tissue + (1 | Indiv) + (1 | Indiv:Amino) + (1 | Indiv:Tissue)
# indicates Model.Ott2 best fits the data and cannot be refined further

# Model.Ott3 <- lmer(d13C ~ Amino * Tissue + (1 | Indiv) + (1 | Indiv:Amino) + (1 | Indiv:Tissue), 
                  # data=AA.Tissue.Ott2, REML=F)
# summary(Model.Ott3, ddf = "Satterthwaite")
# same Results as Model.Ott 2 


#################################
## Proceeding with Model.Ott2 ##
#################################

# tests for normality and homogeneity of variance on Model.Ott2
# sourced from https://ademos.people.uic.edu/Chapter18.html


# checking for linearity of residuals compared to original observations
Plot.Model.Ott2.Linearity<-plot(resid(Model.Ott2),AA.Tissue.Ott2$d13C )
# residual pattern appears random therefore the assumption of linearity has not been violated


# checking for normality of residuals
require("lattice")
qqmath(Model.Ott2, id=0.05)
# the residuals fit the normal line except at the tails
hist(resid(Model.Ott2))
# visual test for normality of residuals 
shapiro.test(resid(Model.Ott2))
# p-value = 3.64e-10 indicates data are not normally distributed 

# test for homogeneity of variance
AA.Tissue.Ott2$Model.Ott.Res<- residuals(Model.Ott2) # extracts the residuals and places them in a new column in our original data table
AA.Tissue.Ott2$Abs.Model.Ott.Res <-abs(AA.Tissue.Ott2$Model.Ott.Res) # creates a new column with the absolute value of the residuals
AA.Tissue.Ott2$Model.Ott.Res2 <- AA.Tissue.Ott2$Abs.Model.Ott.Res^2 # squares the absolute values of the residuals to provide the more robust estimate
Levene.Model.Ott <- lm(Model.Ott.Res2 ~ Indiv, data=AA.Tissue.Ott2) # ANOVA of the squared residuals
anova(Levene.Model.Ott) # displays the results
# p-value is > 0.05 which means the data varies equally among groups
# the assumption of homoscedasticity is met 

# visual analysis of homogeneity of variance 
Plot.Model.Ott2 <- plot(Model.Ott2) # creates a fitted versus residual plot
Plot.Model.Ott2
# trumpet shape indicates there is greater variance in some parts of the data than others 
# this indicates that the variability and our predictor are not totally independent

########################################
## Pull out the results for Table 1 ####
########################################

# pull out pairwise comparisons for all AA and tissues
emmeans(Model.Ott2, ~ Tissue|Amino , adjust = "fdr")
Model.Ott2.results <- emmeans(Model.Ott2, ~ Tissue|Amino , adjust = "fdr")

# transform results into a dataframe and export as csv
pairs(Model.Ott2.results)
export.results.2 <- pairs(Model.Ott2.results)
Model.Ott2.Results = as.data.frame(export.results.2)
write.csv(Model.Ott2.Results, "LinearModelResults.csv")


###########################
## Boxplot - Figure 2 #####
###########################


# AA d13C boxplot - separate tissues within each AA measured
c.t <- ggplot(AA.Tissue.Ott2, aes(x = Amino, y= d13C, fill = Tissue, color = Tissue))
c.t <- c.t + geom_boxplot() + ylim(-30,10)
print(c.t)

# significant differences between tissue types added in Adobe Illustrator by authors - using "LinearModelResults.csv"

#####################################################################
## Linear Regression analysis - Figure 3, Figure S2, Figure S4 ######
#####################################################################


# choose a new folder for PDFs to be exported into 
setwd()

# seperate dataset by tissue type
AA.Tissue.Ott.L <- AA.Tissue.Ott [AA.Tissue.Ott$Tissue == "LIVER",]
AA.Tissue.Ott.M <- AA.Tissue.Ott [AA.Tissue.Ott$Tissue == "MUSCLE",]
AA.Tissue.Ott.W <- AA.Tissue.Ott [AA.Tissue.Ott$Tissue == "WHISKER",]
AA.Tissue.Ott.B <- AA.Tissue.Ott [AA.Tissue.Ott$Tissue == "BONE",]

###########################
##Figure 3 and Figure S2 ##
###########################


# Threonine vs Glycine 
# correlation tests that show equations, R^2, and p-values for each tissue
# cor (x,y) Threonine = x, Glycine =y
cor(AA.Tissue.Ott.B$Gly13C, AA.Tissue.Ott.B$Thr13C)
BoneGlyThr.LM <- lm(AA.Tissue.Ott.B$Gly13C ~ AA.Tissue.Ott.B$Thr13C)
print(BoneGlyThr.LM)
summary(BoneGlyThr.LM)
confint(BoneGlyThr.LM, level = 0.95)

cor(AA.Tissue.Ott.L$Gly13C, AA.Tissue.Ott.L$Thr13C)
LiverGlyThr.LM <- lm(AA.Tissue.Ott.L$Gly13C ~ AA.Tissue.Ott.L$Thr13C)
print(LiverGlyThr.LM)
summary(LiverGlyThr.LM)
confint(LiverGlyThr.LM, level = 0.95)

cor(AA.Tissue.Ott.M$Gly13C, AA.Tissue.Ott.M$Thr13C)
MuscleGlyThr.LM <- lm(AA.Tissue.Ott.M$Gly13C ~ AA.Tissue.Ott.M$Thr13C)
print(MuscleGlyThr.LM)
summary(MuscleGlyThr.LM)
confint(MuscleGlyThr.LM, level = 0.95)

cor(AA.Tissue.Ott.W$Gly13C, AA.Tissue.Ott.W$Thr13C)
WhiskerGlyThr.LM <- lm(AA.Tissue.Ott.W$Gly13C ~ AA.Tissue.Ott.W$Thr13C)
print(WhiskerGlyThr.LM)
summary(WhiskerGlyThr.LM)
confint(WhiskerGlyThr.LM, level = 0.95)

# graph of linear regressions above (no equations shown)
pdf("LMThrGly.pdf")
ggplot(AA.Tissue.Ott, aes(x=Thr13C, y=Gly13C, shape=Tissue, colour = Tissue))+
  geom_smooth(method = "lm", se=FALSE) +
  geom_point(size=3) +
  scale_shape_manual(values = c(15,17,19,8))+
  theme_classic()+
  xlab("Threonine")+
  ylab("Glycine")
dev.off()

# Serine vs Glycine
# correlation tests that show equations, R^2, and p-values for each tissue
# cor (x,y) Glycine = x, Serine =y
cor(AA.Tissue.Ott.B$Gly13C, AA.Tissue.Ott.B$Ser13C)
BoneGlySer.LM <- lm(AA.Tissue.Ott.B$Gly13C ~ AA.Tissue.Ott.B$Ser13C)
print(BoneGlySer.LM)
summary(BoneGlySer.LM)
confint(BoneGlySer.LM, level = 0.95)

cor(AA.Tissue.Ott.L$Gly13C, AA.Tissue.Ott.L$Ser13C)
LiverGlySer.LM <- lm(AA.Tissue.Ott.L$Gly13C ~ AA.Tissue.Ott.L$Ser13C)
print(LiverGlySer.LM)
summary(LiverGlySer.LM)
confint(LiverGlySer.LM, level = 0.95)

cor(AA.Tissue.Ott.M$Gly13C, AA.Tissue.Ott.M$Ser13C)
MuscleGlySer.LM <- lm(AA.Tissue.Ott.M$Gly13C ~ AA.Tissue.Ott.M$Ser13C)
print(MuscleGlySer.LM)
summary(MuscleGlySer.LM)
confint(MuscleGlySer.LM, level = 0.95)

cor(AA.Tissue.Ott.W$Gly13C, AA.Tissue.Ott.W$Ser13C)
WhiskerGlySer.LM <- lm(AA.Tissue.Ott.W$Gly13C ~ AA.Tissue.Ott.W$Ser13C)
print(WhiskerGlySer.LM)
summary(WhiskerGlySer.LM)
confint(WhiskerGlySer.LM, level = 0.95)

# graph of linear regressions above (no equations shown)
pdf("LMSerGly.pdf")
ggplot(AA.Tissue.Ott, aes(x=Ser13C, y=Gly13C, shape=Tissue, colour = Tissue))+
  geom_smooth(method = "lm", se=FALSE) +
  geom_point(size=3) +
  scale_shape_manual(values = c(15,17,19,8))+
  theme_classic()+
  xlab("Serine")+
  ylab("Glycine")
dev.off()

# Serine vs Threonine
# correlation tests that show equations, R^2, and p-values for each tissue
# cor (x,y) Threonine = x, Serine =y
cor(AA.Tissue.Ott.B$Ser13C, AA.Tissue.Ott.B$Thr13C)
BoneSerThr.LM <- lm(AA.Tissue.Ott.B$Ser13C ~ AA.Tissue.Ott.B$Thr13C)
print(BoneSerThr.LM)
summary(BoneSerThr.LM)
confint(BoneSerThr.LM, level = 0.95)

cor(AA.Tissue.Ott.L$Ser13C, AA.Tissue.Ott.L$Thr13C)
LiverSerThr.LM <- lm(AA.Tissue.Ott.L$Ser13C ~ AA.Tissue.Ott.L$Thr13C)
print(LiverSerThr.LM)
summary(LiverSerThr.LM)
confint(LiverSerThr.LM, level = 0.95)

cor(AA.Tissue.Ott.M$Ser13C, AA.Tissue.Ott.M$Thr13C)
MuscleSerThr.LM <- lm(AA.Tissue.Ott.M$Ser13C ~ AA.Tissue.Ott.M$Thr13C)
print(MuscleSerThr.LM)
summary(MuscleSerThr.LM)
confint(MuscleSerThr.LM, level = 0.95)

cor(AA.Tissue.Ott.W$Ser13C, AA.Tissue.Ott.W$Ser13C)
WhiskerSerThr.LM <- lm(AA.Tissue.Ott.W$Ser13C ~ AA.Tissue.Ott.W$Thr13C)
print(WhiskerSerThr.LM)
summary(WhiskerSerThr.LM)
confint(WhiskerSerThr.LM, level = 0.95)

# graph of linear regressions above (no equations shown)
pdf("LMThrSer.pdf")
ggplot(AA.Tissue.Ott, aes(x=Thr13C, y=Ser13C, shape=Tissue, colour = Tissue))+
  geom_smooth(method = "lm", se=FALSE) +
  geom_point(size=3) +
  scale_shape_manual(values = c(15,17,19,8))+
  theme_classic()+
  xlab("Threonine")+
  ylab("Serine")
dev.off()

# Proline vs Glutamic acid
# correlation tests that show equations, R^2, and p-values for each tissue
# cor (x,y) Glutamic acid = x, Proline =y
cor(AA.Tissue.Ott.B$Glu13C,AA.Tissue.Ott.B$Pro13C)
BoneProGlu.LM <- lm(AA.Tissue.Ott.B$Pro13C ~ AA.Tissue.Ott.B$Glu13C)
summary(lm(AA.Tissue.Ott.B$Pro13C ~ AA.Tissue.Ott.B$Glu13C,))
print(BoneProGlu.LM)
summary(BoneProGlu.LM)
confint(BoneProGlu.LM, level = 0.95)

cor(AA.Tissue.Ott.L$Pro13C, AA.Tissue.Ott.L$Glu13C)
LiverProGlu.LM <- lm(AA.Tissue.Ott.L$Pro13C ~ AA.Tissue.Ott.L$Glu13C)
print(LiverProGlu.LM)
summary(LiverProGlu.LM)
confint(LiverProGlu.LM, level = 0.95)

cor(AA.Tissue.Ott.M$Pro13C, AA.Tissue.Ott.M$Glu13C)
MuscleProGlu.LM <- lm(AA.Tissue.Ott.M$Pro13C ~ AA.Tissue.Ott.M$Glu13C)
print(MuscleProGlu.LM)
summary(MuscleProGlu.LM)
confint(MuscleProGlu.LM, level = 0.95)

cor(AA.Tissue.Ott.W$Pro13C, AA.Tissue.Ott.W$Glu13C)
WhiskerProGlu.LM <- lm(AA.Tissue.Ott.W$Pro13C ~ AA.Tissue.Ott.W$Glu13C)
print(WhiskerProGlu.LM)
summary(WhiskerProGlu.LM)
confint(WhiskerProGlu.LM, level = 0.95)

# graph of linear regressions above (no equations shown) 
pdf("LMProGlu.pdf")
ggplot(AA.Tissue.Ott, aes(x=Glu13C, y=Pro13C, shape=Tissue, colour = Tissue))+
  geom_smooth(method = "lm", se=FALSE) +
  geom_point(size=3) +
  scale_shape_manual(values = c(15,17,19,8))+
  theme_classic()+
  xlab("Glutamic acid")+
  ylab("Proline")
dev.off()

# Glutamic acid vs Aspartic acid 
# correlation tests that show equations, R^2, and p-values for each tissue
# cor (x,y) Glutamic acid = x, Aspartic acid =y
cor(AA.Tissue.Ott.B$Glu13C, AA.Tissue.Ott.B$Asp13C)
BoneGluAsp.LM <- lm(AA.Tissue.Ott.B$Glu13C ~ AA.Tissue.Ott.B$Asp13C)
print(BoneGluAsp.LM)
summary(BoneGluAsp.LM)
confint(BoneGluAsp.LM, level = 0.95)

cor(AA.Tissue.Ott.L$Glu13C, AA.Tissue.Ott.L$Asp13C)
LiverGluAsp.LM <- lm(AA.Tissue.Ott.L$Glu13C ~ AA.Tissue.Ott.L$Asp13C)
print(LiverGluAsp.LM)
summary(LiverGluAsp.LM)
confint(LiverGluAsp.LM, level = 0.95)

cor(AA.Tissue.Ott.M$Glu13C, AA.Tissue.Ott.M$Asp13C)
MuscleGluAsp.LM <- lm(AA.Tissue.Ott.M$Glu13C ~ AA.Tissue.Ott.M$Asp13C)
print(MuscleGluAsp.LM)
summary(MuscleGluAsp.LM)
confint(MuscleGluAsp.LM, level = 0.95)

cor(AA.Tissue.Ott.W$Glu13C, AA.Tissue.Ott.W$Asp13C)
WhiskerGluAsp.LM <- lm(AA.Tissue.Ott.W$Glu13C ~ AA.Tissue.Ott.W$Asp13C)
print(WhiskerGluAsp.LM)
summary(WhiskerGluAsp.LM)
confint(WhiskerGluAsp.LM, level = 0.95)

# graph of linear regressions above (no equations shown)
pdf("LMGluAsp.pdf")
ggplot(AA.Tissue.Ott, aes(x=Glu13C, y=Asp13C, shape=Tissue, colour = Tissue))+
  geom_smooth(method = "lm", se=FALSE) +
  geom_point(size=3) +
  scale_shape_manual(values = c(15,17,19,8))+
  theme_classic()+
  xlab("Glutamic acid")+
  ylab("Aspartic acid")
dev.off()

###############
## Figure S4 ##
###############

# please note Figure S4 was made in Microsoft Excel - only regression equations will be shown below


# Threonine vs Glycine 
# inter-tissue 
## cor (x,y) Threonine in liver = x, Glycine in liver, whisker, bone, muscle = y
cor(AA.Tissue.Ott.L$Gly13C, AA.Tissue.Ott.L$Thr13C)
LiverGlyThr.LM <- lm(AA.Tissue.Ott.L$Gly13C ~ AA.Tissue.Ott.L$Thr13C)
print(LiverGlyThr.LM)
summary(LiverGlyThr.LM)
confint(LiverGlyThr.LM, level = 0.95)

cor(AA.Tissue.Ott.W$Gly13C, AA.Tissue.Ott.L$Thr13C)
WLGlyThr.LM <- lm(AA.Tissue.Ott.W$Gly13C ~ AA.Tissue.Ott.L$Thr13C)
print(WLGlyThr.LM)
summary(WLGlyThr.LM)
confint(WLGlyThr.LM, level = 0.95)

cor(AA.Tissue.Ott.B$Gly13C, AA.Tissue.Ott.L$Thr13C)
BLGlyThr.LM <- lm(AA.Tissue.Ott.B$Gly13C ~ AA.Tissue.Ott.L$Thr13C)
print(BLGlyThr.LM)
summary(BLGlyThr.LM)
confint(BLGlyThr.LM, level = 0.95)

cor(AA.Tissue.Ott.M$Gly13C, AA.Tissue.Ott.L$Thr13C)
MLGlyThr.LM <- lm(AA.Tissue.Ott.M$Gly13C ~ AA.Tissue.Ott.L$Thr13C)
print(MLGlyThr.LM)
summary(MLGlyThr.LM)
confint(MLGlyThr.LM, level = 0.95)

# Serine vs Threonine
# inter-tissue
## cor (x,y) Threonine in liver = x, Serine in liver, whisker, bone, muscle = y
cor(AA.Tissue.Ott.L$Ser13C, AA.Tissue.Ott.L$Thr13C)
LiverSerThr.LM <- lm(AA.Tissue.Ott.L$Ser13C ~ AA.Tissue.Ott.L$Thr13C)
print(LiverSerThr.LM)
summary(LiverSerThr.LM)
confint(LiverSerThr.LM, level = 0.95)

cor(AA.Tissue.Ott.W$Ser13C, AA.Tissue.Ott.L$Thr13C)
WLSerThr.LM <- lm(AA.Tissue.Ott.W$Ser13C ~ AA.Tissue.Ott.L$Thr13C)
print(WLSerThr.LM)
summary(WLSerThr.LM)
confint(WLSerThr.LM, level = 0.95)

cor(AA.Tissue.Ott.B$Ser13C, AA.Tissue.Ott.L$Thr13C)
BLSerThr.LM <- lm(AA.Tissue.Ott.B$Ser13C ~ AA.Tissue.Ott.L$Thr13C)
print(BLSerThr.LM)
summary(BLSerThr.LM)
confint(BLSerThr.LM, level = 0.95)

cor(AA.Tissue.Ott.M$Ser13C, AA.Tissue.Ott.L$Thr13C)
MLSerThr.LM <- lm(AA.Tissue.Ott.M$Ser13C ~ AA.Tissue.Ott.L$Thr13C)
print(MLSerThr.LM)
summary(MLSerThr.LM)
confint(MLSerThr.LM, level = 0.95)

# end of inter-tissue regressions

#######################
## LDA and Figure 4  ##
#######################

# pull in source AA data from Elliott Smith, Emma A. et al. "Amino acid d13C dataset for nearshore marine primary producers" (2022)
# https://doi.org/10.5061/dryad.zw3r22897

sources <- read.csv(file="chap4_Final_Dryad.csv", header=T)
str(sources)

# format dataframe structure
sources$Gly13C <- as.numeric(as.character(sources$Gly13C))
sources$Ser13C <- as.numeric(as.character(sources$Ser13C))
sources$Pro13C <- as.numeric(as.character(sources$Pro13C))
sources$Asp13C <- as.numeric(as.character(sources$Asp13C))
sources$Glu13C <- as.numeric(as.character(sources$Glu13C))
sources$Phe13C <- as.numeric(as.character(sources$Phe13C))
sources$Lys13C <- as.numeric(as.character(sources$Lys13C))
sources$Bulk.d13C <- as.numeric(as.character(sources$Bulk.d13C))
sources$Bulk.d15N <- as.numeric(as.character(sources$Bulk.d15N))

# remove samples with missing AA data
sources_clean_EAA <- sources[!(sources$ID == 'KATM-AM-ULV-2012-05' |            # missing Phe
                             sources$ID == 'HUM-ISM6-FITO' |                    # missing Lys
                             sources$ID == 'KATM-AM-LAM-2012-03'),]             # missing Phe
sources_clean_EAA2 <- sources_clean_EAA

#Correcting Phylum designations; see notes column "Analyzed jointly with POM samples due to isotopic similarity" and Elliott Smith et al. 2022 for more details 
colnames(sources_clean_EAA2)[22]= "Phylum2"
sources_clean_EAA2["Phylum2"][sources_clean_EAA2["Phylum2"] == "Unknown_POM"] <- "POM_Phyto"
sources_clean_EAA2["Phylum2"][sources_clean_EAA2["Phylum2"] == "Bacillariophyta"] <- "POM_Phyto"
sources_clean_EAA2["Phylum2"][sources_clean_EAA2["Phylum2"] == "Cyanobacteria"] <- "POM_Phyto"
sources_clean_EAA2["Phylum2"][sources_clean_EAA2["Phylum2"] == "Haptophyta"] <- "POM_Phyto"
sources_clean_EAA2[(sources_clean_EAA2$ID == "EES-SCRIPPS-PORPHYRIDIUM-18"), "Phylum2"] <- "POM_Phyto"

# running a training LDA with these data
sources.lda.nolys <-lda(Phylum2 ~ Ile13C + Leu13C + Phe13C + Thr13C + Val13C, data = sources_clean_EAA2)
sources.lda.cv.nolys <-lda(Phylum2 ~ Ile13C + Leu13C + Phe13C + Thr13C + Val13C, data = sources_clean_EAA2, CV = T)
ct.prod.2.nolys <- table (sources_clean_EAA2$Phylum2, sources.lda.cv.nolys$class)
sum(diag(prop.table(ct.prod.2.nolys)))

# now pull out just the essential AA data from sea otter tissue dataset:
# loading in corrected data for UCM standard - see SI1 for details. 
AA.Tissue.Ott <- read.csv("Corrected_LDA_Data.csv",header=T)
# removing variable Lys data
Ess.nolys <- subset(AA.Tissue.Ott, select = c(Individual, Ile13C, Leu13C, Phe13C, Thr13C, Val13C, Tissue))


# now let's predict the sources sea otters classify with and save this into a dataframe
lda.nolys <- predict(sources.lda.nolys, Ess.nolys)
AA.Tissue.Ott.ess.lda.nolys<- merge(AA.Tissue.Ott, lda.nolys, by.x = 0, by.y = 0) # adding the classifications to the data table - this code merges by the unique ROW ID
write.csv (AA.Tissue.Ott.ess.lda.nolys, "LDA_Classifications.csv")


# Part 2 - graphing the LDA ####
# creating dataframes for plotting
# producer dataframe with LD coordinates
datPred.prod.nolys <- data.frame(Species = sources_clean_EAA2$Phylum2, # create dataframe with producer LD values
                                 predict(sources.lda.nolys)$x, # pulling out producer LD coordinates for plotting
                                 spp = sources_clean_EAA2$Phylum2,  # including producer taxonomic identity/grouping
                                 ID = c('producer'), 
                                 size = c("prod"))   # label for plotting

# consumer dataframe with LD coordinates
datPred.con.nolys <- data.frame(Species = c('Otter'),  # generic label for plotting
                                lda.nolys$x, # pulling out the LD coordinates for consumers
                                spp = Ess.nolys$Individual, # Individual grouping
                                ID = Ess.nolys$Tissue, # separate by tissues
                                size = c("con")) # label for plotting
datPred.all.nolys <- rbind(datPred.prod.nolys, datPred.con.nolys) # merging producer and consumer LD dataframes 

# plotting the fingerprints 
sri.lda.nolys <- ggplot(datPred.all.nolys, aes(x=LD1, y=LD2, col = spp, shape = ID) ) # colors will be individuals (or producer groups) and shape will be tissues (or producer groups)
sri.lda.nolys <- sri.lda.nolys + theme_bw()
sri.lda.nolys <- sri.lda.nolys + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sri.lda.nolys <- sri.lda.nolys + geom_point(aes(size = size)) + scale_size_discrete(range = c(6,3))
sri.lda.nolys <- sri.lda.nolys + guides(color = guide_legend(override.aes = list(size = 4)))
sri.lda.nolys <- sri.lda.nolys + guides(shape = guide_legend(override.aes =list(size = 4)))
sri.lda.nolys <- sri.lda.nolys + scale_colour_manual(values = c("steelblue1", "dodgerblue2" , "darkorchid4", "turquoise2", "hotpink", "deeppink3", "lightpink1", "firebrick2", "coral2","slateblue3",
                                                                "grey36", "black", "grey36", "grey75")) 
sri.lda.nolys <- sri.lda.nolys + scale_shape_manual(values = c(15, 17, 18, 16, 8))
sri.lda.nolys <- sri.lda.nolys + ylim(-5,5)
sri.lda.nolys <- sri.lda.nolys + xlim(-5,6)
print(sri.lda.nolys)
# end LDA

##############################################################################################################################################################################################################





