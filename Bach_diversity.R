#Elizabeth Bach, Ryan Williams
#Code to calculate diversity stats for ITS sequences from environmental samples
#4 June 2015

#rm() removes all designated objects, this pervents confusion with any previous data you've been working with before running this code
#if you are begining a new R session, this is not necessary, but it won't hurt anything to run it
rm(list=ls())

#libarary() will open the packages necessary for the commands we will be using in this code
#if you do not have these packages intalled, you may intall through the menu or by running install.pacakges()
#You only need to install packages once, they need only be opened each time you start a new R session
library(lme4)
library(lmerTest)
library(bbmle)
library(reshape)

library(vegan)
library(ggplot2)
library(plyr)

#Read in the dataset, file.choose() will open a window from which you can select your data file (shold be .csv created in
#Bach_ITS_data_wrangling.R), I typically run this code for both rarified and non-rarified data to see if rarification impacted the
#ecological signal
data.nosing<-read.csv(file.choose())
head(data.nosing[,1:10])
str(data.nosing)
dim(data.nosing)


#calculating richness, shannons, and evenness
#Note richness here is Fisher's alpha, which is a logrithmic estimation of number of OTUs, so will likely differ from actual OTU number
#For raw OTU number per sample, see "Bach_ITS1_taxonomy_merge.R"
#Exclude metadata columns from analysis (in this case, columns 1:3), modify for your dataset
richness<-fisher.alpha(data.nosing[,-c(1:3)], 1)
#Preview
head(richness)
shannons<-diversity(data.nosing[,-c(1:3)])
#Preview
head(shannons)
#note, other diversity indices can be calculated using diversity(), use ?diversity() in console to see other options

evenness<-shannons/log(richness)
#Preview
head(evennes)

#bind diversity stats into new data frame with just the summary statistics and metadata
div_stats<-data.frame(data.nosing[,1:3],richness,shannons,evenness)
#Preview
head(div_stats)
str(div_stats)

#looking at data distribution, this will help determine if transformations would be appropriate
ggplot(div_stats)+geom_histogram(aes(shannons))
ggplot(div_stats)+geom_histogram(aes(richness))
ggplot(div_stats)+geom_histogram(aes(evenness))

#Statistical tests
#your data frame must include all independent factors you are interested in, substitute appropriate factors below (e.g. date, crop, system, depth)
#general linear model ANOVA, testing main effects of all factors on diversity measures
#Using the * between factors will result in a full model testing each main effect and all possible interactions
#Using a + between factors will run a main effects model, ignoring any interactions
#You can customize to include only significant interactions and main effects by mixing and matching + and : 
#(e.g. Factor1+Factor2+Factor3+Factor1:Factor2 will evaluate main effects of factors 1,2,3 and the interaction of factors 1,2 only)
summary(test<-aov(richness~Factor1*Factor2, data=div_stats))
TukeyHSD(test)

summary(test2<-aov(evenness~Factor1*Factor2, data=div_stats))
TukeyHSD(test2)

summary(test3<-aov(shannons~Factor1*Factor2, data=div_stats))
TukeyHSD(test3)

#Mixed models for studies with blocking
#null model (no effect of any factors)
summary(test.null<-lmer(richness~(1|), data=div_stats, REML=FALSE))
#main effects model, with block
summary(test.main<-lmer(richness~Factor1+Factor2+Factor3+(1|block), data=div_stats, REML=FALSE))
#full model (all interactions)
summary(test.full<-lmer(richness~Factor1*Factor2*Factor3+(1|block), data=div_stats, REML=FALSE))
#you can remove non-significant interactions and include on the significant ones, eg
summary(test.part<-lmer(richness~Factor1+Factor2+Factor3+Factor1:Factor2+(1|block), data=div_stats, REML=FALSE))
#compare tests to see which model best fits the data (AIC comaprison), best fitting will have lowest AIC (by at least 5 usually)
AIC(test.null, test.main, test.full, test.part)
#Can run ANOVA on models to see if one fits "significantly" better than another
anova(test.null, test.main, test.full, test.part)
#Once best model is selected, run the model to see output and get the statistics for the data
#This is what you will report in a manuscript
summary(test.main)

#Nested factor
#main effects model, with block, Factor 2 is nested within Factor 1
summary(test.main<-lmer(richness~Factor1+Factor2+Factor3+(1|block)+(1|Factor1/Factor2), data=div_stats, REML=FALSE))

#"quick & dirty" graphs
#Richness
rich<-ddply(div_stats, .(Factor1), summarise, .progress="text", mean=mean(richness), SE=sd(richness)/sqrt(length(richness)-1))
ggplot(rich, aes(Factor1, mean))+geom_pointrange(aes(ymax=mean+SE), ymin=(mean-SE))

#can evaluate interactions between factors:
rich<-ddply(div_stats, .(Factor1, Factor2), summarise, .progress="text", mean=mean(richness), SE=sd(richness)/sqrt(length(richness)-1))
ggplot(rich, aes(Factor1, mean, groups=Factor2))+geom_pointrange(aes(ymax=mean+SE), ymin=(mean-SE))

#Evenness
rich<-ddply(div_stats, .(Factor1), summarise, .progress="text", mean=mean(evenness), SE=sd(evenness)/sqrt(length(evenness)-1))
ggplot(rich, aes(Factor1, mean))+geom_pointrange(aes(ymax=mean+SE), ymin=(mean-SE))

#Shannon's
rich<-ddply(div_stats, .(Factor1), summarise, .progress="text", mean=mean(shannons), SE=sd(shannons)/sqrt(length(shannons)-1))
ggplot(rich, aes(Factor1, mean))+geom_pointrange(aes(ymax=mean+SE), ymin=(mean-SE))