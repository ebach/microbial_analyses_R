#Elizabeth Bach, Ryan Williams
#Relative and Rank abundance of OTUs and taxa for microbial community data
#12 June 2015

#rm() removes all designated objects, this pervents confusion with any previous data you've been working with before running this code
#if you are begining a new R session, this is not necessary, but it won't hurt anything to run it
rm(list=ls())

#libarary() will open the packages necessary for the commands we will be using in this code
#if you do not have these packages intalled, you may intall through the menu or by running install.pacakges()
#You only need to install packages once, they need only be opened each time you start a new R session
library(reshape)
library(lme4)
library(lmerTest)
library(bbmle)
library(ggplot2)
library(plyr)

#Use .csv file generated in "Bach_ITS_taxonomy_merge.R" file (0s NOT removed)
data_taxa2<-read.csv(file.choose())
#There should be a reads column, ensure it is title reads (substitute appropriate column number, 11 is used below)
names(data_taxa2)[11]<-paste("reads")
#If working with the rarified data, this will be the number of reads you rarifid to (e.g. 10,000)
#If read column does not exist, generate using ddply OR cbind rarified read total (e.g. 10,000)
#rarified
reads<-rowSums(data_taxa2[,-c(1:8)])
data_taxa2<-cbind(reads,data_taxa2)
#Ensure OTU column is titled OTU
names(data_taxa2)[2]<-paste("OTU")
head(data_taxa2)
str(data_taxa2)

#melt data to preseve taxonomic levels of interest and consolidate, in the example below, we are looking at taxonomic data at all levels from phylum to family, but consolidating genus and species information at the family level
#Generally genus and sepcies data is not always a strong match, especially for ITS, so summaries at family or order level are often most informative
#taxonomic levels can be added or removed below as desired, the names in id= are the column titles, modify as needed to match your data file
data_melt<-melt(data_taxa2[,-c(12,17,18)],id=c("X.2","X.1","X","OTU","Factor1","Factor2","Factor3","block","Sample","reads","phylum","class","order","family"))
data_melt$CropDate<-factor(paste(data_melt$Crop, data_melt$Date, sep=""))
str(data_melt)
head(data_melt)
#Now order the data in each sample by number of reads
#using "-reads" will rank OTUs/taxa from highest number of reads to lowest, this is what we want for rank abundance as rank 1= most reads
data_ordered<-arrange(data_melt, Sample, -reads)
head(data_ordered)
#Use dimensions to ensure data_ordered matches the OTUs and sample size of the orginal fial (total rows/#samples= total OTUs)
dim(data_ordered)

#Relative abundance of taxa
#generate the relative abundance by dividing #reads by total reads per sample
#In this example, we look at relative abundance of faimlies (reads per family per sample include total reads from all OTUs in each family in each sample), this is set up for a dataset rarified to 10000 reads, modify this number based on number of reads your dataset is rarified to
#modify column names as necessary
data.relative<-ddply(data_ordered, .(Sample,OTU, block, Factor1, Factor2, Factor3, phylum,class,order, family), summarize, .progress="text",
relative_read=(reads/10000))

#check, sum of all proportions in each sample is 1
sum_relative<-ddply(data.relative, .(Sample), summarize, .progress="text", sum=sum(relative_read))
head(sum_relative)

#Preview data.relative file
head(data.relative)
max(data.relative$relative_read)
dim(data.relative)

#remove 0s
data.relative2<-subset(data.relative, data.relative$relative_read>0)
dim(data.relative2)
head(data.relative2)

#If desired, you can arrange data to list taxa relative abundance from most abundant to least abundant within each sample
#this is optional
data.relative_order<-arrange(data.relative2, Sample, -relative_read)
head(data.relative_order)

#Now summarize mean and variance of family abundance among Factor levels, you can look at each factor individually, and interacting factors if desired
family.abund<-ddply(data.relative2, .(family, Factor1), summarize, .progress="text",
mean=mean(relative_read))
family.abund_ordered<-arrange(family.abund, Factor1, -mean)
family.abund_ordered

#To summarize abundance at higher taxonomic levels, relativize abundance at the desire level before summarizing
#Example for phylum
data.relative.p<-ddply(data_ordered, .(Sample,OTU, block, Factor1, Factor2, Factor3, phylum), summarize, .progress="text",
relative_read=(reads/10000))

#summarize phyla abundance among factor levels for Factor 1, you can look at each factor individually, and interacting factors if desired
phylum.abund<-ddply(data.relative2, .(phylum, Factor1), summarize, .progress="text",
mean=mean(relative_read*100), n=length(relative_read), SE=sd(relative_read*100)/(sqrt(n-1)))
head(phylum.abund)
phylum.abund_ordered<-arrange(phylum.abund, Factor1, -mean)
phylum.abund_ordered

#This can be easily repeated/modified for any taxonomic level

#Now run models (ANOVA) to determine statistical diffrences in realtive abundance of taxonomic groups among Factor levels
#This is set up for mixed models, including block,  general linear models can be used as well
#First compare null model (no effects), full model (all interactions), main model (no interactions), and any other interations to determine the best fitting model for the data
model.null<-lmer(relative_read~1+(1|block), data=phylum.relative, REML=FALSE)
model.full<-lmer(relative_read~Factor1*Factor2*Factor3+(1|block), data=phylum.relative, REML=FALSE)
model.main<-lmer(relative_read~Factor1+Factor2+Factor3+(1|block), data=phylum.relative, REML=FALSE)
AICtab(model.null,model.full,model.main)
anova(model.null,model.full,model.main)
#Select model with lowest AIC value (usually by at least 3), if 2 models have simlar AIC values and are not residuals are not statistically different in the anova, you can justify using either one (e.g. if main model is as good of fit as full model, no reason to add extra complexity from full model)
#now run the selected model for results
anova(model.main)

#I often use relative abundance to help me understand my data better, figure out which taxa are driving changes between factor levels
#shifts in relative abundance can be reported in text, tables, or visualized using stacked bar graphs
#code for stacked bar graphs under construction!


#Rank abundance, this is often used for visualizing shifts in dominant and sub-domninant taxa
#We will rank each OTU with 1 being highest number of reads, to do this, you must keep all 0's in the data!  We will prune 0s after ranking
#In this example, number of OTUs is 5057 and sample number is 107, change these numbers to match number of OTUs and smaples in your data
data_ordered$rank<-rep(seq(1,5057,1),107)
dim(data_ordered)
head(data_ordered)
#Preview last 5 lines of the data frame
data_ordered[5055:5060,]
#look at range in rank numbers
range(data_ordered$rank)

#If your data is unblanance (potentially from losing samples in rarification), re-casting the data including fill=0 will fill-in the blanks, it will take some tweaking to ensure you have the correct factor levels in the correct spot in the following code, but this will be different for every dataset, so look at the data, figure out which levels are missing (e.g. if only levels from Factor 3 are missing, the below code will work), this may be an area you will want to ask for personal help to make sure the problem is resolved in the appropraite manner
#indicate the appropraite taxonomic level(s) of interest here
data_ordered2<-data.frame(cast(data_ordered, Factor1 + Factor2 + OTU +class ~ Factor3, fill=0))
head(data_ordered2)
str(data_ordered2)
dim(data_ordered2)
#Now rank the data
data_ordered2$rank<-rep(seq(1,5057,1),107)
dim(data_ordered2)
head(data_ordered2)

#Now you may drop 0s
data_ordered3<-droplevels(subset(data_ordered, data_ordered$rank>0))

#Example questions you can answer with this data
#Which taxa have highest change in abundance (delta abundance) between sampling events (by Factor1)?
range(data_ordered$Factor1)
Factor1.level1<-droplevels(subset(data_ordered, data_ordered$Factor1=="level1"))
range(Factor1.level1$Factor1)
dim(Factor1.level1)
head(Factor1.level1)
str(Factor1.level1)
#look at the higest ranked taxa (top 6)
head(arrange(Factor1.level1, Factor1))

#Rank abundance is most commonly used for visualization
#code for stacked bar graphs under construction!