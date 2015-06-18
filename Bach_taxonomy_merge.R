#Elizabeth Bach
#Merging OTU table with taxonomy
#4 June 2015

rm(list=ls())

#all useful packages, install.packages(...package name...) if you don't have them
library(labdsv)
library(vegan)
library(plyr)
library(reshape)
library(ggplot2)

#Read in the dataset, file.choose() will open a window from which you can select your data file 
#should be .csv created in Bach_ITS_data_wrangling.R
data.nosing<-read.csv(file.choose())
#Preview file to ensure you have the correct one
str(data.nosing)
dim(data.nosing)
head(data.nosing[,1:10])

#"melt" the OTU table so that OTU IDs are now a column, this will generate a data frame with each combination of OTU and sample present as a row
#it will become a very long dataframe.  This is necessary to preserve sample/OTU information and bind taxonomic data
data_melt<-melt(data.nosing, id=c("X","SampleID"))
head(data_melt)
dim(data_melt)

#now read in taxonomy file, if you've performed an open reference OTU picking, the taxonomy file will likely be
#"rep_set_tax_assignments.txt" file in your output folder
#For assign_taxonomy.py commands a similar text file will be created during that step, also likely called something similar
taxonomy<-read.delim(file.choose())
#Preview file
head(taxonomy)
#Assign column  names if necessary
colnames(taxonomy)<-c("OTU","taxonomy","E","OTU2")
head(taxonomy)
#This code will split-up the taxonomy from a single column into seperate columns for kingdom, phylum, class, etc.
taxonomy$kingdom<-as.factor(unlist((lapply(strsplit(as.character(taxonomy$taxonomy), ";"), "[", 1))))
taxonomy$phylum<-as.factor(unlist(lapply(strsplit(as.character(taxonomy$taxonomy), ";"),"[",2)))
taxonomy$class<-as.factor(unlist(lapply(strsplit(as.character(taxonomy$taxonomy), ";"),"[",3)))
taxonomy$order<-as.factor(unlist(lapply(strsplit(as.character(taxonomy$taxonomy), ";"),"[",4)))
taxonomy$family<-as.factor(unlist(lapply(strsplit(as.character(taxonomy$taxonomy), ";"),"[",5)))
taxonomy$genus<-as.factor(unlist(lapply(strsplit(as.character(taxonomy$taxonomy), ";"),"[",6)))
taxonomy$species<-as.factor(unlist(lapply(strsplit(as.character(taxonomy$taxonomy), ";"),"[",7)))
#Preview to se that each taxonomic column is a "factor"
str(taxonomy)
head(taxonomy)
#Double check the melted data frame above to see title of column with OTU ids
head(data_melt)

#Merge the melted OTU dataframe with the taxonomy dataframe but OTU ID, R will match the OTU ID (named "variable" in the data_melt, and "OTU" in taxonomy), it is ok these names are different, as long as the OTU IDs themselves are identical
#make sure R hasn't added anything like an "X" to the OTU IDs in one of these dataframes
#Bonus, as long as OTU IDs are identical, it does not matter if samples or OTUs are listed in the same order in the two data frames
data_taxa<-merge(data_melt,taxonomy,by.x="variable", by.y="OTU")
#Preview to see if taxonomy was added, dimensions should have same number of rows as data_melt, all we did was add columns
head(data_taxa)
dim(data_taxa)
str(data_taxa)

#write .csv for future use, fill in desired file name and avoid spaces
write.csv(data_taxa, file="ITSdata_taxa.csv")


#Looking at data first
str(data_taxa[,1:10])
#129 OTUs
levels(data_taxa$phylum)

#Some examples of quick ways to preview the data and get a sense of taxonomic distribution, this is just for data exploration

#OTUs per sample
#Remove OTUxSample combinations that do not have any reads (e.g. value=0), these exist because melting the dataframe puts every OTU in every sample, even if that OTU only occurs in a few samples
data_taxa2<-droplevels(subset(data_taxa, data_taxa$value>0))
#Summarizes number of OTUs in each sample
# You may also modify this code to a factor level, such as ecosystem, landuse, host species, sample date, etc
sample.richness<-ddply(data_taxa2, .(SampleID), summarise, .progress="text", OTUr=length(levels(droplevels(variable))))
#check dimensions before previewing, just to be sure you can see full dataframe in the console
dim(sample.richness)
sample.richness
#Check for min and max "richness" (number of OTUs)
min(sample.richness$OTUr)
max(sample.richness$OTUr)

#List out OTUs per sample (feasible for small datasets <50 samples, <100 OTUs)
# You may also modify this code to a factor level, such as ecosystem, landuse, host species, sample date, etc
sample.1<-droplevels(subset(data_taxa2, data_taxa2$SampleID=="Sample.1"))
sample.1
sample.2<-droplevels(subset(data_taxa2, data_taxa2$SampleID=="Sample.2"))
sample.2
comp.1<-droplevels(subset(data_taxa2, data_taxa2$SampleID=="Sample.1.comp"))
comp.1
comp.2<-droplevels(subset(data_taxa2, data_taxa2$SampleID=="Sample.2.comp"))
comp.2
R47.fr5.comp<-droplevels(subset(data_taxa2, data_taxa2$SampleID=="R47.fr5.comp"))
R47.fr5.comp
R48.fr5.comp<-droplevels(subset(data_taxa2, data_taxa2$SampleID=="R48.fr5.comp"))
R48.fr5.comp