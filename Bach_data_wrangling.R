#Elizabeth Bach, Ryan Williams
#data wrangling for OTU table generated in QIIME
#This code removes singltons and does rarification, generating a .csv file that can be used for subsequent analyses
#This code needs only to be run once
#4 June 2015

#rm() removes all designated objects, this pervents confusion with any previous data you've been working with before running this code
#if you are begining a new R session, this is not necessary, but it won't hurt anything to run it
rm(list=ls())

#libarary() will open the packages necessary for the commands we will be using in this code
#if you do not have these packages intalled, you may intall through the menu or by running install.pacakges()
#You only need to install packages once, they need only be opened each time you start a new R session
library(labdsv)
library(vegan)
library(plyr)
library(reshape)
library(ggplot2)

#Read in the dataset, file.choose() will open a window from which you can select your data file
#use the tab-delimated text file OTU table generated in QIIME (should be .txt file)
dataset<-read.delim(file.choose(), header=TRUE, skip=1)
#Look at the dimensions of the data frame (rows=samples, columns=OTUs)
dim(dataset)
#Look at the first 6 rows of the dataset, ensure it looks like what you think it should look like
head(dataset)
#check how R is interperting each column (number, character, factor, etc.)
str(dataset)
#you can look at x rows and y columns in our object "dataset" using dataset[x,y], 
#below we look at rows 1-10 and columns 5-9 in object "dataset"
dataset[1:10,5:9]

#transform data frame so it is sample by OTU, not OTU by sample

#creating object "SampleID" with the sample names, which are the column names in the "dataset" object
SampleID<-matrix(colnames(dataset))
#binding "SampleID" with OTU IDs to preserve names and ensure names are assigned as a character rather than number as.character()
#You may need to change the "X.OTU.ID" portion after the $ to the correct name of the OTU name column
OTU<-c("SampleID", as.character(dataset$X.OTU.ID))
#transpose "dataset", excluding column 1 ("-1"), which are the OTU names, R often misinterprets transposed character values
dataset2<-data.frame(t(dataset[,-1]))
#bind the transposed "dataset" (dataset2) with the "SampleID" object created above to include sample names as the first column
dataset3<-cbind(SampleID[-1,], dataset2)
#Set OTU IDs as column headers using the OTU object created above
colnames(dataset3)<-OTU
#preview first 10 row and 10 columns of "dataset3" to ensure samples are rows, OTUs are columns
dataset3[1:10,1:10]
#check the dimensions, this should be the converse of the dimensions of the original "dataset" 
#(note, you may have an additional column due to rebinding of SampleID as a column rather than just rownames)
dim(dataset3)
#Check that R is interpreting datset3 as desired (factors are factors, OTUs are numbers, etc.)
str(dataset3[1:8,1:10])

# First we remove the singletons using the dropspc() function form the labdsv package.  In the line below I bind metadata (columns 1 through 5) 
#you will need to change this based on how many columns of metadata your table may include 
#Note that for the dropspc function I exclude columns 1:5 so only the OTU columns are considered (no metadata) 
#you will need to modify these numbers for your specific dataset
#If you used split_libraries in QIIME, it may alread have eleminated singtons
#Here I'm indicating to exclude OTUs with fewer than 5 reads, this cut-off can be changed as you wish
data.nosing<-cbind(dataset3[1:5],dropspc(dataset3[,-1:5],5))
#Preview data.nosing to see if it looks right and how many (if any) OTUs were excluded
str(data.nosing)
dim(data.nosing)
data.nosing[1:8,1:10]

# Now its time to figure out how many reads to rarefy by...I added a column to our dataset of the total number of reads per sample (row)
#again excluding an metadata columns
reads<-rowSums(data.nosing[,-c(1)])
#preview "reads"
head(reads)
#bind it to the OTU table created above (data.nosing)
data.nosing.reads<-cbind(reads,data.nosing)
#Preview
head(data.nosing.reads[,1:10])
max(data.nosing.reads$reads)
min(data.nosing.reads$reads)

#generate .csv of non-rarified data for analysis and/or reference
#Change the text in file="" to what you wish to name your file, avoid using spaces
write.csv(data.nosing.reads, file="Project_ITS1_data_nonrar.csv")

#Rarification
# lets create rarefaction curves for each sample starting around 10000
rared<-rarefy(subset(data.nosing.reads, reads > 10000)[,-1],sample=c(1,10,25,50,75,100,250,500,700,1250,2500,5000),se=FALSE )
rared_melt<-melt(rared)
names(rared_melt)<-c("sample","sample_size","OTUs")
#in this example we have 6 samples, so each read value should be evaluated 6 times, once for each sample, 
#replace the 6s below with the number of samples in your dataset
rared_melt$sample_size<-c(rep(1,6),rep(10,6),rep(25,6),rep(50,6),rep(75,6),rep(100,6),rep(250,6),rep(500,6),rep(700,6),rep(1250,6),rep(2500,6),rep(5000,6))
head(rared_melt)
#this will generate the plot with the curves
ggplot(rared_melt)+geom_line(aes(x=sample_size,y=OTUs,colour=sample,group=sample))+theme(aspect.ratio=1)+theme_bw()

#You can copy and paste the above code and change the read threshold numer (to 1000, 5000, 100000, whatever makes sense in your dataset)
#You may want to look at a several different thresholds to see where curves are begning to level off (e.g. few new OTUs added with additional reads)
#You can also do additional sensitivity analyses to help make a decision about the best rarification cut-off, see Williams in prep

# In this example, I will decide 5000 reads is a good cut-off
#Again change the column numbers so only OTU columns are rarified and metadata information is preserved
data.nosing.rar<-cbind(subset(data.nosing, reads > 4999)[,1:5],rrarefy(subset(data.nosing,reads > 4999)[,-c(1)],5000))
#Preview file
data.nosing.rar[1:10,1:10]
str(data.nosing.rar[,1:10])
#write .csv  This is the file you will most likely use for subsequent analyses
#Fill in file name of your choice below, avoid using spaces
write.csv(data.nosing.rar, file="ChooseFileName_here.csv")
