#Elizabeth Bach, Ryan Williams
#NMDS with analysis of taxanomic centroids and environmental co-variates
#5 June 2015

#rm() removes all designated objects, this pervents confusion with any previous data you've been working with before running this code
#if you are begining a new R session, this is not necessary, but it won't hurt anything to run it
rm(list=ls())

#libarary() will open the packages necessary for the commands we will be using in this code
#if you do not have these packages intalled, you may intall through the menu or by running install.pacakges()
#You only need to install packages once, they need only be opened each time you start a new R session
library(reshape)
library(grid)
library(ggplot2)
library(lme4)
library(lmerTest)
library(bbmle)
library(vegan)
library(gridExtra)

#To produce the NMDS figures, you will need to run this function at the start of the R session
#code written by R. J. Williams
ggplot.NMDS<-function(XX,ZZ,COLORS){
	library(ggplot2)
MDS1<-data.frame(scores(XX))$NMDS1
MDS2<-data.frame(scores(XX))$NMDS2
Treatment<-ZZ

NMDS<-data.frame(MDS1,MDS2,Treatment)

NMDS.mean=aggregate(NMDS[,1:2],list(group=Treatment),mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }

  df_ell <- data.frame()
  for(g in levels(NMDS$Treatment)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                    veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                    ,group=g))
  }

X1<-ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = Treatment),size=3,alpha=0.75) +
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=2, linetype=5)+theme_bw()+theme(aspect.ratio=1)+scale_color_manual(values=COLORS)+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+theme(legend.title=element_text(size=15),legend.text=element_text(size=15))
X1    
}

#Read in the taxonomic dataset, file.choose() will open a window from which you can select your data file (shold be .csv created in
#Bach_ITS_taxonomy_merge.R)
data_taxa<-read.csv(file.choose())
head(data.nosing[,1:10])
str(data.nosing)
dim(data.nosing)

#Summarize by the taxonomic level you wish to investigate
#You can perform correlations at any taxonomic level, I recommend starting with phylum and working down
#It's likely order or family will be the lowest level with informative relationships worth displaing on NMDS
data_phyla<-data.frame(cast(data_taxa, SampleID~phylum, value="value", fun.aggregate=sum, add.missing=TRUE))
head(data_phyla)
data_order<-data.frame(cast(data_taxa, SampleID~order, value="value", fun.aggregate=sum, add.missing=TRUE))
head(data_order)
data_family<-data.frame(cast(data_taxa, SampleID~family, value="value", fun.aggregate=sum, add.missing=TRUE))
head(data_family)
#The following code will use data_phyla, but you can search and replace with data_order or whichever level you wish to work with

#You will need to read in the OTU table produced in Bach_ITS1_data_wrangling.R to perform the NMDS permutations, recall the taxa file is structured differently (samples and OTUs make up rows in the taxa file, samples are rows and OTUs are columns in the OTU table file)
#This code uses the rarified dataset, but code can always be copy-pasted/modified to look at non-rarified data
data.nosing.rar<-read.csv(file.choose())
head(data.nosing.rar[,1:10])
str(data.nosing.rar)
dim(data.nosing.rar)

#Presence/Absence
#Run the NMDS and select dimensions you wish to use (if you've done this in Bach_NMDS_exploration, you can delete all the rows except the dimensions you want to use)
mds.pa6<-metaMDS(decostand(data.nosing.rar[,-c(1:3)],"pa" ),distance="jaccard", k=6,autotransform=FALSE, na.rm=TRUE)
#now look at best stress score
mds.pa6
mds.pa5<-metaMDS(decostand(data.nosing.rar[,-c(1:3)],"pa" ),distance="jaccard", k=5,autotransform=FALSE, na.rm=TRUE)
mds.pa5
mds.pa4<-metaMDS(decostand(data.nosing.rar[,-c(1:3)],"pa" ),distance="jaccard", k=4,autotransform=FALSE, na.rm=TRUE)
mds.pa4
mds.pa3<-metaMDS(decostand(data.nosing.rar[,-c(1:3)],"pa" ),distance="jaccard", k=3,autotransform=FALSE, na.rm=TRUE)
mds.pa3
mds.pa2<-metaMDS(decostand(data.nosing.rar[,-c(1:3)],"pa" ),distance="jaccard", k=2,autotransform=FALSE, na.rm=TRUE)
mds.pa2
#Preview result you want to use
head(mds.pa2)

#Extract co-ordinates and coalate mds scors into dataframe, remember to use correct mds.pa file
MDS1<-data.frame(scores(mds.pa2))$NMDS1
MDS2<-data.frame(scores(mds.pa2))$NMDS2
SampleID<-data.nosing$SampleID
NMDS.pa<-data.frame(MDS1,MDS2,SampleID)
head(NMDS.pa)

#envfit() will correlate the NMDS co-ordinates with taxanomic data
#use the mds.pa object run above (change number if needed), use the data_phlya (or _order, _family), be sure to exclude metadata columns
TaxVectors1<-envfit(mds.pa2, data_phyla[,2:6], na.rm=TRUE)
#look at all correlations
TaxVectors1
#extract r, p-vals, and co-ordinates for vectors (must extact from envfit object)
vectors<-data.frame(TaxVectors1$vectors[1:4])
vectors
names<-rownames(TaxVectors))
TaxVectors2<-data.frame(names, vectors)
#subset the data frame so only statsitically significant vectors are retained (can change the value to 0.1 if desired)
TaxVectors3<-(subset(TaxVectors2, pvals<0.05))

#Now map taxa "centroids" on NMDS, these points represent the ecological space where each taxon is most common
#remember to ensure the correct objects and factors are being called, change the number in rainbow() to generate correct number of colors
#taxa points will be "grey", any color can be substituted, taxa points will be labeled with taxa name
ggplot.NMDS(mds.pa2, (data.nosing.rar$Factor1), rainbow(3))+geom_point(data=TaxVectors3, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE, hjust=0.5, vjust=0.5)+
geom_text(data=TaxVectors3,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=4)

#the position of taxa names can be adjusted for clarity using hjust and vjust, you can produce a vector with differing values for each taxon by concatonating a list: hjust.list<-c(0.25, 0.5, 0, 1.5) and use hjust=hjust.list in the geom_text code

#Now let's look at environmental variables that may correlate with community types
#Environmental data may be included as metadata in the data.nosing.rar OTU table, or may be read in as a seperate file
#Here we read in an additional file
#Read in environmental data, this will be a .csv you create independently, ex: "Project_taxaEnv.csv", make sure sample_IDs are identical
data.metadata<-read.csv(file.choose())
head(data.metadata[,1:30])
dim(data.metadata)

#now correlate the environmental data with the NMDS coordinates using envfit(), be sure to define the columns you want to run correlations with (here, we're using columns 7-25)
envectors1<-envfit(mds.pa2, data.metadata[,7:25], na.rm=TRUE)
#note, if using environmental data from the data.nosing.rar file, use this line of code instead (remove # to run) this example uses columns 3-7
#envectors1<-envfit(mds.pa2, data.nosing.rar[,3:7], na.rm=TRUE)
head(envectors1)
envectors1
#extract r, p-vals, and co-ordinates for vectors (must extact from envfit object)
names<-rownames(envectors))
envectors2<-data.frame(names, envectors1$vectors[1:4])
#you may subset the data frame to include only statistically significant correlations
envectors3<-subset(envectors2, envectors2$pvals<0.05)
envectors3

#Now map environmental vectors on NMDS, these arrows represent the direction of the positive correlation, all vectors should originate from the plot origin (0,0)
#remember to ensure the correct objects and factors are being called, change the number in rainbow() to generate correct number of colors
#lines and arrows will be "darkgrey", any color can be substituted, arrows will be labeled with name of the variable, here in "bold"
#agian, hjust and vjust can be used to adjust text position
ggplot.NMDS(mds.pa2, (data.nosing.rar$Factor1), rainbow(3))+geom_segment(data=envectors3, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",size=1,inherit_aes=FALSE)+
geom_text(data=envectors3,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),hjust=0,vjust=0.75,size=4,fontface="bold")

#We can produce and NMDS with both the taxonomic centroids and environmental vectors by layering
#assigning the plot as an object allows the plot to be easily incorporated into a multi-panel figure if needed
PA.NMDS<-ggplot.NMDS(mds.pa2, (data.nosing.rar$Factor1), rainbow(3))+geom_point(data=TaxVectors3, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="darkgrey",size=3,inherit_aes=FALSE)+
geom_text(data=TaxVectors3,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),vjust=-1,hjust=0.5,size=4,fontface="bold")+geom_segment(data=envectors3, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",size=1,inherit_aes=FALSE)+
geom_text(data=envectors3,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),hjust=0,vjust=0.75,size=4,fontface="bold")
PA.NMDS

#Figure aesthetics can be manipulated extensively in R, theme() will be most helpful for "cleaning" up the graph for publication purposes
#below is example code of how I manipulated a fiugre for presentation/publication, play around with the code to produce the effect you want
#I prefer to use color palates from iwanthue, which take into account color theory and color blindness etc.
color.1<-rgb(105,128,158, max=255)
color.2<-rgb(129,143,62, max=255)
color.3<-rgb(165,98,191, max=255)
colors.factor

PA.NMDS2<<-ggplot.NMDS(mds.pa2, (data.nosing.rar$Factor1), colors.factor)+geom_point(data=TaxVectors3, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="darkgrey",size=3,inherit_aes=FALSE)+
geom_text(data=TaxVectors3,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),vjust=-1,hjust=0.5,size=4,fontface="bold")+geom_segment(data=envectors3, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",size=1,inherit_aes=FALSE)+
geom_text(data=envectors3,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),hjust=0,vjust=0.75,size=4,fontface="bold")+
theme(axis.line=element_line(size=2), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=2, colour="black"), legend.position=c(0.87, 0.9), legend.background=element_blank(), legend.text=element_text(size=16, face="bold"),legend.key=element_blank(),legend.title=element_blank(),panel.background=element_blank(), axis.text=element_text(size=20, face="bold", colour="black"), axis.title=element_text(size=22, face="bold", colour="black"))
PA.NMDS2

#The same approach can be used for OTU abundance
#This code is the same as above, except we will relatives the abundance of each OTU in each sample using decostand() with "total", agian omit any metadata columns
#we will use Bray-Curtis distance for this NMDS, Bray-Curtis is the default for meatMDS() function, so we do not need to specify it
mds.ab6<-metaMDS(decostand(data.nosing.rar[,-c(1:3)],"total" ), k=6,autotransform=FALSE, na.rm=TRUE)
#now look at best stress score
mds.ab6
mds.ab5<-metaMDS(decostand(data.nosing.rar[,-c(1:3)],"total" ), k=5,autotransform=FALSE, na.rm=TRUE)
mds.ab5
mds.ab4<-metaMDS(decostand(data.nosing.rar[,-c(1:3)],"total" ), k=4,autotransform=FALSE, na.rm=TRUE)
mds.ab4
mds.ab3<-metaMDS(decostand(data.nosing.rar[,-c(1:3)],"total" ), k=3,autotransform=FALSE, na.rm=TRUE)
mds.ab3
mds.ab2<-metaMDS(decostand(data.nosing.rar[,-c(1:3)],"total" ), k=2,autotransform=FALSE, na.rm=TRUE)
mds.ab2
#Preview result you want to use
head(mds.ab2)

#Extract co-ordinates and coalate mds scors into dataframe
MDS1.ab<-data.frame(scores(mds.ab2))$NMDS1
MDS2.ab<-data.frame(scores(mds.ab2))$NMDS2
SampleID<-data.nosing.rar$SampleID
NMDS.ab<-data.frame(MDS1.ab,MDS2.ab,SampleID)
head(NMDS.ab)

#envfit() will correlate the NMDS co-ordinates with taxanomic data
#use the mds.pa object run above (change number if needed), use the data_phlya (or _order, _family), be sure to exclude metadata columns
TaxVectors1<-envfit(mds.ab2, data_phyla[,2:6], na.rm=TRUE)
#look at all correlations
TaxVectors1
#extract r, p-vals, and co-ordinates for vectors (must extact from envfit object)
vectors<-data.frame(TaxVectors1$vectors[1:4])
vectors
names<-rownames(TaxVectors))
TaxVectors2<-data.frame(names, vectors)
#subset the data frame so only statsitically significant vectors are retained (can change the value to 0.1 if desired)
TaxVectors3<-(subset(TaxVectors2, pvals<0.05))

#Now map taxa "centroids" on NMDS, these points represent the ecological space where each taxon is most common
#remember to ensure the correct objects and factors are being called, change the number in rainbow() to generate correct number of colors
#taxa points will be "grey", any color can be substituted, taxa points will be labeled with taxa name
ggplot.NMDS(mds.ab2, (data.nosing.rar$Factor1), rainbow(3))+geom_point(data=TaxVectors3, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="grey",inherit_aes=FALSE, hjust=0.5, vjust=0.5)+
geom_text(data=TaxVectors3,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=4)

#now correlate the environmental data with the NMDS coordinates using envfit(), be sure to define the columns you want to run correlations with (here, we're using columns 7-25)
envectors1<-envfit(mds.ab2, data.metadata[,7:25], na.rm=TRUE)
head(envectors1)
envectors1
#extract r, p-vals, and co-ordinates for vectors (must extact from envfit object)
names<-rownames(envectors))
envectors2<-data.frame(names, envectors1$vectors[1:4])
#you may subset the data frame to include only statistically significant correlations
envectors3<-subset(envectors2, envectors2$pvals<0.05)
envectors3

#Now map environmental vectors on NMDS, these arrows represent the direction of the positive correlation, all vectors should originate from the plot origin (0,0)
#remember to ensure the correct objects and factors are being called, change the number in rainbow() to generate correct number of colors
#lines and arrows will be "darkgrey", any color can be substituted, arrows will be labeled with name of the variable, here in "bold"
#agian, hjust and vjust can be used to adjust text position
ggplot.NMDS(mds.ab2, (data.nosing.rar$Factor1), rainbow(3))+geom_segment(data=envectors3, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",size=1,inherit_aes=FALSE)+
geom_text(data=envectors3,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),hjust=0,vjust=0.75,size=4,fontface="bold")

#We can produce and NMDS with both the taxonomic centroids and environmental vectors by layering
#assigning the plot as an object allows the plot to be easily incorporated into a multi-panel figure if needed
AB.NMDS<-ggplot.NMDS(mds.ab2, (data.nosing.rar$Factor1), rainbow(3))+geom_point(data=TaxVectors3, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="darkgrey",size=3,inherit_aes=FALSE)+
geom_text(data=TaxVectors3,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),vjust=-1,hjust=0.5,size=4,fontface="bold")+geom_segment(data=envectors3, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",size=1,inherit_aes=FALSE)+
geom_text(data=envectors3,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),hjust=0,vjust=0.75,size=4,fontface="bold")
AB.NMDS
