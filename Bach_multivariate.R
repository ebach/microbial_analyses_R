#Elizabeth Bach, Ryan Williams
#12 June 2015
#multivariate stats to analyze ITS communities

#rm() removes all designated objects, this pervents confusion with any previous data you've been working with before running this code
#if you are begining a new R session, this is not necessary, but it won't hurt anything to run it
rm(list=ls())

#libarary() will open the packages necessary for the commands we will be using in this code
#if you do not have these packages intalled, you may intall through the menu or by running install.pacakges()
#You only need to install packages once, they need only be opened each time you start a new R session

library(reshape)
library(vegan)

#Read in the dataset, file.choose() will open a window from which you can select your data file (shold be .csv created in
#Bach_ITS_data_wrangling.R), this code uses the rarified data, but you can run the same analyses with non-rarified
data.nosing.rar<-read.csv(file.choose())

#We'll being with the presence/absence data (stats to match the Jaccard distance NMDS)
#transforming to presence/absence (0 and 1) for multivariate tests.  This uses decostand() from the vegan pacakge
head(data.nosing.rar[,1:8])
data.pa.rar<-cbind(data.nosing.rar[,1:8],decostand(data.nosing.rar[,-c(1:8)],"pa"))

#Use the adonis() function to perform PERMANOVA (non-parametric, multivariate analysis of variance) on the communities of each sample
#Again, exclude any columns of metadata (here, columns 1-8)
#First run the full model, look for any interactions with factors
adonis(data.pa.rar[,-c(1:8)]~data.pa.rar$Factor1*data.pa.rar$Factor2*data.pa.rar$Factor3, permutations=9999)

#Are all interactions significant? If not, you can remove non-significant interactions to a modified model, for example:
adonis(data.pa.rar[,-c(1:8)]~data.pa.rar$Factor1*data.pa.rar$Factor2+data.pa.rar$Factor3, permutations=9999)
#The + adds Factor3 as a main effect only, and the * preserves Factors 1 and 2 as main factors and their interaction
#This eliminates any interaction terms with Factor3 (Factor1*Factor3, Factor2*Factor3, and the 3-way interaction)
#Note, you can copy and past the PERMANOVA table from the console output into this code file with # to preserve it for future reference (without needing to re-run the code)

#If your factor(s) of interest have more than two levels, you can use the functions metaMDSdist() and betadisper() to contrast the dispersion (spatial median) of the samples within each group to determine which groups are different
#exclude metadata comluns, and adjust the dimensions (k) for the number of dimensions you are using in your NMDS
#metaMDSdist() is essentially performing the same NMDS function in the Bach_ITS_NMDS_Exploration code, but summariazing distance between points rather than generating co-ordinates
mds.dist<-metaMDSdist(decostand(data.pa.rar[,1:8], "pa"),k=2, index="jaccard", autotransform=FALSE)
Factor1.groups<-betadisper(mds.dist, data.pa.rar$Factor1, type="median")
#You can also contrast the group centroid by changing type="centroid", you can look at both simultaneously by using type=c("median","centroid"), note "median" is the default
#Now use TukeysHSD to contrast the median dispersio of the groups and find differences
TukeyHSD(Factor1.groups)



#Now looking at the scaled abundance data (stats to match the Bray-Curtis distance NMDS)
#first transform to scaled data within each sample (between 0 and 1) for multivariate tests.  This uses decostand() from the vegan pacakge, exclude any metadata columns (in this case, columns 1-8)
head(data.nosing.rar[,1:10])
data.ab.rar<-cbind(data.nosing.rar[,1:8],decostand(data.nosing.rar[,-c(1:8)],"total"))

#Use the adonis() function to perform PERMANOVA (non-parametric, multivariate analysis of variance) on the communities of each sample
#Again, exclude any columns of metadata (here, columns 1-8)
#First run the full model, look for any interactions with factors
adonis(data.ab.rar[,-c(1:8)]~data.ab.rar$Factor1*data.ab.rar$Factor2*data.ab.rar$Factor3, permutations=9999)

#Are all interactions significant? If not, you can remove non-significant interactions to a modified model, for example:
adonis(data.ab.rar[,-c(1:8)]~data.ab.rar$Factor1*data.ab.rar$Factor2+data.ab.rar$Factor3, permutations=9999)
#The + adds Factor3 as a main effect only, and the * preserves Factors 1 and 2 as main factors and their interaction
#This eliminates any interaction terms with Factor3 (Factor1*Factor3, Factor2*Factor3, and the 3-way interaction)

#Note, you can copy and past the PERMANOVA table from the console output into this code file with # to preserve it for future reference (without needing to re-run the code)

#If your factor(s) of interest have more than two levels, you can use the functions metaMDSdist() and betadisper() to contrast the dispersion (spatial median) of the samples within each group to determine which groups are different
#exclude metadata comluns, and adjust the dimensions (k) for the number of dimensions you are using in your NMDS
#metaMDSdist() is essentially performing the same NMDS function in the Bach_ITS_NMDS_Exploration code, but summariazing distance between points rather than generating co-ordinates
mds.dist<-metaMDSdist(decostand(data.ab.rar[,1:8], "total"),k=2, autotransform=FALSE)
Factor1.groups<-betadisper(mds.dist, data.ab.rar$Factor1, type="median")
#You can also contrast the group centroid by changing type="centroid", you can look at both simultaneously by using type=c("median","centroid"), note "median" is the default
#Now use TukeysHSD to contrast the median dispersio of the groups and find differences
TukeyHSD(Factor1.groups)


#Interactions between factors is more difficult to investigate
#Below is the code to do this for the presence/absence data, code can be copied or modified for the abundance data (change 'data.pa.rar' to 'data.ab.rar' when subsetting data)
#You will need to run an individual PERMANOVA for each interaction combination (of interest) to determine which couplings are different
#For this example, we're untangling an interaction between Factor1 (let's say has 3 levels) and Factor 2 (let's say has 2 levels)
#subset data for level1 of Factor1
data.Fact1level1<-droplevels(subset(data.pa.rar, Factor1=="level1"))
#Preview to ensure subset worked
str(data.Fact1level1[,1:10])

#run PERMANOVA to contrast the 2 levels of Factor 2 within level1 of Factor1
#I'm includeing Factor3 here to account for variation in the data due to Factor3, not because I'm interested in the significance of Factor3 for this contrast
adonis(data.Fact1level1[,1:8]~data.Fact1level1$Factor2+data.Fact1level1$Factor3, permutations=9999)
#Is Factor2 significant in this model?  If yes, then Factor2, level1 and level2 are different from each other at Factor1 level1

#Now subset data to isolate level2 of Factor1
data.Fact1level2<-droplevels(subset(data.pa.rar, Factor1=="level2"))
#Preview to ensure subset worked
str(data.Fact1level2[,1:10])

#run PERMANOVA to contrast the 2 levels of Factor 2 within level2 of Factor1
#I'm includeing Factor3 here to account for variation in the data due to Factor3, not because I'm interested in the significance of Factor3 for this contrast
adonis(data.Fact1level2[,1:8]~data.Fact1level2$Factor2+data.Fact1level2$Factor3, permutations=9999)
#Is Factor2 significant in this model?  If yes, then Factor2, level1 and level2 are different from each other at Factor1 level2

#Now subset data to isolate level3 of Factor1
data.Fact1level3<-droplevels(subset(data.pa.rar, Factor1=="level3"))
#Preview to ensure subset worked
str(data.Fact1level3[,1:10])

#run PERMANOVA to contrast the 2 levels of Factor 2 within level2 of Factor1
#I'm includeing Factor3 here to account for variation in the data due to Factor3, not because I'm interested in the significance of Factor3 for this contrast
adonis(data.Fact1level3[,1:8]~data.Fact1level3$Factor3+data.Fact1level3$Factor3, permutations=9999)
#Is Factor2 significant in this model?  If yes, then Factor2, level1 and level2 are different from each other at Factor1 level3

#Now isolate the 2 levels of Factor2 to see if Factor1 differs
#subset data for level1 of Factor2
data.Fact2level1<-droplevels(subset(data.pa.rar, Factor2=="level1"))
#Preview to ensure subset worked
str(data.Fact2level1[,1:10])

#Run PERMANOVA to see if Factor1 is significant at Factor2, level1
#exclude metadata (columns 1-8), include Factor3 to account for variation caused by that factor
adonis(data.Fact2level1[,1:8],~data.Fact2level1$Factor1+data.Fact2level1$Factor3, permutations=9999)

#If Factor1 is significant, we still can't tell which of the three levels are diffrent from one another, so we isolate each pairing and run a new PERMANOVA
#If Factor1 is not significant in the above model, there is no effect and there is no reason to do the contrasts for Factor1 levels
#We'll isolate levels 1 and 2 from Factor1 to see if they differ at Factor2 level1
data.Fact2level1a<-droplevels(subset(data.Fact2level1, (Factor1=="level1"|"level2")))
#PERMANOVA contrast of levels 1 and 2, retaining Factor3 to account for that variation only
adonis(data.Fact2level1a[,1:8],~data.Fact2level1a$Factor1+data.Fact2level1a$Factor3, permutations=9999)

#Now contrast levels 1 and 3 from Factor1
data.Fact2level1b<-droplevels(subset(data.Fact2level1, (Factor1=="level1"|"level3")))
adonis(data.Fact2level1b[,1:8],~data.Fact2level1b$Factor1+data.Fact2level1b$Factor3, permutations=9999)

#Now contrast levels 2 and 3 from Factor1
data.Fact2level1c<-droplevels(subset(data.Fact2level1, (Factor1=="level2"|"level3")))
adonis(data.Fact2level1c[,1:8],~data.Fact2level1c$Factor1+data.Fact2level1c$Factor3, permutations=9999)

#If the Factor1 contrast is significant in each of these models, then levels 1, 2, and 3 are all different from each other
#If only the 3rd iteration has a significant effect for Factor1, then levels 2 and 3 are different communities, but levels 1 and 2 are not

#Now repeat for Factor2 level 2
#subset data for level2 of Factor2
data.Fact2level2<-droplevels(subset(data.pa.rar, Factor2=="level2"))
#Preview to ensure subset worked
str(data.Fact2level2[,1:10])

#Run PERMANOVA to see if Factor1 is significant at Factor2, level2
#exclude metadata (columns 1-8), include Factor3 to account for variation caused by that factor
adonis(data.Fact2level2[,1:8],~data.Fact2level2$Factor1+data.Fact2level2$Factor3, permutations=9999)

#If Factor1 is significant, we still can't tell which of the three levels are diffrent from one another, so we isolate each pairing and run a new PERMANOVA
#If Factor1 is not significant in the above model, there is no effect and there is no reason to do the contrasts for Factor1 levels
#We'll isolate levels 1 and 2 from Factor1 to see if they differ at Factor2 level2
data.Fact2level2a<-droplevels(subset(data.Fact2level2, (Factor1=="level1"|"level2")))
#PERMANOVA contrast of levels 1 and 2, retaining Factor3 to account for that variation only
adonis(data.Fact2level2a[,1:8],~data.Fact2level2a$Factor1+data.Fact2level2a$Factor3, permutations=9999)

#Now contrast levels 1 and 3 from Factor1
data.Fact2level2b<-droplevels(subset(data.Fact2level2, (Factor1=="level1"|"level3")))
adonis(data.Fact2level2b[,1:8],~data.Fact2level2b$Factor1+data.Fact2level2b$Factor3, permutations=9999)

#Now contrast levels 2 and 3 from Factor1
data.Fact2level2c<-droplevels(subset(data.Fact2level2, (Factor1=="level2"|"level3")))
adonis(data.Fact2level2c[,1:8],~data.Fact2level2c$Factor1+data.Fact2level2c$Factor3, permutations=9999)

#If the Factor1 contrast is significant in each of these models, then levels 1, 2, and 3 are all different from each other
#If only the 3rd iteration has a significant effect for Factor1, then levels 2 and 3 are different communities, but levels 1 and 2 are not

#You can look at all of the above contrasts and see which ones were significant and which weren't to piece together where the interaction(s) occurr (e.g. in Factor1, level3, Factor2 was statistically different, but not at Factor1, levels 1 and 2)
#This will take some careful thought!

#Bonus coding trick:  You can combine the subsetting of factors into a single step:
data.Fact2level2a<-droplevels(subset(data.Fact2level2, Factor2=="level2" & (Factor1=="level1"|"level2")))