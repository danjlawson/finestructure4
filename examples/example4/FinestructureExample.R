##################################################################
## Finestructure R Example
## Author: Daniel Lawson (dan.lawson@bristol.ac.uk)
## For more details see www.paintmychromosomes.com ("R Library" page)
## Date: 13/04/2016
## Notes:
##    These functions are provided for help working with fineSTRUCTURE output files
## but are not a fully fledged R package for a reason: they are not robust
## and may be expected to work only in some very specific cases! USE WITH CAUTION!
## SEE FinestrictureLibrary.R FOR DETAILS OF THE FUNCTIONS
##
## Licence: GPL V3
## 
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.

##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.

##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.

##########################################
## STEP 1: GET THE DATA:
## On Linux:
system("wget http://www.maths.bristol.ac.uk/~madjl/finestructure/FinestructureRexampledata.zip")
## On a mac:
system("curl -O http://www.maths.bristol.ac.uk/~madjl/finestructure/FinestructureRexampledata.zip")

system("unzip FinestructureRexampledata.zip")

##########################################
source("FinestructureLibrary.R") # read in the R functions, which also calls the needed packages

## make some colours
some.colors<-MakeColorYRP() # these are yellow-red-purple
some.colorsEnd<-MakeColorYRP(final=c(0.2,0.2,0.2)) # as above, but with a dark grey final for capped values

### Define our input files
chunkfile<-"EastAsiaSimple.EMlinked.chunkcounts.out" ## chromopainter chunkcounts file
mcmcfile<-"EastAsiaSimple.EMlinked.mcmc.xml" ## finestructure mcmc file
treefile<-"EastAsiaSimple.EMlinked.tree.xml" ## finestructure tree file

## Additional files that you can extract from finestructure
mappopchunkfile<-"EastAsiaSimple.EMlinked.mapstate.csv" # population-by-population chunkcount file for the populations used in the MAP (i.e tree)
system( paste("fs fs -X -Y -e X2",chunkfile,treefile,mappopchunkfile) )
meancoincidencefile<-"EastAsiaSimple.EMlinked.meancoincidence.csv" # pairwise coincidence, .i.e. proportion of MCMC files where individuals are found in the same 
system( paste("fs fs -X -Y -e meancoincidence",chunkfile,mcmcfile,meancoincidencefile) )
## there are ways of generating these within R but are either slower or more annoying - its your call how you do it

###### READ IN THE CHUNKCOUNT FILE
dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=1)) # read in the pairwise coincidence 

###### READ IN THE MCMC FILES
mcmcxml<-xmlTreeParse(mcmcfile) ## read into xml format
mcmcdata<-as.data.frame.myres(mcmcxml) ## convert this into a data frame

###### READ IN THE TREE FILES

treexml<-xmlTreeParse(treefile) ## read the tree as xml format
ttree<-extractTree(treexml) ## extract the tree into ape's phylo format
## If you dont want to plot internal node labels (i.e. MCMC posterior assignment probabilities)
## now is a good time to remove them via:
#     ttree$node.label<-NULL
## Will will instead remove "perfect" node labels
ttree$node.label[ttree$node.label=="1"] <-""
## And reduce the amount of significant digits printed:
ttree$node.label[ttree$node.label!=""] <-format(as.numeric(ttree$node.label[ttree$node.label!=""]),digits=2)

tdend<-myapetodend(ttree,factor=1) # convert to dendrogram format

####################################
## PLOT 1: RAW DENDROGRAM PLOT
pdf(file="EastAsiaSimpleFullDendrogram.pdf",height=6,width=14)
par(mar=c(6,0,2,0),mfrow=c(1,1))
fs.plot.dendrogram(tdend,horiz=FALSE,nodePar=list(cex=0,lab.cex=0.6),edgePar=list(p.lwd=0,t.srt=90,t.off=-0.5),axes=F)
dev.off()

## Now we work on the MAP state
mapstate<-extractValue(treexml,"Pop") # map state as a finestructure clustering
mapstatelist<-popAsList(mapstate) # .. and as a list of individuals in populations

popnames<-lapply(mapstatelist,NameSummary) # population names IN A REVERSIBLE FORMAT (I.E LOSSLESS)
## NOTE: if your population labels don't correspond to the format we used (NAME<number>) YOU MAY HAVE TROUBLE HERE. YOU MAY NEED TO RENAME THEM INTO THIS FORM AND DEFINE YOUR POPULATION NAMES IN popnamesplot BELOW
popnamesplot<-lapply(mapstatelist,NameMoreSummary) # a nicer summary of the populations
names(popnames)<-popnamesplot # for nicety only
names(popnamesplot)<-popnamesplot # for nicety only


popdend<-makemydend(tdend,mapstatelist) # use NameSummary to make popdend

popdendclear<-makemydend(tdend,mapstatelist,"NameMoreSummary")# use NameMoreSummary to make popdend

########################
## PLOT 2: population tree
pdf(file="EastAsiaSimplePopulationDendrogram.pdf",height=6,width=12)
par(mar=c(8,2,2,2),cex=0.8)
fs.plot.dendrogram(popdendclear,horiz=FALSE,nodePar=list(cex=0,lab.cex=1.2,las=2),edgePar=list(p.lwd=0,t.srt=0,t.off=0.3),yaxt="n",height=0.5,dLeaf=0.2)
dev.off()
	
########################
## PAIRWISE COINCIDENCES

fullorder<-labels(tdend) # the order according to the tree
mcmcmatrixraw<-as.matrix(read.csv(meancoincidencefile,row.names=1)) # read in the pairwise coincidence file we created earlier
mcmcmatrix<-mcmcmatrixraw[fullorder,fullorder] 
mapstatematrix<-groupingAsMatrix(mapstatelist)[fullorder,fullorder] # map state for reference

#########################
## PLOT 3: Pairwise coincidence, showing the MAP state

source("FinestructureLibrary.R")
pdf(file="EastAsiaSimplePairwiseCoincidence.pdf",height=12,width=12)
plotFinestructure(mcmcmatrix,dimnames(mcmcmatrix)[[1]],dend=tdend,optpts=mapstatematrix,cex.axis=0.6,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
dev.off()
	
########################
## COANCESTRY MATRIX

datamatrix<-dataraw[fullorder,fullorder] # reorder the data matrix

tmatmax<-500 # cap the heatmap
tmpmat<-datamatrix 
tmpmat[tmpmat>tmatmax]<-tmatmax # 
pdf(file="EastAsiaSimpleCoancestry.pdf",height=12,width=12)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=0.6,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
dev.off()

## Population averages
popmeanmatrix<-getPopMeanMatrix(datamatrix,mapstatelist)

tmatmax<-500 # cap the heatmap
tmpmat<-popmeanmatrix
tmpmat[tmpmat>tmatmax]<-tmatmax # 
pdf(file="EastAsiaPopAveragedCoancestry.pdf",height=12,width=12)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=0.6,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
dev.off()


### Useful tricks with labels

mappopcorrectorder<-NameExpand(labels(popdend))
mappopsizes<-sapply(mappopcorrectorder,length)
labellocs<-PopCenters(mappopsizes)
labelcols<-c(2,2,3,3,4,4,1,1,1,5,6,6,6,7,8,8,2) # different label colours allow clearer identification of individuals, too
labelcrt=45

pdf(file="EastAsiaSimpleCoancestry2.pdf",height=12,width=12)
plotFinestructure(tmpmat,labelsx=labels(popdendclear),labelsatx=labellocs,crt=labelcrt,dend=tdend,text.col=labelcols,cex.axis=1.0,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
dev.off()


####################################
## PCA Principal Components Analysis
pcares<-mypca(dataraw)
# For figuring out how many PCs are important; see Lawson & Falush 2012
# You need packages GPArotation and paran, and psych

## If you don't already have them:
## install.packages("psych")
## install.packages("paran")
## install.packages("GPArotation")
## install.packages("mclust")

tmap<-optimalMap(dataraw)
thorn<-optimalHorn(dataraw)
c(tmap,thorn) # 11 and 5. Horn typically underestimates, Map is usually better
pcapops<-getPopIndices(rownames(dataraw),mapstatelist)
pcanames<-rownames(dataraw)
rcols<-rainbow(max(pcapops))

pdf("EasAsiaSimpleExample_PCA.pdf",height=16,width=12)
par(mfrow=c(4,3))
for(i in 1:4) for(j in (i+1):5) {
  plot(pcares$vectors[,i],pcares$vectors[,j],col=rcols[pcapops],xlab=paste("PC",i),ylab=paste("PC",j),main=paste("PC",i,"vs",j),pch=rcols)
  text(pcares$vectors[,i],pcares$vectors[,j],labels=pcanames,col=rcols[pcapops],cex=0.5,pos=1)
}
dev.off()

###################################
## Clustering of the PCs
library("mclust")

## Using the MAP estimate of the number of PCs to retain:
pcastatemap<-pcaMclust(pcares$vectors[,1:tmap],40)$state
colnames(pcastatemap)<-paste0("Pop",1:dim(pcastatemap)[2])
rownames(pcastatemap)<-pcanames
pcapopmap<-matrixAsPopList(pcastatemap)
length(pcapopmap) # I get 24 populations, compared to the 17 from finestructure

# Using the Horn estimate of the number of PCs to retain:
pcastatehorn<-pcaMclust(pcares$vectors[,1:thorn],40)$state
colnames(pcastatehorn)<-paste0("Pop",1:dim(pcastatehorn)[2])
rownames(pcastatehorn)<-pcanames
pcapophorn<-matrixAsPopList(pcastatehorn)
length(pcapophorn) # I now get 16 populations

stateCor(pcapophorn,mapstatelist) # High; 0.81
stateCor(pcapopmap,mapstatelist) # Low; 0.36
stateCor(pcapopmap,pcapophorn) # Low; 0.33
# The map solution in this case appears to have oversplit; it created lots of singleton populations

#########################
## CHROMOPAINTER
## NOTE: Requires downloading "ChromoPainterExampleHGDPdata.zip" containing the additional Chromosome 1 example files

copyprobsfile<-"EastAsiaSimple.chrom1.linked.hap1.copyprobsperlocus.out.gz"
## file contains only haplotypes from individual 1 (it was run with "chromopainter -a 1 1 -b -in -iM -i 10 -g EastAsiaSimple.chrom1.phase -r EastAsiaSimple.chrom1.trecombfile -o EastAsiaSimple.chrom1.linked.hap1")
myhap<-getHap(1,copyprobsfile,verbose=TRUE) # read in first haplotpe (it takes a minute or two)
myhap2<-getHap(2,copyprobsfile,nlines=length(myhap$snps),verbose=TRUE) # second haplotype (takes half the time as provided with SNP count)
simplecollist<-MakeColorYRP(0.1) # construct a list of colours
cpdensityplot(myhap$snps[1:1000],myhap$probs[1:1000,],simplecollist) # plot the first 1000 SNPs
###########
# Now we will use the finestructure run to cluster the individuals

dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=1)) # read in the pairwise coincidence 
treexml<-xmlTreeParse(treefile) ## read the tree as xml format
mapstate<-extractValue(treexml,"Pop") # map state as a finestructure clustering
mapstatelist<-popAsList(mapstate) # .. and as a list of individuals in populations
## note: of course, you can just cluster by labels here if you prefer

names(mapstatelist)<-sapply(mapstatelist,NameMoreSummary) # choose how we will name populations
colnames(myhap$probs)<-dimnames(dataraw)[[1]] # name the probability matrix (not named otherwise, and needed for summing)
colnames(myhap2$probs)<-dimnames(dataraw)[[1]] # name the probability matrix (not named otherwise, and needed for summing)

popsnpmathap1<-matColSums(myhap$probs,mapstatelist) # construct a population level SNP matrix
popsnpmathap2<-matColSums(myhap2$probs,mapstatelist) # construct a population level SNP matrix
## CAREFUL with matColSums: if the names don't match, you'll miss individuals
collist2<-c(MakeColorYRP(0.2),rgb(0.3,0.3,0.3)) # contstruct a list of colours the correct length
collist2<-collist2[c(1,4,7,10,13,16,2,5,8,11,14,17,3,6,9,12,15)] # reorder for higher contrast
## FINALLY MAKE THE PLOT
png(file="EastAsiaExampleHaplotype.png",height=768,width=1024)
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
layout(matrix(1:2,nrow=2,ncol=1),height=c(3,1))
cpdensityplot(myhap2$snps[1:1000],popsnpmathap2[1:1000,],collist2) # plot the first 1000 SNPs
## Make a legend
par(mar=c(0,4,2,2)+0.1)
plot(c(0,1),c(0,1),type="n",axes=F,xlab="",ylab="")
legend("topleft",legend=names(mapstatelist)[1:6],col=collist2[1:6],lty=1,lwd=2,bty="n")
legend("top",legend=names(mapstatelist)[7:12],col=collist2[7:12],lty=1,lwd=2,bty="n")
legend("topright",legend=names(mapstatelist)[13:17],col=collist2[13:17],lty=1,lwd=2,bty="n")
dev.off()
