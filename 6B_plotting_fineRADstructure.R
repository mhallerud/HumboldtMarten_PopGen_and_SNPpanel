##################################################################
## A simple R example for plotting fineRADstructure output
## Author: Milan Malinsky (millanek@gmail.com), adapted from a Finestructure R Example by Daniel Lawson (dan.lawson@bristol.ac.uk) and using his library of R functions
## Date: 04/04/2016
## Notes:
##    These functions are provided for help working with fineSTRUCTURE output files
## but are not a fully fledged R package for a reason: they are not robust
## and may be expected to work only in some specific cases - often they may require 
## at least minor modifications! USE WITH CAUTION!
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


### 1) EDIT THE FOLLOWING THREE LINES TO PROVIDE PATHS TO THE fineRADstructure OUTPUT 
setwd('/Users/maggiehallerud/Desktop/Marten_Fisher_Population_Genomics_Results/Marten/RADseq_Plate3_Results/fineradstructure/') ## The directory where the files are located
prefix <- "LowMissingnessRAD"
chunkfile <- paste0(prefix,".p.haps_chunks.out")
mcmcfile <- paste0(prefix, ".p.haps_chunks.mcmc.xml")
treefile <- paste0(prefix, ".p.haps_chunks.mcmcTree.xml")


### 2) EDIT THIS PATH TO WHERE YOU WANT THE PLOTS:
plotsFolder <- getwd()

### 3) SET VALUES FOR THESE VARIABLES: "analysisName" will be included in output plots
analysisName <- prefix

### 4) EDIT THE PATH TO YOUR COPY of FinestructureLibrary.R
#source("/Users/maggiehallerud/Downloads/FinestructureRcode-2/FinestructureLibrary.R", chdir = TRUE) # read in the R functions, which also calls the needed packages
source("/Users/maggiehallerud/Desktop/Mustelid_Genomics/R/FinestructureLibrary.R")
source("/Users/maggiehallerud/Desktop/Mustelid_Genomics/R/FinestructureDendrogram.R")


### 5) EXECUTE THE CODE ABOVE AND THE REST OF THE CODE BELOW
## make some colours
some.colors<-MakeColorYRP() # these are yellow-red-purple
some.colorsEnd<-MakeColorYRP(final=c(0.2,0.2,0.2)) # as above, but with a dark grey final for capped values
###### READ IN THE CHUNKCOUNT FILE
dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=1)) # read in the pairwise coincidence 
hist(datamatrix)
maxIndv <- max(datamatrix); maxPop<- max(datamatrix)


###### READ IN THE MCMC FILES
mcmcxml<-xmlTreeParse(mcmcfile, asText=FALSE) ## read into xml format
mcmcdata<-as.data.frame.myres(mcmcxml) ## convert this into a data frame
###### READ IN THE TREE FILES
treexml<-xmlTreeParse(treefile) ## read the tree as xml format
ttree<-extractTree(treexml) ## extract the tree into ape's phylo format

## Reduce the amount of significant digits printed in the posteror assignment probabilities (numbers shown in the tree):
ttree$node.label[ttree$node.label!=""] <-format(as.numeric(ttree$node.label[ttree$node.label!=""]),digits=2)
# convert to dendrogram format
tdend<-myapetodend(ttree,factor=1)
## Now we work on the MAP state
mapstate<-extractValue(treexml,"Pop") # map state as a finestructure clustering
mapstatelist<-popAsList(mapstate) # .. and as a list of individuals in populations
popnames<-lapply(mapstatelist,NameSummary) # population names IN A REVERSIBLE FORMAT (I.E LOSSLESS)
## NOTE: if your population labels don't correspond to the format we used (NAME<number>) YOU MAY HAVE TROUBLE HERE. YOU MAY NEED TO RENAME THEM INTO THIS FORM AND DEFINE YOUR POPULATION NAMES IN popnamesplot BELOW
popnamesplot<-lapply(mapstatelist,NameMoreSummary) # a nicer summary of the populations
mapstatelist <- popAsList(mapstate)
names(popnames)<-popnamesplot # for nicety only
names(popnamesplot)<-popnamesplot # for nicety only
popdend<-makemydend(tdend,mapstatelist) # use NameSummary to make popdend
popdend<-fixMidpointMembers(popdend) # needed for obscure dendrogram reasons
popdendclear<-makemydend(tdend,mapstatelist,"NameMoreSummary")# use NameMoreSummary to make popdend
popdendclear<-fixMidpointMembers(popdendclear) # needed for obscure dendrogram reasons


########################
## Plot 1: COANCESTRY MATRIX
# FOR PLOTTING!
fullorder<-labels(tdend) # the order according to the tree
which(!colnames(dataraw)%in%fullorder)
datamatrix<-dataraw[fullorder,fullorder] # reorder the data matrix

tmpmat<-datamatrix 

##FIX NAMES
dimnames(tmpmat)[[1]]
ALLMARTEN_NAMES <- c("NDunes_F07","NDunes_F02","NDunes_F06","NDunes_M02","NDunes_F31",
                     "SDunes_M01","SDunes_Mrdkl1","SDunes_Mrdkl3", "SDunes_M04","SDunes_Mrdkl2",
                     "NCali_M05","Trinidad_F01","Trinidad_M21","Trinidad_M12","SOreg_F05",
                     "SOreg_M15","NCali_Mrdkl3","NCali_Mrdkl2","NCali_F03","NCali_M02",
                     "SOreg_F07","SOreg_F08","SOreg_M22","SOreg_F04","SOreg_M09","SOreg_M10",
                     "SOreg_M13","SOreg_M16","SOreg_M07","SOreg_M08","SOreg_M17","SOreg_F10",
                     "SOreg_F11","SDunes_F04","ORcasc_M258","ORcasc_M252","ORcasc_F8943",
                     "Lssn_M65","Lssn_M66","Lssn_M61","Lssn_M64","Lssn_M67","NCali_M03",
                     "WY_F31","nCO_M21256","nCO_M21921","nCO_M21257","nCO_M21922","nCO_F21258",
                     "WAcasc_F215","BMT_03","BMT_M01","BMT_M06","BMT_M04","BMT_M11","BMT_M07","BMT_F08","BMT_F10")
cbind(dimnames(tmpmat)[[1]], ALLMARTEN_NAMES)#check
colnames(tmpmat) <- rownames(tmpmat) <- ALLMARTEN_NAMES

tmpmat[tmpmat>maxIndv] <- maxIndv #  # cap the heatmap
pdf(file=paste(plotsFolder,analysisName,"-SimpleCoancestry.pdf",sep=""),height=25,width=25)
plotFinestructure(tmpmat, colnames(tmpmat), dend=tdend,cols=some.colorsEnd,cex.axis=1.1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.1))
dev.off()



########################
## Plot 2: POPULATIONS AND COANCESTRY AVERAGES
#######################
## FIX MAPSTATELIST
mapstatelist
mapstatelist <- list(NDunes_F07="NDunes_F07", SDunes_F04="SDunes_F04", WY_F31="WY_F31",
                     WAcasc_F215="WAcasc_F215", NCali_M05="NCali_M05", NCali_M03="NCali_M03",
                     `SOreg_M22,SOreg_F08,SOreg_F07`=c("SOreg_M22","SOreg_F08","SOreg_F07"),
                     `Trinidad_F01,Trinidad_M21,Trinidad_M12`=c("Trinidad_F01","Trinidad_M21","Trinidad_M12"),
                     `BMT_M01,BMT_03`=c("BMT_M01","BMT_03"), `ORcasc_F8943,ORcasc_M252,ORcasc_M258`=c("ORcasc_F8943","ORcasc_M252","ORcasc_M258"),
                     BMT_F10="BMT_F10",BMT_M06="BMT_M06", 
                     `NDunes_F02,NDunes_F06,NDunes_M02,NDunes_F31`=c("NDunes_F02","NDunes_F06","NDunes_M02","NDunes_F31"),
                     `nCO_M21256,nCO_M21921`=c("nCO_M21256", "nCO_M21921"), `Lssn_M66,Lssn_M65`=c("Lssn_M66","Lssn_M65"),
                     `BMT_M07,BMT_M11,BMT_F08,BMT_M04`=c("BMT_M07","BMT_M11","BMT_F08","BMT_M04"),
                     `Lssn_M64,Lssn_M61,Lssn_M67`=c("Lssn_M64","Lssn_M61","Lssn_M67"), 
                     `SOreg_M10,SOreg_M09,SOreg_M13,SOreg_F04,SOreg_M16`=c("SOreg_M10","SOreg_M09","SOreg_M13","SOreg_F04","SOreg_M16"),
                     `nCO_M21922,nCO_M21257`=c("nCO_M21922","nCO_M21257"), 
                     `SOreg_M08,SOreg_M17,SOreg_F10,SOreg_M07,SOreg_F11`=c("SOreg_M08","SOreg_M17","SOreg_F10","SOreg_M07","SOreg_F11"),
                     `SOreg_M15,NCali_Mrdkl3,NCali_Rdkl2,SOreg_F05`=c("SOreg_M15","NCali_Mrdkl3","NCali_Mrdkl2","SOreg_F05"),
                     nCO_F21258="nCO_F21258", `NCali_M02,NCali_F03`=c("NCali_M02","NCali_F03"),
                     `SDunes_M01,SDunes_Mrdkl1,SDunes_Mrdkl3,SDunes_M04,SDunes_Mrdkl2`=c("SDunes_M01","SDunes_Mrdkl1","SDunes_Mrdkl3","SDunes_M04","SDunes_Mrdkl2")
)

rownames(datamatrix) <- colnames(datamatrix) <- ALLMARTEN_NAMES
popmeanmatrix <- getPopMeanMatrix(datamatrix, mapstatelist)
tmpmat <- popmeanmatrix
tmpmat[tmpmat>maxPop] <- maxPop # cap the heatmap
pdf(file=paste(plotsFolder,"/",analysisName,"-PopAveragedCoancestry.pdf",sep=""),height=20,width=20)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=1.1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.01))
dev.off()


########################
## Plot 3: POPULATIONS AND COANCESTRY AVERAGES WITH PERHAPS MORE INFORMATIVE LABELS
#fix duplicate label
labels <- labels(popdend)

mappopcorrectorder<-NameExpand(labels)
mappopsizes<-sapply(mappopcorrectorder,length)
labellocs<-PopCenters(mappopsizes)
xcrt=0
ycrt=45

pdf(file=paste(plotsFolder,analysisName,"-PopAveragedCoancestry2.pdf",sep=""),height=25,width=25)
plotFinestructure(tmpmat, dimnames(tmpmat)[[1]], labelsx=labels(popdendclear), labelsatx=labellocs,
                  cols=some.colorsEnd, dend=tdend, cex.axis=1.1, edgePar=list(lwd=1.5, p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0), hmmar=c(3,0,0,1))
dev.off()
  

###############################################
### PLOT 4: MIXED SIMPLE + POP-AVERAGED PLOT ##
###############################################
# upper=simple coancestry, lower+diagonal=pop means
tmpmat=datamatrix
lwrpopmeans = upper.tri(popmeanmatrix, diag=F)
tmpmat[lwrpopmeans] <- popmeanmatrix[lwrpopmeans]

#fix duplicate labels
labels <- labels(popdend)
mappopcorrectorder<-NameExpand(labels)
mappopsizes<-sapply(mappopcorrectorder,length)
labellocs<-PopCenters(mappopsizes)
xcrt=0
ycrt=45

#fix pop labels
labels(popdend)
labels <- c("NDunes_F07","NDunes","SDunes","NCali_M05","Trinidad","NCali/SOreg",
            "NCali","1SOreg","2SOreg","3SOreg","SDunes_F04","ORcasc","1Lassen","2Lassen","NCali_M03",
            "WY_F31","1nCO","2nCO","nCO_F21258","WAcasc_F215","1BMT","2BMT","3BMT","4BMT")


pdf(file=paste(plotsFolder,analysisName,"-SimpleAndPopCoancestry.pdf",sep=""),height=25,width=25)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],
                  labelsx=rep("",length(labellocs)),#labels,
                  labelsatx=labellocs,
                  cols=some.colorsEnd,
                  dend=tdend,
                  cex.axis=1.2,
                  edgePar=list(lwd=2, p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.01), hmmar=c(0,0,0,1))
dev.off()
