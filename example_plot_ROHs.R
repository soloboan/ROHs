
### Input files (5 argument used)
# 1. rohbinaryfile  ==== name of binary roh file
# 2. mapfile    ==== mapfile used to generate the ROH files
# 3. chrNumber  ==== just a single chromosome number is allowed
# 4. anim   ==== the row numbers of the animals you want to plot
# 5. colourRohs ==== do you want to colour the ROHs or the non-ROHs
#                    if YES, the rohs are coloured, if NO then non-ROHs are coloured
# 6. title ==== title of the graph
# 7. outname ==== output name of the graph

########################################################################################
setwd('~/plot_ROHs/')
### This is the whole script
source("plot_ROHs.R")

### example to plot for chromosome 12
## see "anim" == am selecting the animals i want to plot (70 animals in total)
## see "colourRohs" == i dont want to colour ROHs so i chose NO (instead i want to colour the non-ROHs)
source("plot_ROHs.R")
ROH_plot(rohbinaryfile="chr12_binary.roh",
         mapfile='HD_imp.bim',
         chrNumber=12,
         anim=c(1:5,25:30,100:120,220:250,385:395,455:470),
         colourRohs=c("NO","darkgreen"),
         title="chromosome 12 ",
         outname="chr_12_nocolROH")


## see "colourRohs" == i want to colour ROHs so i chose YES (instead i want to colour the non-ROHs)
source("plot_ROHs.R")
ROH_plot(rohbinaryfile="chr12_binary.roh",
         mapfile='HD_imp.bim',
         chrNumber=12,
         anim=c(1:5,25:30,100:120,220:250,385:395,455:470),
         colourRohs=c("YES","black"),
         title="chromosome 12 ",
         outname="chr_12_colROH")


#### if i want to plot for multiple chromosomes
#### like chromosome 12,25,29
source("plot_ROHs.R")
for (i in c(25,29,12)){
  ROH_plot(rohbinaryfile=paste("chr",i,"_binary.roh",sep=""),
           mapfile='HD_imp.bim',
           chrNumber=i,
           anim=c(1:5,25:50,100:155,250:270,361:402),
           colourRohs=c("NO","darkgreen"),
           title=paste("chromosome",i),
           outname=paste("chr_",i,"_nocolROH",sep=""))
  
  cat(".... Chromosome ",i," .....\n")
}


source("plot_ROHs.R")
for (i in c(25,29,12)){
  ROH_plot(rohbinaryfile=paste("chr",i,"_binary.roh",sep=""),
           mapfile='HD_imp.bim',
           chrNumber=i,
           anim=c(1:5,25:50,100:155,250:270,361:402),
           colourRohs=c("YES","black"),
           title=paste("chromosome",i),
           outname=paste("chr_",i,"_colROH",sep=""))
  
  cat(".... Chromosome ",i," .....\n")
}


