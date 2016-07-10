setwd('haplotype_ROHS/')


## STEP by STEP
## suppose after PLINK ROH analysis we got the file "roh.hom.indiv"; "roh.hom.summary"; "roh.hom"

####### We create a ROH binary file as follows
## list of animals
fam <- read.table('roh.hom.indiv',header=T)
write.table(fam[,2],'animals.txt',quote=F,row.names=F,col.names=F)
rm(fam)

### SNP/MAP information
map <- read.table('roh.hom.summary',header=T); map$gp <- 0
write.table(map[,c('CHR','SNP','gp','BP')],'HD_snpmap.txt',quote=F,row.names=F,col.names=F)
rm(map)

### want to analyze chromosome 24
chr=24
## run script to generate the ROH binary file (same as what we have been doing before)
## e.g. ROH file "roh.hom"
source('ROH_binary.R')
roh_binary(rohsummary='roh.hom',mapfile='HD_snpmap.txt',
           animlist='animals.txt',chrNumber=chr,outputname='NRD')


####### extracting the region that needs to be analyzed
chr=24
binarROH_filename=paste('NRD_chr',chr,'_binary.roh',sep='')
snpslist=paste('NRD_chr',chr,'_binary.snpslist',sep='')

#### importing map, binary and snp files
maporig <- read.table('HD_snpmap.txt')
data <- read.table(paste(binarROH_filename),colClasses='numeric')
snpmap <- read.table(paste(snpslist,sep=''))
maporigSNP <- maporig[which(maporig[,1]==chr),]
maporigSNP <- merge.data.frame(snpmap,maporigSNP,by.x=1,by.y=2,sort=F)[,-c(3,5:6)]  ## dont worry about warning message here
maporigSNP$nr <- 1:nrow(maporigSNP)

#### calculating proportion of animals in a ROH for each SNP
hotspot <- data.frame(mean=colMeans(data))
plot(hotspot$mean,pch=20) ### visual look

#### SNPs with ROH proportion lower than the (maximum ROH propotion - lowercut) are discarded 
lowercut=0.01
maxhotspot <- max(hotspot)-lowercut
hotspot <- cbind(maporigSNP,hotspot)
snpmaphot <- hotspot[hotspot$mean>maxhotspot,]
diff(snpmaphot$nr)  ### the values should all be 1's 
colsnphot <- which(hotspot$mean>maxhotspot)

### export the data 
## SNP information for only the regions of interest
write.table(snpmaphot[,1],'snplist.txt',quote=F,row.names=F,col.names=F)
write.table(snpmaphot,paste('NRD_chr',chr,'.haploSNPs',sep=''),quote=F,row.names=F,col.names=F)

### extract genotype data
system(paste('plink.exe --cow --bfile NRD --extract snplist.txt --recode --out NRD_chr',chr,sep=''))

### export the new binary data for only that region of interest
write.table(data[,colsnphot],paste('NRD_chr',chr,'.rohbin',sep=''),quote=F,row.names=F,col.names=F)


######################################################################
##### report haplotypes and frequencies
### Name of PEDFILE for only the short segment (the file is generated above)
### islanddata=paste('NRD_chr',chr,'.ped',sep='')

### Name of ROH binary file for only the short segment (the file is generated above)
### ROHbinary=paste('NRD_chr',chr,'.rohbin',sep='')

## changemisshet=T
##  1. This is to change heteroyzyotes or not
##  if 'T' that means chage heterozygote
##  if 'F' that means dont change the heterozgyote

### animids="animals.txt"
### list of animals 

####################################################################
source('haplo_ROHs.R')
chr=24
haplo_ROHs(islanddata=paste('NRD_chr',chr,'.ped',sep=''),
           ROHbinary=paste('NRD_chr',chr,'.rohbin',sep=''),
           changemisshet=T,
           animids="animals.txt",
           outname=paste("chr",chr,"_region1",sep=""))

haplo_ROHs(islanddata=paste('NRD_chr',chr,'.ped',sep=''),
           ROHbinary=paste('NRD_chr',chr,'.rohbin',sep=''),
           changemisshet=F,
           animids="animals.txt",
           outname=paste("chr",chr,"_region1_misshet",sep=""))
