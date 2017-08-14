## Generating ROH binary file from PLINK output

### Input files (5 argument used)
# 1. rohsummary ==== name of PLINK '.hom' file
# 2. mapfile    ==== mapfile used to generate the ROH files
# 3. animlist   ==== file containing all animals used in generating ROH files
# 4. chrNumber  ==== just a single chromosome number is allowed
# 5. outputname ==== outputname

#### Output files (2 files are produced)
# 1. "*_binary.roh"  ==== This contains the ROH binary file
# 2. ".snplist"      ==== SNP list


########################################################################################
### This is the whole script
source('ROH_binary.R')
## run example file
roh_binary(rohsummary='dataset.hom',mapfile='dataset.map',
           animlist='animals.txt',chrNumber=25,outputname='ResultDataset')

### The output files (ResultDataset_chr25_binary.roh and ResultDataset_25_binary.snpslist) 
########################################################################################



### checking if what was done is correct
### comparing the columns sums of the created file to the PLINK (.hom.summary)

### reading the output file from the script
rohbin<-read.table('ResultDataset_25_binary.roh',header=F)
res.roh <- data.frame(ROHs=colSums(rohbin))

### reading the PLINK output file (hom.summary)
dat<-read.table('dataset.hom.summary',header=T)
dat <- dat[which(dat$CHR==25),]

### if the diff is zero then it matchjes perfectly
checresut <- res.roh$ROHs-dat$UNAFF
unique(checresut)
