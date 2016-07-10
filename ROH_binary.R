roh_binary <- function(rohsummary,mapfile,animlist,chrNumber,outputname){
  ##### reading all the files
  cat('@******************************************************@\n')
  cat('@    making binary ROH files from PLINK (ROH) files    @\n')
  cat('@                                                      @\n')
  cat('@                                      by: S.A. Boison @\n')
  cat('@******************************************************@\n')
  cat('\n')
  cat('... reading PLINK .hom file ...\n')
  
  ROHs <- read.table(rohsummary,header=T,stringsAsFactors=F)
  ROHs <- ROHs[ROHs[,"CHR"]==chrNumber,]
  if(nrow(ROHs)==0)
    stop("... There are no ROHs on this Chromosome .......")
  cat('... finished reading .hom file ...\n')
  
  cat('... reading map file and the list of animals ...\n')
  map <- read.table(mapfile)
  map <- map[map[,1]==chrNumber,] 
  fam <- read.table(animlist,stringsAsFactors=F)
  cat('... finished reading map  file and list of animals ...\n')
  
  #### number of animal and SNPs
  nanim <- nrow(fam)
  nsnp <- nrow(map)
  
  ### creating a ROH binary file
  ROH_binary <-matrix(0,nrow=nanim,ncol=nsnp)
  colnames(ROH_binary) <- map[,2]
  rownames(ROH_binary) <- fam[,1]
  
  #### which animals have a ROH and how many are they
  animrohs <- unique(ROHs[,"IID"])
  Nanimrohs <- length(animrohs)
  
  ### printing info  
  iterchecks.anim <- round(Nanimrohs/5,digits=0)
  
  #### start filling the matrix with 1's represent 
  ## presence of a ROH for a region
  cat('... building the binary ROH file ...\n')
  
  for (i in 1:Nanimrohs){
    anim <- animrohs[i]
    rownameanim <- which(rownames(ROH_binary)==anim)
    roh_anim <- ROHs[ROHs[,"IID"] %in% anim,]
    
    for (m in 1:nrow(roh_anim)){
      snp1 <-which(colnames(ROH_binary)==roh_anim[m,"SNP1"])
      snp2 <-which(colnames(ROH_binary)==roh_anim[m,"SNP2"])
      ROH_binary[rownameanim,snp1:snp2] <- 1
    }
    if(i %% iterchecks.anim==0){
      cat(paste("...",i,"... out of ",Nanimrohs,
                     " (total animals ==",nanim,") with ROHs ...",sep=""),"\n")}
  }
  cat('... writing out binray ROH file ...\n')
  write.table(ROH_binary,paste(outputname,"_chr",chrNumber,"_binary.roh",sep=''),
              quote=F,row.names=F,col.names=F)
  write.table(colnames(ROH_binary),paste(outputname,"_chr",chrNumber,"_binary.snpslist",sep=''),
              quote=F,row.names=F,col.names=F)
  cat('... processed finished ...\n')
}
