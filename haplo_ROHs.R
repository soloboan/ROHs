haplo_ROHs <- function(islanddata,ROHbinary,animids,changemisshet,outname){
  require(psych)
  genodata <- read.table(paste(islanddata,sep=''),header=F)
  geno <- genodata[,1:6]
  for(m in seq(7,ncol(genodata),2)){
    dat <- data.frame(gen=paste(genodata[,m],genodata[,m+1],sep=''))
    geno <- cbind.data.frame(geno,dat)
    rm(dat)
  }
  rm(genodata)
  ROHstatus <- read.table(paste(ROHbinary,sep=''),header=F)
  anim <- read.table(paste(animids,sep=''),header=F)
  ROHstatus <- cbind.data.frame(anim,ROHstatus)
  ROHstatus <- merge(geno[,2],ROHstatus,by=1,sort=F)
  geno <- merge(ROHstatus[,1],geno,by.x=1,by.y=2,sort=F)
  
  nROHSNPs <- ncol(ROHstatus)-1
  animinROH <- which(rowSums(ROHstatus[,-1])>=(nROHSNPs))
  ROHstatus <- ROHstatus[animinROH,]
  geno <- geno[animinROH,]
  
  if(changemisshet==T){
    for (i in 1:nROHSNPs){
      recode012 <- function(x){
        if(x=='AA'){x=0
        } else if(x=='BB'){x=2
        } else if(x=='AB' | x=='BA' | x=='00'){x=1}
      }
      
      res <- data.frame(geno=geno[,6+i])
      res <- t(t(apply(res,1,recode012)))
      if(i==1){geno012 <- res} else {geno012 <- cbind.data.frame(geno012,res)}
      rm(res)
    }
      for (s in 1:nROHSNPs){
        res <- data.frame(geno=geno012[,s])
        frqhap <- data.frame(table(res))
        frqhap <- data.frame(frqhap[which(frqhap[,1]!=1),])
        frqhap <- as.numeric(as.vector(frqhap[frqhap[,2]==max(frqhap[,2]),1]))
        res <- data.frame(geno=ifelse(as.vector(res[,1])==1,frqhap,as.vector(res[,1])))
        if(s==1){replace_frq <- res
        } else {replace_frq <- cbind.data.frame(replace_frq,res)}
        rm(res,frqhap) 
    }
  } else if(changemisshet!=T){
    for (i in 1:nROHSNPs){
    recode012 <- function(x){
      if(x=='AA'){x=0
      } else if(x=='BB'){x=2
      } else if(x=='AB' | x=='BA'){x=1
      } else if(x=='00'){x=5}
    }
    res <- data.frame(geno=geno[,6+i])
    res <- t(t(apply(res,1,recode012)))
    if(i==1){geno012 <- res} else {geno012 <- cbind.data.frame(geno012,res)}
    rm(res)
    }
    for (s in 1:nROHSNPs){
      res <- data.frame(geno=geno012[,s])
      frqhap <- data.frame(table(res))
      frqhap <- data.frame(frqhap[which(frqhap[,1]!=1 | frqhap[,1]!=5),])
      frqhap <- as.numeric(as.vector(frqhap[frqhap[,2]==max(frqhap[,2]),1]))
      res <- data.frame(geno=ifelse(as.vector(res[,1])==5,frqhap,as.vector(res[,1])))
      if(s==1){replace_frq <- res
      } else {replace_frq <- cbind.data.frame(replace_frq,res)}
      rm(res)
    }      
  }
  genohaplo <- data.frame(haplo=apply(replace_frq,1,function(x) paste(x,collapse='')))
  freqhaplotye <- data.frame(table(genohaplo))
  freqhaplotye$percent <- (freqhaplotye$Freq/sum(freqhaplotye$Freq))*100
  freqhaplotye$hapNO <- as.numeric(freqhaplotye$genohaplo)
  colnames(freqhaplotye) <- c('haplotype','Num_animals','haplofrequency','haplotypeNumber')
  freqhaplotye$totanim <- nrow(ROHstatus)
  freqhaplotye$NSNPs <- nROHSNPs
  freqhaplotye <- freqhaplotye[order(freqhaplotye$haplofrequency,decreasing=T),]
  write.table(freqhaplotye,paste(outname,".csv",sep=''),quote=F,row.names=F,sep=',')
  
  plothap <- freqhaplotye[order(freqhaplotye$haplofrequency,decreasing=T),]
  plothap <- plothap[1:5,]
  dpi=300;width.cm<-15;height.cm<-13;width.in<-width.cm/2.54;height.in<-height.cm/2.54
  png(paste("haplotype_",outname,".png",sep=""),width=width.in*dpi,height=height.in*dpi,pointsize=8,units="px",res=dpi)
  barplot(as.vector(plothap$haplofrequency),names.arg=as.vector(plothap$haplotypeNumber),
          col=2:(nrow(plothap)+1),ylim=c(0,100),
          xlab="Haplotypes (represeted by Numbers)",ylab="Haplotype Frequency")
  dev.off()
}