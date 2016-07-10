ROH_plot <- function(rohbinaryfile,mapfile,chrNumber,anim,
                     colourRohs=c("YES","darkgrey"),title,outname){
  ####### importing and initial editing
  del <- read.table(rohbinaryfile,header=F,nrows=3)
  classes <- sapply(del,class); rm(del)
  rohbinary <- read.table(rohbinaryfile,header=F,colClasses=classes)
  cat('\n....... binary ROH file imported.......\n')
  
  del <- read.table(mapfile,header=F,nrows=10)
  classes <- sapply(del,class); rm(del)
  map <- read.table(mapfile,header=F,colClasses=classes)
  map <- map[which(map[,1]==chrNumber),]
  cat('....... MAP file imported.......\n')
  
  ## select animals you want to plot
  rohbinary <- rohbinary[anim,]
  #newanim <- which(rowMeans(rohbinary)>0)
  #rohbinary <- rohbinary[newanim,]
  ####
  nsnps <- 1:ncol(rohbinary)
  nanim <- 1:nrow(rohbinary)
  max.map <- max(map[,4]/1000000)
  
  #   if(max.map<50){
  #     ticks<-round(length(nsnps)/5); ticks<-seq(from=1,to=length(nsnps),ticks)
  #     axis.map <- round(map[ticks,4]/1000000,1); ticks[1] <- 0;axis.map[1]<-0    
  #   } else if(max.map>50 & max.map<=100){
  #     ticks<-round(length(nsnps)/10);ticks<-seq(from=1,to=length(nsnps),ticks) 
  #     axis.map <- round(map[ticks,4]/1000000,1);ticks[1] <- 0;axis.map[1]<-0    
  #   } else if(max.map>100){
  #     ticks<-round(length(nsnps)/12);ticks<-seq(from=1,to=length(nsnps),ticks) 
  #     axis.map <- round(map[ticks,4]/1000000,1);ticks[1] <- 0;axis.map[1]<-0    
  #   } 
  
  cat('....... Plotting of ROHs started .......\n')
  dpi=300;width.cm<-15;height.cm<-12;width.in<-width.cm/2.54;height.in<-height.cm/2.54
  tiff(file=paste(outname,'.tiff'),width=width.in*dpi,height=height.in*dpi,pointsize=8,
       units="px",res=dpi,compress="lzw")
  plot(x="",y="",xlab="Genomic Position (Mb)",ylab="Animals",cex=2,
       xlim=c(1,(round(max.map)+1)),xaxt="n",ylim=c(1,(length(nanim)*2)+2),yaxt="n",
       main=title)
  axis(side=1)
  
  #   plot(x="",y="",xlab="Genomic Position (Mb)",ylab="Animals",cex=2,
  #        xlim=c(1,max(nsnps)),xaxt="n",ylim=c(1,(length(nanim)*1)+2),yaxt="n",
  #        main=title)
  #   
  #   axis(side=1,at=ticks,labels=axis.map)
  
  if(colourRohs[1]=="YES"){
    for(i in nanim){
      datroh <- cbind.data.frame(i*2,map[,4]/1000000,t(rohbinary[i,]),1:nrow(map))
      datroh <- datroh[which(datroh[,3]==1),]
      if(nrow(datroh)>1){
        datroh$diff <- c(1,diff(datroh[,4]))
        segcut <- c(1,which(datroh$diff!=1),nrow(datroh))
        homanseg <- (length(segcut)-1)
        for (m in 1:homanseg){
          mincut=segcut[m];maxcut=(segcut[(m+1)])-1
          if(m==homanseg){maxcut=(segcut[(m+1)])}
          segments(datroh[mincut,2],datroh[mincut,1],datroh[maxcut,2],
                   datroh[maxcut,1],lwd=2,col=colourRohs[2],lty=1) 
        }          
      }
      #datroh <- cbind.data.frame(i*1,nsnps,t(rohbinary[i,]))
      #datroh <- datroh[which(datroh[,3]==1),]
      #points(x=datroh[,2],y=datroh[,1],pch="|",cex=0.20,col=colourRohs[2])
      
      if(i %% round(length(nanim)/10,0)==0){
        cat(" ...",i,"... out of ",length(nanim)," done !!! ...\n")}
    } 
  } else if (colourRohs[1]=="NO"){
    for(i in nanim){
      datroh <- cbind.data.frame(i*2,map[,4]/1000000,t(rohbinary[i,]),1:nrow(map))
      datroh <- datroh[which(datroh[,3]==0),]
      if(nrow(datroh)>1){
        datroh$diff <- c(1,diff(datroh[,4]))
        segcut <- c(1,which(datroh$diff!=1),nrow(datroh))
        homanseg <- (length(segcut)-1)
        for (m in 1:homanseg){
          mincut=segcut[m];maxcut=(segcut[(m+1)])-1
          if(m==homanseg){maxcut=(segcut[(m+1)])}
          segments(datroh[mincut,2],datroh[mincut,1],datroh[maxcut,2],
                   datroh[maxcut,1],lwd=2,col=colourRohs[2],lty=1) 
        }
      }
      #datroh <- cbind.data.frame(i*1,nsnps,t(rohbinary[i,]))
      #datroh <- datroh[which(datroh[,3]==0),]
      #points(x=datroh[,2],y=datroh[,1],pch="|",cex=0.20,col=colourRohs[2])
      if(i %% round(length(nanim)/10,0)==0){
        cat(" ...",i,"... out of ",length(nanim)," done !!! ...\n")}
    }
  }
  dev.off()
  cat('....... !!! Finished !!!! .......\n')    
}
