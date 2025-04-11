require("preprocessCore")
require("geneplotter")
library("grDevices")
library("Rsamtools")
library("IRanges")
library("rtracklayer")
library("mclust")
library("geneplotter")
library("colorRamps")
require("limma")
require("topGO")
require("GO.db")
require("org.Hs.eg.db")
require("gplots")

#Import ranges from introns
load("./GR_expressed_genes_VCP_study_April2017.RData")#"GR.pairs.time1","GR.pairs.time2","GR.pairs.strain","myGR.event","myCoordinates.focus","myCoordinates.full"; created in import_VARTOOLS_v2.R

#Import totdat
totdat                            <-  read.table("./INCLUSION_LEVELS_FULL-Hsa31-hg19.tab",header=T,fill=T,sep="\t")
myCoordinates.full$GENE           <- totdat$GENE[match(myCoordinates.full$id,totdat$EVENT)]
myCoordinates.IR                  <- myCoordinates.full[myCoordinates.full$type%in%c("IR-C","IR-S")&myCoordinates.full$GENE%in%sel_vcp$GENE,]
colnames(myCoordinates.IR)[c(3,4)]<-c("start","end")
myCoordinates.IR$start        <- as.numeric(as.character(myCoordinates.IR$start))
myCoordinates.IR$end          <- as.numeric(as.character(myCoordinates.IR$end))
myGR.IR                       <- makeGRangesFromDataFrame(myCoordinates.IR,keep.extra.columns=TRUE,ignore.strand=FALSE, seqinfo=NULL,seqnames.field="chr", start.field="start", end.field="end", strand.field="strand",starts.in.df.are.0based=FALSE)



load("./fe_clip_IR.RData")


pdf("enrichment_FUS_CLIP.pdf")
par(mfrow=c(1,1),mar=c(4,6,3,3))
seltop<- c(308:338)
temp_dat <- c(fe.v1.iclip[,"FUS"],fe.v1.HepG2[,"FUS"],fe.v1.K652[,"FUS"])
mycols=rep("grey",length(temp_dat))
mycols[seltop]<- "red"
barplot(sort(temp_dat),horiz=TRUE,las=1,col=rev(mycols))
mtext(side=1,line=2,"CLIP FE enrichment in FUS retained intron",cex=0.7)
barplot(sort(temp_dat)[seltop],horiz=TRUE,las=1,col="black")
dev.off()


#Compute enrichment
all_events   <- unlist(lapply(myEventsOI,function(Z)return(as.character(Z$EVENT))))
myGR.hg19    <- import.gff("./myGR.introns_hg19.gtf")
totdat.hg19  <-  read.table("./INCLUSION_LEVELS_FULL-Hsa31-hg19.tab",header=T,fill=T,sep="\t")
myGR.hg19$GENE <- as.character(totdat.hg19$GENE)[match(myGR.hg19$ID,totdat.hg19$EVENT)]
all_events     <- all_events[which(all_events%in%myGR.hg19$ID)]
GENES          <- as.character(totdat.hg19$GENE)[match(all_events,totdat.hg19$EVENT)]
myGR.IR        <- myGR.hg19[myGR.hg19$GENE%in%GENES]

#eCLIP -- HepG2
foi              <- list.files("./HepG2/")
foi              <- foi[grep(foi,pattern=".sum.gff")]
MOTS             <- unlist(lapply(foi,function(Z)return(unlist(strsplit(Z,split="_"))[1])))
REP              <- unlist(lapply(foi,function(Z)return(unlist(strsplit(Z,split="_"))[2])))

#FE 
temp <- lapply(foi,function(Z){
  sfpq.1           <- import.gff(paste("./HepG2/",Z,sep=""))
  mygenes          <- as.character(unique(myGR.IR$GENE))
  mymat            <- matrix(0,ncol=2,nrow=length(mygenes))
  for(IX in c(1:length(mygenes))){
    myGR.IR.oi       <- myGR.IR[myGR.IR$GENE==mygenes[IX],]
    gOver            <- data.frame(findOverlaps(query=myGR.IR.oi,subject=sfpq.1,ignore.strand=FALSE))
    wL               <- width(myGR.IR.oi)
    mysum            <- tapply(sfpq.1$score[gOver$subjectHits],INDEX=factor(gOver$queryHits,levels=c(1:length(myGR.IR.oi))),FUN=sum)
    mysum[is.na(mysum)]<-0
    oi          <- match(all_events[which(GENES==mygenes[IX])],as.character(myGR.IR.oi$ID))
    fg               <- sum(mysum[oi])/sum(wL[oi])
    bg               <- sum(mysum[-oi])/sum(wL[-oi])
    mymat[IX,1]      <- fg/bg
    mymat[IX,2]      <- (fg+1)/(bg+1)
    
  }
  colnames(mymat)<- c("FE.1","FE.2")
  rownames(mymat)<- mygenes
  print(Z)
  return(mymat)
})


fe.v1.HepG2            <- do.call(what=rbind,args=lapply(temp,function(Z)return(Z[,1])))
colnames(fe.v1.HepG2)  <- as.character(unique(myGR.IR$GENE))
rownames(fe.v1.HepG2)  <- paste(MOTS,REP,sep=".")

fe.v1.HepG2[is.nan(fe.v1.HepG2)] <- 0
for(i in c(1:ncol(fe.v1.HepG2))){
  fe.v1.HepG2[is.infinite(fe.v1.HepG2[,i]),i]<- max(fe.v1.HepG2[!is.infinite(fe.v1.HepG2[,i]),i])
}



#iCLIP
foi              <- list.files("./iCLIP/")
foi              <- foi[grep(foi,pattern=".sorted.gff")]
foi              <- foi[-grep(foi,pattern=".sorted.gff.idx")]
anno.files       <- lapply(foi,function(Z)return(unlist(strsplit(Z,split="[_-]"))))
anno.files       <- lapply(anno.files,function(Z)return(Z[Z!="all"]))
MOTS             <- unlist(lapply(anno.files,function(Z)return(Z[3])))



#FE 
temp <- lapply(foi,function(Z){
  sfpq.1           <- import.gff(paste("./iCLIP/",Z,sep=""))
  mygenes          <- as.character(unique(myGR.IR$GENE))
  mymat            <- matrix(0,ncol=2,nrow=length(mygenes))
  for(IX in c(1:length(mygenes))){
    myGR.IR.oi       <- myGR.IR[myGR.IR$GENE==mygenes[IX],]
    gOver            <- data.frame(findOverlaps(query=myGR.IR.oi,subject=sfpq.1,ignore.strand=FALSE))
    wL               <- width(myGR.IR.oi)
    mysum            <- tapply(sfpq.1$score[gOver$subjectHits],INDEX=factor(gOver$queryHits,levels=c(1:length(myGR.IR.oi))),FUN=sum)
    mysum[is.na(mysum)]<-0
    oi          <- match(all_events[which(GENES==mygenes[IX])],as.character(myGR.IR.oi$ID))
    fg               <- sum(mysum[oi])/sum(wL[oi])
    bg               <- sum(mysum[-oi])/sum(wL[-oi])
    mymat[IX,1]      <- fg/bg
    mymat[IX,2]      <- (fg+1)/(bg+1)
    
  }
  colnames(mymat)<- c("FE.1","FE.2")
  rownames(mymat)<- mygenes
  print(Z)
  return(mymat)
})


fe.v1.iclip            <- do.call(what=rbind,args=lapply(temp,function(Z)return(Z[,1])))
colnames(fe.v1.iclip)  <- as.character(unique(myGR.IR$GENE))
rownames(fe.v1.iclip)  <- MOTS



fe.v1.iclip[is.nan(fe.v1.iclip)] <- 0
for(i in c(1:ncol(fe.v1.iclip))){
  fe.v1.iclip[is.infinite(fe.v1.iclip[,i]),i]<- max(fe.v1.iclip[!is.infinite(fe.v1.iclip[,i]),i])
}


#Z
temp <- lapply(foi,function(Z){
  sfpq.1           <- import.gff(paste("./iCLIP/",Z,sep=""))
  mygenes          <- as.character(unique(myGR.IR$GENE))
  mymat            <- matrix(0,ncol=2,nrow=length(mygenes))
  for(IX in c(1:length(mygenes))){
    myGR.IR.oi       <- myGR.IR[myGR.IR$GENE==mygenes[IX],]
    gOver            <- data.frame(findOverlaps(query=myGR.IR.oi,subject=sfpq.1,ignore.strand=FALSE))
    wL               <- width(myGR.IR.oi)
    mysum            <- tapply(sfpq.1$score[gOver$subjectHits],INDEX=factor(gOver$queryHits,levels=c(1:length(myGR.IR.oi))),FUN=sum)
    mysum[is.na(mysum)]<-0
    mydensity   <- 1000*mysum/wL
    oi          <- match(all_events[which(GENES==mygenes[IX])],as.character(myGR.IR.oi$ID))
    mymat[IX,1] <-  mean(mydensity[oi])/mean(mydensity[-oi])
    mymat[IX,2] <- (mean(mydensity[oi])-mean(mydensity[-oi]))/sd(mydensity[-oi])
  }
  colnames(mymat)<- c("FE","Zscore")
  rownames(mymat)<- mygenes
  print(Z)
  return(mymat)
})

fe.iCLIP          <- do.call(what=rbind,args=lapply(temp,function(Z)return(Z[,1])))
colnames(fe.iCLIP)<- as.character(unique(myGR.IR$GENE))
rownames(fe.iCLIP)<- MOTS




#eCLIP -- K652
foi              <- list.files("./K652/")
foi              <- foi[grep(foi,pattern=".sum.gff")]
MOTS             <- unlist(lapply(foi,function(Z)return(unlist(strsplit(Z,split="_"))[1])))
REP              <- unlist(lapply(foi,function(Z)return(unlist(strsplit(Z,split="_"))[2])))

#FE 
temp <- lapply(foi,function(Z){
  sfpq.1           <- import.gff(paste("./K652/",Z,sep=""))
  mygenes          <- as.character(unique(myGR.IR$GENE))
  mymat            <- matrix(0,ncol=2,nrow=length(mygenes))
  for(IX in c(1:length(mygenes))){
    myGR.IR.oi       <- myGR.IR[myGR.IR$GENE==mygenes[IX],]
    gOver            <- data.frame(findOverlaps(query=myGR.IR.oi,subject=sfpq.1,ignore.strand=FALSE))
    wL               <- width(myGR.IR.oi)
    mysum            <- tapply(sfpq.1$score[gOver$subjectHits],INDEX=factor(gOver$queryHits,levels=c(1:length(myGR.IR.oi))),FUN=sum)
    mysum[is.na(mysum)]<-0
    oi          <- match(all_events[which(GENES==mygenes[IX])],as.character(myGR.IR.oi$ID))
    fg               <- sum(mysum[oi])/sum(wL[oi])
    bg               <- sum(mysum[-oi])/sum(wL[-oi])
    mymat[IX,1]      <- fg/bg
    mymat[IX,2]      <- (fg+1)/(bg+1)
    
  }
  colnames(mymat)<- c("FE.1","FE.2")
  rownames(mymat)<- mygenes
  print(Z)
  return(mymat)
})


fe.v1.K652            <- do.call(what=rbind,args=lapply(temp,function(Z)return(Z[,1])))
colnames(fe.v1.K652)  <- as.character(unique(myGR.IR$GENE))
rownames(fe.v1.K652)  <- paste(MOTS,REP,sep=".")

fe.v1.K652[is.nan(fe.v1.K652)] <- 0
for(i in c(1:ncol(fe.v1.K652))){
  fe.v1.K652[is.infinite(fe.v1.K652[,i]),i]<- max(fe.v1.K652[!is.infinite(fe.v1.K652[,i]),i])
}




#FE 
temp <- lapply(foi,function(Z){
  sfpq.1           <- import.gff(paste("./K652/",Z,sep=""))
  mygenes          <- as.character(unique(myGR.IR$GENE))
  mymat            <- matrix(0,ncol=2,nrow=length(mygenes))
  for(IX in c(1:length(mygenes))){
    
    myGR.IR.oi       <- myGR.IR[myGR.IR$GENE==mygenes[IX],]
    gOver            <- data.frame(findOverlaps(query=myGR.IR.oi,subject=sfpq.1,ignore.strand=FALSE))
    wL               <- width(myGR.IR.oi)
    mysum            <- tapply(sfpq.1$score[gOver$subjectHits],INDEX=factor(gOver$queryHits,levels=c(1:length(myGR.IR.oi))),FUN=sum)
    mysum[is.na(mysum)]<-0
    mydensity        <- mysum/wL
    oi          <- match(all_events[which(GENES==mygenes[IX])],as.character(myGR.IR.oi$ID))
    fg               <- sum(mysum[oi])/sum(wL[oi])
    bg               <- sum(mysum[-oi])/sum(wL[-oi])
    mymat[IX,1]      <- fg/bg
    mymat[IX,2]      <- (fg+1)/(bg+1)
  }
  colnames(mymat)<- c("FE","FEm")
  rownames(mymat)<- mygenes
  print(Z)
  return(mymat)
})

fe.v1.K652            <- do.call(what=rbind,args=lapply(temp,function(Z)return(Z[,1])))
colnames(fe.v1.K652)  <- as.character(unique(myGR.IR$GENE))
rownames(fe.v1.K652)  <- paste(MOTS,REP,sep=".")

fe.v1.K652[is.nan(fe.v1.K652)] <- 0
for(i in c(1:ncol(fe.v1.K652))){
  fe.v1.K652[is.infinite(fe.v1.K652[,i]),i]<- max(fe.v1.K652[!is.infinite(fe.v1.K652[,i]),i])
}



save(list=c("fe.v1.HepG2","fe.v1.K652","fe.v1.iclip"),file="./fe_clip_IR.RData")


