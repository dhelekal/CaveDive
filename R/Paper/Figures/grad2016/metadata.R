#' readResist
#' @description This function reads the metadata file and returns resistance status according to EUCAST guidelines
#' @export
readResist=function() {
  metadata=read.table(system.file("extdata", "grad2016.tab", package = "gonophylo"),sep='\t',comment.char='',header=T,as.is=T)
  abx=3:9  #PEN,TET, SPC,CFX,  CRO,  CIP, AZI
  cutoff1=c(1   ,1  ,64 ,0.125,0.125,0.06,0.5 )#EUCAST higher breakpoints - resistant only if >  cutoff
  cutoff2=c(0.06,0.5,64 ,0.125,0.125,0.03,0.25)#EUCAST lower  breakpoints - sensitive only if <= cutoff
  cutoff1=c(1   ,1  ,64 ,0.125,0.06 ,0.5 ,1   )#used   higher breakpoints - resistant only if >  cutoff
  abxnames=colnames(metadata)[abx]
  resist=matrix(1,nrow(metadata),length(abx))
  for (i in 1:length(abx)) {
    resist[which(is.na(metadata[,abx[i]])),i]=NA
    resist[which(metadata[,abx[i]]> cutoff1[i]),i]=2#Resistant
    resist[which(metadata[,abx[i]]<=cutoff2[i]),i]=0#Sensitive
  }
  colnames(resist)=abxnames
  rownames(resist)=metadata[,1]
  return(resist)
}

#' readDates
#' @description Read and return dates from metadata file
#' @export
readDates=function() {
  metadata=read.table(system.file("extdata", "grad2016.tab", package = "gonophylo"),sep='\t',comment.char='',header=T,as.is=T)
  dates=metadata[,'Year']
  names(dates)=metadata[,1]
  return(dates)
}

#' plotResist
#' @description Plot resistance data against dates
#' @export
plotResist=function(pdfout=NA) {
  resist=readResist()
  dates=readDates()
  r=range(dates)
  r=r[1]:r[2]
  abxnames=colnames(resist)
  if (!is.na(pdfout)) pdf(pdfout,14,7)
  par(mfrow=c(2,4))
  cols=c('blue','orange','red')
  for (i in 1:ncol(resist)) {  
    ta=table(c(dates,r,c(NA,NA,NA)),c(resist[,i],rep(NA,length(r)),c(0,1,2)))
    barplot(t(ta),col=cols,legend.text=c('Sensitive','Intermediate',sprintf('Resistant to %s',abxnames[i])),args.legend=list(x = "topleft"))
  }
  if (!is.na(pdfout)) dev.off()
}

#' readUsage
#' @description This function reads the usage file and returns as a convenient matrix
#' @param smooth Factor of smoothing (default is 1=no smoothing)
#' @export
readUsage=function(smooth=1) {
  usage=read.table(system.file("extdata", "usage-clean.csv", package = "gonophylo"),sep=',',comment.char='',header=T,as.is=T,row.names = 1)
  xout=seq(1988,2016,1/smooth)
  if (smooth>1) {
    oldusage=usage
    usage=matrix(NA,nrow(oldusage),length(xout))
    for (i in 1:nrow(oldusage)) usage[i,]=spline(1988:2016,oldusage[i,],xout=xout,method='natural')$y
    usage=as.data.frame(usage)
  }
  colnames(usage)<-xout
  rownames(usage)<-c('NON','CRO2','CRO1','CIP','CFX','PEN','AZI','GEN','OTH','OCE','SPC','OFL','TET')
  return(usage)
}
