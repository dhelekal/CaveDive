library(ape)
library(BactDating)
rm(list=ls())
tree=read.tree('pathogenwatch.nwk')
metadata=read.table('metadata.csv',sep=',',header=T)

months=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
date=c()
for (i in 1:nrow(metadata)) {
  l=nchar(metadata[i,2])
  date[i]=as.numeric(substr(metadata[i,2],l-3,l))+0.5
  for (j in 1:length(months))
    if (length(grep(months[j],metadata[i,2]))>0) date[i]=date[i]-0.5+(j-0.5)/12
}
names(date)<-metadata[,1]

#rate=3e-6*2872769/3 # for USA300

#roottotip(tree,date)
res=bactdate(tree,date,showProgress = T,updateRoot = 'branch')
plot(res,'trace')
write.tree(res$tree,'dated.nwk')