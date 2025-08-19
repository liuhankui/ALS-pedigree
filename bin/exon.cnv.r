args<-commandArgs()
dfa<-read.table(gzfile('AA.regions.bed.gz'))
dfb<-read.table(gzfile('BB.regions.bed.gz'))
names(dfa)<-c('chromosome','start','end','gene','depth')
names(dfb)<-c('chromosome','start','end','gene','depth')
dfa$log2<-round(log2(dfa$depth),5)
dfb$log2<-round(log2(dfb$depth),5)
dfa$log2[is.infinite(dfa$log2)]<- -20
dfb$log2[is.infinite(dfb$log2)]<- -20
write.table(dfa,file=paste0(args[6],'.targetcoverage.cnn'),quote=F,sep='\t',col.names=T,row.names=F)
write.table(dfb,file=paste0(args[6],'.antitargetcoverage.cnn'),quote=F,sep='\t',col.names=T,row.names=F)
