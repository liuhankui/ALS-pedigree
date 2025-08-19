library(paramlink2)
library(ascii)

map<-read.table('merlin.map',head=T)
names(map)<-c('chrom','marker','cm')
ped<-read.table('merlin.ped')
aff<-as.numeric(ped$V6)
ped0<-ped[,1:5]
nmax<-ncol(ped)
ped1<-ped[,seq(7,nmax,by=2)]
ped2<-ped[,seq(8,nmax,by=2)]
ped3<-paste.matrix(ped1,ped2,sep="/")
ped<-cbind(ped0,ped3)
ped[ped=='0/0']<-'-/-'
colnames(ped)<-c('famid','id','fid','mid','sex',map$marker)
ped<-as.ped(ped)
modAD = diseaseModel("AD")
lods = lod(ped, aff = aff, model = modAD)##,rho=0.1)
write.table(lods,file='Lods.txt',row.names=F,col.names=T,quote=F,sep=' ')
