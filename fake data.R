rm(list=ls(all=TRUE))
library('MCMCpack')
set.seed(49)

#basic settings
ntsegm.ind=20 #number of time segments per individual
nind=40
#be careful with large nbehavior because it is easy to:
#- generate phi1 and phi2 in such a way that things are not clearly separable
#- generate theta in such a way that we don't have approx. pure time segments
nbehavior=6
b1=6
b2=10
nobs.tsegm=100 #number of observations per time segment

#get parameters
phi1=rdirichlet(nbehavior,alpha=rep(0.1,b1))
phi2=rdirichlet(nbehavior,alpha=rep(0.1,b2))
theta.true=theta=rdirichlet(ntsegm.ind*nind,alpha=rep(0.1,nbehavior))

#make sure that behaviors are not too similar
# for (i in 1:(nbehavior-1)){
#   for (j in (i+1):nbehavior){
#     #phi1
#     tmp=phi1[c(i,j),]
#     tmp1=cor(t(tmp))[1,2]
#     if (tmp1>0.9) phi1[j,]=rdirichlet(1,alpha=rep(0.5,b1))
#     
#     #phi2
#     tmp=phi2[c(i,j),]
#     tmp1=cor(t(tmp))[1,2]
#     if (tmp1>0.9) phi2[j,]=rdirichlet(1,alpha=rep(0.5,b2))
#   }
# }
phi1.true=phi1
phi2.true=phi2

#look at these parameters
par(mfrow=c(ceiling(nbehavior/2),2),mar=rep(1,4))
for (i in 1:nbehavior) plot(phi1.true[i,],type='h',main=i)
par(mfrow=c(ceiling(nbehavior/2),2),mar=rep(1,4))
for (i in 1:nbehavior) plot(phi2.true[i,],type='h',main=i)
boxplot(theta)
apply(theta>0.95,2,mean)

#generate data
oo=1
y1=matrix(0,ntsegm.ind*nind,b1)
y2=matrix(0,ntsegm.ind*nind,b2)
info=matrix(0,ntsegm.ind*nind,2)
z2.true=z1.true=matrix(0,ntsegm.ind*nind,nbehavior)
z1.agg.true=array(NA,dim=c(ntsegm.ind*nind,b1,nbehavior))
z2.agg.true=array(NA,dim=c(ntsegm.ind*nind,b2,nbehavior))

for (i in 1:nind){
  for (t in 1:ntsegm.ind){
    z1=rmultinom(1,size=nobs.tsegm,prob=theta[oo,])
    z1.true[oo,]=z1
    z2=rmultinom(1,size=nobs.tsegm,prob=theta[oo,])
    z2.true[oo,]=z1
    
    for (k in 1:nbehavior){
      z1.agg.true[oo,,k]=rmultinom(1,size=z1[k],prob=phi1[k,])
      z2.agg.true[oo,,k]=rmultinom(1,size=z2[k],prob=phi2[k,])
    }
    info[oo,]=c(i,t)
    oo=oo+1
  }
}
y1=apply(z1.agg.true,1:2,sum)
y2=apply(z2.agg.true,1:2,sum)

colnames(y1)=paste0('y1','_',1:b1)
colnames(y2)=paste0('y2','_',1:b2)
colnames(info)=c('ind.id','tsegm')
fim=cbind(info,y1,y2)
image(cbind(y1,y2))

#export fake data
setwd('U:\\GIT_models\\git_LDA_behavior')
write.csv(fim,'fake data.csv',row.names=F)
