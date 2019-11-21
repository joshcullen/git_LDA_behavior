rm(list=ls(all=TRUE))
library('MCMCpack')
set.seed(3)

#basic settings
ntsegm.ind=20 #number of time segments per individual
nind=40
nbehavior=6
b1=6
b2=10
nobs.tsegm=100 #number of observations per time segment

#get parameters
phi1.true=phi1=rdirichlet(nbehavior,alpha=rep(0.1,b1))
phi2.true=phi2=rdirichlet(nbehavior,alpha=rep(0.1,b2))
theta.true=theta=rdirichlet(ntsegm.ind*nind,alpha=rep(0.1,nbehavior))

#look at these parameters
for (i in 1:nbehavior) plot(phi1.true[i,],type='h',main=i)
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
