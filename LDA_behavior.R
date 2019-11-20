# rm(list=ls(all=TRUE))
library('MCMCpack')
set.seed(2)

setwd('U:\\GIT_models\\git_LDA_behavior')
source('LDA_behavior_function.R')
dat=read.csv('fake data.csv',as.is=T)
nobs=nrow(dat)

#separate variables
ind=grep('y1',colnames(dat))
y1=data.matrix(dat[,ind])
b1=length(ind)
ind=grep('y2',colnames(dat))
y2=data.matrix(dat[,ind])
b2=length(ind)

#prior
gamma1=0.1
alpha=0.1

#initial values
nbehav=3
phi1=matrix(1/b1,nbehav,b1)
phi2=matrix(1/b2,nbehav,b2)
theta=matrix(1/nbehav,nobs,nbehav)
z1.agg=array(NA,dim=c(nobs,b1,nbehav))
z2.agg=array(NA,dim=c(nobs,b2,nbehav))
for (i in 1:nobs){
  for (j in 1:b1){
    z1.agg[i,j,]=rmultinom(1,size=y1[i,j],prob=rep(1/nbehav,nbehav))
  }
  for (j in 1:b2){
    z2.agg[i,j,]=rmultinom(1,size=y2[i,j],prob=rep(1/nbehav,nbehav))
  }
}

#prepare for gibbs
ngibbs=1000
store.phi1=matrix(NA,ngibbs,nbehav*b1)
store.phi2=matrix(NA,ngibbs,nbehav*b2)
store.theta=matrix(NA,ngibbs,nobs*nbehav)
store.loglikel=rep(NA,1)
for (i in 1:ngibbs){
  print(i)
  
  #sample from FCD's
  z1.agg=sample.z1.agg(lphi1=log(phi1),ltheta=log(theta),y1=y1,
                       nobs=nobs,b1=b1)
  z2.agg=sample.z2.agg(lphi2=log(phi2),ltheta=log(theta),y2=y2,
                       nobs=nobs,b2=b2)
  # z1.agg=z1.agg.true
  # z2.agg=z2.agg.true
  
  v=sample.v(z1.agg=z1.agg,z2.agg=z2.agg,gamma1=gamma1,
             nobs=nobs,nbehav=nbehav)
  theta=get.theta(v=v,nbehav=nbehav,nobs=nobs)
  # theta=theta.true
  
  phi1=sample.phi1(z1.agg=z1.agg,alpha=alpha,nbehav=nbehav,b1=b1)
  # phi1=phi1.true
  phi2=sample.phi2(z2.agg=z2.agg,alpha=alpha,nbehav=nbehav,b2=b2)
  # phi2=phi2.true
  
  #calculate log-likelihood
  prob1=theta%*%phi1
  p1=sum(y1*log(prob1))
  prob2=theta%*%phi2
  p2=sum(y2*log(prob2))
  
  #store results
  store.phi1[i,]=phi1
  store.phi2[i,]=phi2
  store.theta[i,]=theta
  store.loglikel[i]=p1+p2
}
plot(store.loglikel,type='l')
