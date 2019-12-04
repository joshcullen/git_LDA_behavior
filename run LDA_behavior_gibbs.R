# rm(list=ls(all=TRUE))
library('MCMCpack')
library('Rcpp')
set.seed(2)

source('LDA_behavior_function.R')
source('LDA_behavior_gibbs.R')
sourceCpp('aux1.cpp')
dat=read.csv('fake data.csv',as.is=T)

#prior
gamma1=0.1
alpha=0.1

#prepare for gibbs
ngibbs=1000
nburn=ngibbs/2
ind1=grep('y1',colnames(dat))
ind2=grep('y2',colnames(dat))
nmaxclust=10#max(length(ind1),length(ind2))-1
res=LDA_behavior_gibbs(dat=dat,gamma1=gamma1,alpha=alpha,
                       ngibbs=ngibbs,nmaxclust=nmaxclust,
                       nburn=nburn)
