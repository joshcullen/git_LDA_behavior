LDA_behavior_gibbs=function(dat,gamma1,alpha,ngibbs,nmaxclust,nburn){
  nobs=nrow(dat)
  
  #separate variables
  ind=grep('y1',colnames(dat))
  y1=data.matrix(dat[,ind])
  b1=length(ind)
  ind=grep('y2',colnames(dat))
  y2=data.matrix(dat[,ind])
  b2=length(ind)

  #initial values
  phi1=matrix(1/b1,nmaxclust,b1)
  phi2=matrix(1/b2,nmaxclust,b2)
  theta=matrix(1/nmaxclust,nobs,nmaxclust)
  z1.agg=array(NA,dim=c(nobs,b1,nmaxclust))
  z2.agg=array(NA,dim=c(nobs,b2,nmaxclust))
  for (i in 1:nobs){
    for (j in 1:b1){
      z1.agg[i,j,]=rmultinom(1,size=y1[i,j],prob=rep(1/nmaxclust,nmaxclust))
    }
    for (j in 1:b2){
      z2.agg[i,j,]=rmultinom(1,size=y2[i,j],prob=rep(1/nmaxclust,nmaxclust))
    }
  }

  #prepare for gibbs
  store.phi1=matrix(NA,ngibbs,nmaxclust*b1)
  store.phi2=matrix(NA,ngibbs,nmaxclust*b2)
  store.theta=matrix(NA,ngibbs,nobs*nmaxclust)
  store.loglikel=rep(NA,1)
  for (i in 1:ngibbs){
    print(i)
    
    #re-order clusters
    if (i < nburn & i%%50==0){
      med=apply(theta,2,mean)
      ordem=order(med,decreasing=T)
      theta=theta[,ordem]
      z1.agg=z1.agg[,,ordem]
      z2.agg=z2.agg[,,ordem]
      phi1=phi1[ordem,]
      phi2=phi2[ordem,]
    }
    
    #sample from FCD's
    z1.agg=sample.z1.agg(lphi1=log(phi1),ltheta=log(theta),y1=y1,
                         nobs=nobs,b1=b1,nbehav=nmaxclust)
    z2.agg=sample.z2.agg(lphi2=log(phi2),ltheta=log(theta),y2=y2,
                         nobs=nobs,b2=b2,nbehav=nmaxclust)
    # z1.agg=z1.agg.true
    # z2.agg=z2.agg.true
    
    v=sample.v(z1.agg=z1.agg,z2.agg=z2.agg,gamma1=gamma1,
               nobs=nobs,nbehav=nmaxclust)
    theta=get.theta(v=v,nbehav=nmaxclust,nobs=nobs)
    # theta=theta.true
    
    phi1=sample.phi1(z1.agg=z1.agg,alpha=alpha,nbehav=nmaxclust,b1=b1)
    # phi1=phi1.true
    phi2=sample.phi2(z2.agg=z2.agg,alpha=alpha,nbehav=nmaxclust,b2=b2)
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
  
  list(phi1=store.phi1,phi2=store.phi2,theta=store.theta,
       loglikel=store.loglikel,z1.agg=z1.agg,z2.agg=z2.agg)  
}



