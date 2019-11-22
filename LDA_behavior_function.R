sample.z1.agg=function(lphi1,ltheta,y1,nobs,b1,nbehav){
  z1.agg=array(0,dim=c(nobs,b1,nbehav))
  for(i in 1:nobs){
    for (j in 1:b1){
      if (y1[i,j]>0){
        lprob=lphi1[,j]+ltheta[i,]
        lprob=lprob-max(lprob)
        prob=exp(lprob)
        prob=prob/sum(prob)
        z1.agg[i,j,]=rmultinom2(prob,y1[i,j],runif(y1[i,j]),nmaxclust)#rmultinom(1,size=y1[i,j],prob=prob)
      }
    }
  }
  z1.agg
}
#-----------------------------------
sample.z2.agg=function(lphi2,ltheta,y2,nobs,b2,nbehav){
  z2.agg=array(0,dim=c(nobs,b2,nbehav))
  for(i in 1:nobs){
    for (j in 1:b2){
      if (y2[i,j]>0){
        lprob=lphi2[,j]+ltheta[i,]
        lprob=lprob-max(lprob)
        prob=exp(lprob)
        prob=prob/sum(prob)
        z2.agg[i,j,]=rmultinom2(prob,y2[i,j],runif(y2[i,j]),nmaxclust)#rmultinom(1,size=y2[i,j],prob=prob)
      }
    }
  }
  z2.agg
}
#-----------------------------------
sample.v=function(z1.agg,z2.agg,gamma1,nobs,nbehav){
  soma1=apply(z1.agg,c(1,3),sum)
  soma2=apply(z2.agg,c(1,3),sum)
  soma.fim=soma1+soma2
  
  cumsum1=t(apply(soma1[,nbehav:1],1,cumsum))
  cumsum1=cumsum1[,nbehav:1]
  cumsum1=cumsum1[,-1]
  cumsum2=t(apply(soma2[,nbehav:1],1,cumsum))
  cumsum2=cumsum2[,nbehav:1]
  cumsum2=cumsum2[,-1]
  cumsum.fim=cumsum1+cumsum2
  
  v=matrix(NA,nobs,nbehav-1)
  for (i in 1:(nbehav-1)){
    v[,i]=rbeta(nobs,soma.fim[,i]+1,cumsum.fim[,i]+gamma1)
  }
  cbind(v,1)
}
#-----------------------------------
get.theta=function(v,nbehav,nobs){
  theta=matrix(NA,nobs,nbehav)
  theta[,1]=v[,1]
  prod=rep(1,nobs)
  for (i in 2:nbehav){
    prod=prod*(1-v[,i-1])
    theta[,i]=v[,i]*prod
  }
  theta
}
#-----------------------------------
sample.phi1=function(z1.agg,alpha,nbehav,b1){
  soma=t(apply(z1.agg,2:3,sum))
  phi1=matrix(NA,nbehav,b1)
  for (i in 1:nbehav){
    phi1[i,]=rdirichlet(1,soma[i,]+alpha)
  }
  phi1
}
#-----------------------------------
sample.phi2=function(z2.agg,alpha,nbehav,b2){
  soma=t(apply(z2.agg,2:3,sum))
  phi2=matrix(NA,nbehav,b2)
  for (i in 1:nbehav){
    phi2[i,]=rdirichlet(1,soma[i,]+alpha)
  }
  phi2
}

