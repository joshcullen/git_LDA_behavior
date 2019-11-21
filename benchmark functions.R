library(microbenchmark)

microbenchmark(
  z1.agg=sample.z1.agg(lphi1=log(phi1),ltheta=log(theta),y1=y1,
                       nobs=nobs,b1=b1,nbehav=nmaxclust),
  z2.agg=sample.z2.agg(lphi2=log(phi2),ltheta=log(theta),y2=y2,
                       nobs=nobs,b2=b2,nbehav=nmaxclust),
  v=sample.v(z1.agg=z1.agg,z2.agg=z2.agg,gamma1=gamma1,
             nobs=nobs,nbehav=nmaxclust),
  theta=get.theta(v=v,nbehav=nmaxclust,nobs=nobs),
  phi1=sample.phi1(z1.agg=z1.agg,alpha=alpha,nbehav=nmaxclust,b1=b1),
  phi2=sample.phi2(z2.agg=z2.agg,alpha=alpha,nbehav=nmaxclust,b2=b2)
)
