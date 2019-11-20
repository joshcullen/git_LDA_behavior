compare=function(true1,estim1){
  rango=range(c(true1,estim1))
  plot(true1,estim1,ylim=rango,xlim=rango)
  lines(rango,rango,col='red')
}

z1.tmp=apply(z1.agg,c(1,3),sum); head(cbind(z1.tmp,z1.true))
ordem=c(2,1,3)
head(cbind(z1.tmp[,ordem],z1.true))

z2.tmp=apply(z2.agg,c(1,3),sum)
head(cbind(z2.tmp[,ordem],z2.true))

compare(z1.true,z1.tmp[,ordem])
compare(z2.true,z2.tmp[,ordem])

compare(theta.true,theta[,ordem])

compare(phi1.true,phi1[ordem,])
compare(phi2.true,phi2[ordem,])
