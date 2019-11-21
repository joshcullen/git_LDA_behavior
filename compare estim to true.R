compare=function(true1,estim1){
  rango=range(c(true1,estim1))
  plot(true1,estim1,ylim=rango,xlim=rango)
  lines(rango,rango,col='red')
}

boxplot(theta)

z1.tmp=apply(z1.agg,c(1,3),sum)[,1:4] 
z2.tmp=apply(z2.agg,c(1,3),sum)[,1:4]

#find best order
ordem=numeric()
for (i in 1:ncol(z1.true)){
  res=rep(NA,ncol(z1.true))
  for (j in 1:ncol(z1.tmp)){
    res[j]=cor(cbind(z1.tmp[,j],z1.true[,i]))[1,2]
  }
  ind=which(res==max(res))
  ordem=c(ordem,ind)
}

# head(cbind(z1.tmp[,ordem],z1.true))
# head(cbind(z2.tmp[,ordem],z2.true))

compare(z1.true,z1.tmp[,ordem])
compare(z2.true,z2.tmp[,ordem])

compare(theta.true,theta[,ordem])

compare(phi1.true,phi1[ordem,])
compare(phi2.true,phi2[ordem,])
