plot(res$loglikel[100:ngibbs],type='l')

compare=function(true1,estim1){
  rango=range(c(true1,estim1))
  plot(true1,estim1,ylim=rango,xlim=rango)
  lines(rango,rango,col='red')
}

theta.estim=matrix(res$theta[ngibbs,],nrow(dat),nmaxclust)
boxplot(theta.estim)

z1.tmp=apply(res$z1.agg,c(1,3),sum)[,1:6] 
z2.tmp=apply(res$z2.agg,c(1,3),sum)[,1:6]

#find best order
ordem=numeric()
for (i in 1:ncol(z1.true)){
  tmp=rep(NA,ncol(z1.true))
  for (j in 1:ncol(z1.tmp)){
    tmp[j]=cor(cbind(z1.tmp[,j],z1.true[,i]))[1,2]
  }
  ind=which(tmp==max(tmp))
  ordem=c(ordem,ind)
}

# head(z1.tmp[,ordem])
# head(z1.true)
# head(cbind(z2.tmp[,ordem],z2.true))

compare(z1.true,z1.tmp[,ordem])
compare(z2.true,z2.tmp[,ordem])

compare(theta.true,theta.estim[,ordem])

ind1=grep('y1',colnames(dat))
phi1.estim=matrix(res$phi1[ngibbs,],nmaxclust,length(ind1))
compare(phi1.true,phi1.estim[ordem,])

ind2=grep('y2',colnames(dat))
phi2.estim=matrix(res$phi2[ngibbs,],nmaxclust,length(ind2))
compare(phi2.true,phi2.estim[ordem,])

# par(mfrow=c(2,1))
# phi1.estim1=phi1.estim[ordem,]
# for (i in 1:5){
#   barplot(phi1.true[i,])
#   barplot(phi1.estim1[i,])
# }
# 
# par(mfrow=c(2,1))
# phi2.estim1=phi2.estim[ordem,]
# for (i in 1:5){
#   barplot(phi2.true[i,])
#   barplot(phi2.estim1[i,])
# }
