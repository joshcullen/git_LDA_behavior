df.to.list=function(dat) {  #only for id as col in dat
  id<- unique(dat$id)
  n=length(id)
  dat.list<- vector("list", n)
  names(dat.list)<- id
  
  for (i in 1:length(id)) {
    dat.list[[i]]<- dat[dat$id==id[i],]
  }
  dat.list
}
#------------------------------------------------
get.summary.stats_behav=function(dat){  #dat must have time.seg assigned; for all IDs
  
  #create list of input and to store output
  dat.list<- df.to.list(dat = dat)
  id<- unique(dat$id)
  n<- length(id)
  obs.list<- vector("list", n)
  names(obs.list)<- id
  
  #calculate # of obs in each bin (per move param) by behav.seg
  for (i in 1:length(dat.list)) {
    dat.ind=dat.list[[i]]
    ntseg=max(dat.ind$behav.seg)
    
    
    #TA
    TA<- matrix(0, ntseg, max(dat.list[[i]]$TA, na.rm = T))
    colnames(TA)<- paste0("y1.",1:max(dat.list[[i]]$TA, na.rm = T))
    for (j in 1:ntseg){
      tmp<- dat.list[[i]] %>% filter(behav.seg == j) %>% dplyr::select(TA) %>% table()
      TA[j,as.numeric(names(tmp))]<- tmp
    }
    
    
    #SL
    SL<- matrix(0, ntseg, max(dat.ind$SL, na.rm = T))
    colnames(SL)<- paste0("y2.",1:max(dat.ind$SL, na.rm = T))
    for (j in 1:ntseg){
      tmp<- dat.ind %>% filter(behav.seg == j) %>% dplyr::select(SL) %>% table()
      SL[j,as.numeric(names(tmp))]<- tmp
    }
    
    
    id<- rep(unique(dat.ind$id), ntseg)
    tseg<- 1:ntseg
    behav.res<- cbind(id, tseg, TA, SL) %>% data.frame()
    obs.list[[i]]<- behav.res
  }
  #obs<- do.call(rbind.data.frame, obs.list)
  obs<- map_dfr(obs.list, `[`)
  obs
}
#------------------------------------------------
get_behav_hist=function(res) {  #generate DF of bin counts for histogram per behavior
  
  #summarize TA results by frequency and proportion
  behav.res.TA<- matrix(0, dim(res$z1.agg)[2]*dim(res$z1.agg)[3], 4)
  oo = 1
  for (i in 1:dim(res$z1.agg)[3]) {
    tmp<-  apply(res$z1.agg[,,i], 2, sum) %>% data.frame(count = ., prop = ./sum(.)) %>%
      mutate(bin = 1:dim(res$z1.agg)[2]) %>% mutate(behav = i) %>% as.matrix()
    
    behav.res.TA[oo:(oo+dim(res$z1.agg)[2] - 1),]<- tmp
    oo = oo+dim(res$z1.agg)[2]
  }
  
  behav.res.TA<- data.frame(behav.res.TA)
  names(behav.res.TA)<- c("count","prop","bin","behav")
  behav.res.TA$param<- "TA"
  
  
  #summarize SL results by frequency and proportion
  behav.res.SL<- matrix(0, dim(res$z2.agg)[2]*dim(res$z2.agg)[3], 4)
  oo = 1
  for (i in 1:dim(res$z2.agg)[3]) {
    tmp<-  apply(res$z2.agg[,,i], 2, sum) %>% data.frame(count = ., prop = ./sum(.)) %>%
      mutate(bin = 1:dim(res$z2.agg)[2]) %>% mutate(behav = i) %>% as.matrix()
    
    behav.res.SL[oo:(oo+dim(res$z2.agg)[2] - 1),]<- tmp
    oo = oo+dim(res$z2.agg)[2]
  }
  
  behav.res.SL<- data.frame(behav.res.SL)
  names(behav.res.SL)<- c("count","prop","bin","behav")
  behav.res.SL$param<- "SL"
  
  #combine params
  behav.res<- rbind(behav.res.TA, behav.res.SL)
  
  behav.res
}
#------------------------------------------------
obs.per.tseg=function(dat.list) {  #count the number of observations w/in time segments
  tmp<- dat.list %>% group_by(behav.seg) %>% tally()
  tmp$id<- unique(dat.list$id)
  
  tmp
}
#------------------------------------------------
aug_behav_df=function(dat, theta.estim, nobs) {  #augment from time segments to observations
  
  for (i in 1:nrow(theta.estim)) {
    ind<- which(dat$id == theta.estim$id[i] & dat$behav.seg == theta.estim$tseg[i])
    
    if (i == 1) {
      theta.estim2<- rep(theta.estim[i,], nobs$n[i]) %>%
        unlist() %>%
        matrix(nrow = nobs$n[i], ncol = ncol(theta.estim), byrow = TRUE)
    } else {
      tmp<- rep(theta.estim[i,], nobs$n[i]) %>%
        unlist() %>%
        matrix(nrow = nobs$n[i], ncol = ncol(theta.estim), byrow = TRUE)
      
      theta.estim2<- rbind(theta.estim2, tmp)
    }
  }
  
  colnames(theta.estim2)<- names(theta.estim)
  theta.estim2<- data.frame(theta.estim2, time1 = dat$time1, ESTtime = dat$ESTtime)
  
  theta.estim2
}
#------------------------------------------------
assign_behav=function(dat.list, theta.estim2) {  #assign dominant behavior to observations
  
  for (i in 1:length(dat.list)) {
    tmp<- matrix(0, nrow(dat.list[[i]]), 2)
    sub<- theta.estim2[theta.estim2$id == unique(dat.list[[i]]$id),]
    
    for (j in 1:nrow(sub)) {
      ind<- which.max(sub[j,3:9])
      tmp[j,1]<- names(ind) %>% as.character()
      tmp[j,2]<- round(sub[j,(2+ind)], 3) %>% as.numeric()
    }
    colnames(tmp)<- c("behav","prop")
    dat.list[[i]]<- cbind(dat.list[[i]], tmp)
    dat.list[[i]]$behav<- dat.list[[i]]$behav %>% as.character()
    dat.list[[i]]$prop<- dat.list[[i]]$prop %>% as.character()
  }
  
  #Convert to DF
  dat2<- map_dfr(dat.list, `[`)
  dat2$prop<- dat2$prop %>% as.numeric() 
  
  dat2
}