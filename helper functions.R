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
get_behav_hist=function(res) {
  
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