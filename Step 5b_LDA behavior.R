set.seed(2)

library('MCMCpack')
library('Rcpp')
library(progress)
library(tidyverse)

source('LDA_behavior_function.R')
source('gibbs sampler.R')
source('helper functions.R')
sourceCpp('aux1.cpp')




############################
#### Load and Prep Data ####
############################

#get data
dat<- read.csv('Snail Kite Gridded Data_behav.csv', header = T, sep = ',')
dat.list<- df.to.list(dat)
obs<- get.summary.stats_behav(dat)


#prepare for Gibbs sampler
ngibbs=1000
nburn=ngibbs/2
ind1=grep('y1',colnames(obs))
ind2=grep('y2',colnames(obs))
nmaxclust=max(length(ind1),length(ind2))-1


#####################################################
#### Run Gibbs Sampler on All IDs Simultaneously ####
#####################################################

res=LDA_behavior_gibbs(dat=obs,gamma1=gamma1,alpha=alpha,
                       ngibbs=ngibbs,nmaxclust=nmaxclust,
                       nburn=nburn)


plot(res$loglikel,type='l')

theta.estim=matrix(res$theta[ngibbs,],nrow(obs),nmaxclust)
boxplot(theta.estim)


## Determine proportion of behavior assignment (across all time segments)
apply(theta.estim, 2, sum)/nrow(theta.estim)


behav.res<- get_behav_hist(res)

#Plot histograms of frequency data
ggplot(behav.res, aes(x = bin, y = count, fill = behav)) +
  geom_bar(stat = 'identity') +
  scale_fill_viridis_c(guide = F, direction = -1) +
  labs(x = "\nBin", y = "Frequency\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14), strip.text.x = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90)) +
  facet_grid(param ~ behav, scales = "free_y")



#Plot histograms of proportion data
ggplot(behav.res, aes(x = bin, y = prop, fill = behav)) +
  geom_bar(stat = 'identity') +
  scale_fill_viridis_c(guide = F, direction = -1) +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14), strip.text.x = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90)) +
  facet_grid(param ~ behav, scales = "fixed")
