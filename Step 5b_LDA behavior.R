set.seed(2)

library('MCMCpack')
library('Rcpp')
library(progress)
library(tidyverse)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)


source('LDA_behavior_function.R')
source('gibbs sampler.R')
source('helper functions.R')
sourceCpp('aux1.cpp')




############################
#### Load and Prep Data ####
############################

#get data
dat<- read.csv('Snail Kite Gridded Data_behav.csv', header = T, sep = ',')
dat$ESTtime<- dat$ESTtime %>% as_datetime()
dat.list<- df.to.list(dat)
obs<- get.summary.stats_behav(dat)


#prepare for Gibbs sampler
ngibbs=1000
nburn=ngibbs/2
ind1=grep('y1',colnames(obs))
ind2=grep('y2',colnames(obs))
nmaxclust=max(length(ind1),length(ind2))-1  #max possible is 1 fewer than largest number of bins


#####################################################
#### Run Gibbs Sampler on All IDs Simultaneously ####
#####################################################

res=LDA_behavior_gibbs(dat=obs,gamma1=gamma1,alpha=alpha,
                       ngibbs=ngibbs,nmaxclust=nmaxclust,
                       nburn=nburn)

#Check traceplot of log marginal likelihood
plot(res$loglikel,type='l')

#Extract and plots proportions of behaviors per time segment
theta.post<- res$theta[(nburn+1):ngibbs,] %>% apply(2, mean)  #calc mean of posterior
theta.estim=matrix(theta.post, nrow(obs), nmaxclust)
boxplot(theta.estim)


#Determine proportion of behaviors (across all time segments)
#Possibly set threshold below which behaviors are excluded
round(apply(theta.estim, 2, sum)/nrow(theta.estim), digits = 3)



#################################################################
#### Visualize Histograms of Movement Parameters by Behavior ####
#################################################################

behav.res<- get_behav_hist(res)

#Plot histograms of frequency data
ggplot(behav.res, aes(x = bin, y = count, fill = behav)) +
  geom_bar(stat = 'identity') +
  scale_fill_viridis_c(guide = F, direction = -1) +
  labs(x = "\nBin", y = "Frequency\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14), strip.text.x = element_text(face = "bold")) +
  facet_grid(param ~ behav, scales = "free_y")



#Plot histograms of proportion data
ggplot(behav.res, aes(x = bin, y = prop, fill = behav)) +
  geom_bar(stat = 'identity') +
  scale_fill_viridis_c(guide = F, direction = -1) +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14), strip.text.x = element_text(face = "bold")) +
  facet_grid(param ~ behav, scales = "fixed")



################################################
#### Visualize Behavior Estimates Over Time ####
################################################

#Assign behaviors (via theta) to each time segment
theta.estim<- data.frame(id = obs$id, tseg = obs$tseg, theta.estim)
names(theta.estim)<- c("id", "tseg", 1:nmaxclust)  #define behaviors
nobs<- purrr::map_dfr(dat.list, obs.per.tseg)  #calc obs per tseg

#Create augmented matrix by replicating rows (tsegs) according to obs per tseg
theta.estim2<- aug_behav_df(dat = dat, theta.estim = theta.estim, nobs = nobs)

#Change into long format
theta.estim.long<- theta.estim2 %>% gather(key, value, -id, -tseg, -time1, -ESTtime)
theta.estim.long$ESTtime<- theta.estim.long$ESTtime %>% as_datetime()
names(theta.estim.long)[5:6]<- c("behavior","prop")




### Aligned by first observation

#lines
ggplot(theta.estim.long) +
  geom_path(aes(x=time1, y=prop, color = behavior)) +
  labs(x = "\nObservation", y = "State Probability\n") +
  scale_color_viridis_d("Behavior", direction = -1) +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id, scales = "free_x")

#stacked area
ggplot(theta.estim.long) +
  geom_area(aes(x=time1, y=prop, fill = behavior), color = "black", size = 0.25,
            position = "fill") +
  labs(x = "\nObservation", y = "State Probability\n") +
  scale_fill_viridis_d("Behavior", direction = -1) +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id, scales = "free_x")

### Aligned by date


#Window of peak breeding (March 1 - June 30)
breed<- data.frame(xmin = as_datetime(c("2016-03-01 00:00:00","2017-03-01 00:00:00",
                                        "2018-03-01 00:00:00","2019-03-01 00:00:00")),
                   xmax = as_datetime(c("2016-06-30 23:59:59","2017-06-30 23:59:59",
                                        "2018-06-30 23:59:59","2019-06-30 23:59:59")),
                   ymin = -Inf, ymax = Inf)


#lines
ggplot() +
  geom_rect(data = breed, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.5) +
  geom_path(data = theta.estim.long, aes(x=ESTtime, y=prop, color = behavior)) +
  labs(x = "\nTime", y = "State Probability\n") +
  scale_color_viridis_d("Behavior", direction = -1) +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id)

#stacked area
ggplot(theta.estim.long) +
  geom_area(aes(x=ESTtime, y=prop, fill = behavior), color = "black", size = 0.25,
            position = "fill") +
  geom_rect(data = breed, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.25) +
  labs(x = "\nTime", y = "State Probability\n") +
  scale_fill_viridis_d("Behavior", direction = -1) +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id)




########################################
#### Map Dominant Behavioral States ####
########################################


#Add cluster assignments to original data; one column for dominant behavior and another for prop/prob to use for alpha of points
dat2<- assign_behav(dat.list = dat.list, theta.estim2 = theta.estim2)


#load map data
usa <- ne_states(country = "United States of America", returnclass = "sf")
fl<- usa %>% filter(name == "Florida")
fl<- st_transform(fl, crs = "+init=epsg:32617") #change projection to UTM 17N


# Facet plot of maps
ggplot() +
  geom_sf(data = fl) +
  coord_sf(xlim = c(min(dat$utmlong-120000), max(dat$utmlong+40000)),
           ylim = c(min(dat$utmlat-20000), max(dat$utmlat+20000)), expand = FALSE) +
  geom_path(data = dat2, aes(x=utmlong, y=utmlat), color="gray60", size=0.25) +
  geom_point(data = dat2, aes(utmlong, utmlat, fill=behav), size=2.5, pch=21, alpha=dat2$prop) +
  scale_fill_viridis_d("Behavior", direction = -1) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(~id)
