library(lme4)
library(dplyr)
library(ggplot2)

d <- read.csv('Downloads/final.csv')

d %>%
  mutate(hrs_since_4 = (Minutes_since_4am/60)/60,
         grayvol_std = scale(TotalGrayVol)[,1], #this variable is constant within subject.
         cortvol_std = scale(rhCortexVol)[,1]) %>%
  select(cortvol_std, grayvol_std, hrs_since_4, AgeScan, subject) -> mod.dat

m <- lmer(cortvol_std ~ hrs_since_4 + AgeScan + (1|subject), data=mod.dat) # fits, but there's something weird about using hrs_since_4 in this model...

m <- lm(grayvol_std ~ minutes_std, data=mod.dat)

d %>%
  select(rhCortexVol, subject, Minutes_since_4am, AgeScan) %>%
  filter(subject!='') %>%
  group_by(subject) %>%
  mutate(n_obs = n()) %>%
  ungroup() %>%
  filter(n_obs>5) %>%
  arrange(subject, AgeScan) #just to take a look at some of those who contributed a lot of data (relatively speaking)
  
ids <- as.character(unique(d$subject))
ids <- ids[ids!=''] #all unique ids


d %>%
  select(rhCortexVol, subject, Minutes_since_4am, AgeScan) %>%
  filter(subject!='') %>%
  group_by(subject) %>%
  filter(AgeScan == min(AgeScan, na.rm=T)) %>%
  mutate(first_vol = rhCortexVol) %>%
  select(subject, first_vol) %>%
  left_join(d) %>%
  select(rhCortexVol, subject, Minutes_since_4am, AgeScan, first_vol) %>%
  mutate(perc_change = ((rhCortexVol/first_vol)*100)-100,
         hrs_since_4 = (Minutes_since_4am/60)/60) -> mod.dat #adam wants models in terms of % change

samp_ids <- sample(ids, 25) #look at a sample of 25 people
mod.dat %>%
  filter(subject %in% samp_ids) %>%
  ggplot(aes(x=AgeScan, y=hrs_since_4)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~subject)

mod.dat %>%
  filter(subject %in% samp_ids) %>%
  ggplot(aes(x=AgeScan, y=perc_change)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~subject) + 
  geom_hline(yintercept = 0)

mod.dat %>%
  filter(subject %in% samp_ids) %>%
  ggplot(aes(x=hrs_since_4, y=perc_change)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~subject) + 
  geom_hline(yintercept = 0) #illustrates the issue with using first scan (in terms of year) as baseline to estimate % change: that first scan could have happened at any time of day!

m1 <- lmer(perc_change ~ hrs_since_4 + (1|subject), data=mod.dat) #bad idea
m2 <- lmer(perc_change ~ AgeScan + (1|subject), data=mod.dat) #fine
m3 <- lmer(perc_change ~ AgeScan + hrs_since_4 + (1|subject), data=mod.dat) #also probably a bad idea
m4 <- lmer(perc_change ~ hrs_since_4 + (hrs_since_4|subject), data=mod.dat) #not identifiable
m5 <- lmer(perc_change ~ AgeScan + (AgeScan|subject), data=mod.dat) #not identifiable
m5 <- lmer(perc_change ~ AgeScan + hrs_since_4 + (AgeScan + hrs_since_4|subject), data=mod.dat) #not identifiable
