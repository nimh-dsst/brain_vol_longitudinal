library(lme4)
library(dplyr)
library(tidyr)
library(ggplot2)
library(brms)
library(purrr)
library(broom)
options(mc.cores=parallel::detectCores())
`%ni%` <- Negate('%in%')

d <- read.csv('Downloads/final.csv')

d %>%
  filter(subject!='') %>%
  select(subject, visit, ScanID, AgeScan, Left.Lateral.Ventricle:EstimatedTotalIntraCranialVol) %>%
  gather(structure, volume, Left.Lateral.Ventricle:EstimatedTotalIntraCranialVol) -> mod.dat

ggplot(mod.dat, aes(x=AgeScan, y=volume)) + 
  geom_point() +
  facet_wrap(~structure, scales = 'free_y') #outliers? #different scales? #constants?


mod.dat %>%
  group_by(subject, structure) %>%
  filter(AgeScan == min(AgeScan, na.rm=T)) %>%
  mutate(first_vol = volume) %>%
  select(subject, structure, first_vol) %>%
  left_join(mod.dat) %>%
  mutate(perc_change = ((volume/first_vol)*100)-100) -> mod.dat
  
ggplot(mod.dat, aes(x=AgeScan, y=perc_change)) + 
  geom_point() +
  facet_wrap(~structure, scales = 'free_y') #outliers? x5th ventricle? right inf lat vent? some no changesr?

ids <- unique(mod.dat$subject)
samp_ids <- sample(ids, 5) #look at a sample of 5 people
mod.dat %>%
  filter(subject %in% samp_ids) %>%
  ggplot(aes(x=AgeScan, y=perc_change, group=subject, color=subject)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~structure, scales='free_y')

m <- lmer(volume~AgeScan + (-1+AgeScan|subject) + (AgeScan|structure), data=mod.dat) #complains

mod.dat %>%
  ungroup() %>%
  group_by(structure) %>%
  mutate(vol_std = scale(volume)[,1]) -> mod.dat

m <- lmer(vol_std~AgeScan + (AgeScan|subject) + (AgeScan|structure), data=mod.dat) #also complains
m <- lmer(vol_std~AgeScan + (1|subject) + (AgeScan|structure), data=mod.dat) #also complains

#structs that show problems:
# Right.inf.lat.vent, left.inf.lat.vent, right.non.WM.hypointensities, right.wm.hypointensities, x5th.ventricle

mod.dat %>%
  filter(structure %ni% c('EstimatedTotalIntraCranialVol', 'Left.non.WM.hypointensities', 'Left.WM.hypointensities',
                          'Right.Inf.Lat.Vent', 'Left.Inf.Lat.Vent', 'Right.non.WM.hypointensities', 'Right.WM.hypointensities', 'X5th.Ventricle')) -> sub.mod.dat

m <- lmer(vol_std~AgeScan + (AgeScan|subject) + (AgeScan|structure), data=sub.mod.dat) #also complains
m <- lmer(vol_std~AgeScan + (1|subject) + (AgeScan|structure), data=sub.mod.dat) #also complains

sub.mod.dat %>%
  group_by(subject) %>%
  mutate(n_timepoints = n_distinct(AgeScan)) %>%
  mutate(AgeScan_minus_15 = AgeScan-15) %>%
  filter(n_timepoints>1) -> long_part

ggplot(long_part, aes(x=AgeScan_minus_15, y=vol_std)) + geom_point() + facet_wrap(~structure, scales='free_y')

m <- lmer(vol_std~AgeScan + (AgeScan|subject) + (AgeScan|structure), data=long_part) #also complains
m <- lmer(vol_std~AgeScan + (1|subject) + (AgeScan|structure), data=long_part) #also complains
m <- lmer(vol_std~AgeScan_minus_15 + (AgeScan_minus_15|subject) + (AgeScan_minus_15|structure), data=long_part) #also complains
m <- lmer(vol_std~AgeScan_minus_15 + (1|subject) + (AgeScan_minus_15|structure), data=long_part) #also complains


#these have problems
priors <- set_prior('normal(0,2.5)', class='b')

m1 <- brm(vol_std~AgeScan + (AgeScan|subject) + (AgeScan|structure), 
          data=mod.dat, prior=priors, chains=4)
m2 <- brm(vol_std~AgeScan + (1|subject) + (AgeScan|structure), 
          data=mod.dat, prior=priors, chains=4)
m3 <- brm(vol_std~AgeScan + (AgeScan|subject) + (AgeScan|structure), 
          data=long_part, prior=priors, chains=4)
m4 <- brm(vol_std~AgeScan + (1|subject) + (AgeScan|structure), 
          data=long_part, prior=priors, chains=4)
m5 <- brm(vol_std~AgeScan_minus_15 + (AgeScan_minus_15|subject) + (AgeScan|structure), 
          data=long_part, prior=priors, chains=4)
m6 <- brm(vol_std~AgeScan_minus_15 + (1|subject) + (AgeScan|structure), 
          data=long_part, prior=priors, chains=4)

model_vol <- function(df) {
  lmer(volume ~ AgeScan + (1|subject), data=df)
}
model_vol_std <- function(df){
  lmer(vol_std ~ AgeScan + (1|subject), data=df)
}
model_perc_chng <- function(df){
  df <- df[which(!is.nan(df$perc_change)),]
  df <- df[which(is.finite(df$perc_change)),]
  lmer(perc_change ~ AgeScan + (1|subject), data=df)
}

sub.mod.dat %>%
  group_by(structure) %>%
  nest() %>%
  mutate(vol_mod = map(data, model_vol),
         vol_std_model = map(data, model_vol_std),
         perc_ch_model = map(data, model_perc_chng)) %>%
  mutate(tidied_vol = map(vol_mod, tidy),
         tidied_vol_std = map(vol_std_model, tidy),
         tidied_perc_ch = map(perc_ch_model, tidy)) %>%
  unnest(tidied_vol, tidied_vol_std, tidied_perc_ch) -> fig.dat

fig.dat %>%
  filter(term=='AgeScan') %>%
  mutate(upper = estimate2 + std.error2,
         lower = estimate2 - std.error2) %>%
  ggplot(aes(x=reorder(structure, -estimate2), y=estimate2)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lower, ymax=upper)) + 
  coord_flip() + 
  geom_hline(yintercept=0) +
  ylab('Estimated % change per year') +
  xlab('')

fig.dat %>%
  filter(term=='AgeScan') %>%
  mutate(upper = estimate1 + std.error1,
         lower = estimate1 - std.error1) %>%
  ggplot(aes(x=reorder(structure, -estimate1), y=estimate1)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lower, ymax=upper)) + 
  coord_flip() + 
  geom_hline(yintercept=0) +
  ylab('Standard Deviation Change per year') +
  xlab('')
  
fig.dat %>%
  filter(term=='AgeScan') %>%
  ggplot(aes(x=estimate1, y=estimate2)) + geom_point()
