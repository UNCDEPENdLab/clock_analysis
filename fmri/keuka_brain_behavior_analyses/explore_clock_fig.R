# make IEV/DEV plot for EXPLORE A1
library(tidyverse)
library(ggpubr)

# source('~/code/Rhelpers/')
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
# load('trial_df_and_vhdkfpe_clusters.Rdata')
load('explore_clock_rev_vba_out.rdata')
df <- vba_output
df <- df %>% group_by(id) %>% mutate(total_earnings = sum(score_csv))
median <- median(df$total_earnings)
df$performance <- 'NA'
df$performance[df$total_earnings>median] = 'Good performance'
df$performance[df$total_earnings<=median] = 'Poor performance'
# d <- df %>% filter((rewFunc=='IEV' | rewFunc=='DEV') & rt_csv < 4000 & emotion == "scram")
d <- df %>% filter(rt_csv < 5000 )

setwd('~/OneDrive/grants/explore_renewal_drafts/A1/figs/')
d$Contingency <- d$rewFunc
d$rts <- d$rt_csv/1000
df <- df %>%
  group_by(id, run) %>%  dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))),
                                       rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
                                       rt_cs = scale(rt_csv),
                                       rt_lag_cs = lag(rt_cs),
                                       rt_lag = lag(rt_csv) ,
                                       rt_swing_lag = lag(rt_swing),
                                       omission_lag = lag(score_csv==0),
                                       rt_vmax_lag = lag(rt_vmax),
                                       v_max_wi = scale(v_max),
                                       v_max_wi_lag = lag(v_max_wi),
                                       v_entropy_wi = scale(v_entropy),
                                       v_max_b = mean(na.omit(v_max)),
                                       v_entropy_b = mean(na.omit(v_entropy)),
                                       rt_change = rt_csv - rt_lag,
                                       pe_max_lag = lag(pe_max), 
                                       abs_pe_max_lag = abs(pe_max_lag), 
                                       rt_vmax_change = rt_vmax - rt_vmax_lag) %>% ungroup() %>% 
  dplyr::mutate(block=ceiling(trial/40),
                rev_trial = trial - (block-1)*40, rev_trial_ctr = rev_trial - 20) %>% ungroup()

s <- explore_demo
s$id <- s$registration_redcapid
s$group <- factor(s$registration_group, levels = c("ATT", "HC", "DEP", "IDE"))
s <- s %>% select(id,group)
dfs <- merge(s,df)
save(dfs,file = 'explore_prelim_behav_data.Rdata')


m1 <- lmer(rt_csv ~ rev_trial_ctr*rewFunc + (1|id), dfs)
m2 <- lmer(rt_csv ~ I(-10/(rev_trial_ctr+.5))*rewFunc*group + (1|id), dfs)
summary(m2)
car::Anova(m2)
anova(m1,m2)
# raw rts
ggplot(dfs, aes(rev_trial, rt_csv, color = group)) + geom_smooth(method = 'gam') + facet_wrap(~rewFunc)
library(emmeans)
plot(em2 <- emtrends(m2, var = "rev_trial_ctr", specs =  c("rewFunc", "group")))
plot(em2 <- emmeans(m2,  specs =  c("rewFunc", "group")))

# old code:

rt3<-lme4::lmer(rt_csv ~ rev_trial_ctr * rewFunc * group
                + (1|id), 
                data = dfs,verbose = T) 
summary(rt3)    
car::Anova(rt3)
# this is my favorite model before grant submission -- plot for preliminary data section
plot_model(rt3,show.p = TRUE, show.values = T,  vline.color = "slategray3", value.offset = 0.4,axis.lim = c(-.8,.8),
           axis.title = "Faster  < - >  Slower")
plot_model(rt3)

# plot only the high-order interaction for the grant
s <- summary(rt3)
coef <- s$coefficients
terms1 <- labels(coef)[[1]]
terms1 <- terms1[c(7,14:16)]
labels1 <- rev(c("Trial * Contingency",  "Trial * Contingency * Controls (n=26) vs. attempters (n=26)",
                 "Trial * Contingency * Depressed (n=14) vs. attempters",
                 "Trial * Contingency * Ideators (n=18) vs. attempters" ));

p <- plot_model(rt3,show.p = TRUE, terms = terms1,axis.labels = labels1, show.values = T,  vline.color = "slategray3", value.offset = 0.3,axis.lim = c(-.1,.5),
                axis.title = "Worse <-> Better learning")
p <- p + scale_y_continuous(limits = c(-.05,.25)) + ggtitle("Response time")
ggsave("clock_reversal_learning_attempters.pdf", width = 3.5,height = 2.5)





p <- ggplot(d, aes(rts,probability, color = Contingency)) + geom_smooth() + theme_bw() + 
  theme(legend.position = 'none') + xlab('Response time, s')
m <- ggplot(d, aes(rts,magnitude, color = Contingency)) + geom_smooth() + theme_bw() +
  theme(legend.position = c(0.4,.75)) + xlab('Response time, s')
v  <- ggplot(d, aes(rts,ev, color = Contingency)) + geom_smooth() + theme_bw() + 
  theme(legend.position = 'none') + xlab('Response time, s') + ylab("expected value")
h  <- ggplot(d, aes(ev,v_entropy, color = Contingency)) + geom_smooth(method = 'gam') + theme_bw() + 
  theme(legend.position = c(0.3,.2)) + xlab('expected value') + ylab("entropy")
hrt  <- ggplot(d, aes(rt_csv,v_entropy, color = Contingency)) + geom_smooth(method = 'gam') + theme_bw() + 
  theme(legend.position = c(0.2,.2)) + xlab('expected value') + ylab("Response time, s")

pv  <- ggplot(d, aes(ev,probability, color = Contingency)) + geom_smooth(method = 'gam') + theme_bw() + 
  theme(legend.position = c(0.2,.2)) + xlab('expected value') + ylab("probability")
mv  <- ggplot(d, aes(probability,magnitude, color = Contingency)) + geom_smooth(method = 'gam') + theme_bw() + 
  theme(legend.position = c(0.7,.25))


ggarrange(p,m,v, nrow = 1, ncol = 3)


