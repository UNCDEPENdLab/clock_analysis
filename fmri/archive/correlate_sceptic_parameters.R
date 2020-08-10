df <- read.csv("/Users/mnh5174/Data_Analysis/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_fixedparams_ffx_trial_statistics.csv.gz")

#str(df)
library(dplyr)
#library(dlookr)
library(broom)
library(Hmisc)

df <- df %>% mutate(pe_max_abs=abs(pe_max)) %>% group_by(id, run) %>% arrange(asc_trial) %>%
  mutate(pe_max_abs_lag = lag(pe_max_abs, 1),
         pe_max_lag = lag(pe_max, 1),
         rt_vmax_lag = lag(rt_vmax, 1),
         rt_vmax_change = abs(rt_vmax - rt_vmax_lag),
         v_entropy_lag = lag(v_entropy, 1),
         rt_lag = lag(rt_csv, 1),
         rt_swing = abs(rt_csv - rt_lag),
         v_entropy_change_abs = abs(v_entropy - v_entropy_lag),
         v_entropy_change = v_entropy - v_entropy_lag,
         v_entropy_change_pos = v_entropy_change*(v_entropy_change > 0),
         v_entropy_change_neg = abs(v_entropy_change*(v_entropy_change < 0))
         ) %>% ungroup()

df %>% arrange(id, asc_trial) %>% select(id, trial, rt_vmax, rt_vmax_lag, rt_vmax_change, v_entropy, v_entropy_change, v_entropy_change_pos, v_entropy_change_neg) %>% print(n=100)

table(df$v_entropy_change > 0)

xx <- df %>% group_by(id, run) %>% do({
  df <- .
  #cmat <- df %>% select(v_max, d_auc, v_entropy, mean_kld, intrinsic_discrepancy, pe_max, pe_max_abs) %>% cor(use="pairwise") %>% reshape2::melt() %>% filter(Var1 != Var2) 
  cmat <- df %>% select(v_max, d_auc, v_entropy, mean_kld, intrinsic_discrepancy, pe_max_lag, pe_max_abs_lag, rt_vmax_change, v_entropy_change, rt_swing) %>% cor(use="pairwise") %>% reshape2::melt() %>% filter(Var1 != Var2) 
  #cmat <- df %>% select(v_max, d_auc, v_entropy, mean_kld, intrinsic_discrepancy, pe_max, pe_max_lag, pe_max_abs, pe_max_abs_lag, kld_newlearn, kld_forget, rt_vmax_lag, v_entropy_change) %>% as.matrix() %>% rcorr(type="pearson") %>% tidy()# %>% filter(Var1 != Var2) 
}) %>% ungroup()

library(ggplot2)
#ggplot(xx, aes(x=estimate)) + geom_density() + facet_grid(column1 ~ column2) + geom_vline(xintercept=0)
ggplot(xx, aes(x=value)) + geom_density() + facet_grid(Var1 ~ Var2) + geom_vline(xintercept=0)

library(lme4)

#yes, KL decomposable into abs entropy change (stochastic, vertical shift) and change in rt vmax (horizontal shift)
mm <- lmer(mean_kld ~ v_entropy_change + rt_vmax_change + (1 | id), df)
vif.lme(mm)


mm1 <- lmer(sqrt(rt_swing) ~ v_entropy + rt_vmax_change + (1 | id), df)
summary(mm1)
vif.lme(mm1)

mm2 <- lmer(sqrt(rt_swing) ~ mean_kld + v_entropy + (1 | id/run), df)
summary(mm2)

mm3 <- lmer(sqrt(rt_swing) ~ pe_max_abs_lag + v_entropy + (1 | id/run), df)
summary(mm3)

mm4 <- lmer(v_entropy_change ~ pe_max_abs_lag + (1 | id/run), df)
summary(mm4)

library(MuMIn)
MuMIn::r.squaredGLMM(mm4)

mm5 <- lmer(v_entropy_change ~ 1 + (1 | id/run), df)
summary(mm5)


xx %>% group_by(column1, column2) %>% dplyr::summarize(mean(estimate)) %>% print(n=1000)
library(cowplot)



by_id <- split(df, df$id)
pdf("kld_trajectories.pdf", width=20, height=7)
lapply(by_id, function(sub) {
  to_plot <- sub #%>% select(id, run, trial, mean_kld, v_entropy)
  kl_all <- to_plot %>% gather(mean_kld, kld_newlearn, kld_forget, key="kl_measure", value="kl_value")
  to_plot <- to_plot %>% group_by(run) %>%
    mutate(
      pe_max_abs_lag = lag(pe_max_abs, 1),
      med_split_pe = factor(pe_max_abs_lag >  median(pe_max_abs_lag, na.rm=TRUE), levels=c(TRUE, FALSE), labels=c("Big PE", "Small PE"))) %>%
    ungroup()
  g1 <- ggplot(to_plot, aes(x=trial, y=mean_kld)) + geom_line() + geom_point(aes(color=med_split_pe)) + geom_vline(xintercept=c(0, 50, 100, 150, 200, 250, 300, 350, 400), color="blue") + ggtitle(paste0("Subject: ", to_plot$id[1])) #facet_wrap(~run) +
  #g1 <- ggplot(kl_all, aes(x=trial, y=kl_value, color=kl_measure)) + geom_line() + geom_vline(xintercept=c(0, 50, 100, 150, 200, 250, 300, 350, 400), color="blue") + ggtitle(paste0("Subject: ", to_plot$id[1])) + guides(color=FALSE) #facet_wrap(~run) +
  g2 <- ggplot(to_plot, aes(x=trial, y=v_entropy)) + geom_line() + geom_vline(xintercept=c(0, 50, 100, 150, 200, 250, 300, 350, 400), color="blue")
  g <- plot_grid(g1, g2, nrow=2, align="hv")  
  plot(g)
})
dev.off()


df <- df %>% group_by(id) %>% mutate(tot_score=sum(score_vba, na.rm=T)) %>% ungroup() %>% mutate(good_score = tot_score > median(tot_score, na.rm=T))

g1 <- ggplot(df, aes(x=rep(1:50, 600), y=mean_kld, color=good_score)) + stat_smooth() + facet_wrap(~rewFunc) # + geom_vline(xintercept=c(0, 50, 100, 150, 200, 250, 300, 350, 400), color="blue")
g2 <- ggplot(df, aes(x=rep(1:50, 600), y=v_entropy, color=good_score)) + stat_smooth() + facet_wrap(~rewFunc) # + geom_vline(xintercept=c(0, 50, 100, 150, 200, 250, 300, 350, 400), color="blue")
plot_grid(g1, g2, nrow=2, align='hv')


