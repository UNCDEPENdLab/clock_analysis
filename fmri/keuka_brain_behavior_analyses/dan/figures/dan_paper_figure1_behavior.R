##DAN Figure 1 behavior plots
library(ggplot2)
library(dplyr)
library(wesanderson)
library(readr)

pal <- wes_palette(4, name = "IsleofDogs1", type = "discrete")[c(4,2,1,3)] %>%
  setNames(c("IEV", "DEV", "CEV", "CEVR"))

##true underlying IEV contingency.
setwd("~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/figures")
source("~/Data_Analysis/clock_analysis/td_model/getrew.R")
fm <- getMagFreq(0:4000, "IEV")
f <- fm$Freq
m <- fm$Mag
ev <- f*m

library(ggplot2)

df <- c()
#take CEVR out since CEV and CEVR are identical wrt EV (and we are not showing prob + freq)
for (cont in c("IEV", "DEV", "CEV", "CEVR")) { #, "CEVR"
  fm <- getMagFreq(0:4000, cont)
  #if (cont=="CEV") { cont="CEV/CEVR" } #for plot name
  df <- rbind(df, data.frame(contingency=cont, time=0:4000, mag=fm$Mag, freq=fm$Freq, ev=fm$Mag*fm$Freq))
}

#Figure: plot of EV in clock task
pdf("Clock contingencies.pdf", width=5, height=3.4)
g <- ggplot(df, aes(x=time/1000, y=ev, color=contingency)) + geom_line(linewidth=2) + 
  ylab("Expected value (points)") + xlab("Time (seconds)") + 
  scale_color_manual("Contingency", values=pal) +
    theme_bw(base_size=18) + theme(axis.title.x=element_text(margin = margin(t = 10)), 
                                   axis.title.y=element_text(margin = margin(r = 10)), legend.spacing = unit(0.15, "cm"), 
        plot.margin=margin(r=3, l=3, t=10, b=10))
plot(g)
dev.off()

#freq, prob, and ev
library(cowplot)
gcommon <- list(geom_line(linewidth=2), xlab("Time (seconds)"), scale_color_manual("Contingency", values=pal), theme_bw(base_size=18), 
  theme(axis.title.x=element_text(margin = margin(t = 10)), axis.title.y=element_text(margin = margin(r = 8)), 
        legend.spacing = unit(0.15, "cm"), 
        plot.margin=margin(r=10, l=10, t=10, b=5)))


g1 <- ggplot(df, aes(x=time/1000, y=mag, color=contingency)) + gcommon + ylab("Reward magnitude (points)") + 
  theme(legend.position="none", plot.margin=margin(r=10, l=0, t=10, b=5))

g2 <- ggplot(df, aes(x=time/1000, y=freq, color=contingency)) + gcommon + ylab("Reward probability") + theme(legend.position="none")

df$ev[df$contingency=="CEVR"] <- df$ev[df$contingency=="CEVR"] + 0.5 #offset for plotting
df$ev[df$contingency=="CEV"] <- df$ev[df$contingency=="CEV"] - 0.5 #offset for plotting
g3 <- ggplot(df, aes(x=time/1000, y=ev, color=contingency)) + gcommon + ylab("Expected value (points)") + theme(legend.position="none")

pdf("Fig1c_contingencies.pdf", width=8.5, height=4)
pg <- plot_grid(g1, g2, g3, nrow=1)
plot(pg)
dev.off()

pdf("Clock contingencies legend.pdf", width=2, height=2)
legend_b <- get_legend(g1 + theme(legend.position="right") + theme_bw(base_size=25) + guides(color = guide_legend(keywidth=2, keyheight=2)))

p <- plot_grid(legend_b, nrow=1)
plot(p)
dev.off()


# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
# p <- plot_grid(pg, legend_b, nrow=1, rel_widths = c(.9, .15))
# plot(p)

source("~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R")

fdf <- get_trial_data(repo_directory = "/Users/hallquist/Data_Analysis/clock_analysis", dataset = "mmclock_fmri") %>%
  mutate(id=as.integer(id))
mdf <- get_trial_data(repo_directory = "/Users/hallquist/Data_Analysis/clock_analysis", dataset = "mmclock_meg") %>%
  mutate(id=as.integer(id))

#DAN paper: Look up ids used in fMRI and MEG analyses to ensure that the plot matches other data
mr_ids <- data.table::fread("~/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofPittsburgh/DNPLskinner - SCEPTIC_fMRI/wholebrain_betas/L1m-echange/Schaefer_444_final_2009c_2.3mm_cope_l1.csv.gz") %>%
  pull(id) %>% unique()

fdf <- fdf %>% filter(id %in% !!mr_ids)

meg_ids <- readRDS("~/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofPittsburgh/DNPLskinner - SCEPTIC_fMRI/MEG/Time_RT_20Hz_n63/MEG0111_20Hz.rds") %>%
  dplyr::rename(id=Subject) %>% pull(id) %>% unique()

mdf <- mdf %>% dplyr::filter(id %in% !!meg_ids)

fdf <- fdf %>% group_by(id) %>% arrange(id, run, trial) %>% 
  mutate(totreward=sum(score_csv), cumreward=cumsum(score_csv)) %>% ungroup() %>%
  mutate(medreward=median(totreward), #between subjects
         msplit=factor(as.numeric(totreward > medreward), levels=c(0,1), labels=c("Total~earnings<median", "Total~earnings>=median")))

fdf$rewFunc <- factor(fdf$rewFunc, levels=c("IEV", "DEV", "CEV", "CEVR")) #to match contingency plot


mdf <- mdf %>% group_by(id) %>% arrange(id, run, trial) %>% 
  mutate(totreward=sum(score_csv), cumreward=cumsum(score_csv)) %>% ungroup() %>%
  mutate(medreward=median(totreward), #between subjects
         msplit=factor(as.numeric(totreward > medreward), levels=c(0,1), labels=c("Total~earnings<median", "Total~earnings>=median")))

mdf$rewFunc <- factor(mdf$rewFunc, levels=c("IEV", "DEV", "CEV", "CEVR")) #to match contingency plot

both <- dplyr::bind_rows(fdf, mdf)

# Figure 1b

pdf("Fig1b_earnings.pdf", width = 8.5, height = 3.75)
g <- ggplot(both %>% filter(run_trial <= 50), aes(x=run_trial, y=rt_csv, color = rewFunc)) + stat_smooth(method="loess", linewidth = 2) +
  #scale_color_brewer("Contingency", palette="Set2") +
  scale_color_manual("Contingency", values=pal) +
  theme_bw(base_size=18) + facet_wrap(~msplit, labeller=label_parsed) + xlab("Trial") +
  ylab("Response time (seconds)") + scale_y_continuous(breaks=c(1.0, 1.5, 2, 2.5)) +
  theme(axis.title.x=element_text(margin = margin(t = 12)), axis.title.y=element_text(margin = margin(r = 12)), 
        legend.margin = margin(t=0, r=2, b=0, l=5), plot.margin=margin(r=10, l=10, t=10, b=5)) + theme(legend.position="none")
  
plot(g)
dev.off()

# 1d
# pdf("Fig_1d.pdf", width = 10, height = 4)
# ggplot(subset(bdf), aes(x=trial, y=abstschange*100, color = rewFunc)) + stat_smooth(method="loess", size = 2) + theme_bw(base_size=25) + facet_wrap(~msplit) + ylab("RT swings, ms") + labs(colour = "Contingency") #facet_wrap(~msplit) #geom_jitter(alpha=0.2) +
# dev.off()

#pdf("Fig_1d.pdf", width = 8.5, height = 3.75)
pdf("Fig_1d_swings.pdf", width = 8.5, height = 3.75)
g <- ggplot(both %>% filter(run_trial <= 50 & rt_swing > 0 & rt_swing < 4), aes(x=run_trial, y=rt_swing, color = rewFunc)) + stat_smooth(method="loess", linewidth = 2) +
  #scale_color_brewer("Contingency", palette="Dark2") +
  scale_color_manual("Contingency", values=pal) +
  theme_bw(base_size=18) + facet_wrap(~msplit, labeller=label_parsed) + xlab("Trial") +
  ylab("Change in RT (seconds)") +
  theme(axis.title.x=element_text(margin = margin(t = 12)), axis.title.y=element_text(margin = margin(r = 12)), 
        legend.margin = margin(t=0, r=2, b=0, l=5), plot.margin=margin(r=10, l=10, t=10, b=5)) + theme(legend.position="none")

plot(g)
dev.off()
