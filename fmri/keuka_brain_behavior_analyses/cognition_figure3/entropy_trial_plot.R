#fixed versus selective entropy evolution
library(R.matlab)
library(tidyverse)
fig_folder <- "~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/cognition_figure3"
setwd(fig_folder)
fixed <- readMat("fixed_0_multi.mat")
decay <- readMat("fixed_decay_0_multi.mat")
centers <- as.vector(readMat("centers_from_rbf.mat")$c / 10)

getTrialHeights <- function(mat, trial, normalize=TRUE) {
  tvalue <- mat[[1]][[1]][,trial] #funny storage structure from MATLAB import
  if (normalize) { tvalue <- tvalue/max(tvalue) }
  return(tvalue)
}

#trial 3 data
t3 <- bind_rows(data.frame(height=getTrialHeights(fixed, 203), rbf=1:24, time=centers, model="Selective Maintenance", stringsAsFactors = FALSE), #there is an odd flip in t3 for full/selective (selective was higher)
                data.frame(height=getTrialHeights(decay, 203), rbf=1:24, time=centers, model="Full Maintenance", stringsAsFactors = FALSE))
t3$trial <- "Trial 3"

#trial 50 data
t50 <- bind_rows(data.frame(height=getTrialHeights(fixed, 250), rbf=1:24, time=centers, model="Full Maintenance", stringsAsFactors = FALSE),
                data.frame(height=getTrialHeights(decay, 250), rbf=1:24, time=centers, model="Selective Maintenance", stringsAsFactors = FALSE))
t50$trial <- "Trial 50"

both <- bind_rows(t3, t50)

#pdf("value_vector_v1.pdf", width=13, height=6)
g <- ggplot(both, aes(x=time, y=height, color=model)) +
  geom_point(position=position_dodge(width=0.13), size=2, shape = 21, stroke=1.3) +
  geom_line(position=position_dodge(width=0.13), size=1.1) +
  facet_wrap(~trial) + theme_bw(base_size=20) +
  scale_color_brewer("", palette="Dark2") +
  xlab("Time (s)") + ylab("Relative value (TBF weights; a.u.)") +
  theme(axis.title.y = element_text(margin=margin(r=20)),
        axis.title.x = element_text(margin=margin(t=20)),
        panel.spacing = unit(20, "pt"),
        legend.position="top"
  ) + xlim(c(-0.1,4.1))
ggsave(filename="value_vector_v1.pdf", width=10, height=5.5, plot=g, useDingbats=FALSE) #turn off dingbats, which causes odd problem in Illustrator import
#dev.off()

pdf("value_vector_v2.pdf", width=10, height=7)
ggplot(both, aes(x=time, y=height)) +
  geom_point(position=position_dodge(width=0.13), size=2, shape = 21) +
  geom_line(position=position_dodge(width=0.13), size=1.1) +
  facet_grid(model ~ trial) + theme_bw(base_size=20) +
  xlab("Time (s)") + ylab("Relative value (normalized TBF weights)") +
  theme(axis.title.y = element_text(margin=margin(r=20)),
        axis.title.x = element_text(margin=margin(t=20)),
        panel.spacing = unit(20, "pt")
  )
dev.off()

#entropy plot
entropy_df <- data.frame(
  entropy_fixed = fixed$fixed.0.multi[[2]][203:250],
  entropy_decay = decay$fixed.decay.0.multi[[2]][203:250],
  trial = 3:50) %>%
  gather(key="model", value="entropy", -trial) %>%
  mutate(model=recode(model, entropy_fixed="Full Maintenance", entropy_decay="Selective Maintenance"))

pdf("entropy_unfolding.pdf", width=6, height=5)
ggplot(entropy_df, aes(x=trial, y=entropy, color=model)) + geom_line(size=1.3, show.legend = FALSE) +
  theme_bw(base_size=20) + xlab("Trial") + ylab("Shannon's Entropy (nats)") +
  scale_color_brewer("", palette="Dark2") +
  theme(axis.title.y = element_text(margin=margin(r=20)),
        axis.title.x = element_text(margin=margin(t=20))
  )
dev.off()
