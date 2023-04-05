##just a quick plot of the basis
library(ggplot2)
library(dplyr)
library(cowplot)
library(extrafont)
library(ggforce)
library(tidyverse)

plot_dir <- "~/OneDrive/collected_letters/papers/meg/figures/conceptual"
#font_import()
fonts()
loadfonts()
#extrafont::choose_font()
nbasis <- 24
tvec <- 0:4000
margin_offset <- .10
margin_offset = (max(tvec) - min(tvec))*margin_offset; #convert margin_offset into time scale of tvec
basis_overlap <- 1.52 #50%

tmin = min(tvec) - margin_offset
tmax=max(tvec) + margin_offset

#new basis plot for showing entropy etc.
centers <- seq(tmin, tmax, by = (tmax-tmin)/(nbasis-1))

gaussmat <- sapply(centers, function(x) {
  dvec <- dnorm(x=tvec, mean=x, sd=sig)
  dvec <- dvec/max(dvec) #renormalize to max=1
})
sig = (centers[2] - centers[1])/basis_overlap; #SD of the basis functions themselves


###########################
# basis examples
weights <- runif(nbasis, min=10, max=10)
#initial values for basis elements
v <- t(sapply(1:nrow(gaussmat), function(r) {
  gaussmat[r,]*weights
}))
weights_df <- data.frame(time=centers/1000, weight=weights) %>% filter(time > 0 & time <= 4)
vm <- reshape2::melt(v, varnames=c("time", "basis")) %>% mutate(basis=factor(basis), time=time/1000) #seconds
vfunc <- apply(v, 1, sum) %>% as.data.frame() %>% setNames("value") %>% mutate(time=1:length(value)/1000)

# polar version
fsize = 6
basis_range <- c(1:24) # which elements to plot
pdf(paste0("basis", length(basis_range), "_polar.pdf"),  height = 3, width = 3)
{if (length(basis_range) > 2) ggplot(vm, aes(x=time, y=value, group=basis, color = basis)) + scale_color_viridis(option = "magma", discrete = T) +
    geom_rect(data=weights_df, aes(xmin=time-.04, xmax=time+.04, ymin=0, ymax=weight, y=NULL, group=NULL, color=NULL),
              fill="grey70", show.legend = FALSE)} #weights
{if (length(basis_range) < 3) ggplot(vm %>% filter(basis %in% basis_range), aes(x=time, y=value, group=basis)) + 
    geom_rect(data=weights_df[basis_range - 2,], aes(xmin=time-.04, xmax=time+.04, ymin=0, ymax=weight, y=NULL, group=NULL, color=NULL),
              fill="grey70", show.legend = FALSE)} + #weights
  geom_line(show.legend = FALSE, size=1.2) + #basis elements
  cowplot::theme_cowplot(font_size = 20) + #ylab("Location value") +
  geom_hline(yintercept = 0, size = 1.2, color = "black") +
  annotate(geom="text", x=1, y= -5, label="1 s", size=fsize,  family=plot_font) +
  annotate(geom="text", x=2, y=- 4, label="2 s", size=fsize, family=plot_font) +
  annotate(geom="text", x=3, y=-5, label="3 s", size=fsize,  family=plot_font) +
  annotate(geom="text", x=3.8, y=-3, label="4 s", size=fsize,  family=plot_font) +
  
  annotate("segment", x = 0, y = -10, xend = 1.9, yend = -10,
           arrow = arrow(length = unit(0.5, "cm")), color = "black") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  ylim(c(-20,10)) + 
  coord_polar(theta = "x")
dev.off()
# linear version
pdf("basis_linear.pdf", height = 2, width = 5)
ggplot(vm, aes(x=time, y=value, group=basis, color = basis)) + scale_color_viridis(option = "magma", discrete = T) +
  geom_line(show.legend = FALSE, size=1.2) + #basis elements
  cowplot::theme_cowplot(font_size = 24) + ylab("Location value") +
  geom_hline(yintercept = 0, size = 1.5, color = "white") +
  xlab("Time, seconds") +
  # theme_minimal() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank())
dev.off()

# draw selections and outcomes

# just the selection
fsize = 6
pdf("rt_example_1s_110_points_line.pdf", height = 4, width = 4)
ggplot(vm , aes(x=time, y=value)) + 
  cowplot::theme_cowplot(font_size = 10) + #ylab("Location value") +
  geom_hline(yintercept = 0, size = 1.2, color = "black") +
  # annotate(geom="point", x=1, y= 3, size=fsize + 2,  shape = 16, color = "darkgreen", fill = "darkgreen") +
  # annotate(geom="point", x=.05, y= 3, size=fsize + 2,  shape = 16, color = "darkgreen", fill = "darkgreen") +
  annotate(geom="point", x=.05, y= 3, size=fsize + 2,  shape = 16, color = "darkgreen", fill = "darkgreen", alpha = .1) +
  annotate(geom="point", x=.2, y= 3, size=fsize + 2,  shape = 16, color = "darkgreen", fill = "darkgreen", alpha = .2) +
  annotate(geom="point", x=.4, y= 3, size=fsize + 2,  shape = 16, color = "darkgreen", fill = "darkgreen", alpha = .3) +
  annotate(geom="point", x=.6, y= 3, size=fsize + 2,  shape = 16, color = "darkgreen", fill = "darkgreen", alpha = .4) +
  annotate(geom="point", x=.8, y= 3, size=fsize + 2,  shape = 16, color = "darkgreen", fill = "darkgreen", alpha = .5) +
  annotate(geom="point", x=1, y= 3, size=fsize + 2,  shape = 16, color = "darkgreen", fill = "darkgreen", alpha = 1) +
  annotate(geom="text", x=1, y= -5, label="1 s", size=fsize,  family=plot_font) +
  # annotate(geom="text", x=1.1, y= 9, label="RT = 1 s", size=fsize,  family=plot_font) +
  annotate(geom="text", x=1, y= -20, label="You won\n110 points", size=fsize,  family=plot_font) +
  annotate(geom="text", x=2, y=- 4, label="2 s", size=fsize, family=plot_font, alpha = 1) +
  annotate(geom="text", x=2.95, y=-5, label="3 s", size=fsize,  family=plot_font, alpha = 1) +
  annotate(geom="text", x=3.9, y=-4, label="4 s      ", size=fsize,  family=plot_font, alpha = 1) +
  annotate("segment", x = 0.05, y = 8, xend = .9, yend = 8,
           arrow = arrow(length = unit(0.5, "cm")), color = "darkgreen", alpha = .8) +
  annotate(geom = "segment", x = 0, xend = 0, y = -1, yend = 7, size = 1.2, color="black", lineend="butt") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  ylim(c(-20,10)) + 
  coord_polar(theta = "x")
dev.off()

# basis update
time = "post_update"

weights <- runif(nbasis, min=10, max=10)

if (time == "post_update") {
  weights[6] <- 13
  weights[7] <- 24
  weights[8] <- 30
  weights[9] <- 17
  weights[10] <- 12
}
#initial values for basis elements
v <- t(sapply(1:nrow(gaussmat), function(r) {
  gaussmat[r,]*weights
}))
weights_df <- data.frame(time=centers/1000, weight=weights) %>% filter(time > 0 & time <= 4)
vm <- reshape2::melt(v, varnames=c("time", "basis")) %>% mutate(basis=factor(basis), time=time/1000) #seconds
vfunc <- apply(v, 1, sum) %>% as.data.frame() %>% setNames("value") %>% mutate(time=1:length(value)/1000)
fsize = 5
pdf(paste0("basis_", time, "_example_1s_110points.pdf"), height = 4, width = 4)
ggplot(vm, aes(x=time, y=value, group=basis, color = basis)) + scale_color_viridis(option = "magma", discrete = T) +
  geom_line(show.legend = FALSE, size=1.2) + #basis elements
  geom_rect(data=weights_df, aes(xmin=time-.04, xmax=time+.04, ymin=0, ymax=weight, y=NULL, group=NULL, color=NULL),
            fill="grey70", show.legend = FALSE) + #weights
  # geom_line(data=vfunc, aes(group=NULL), color="black", size=1.5) + #integrated value
  cowplot::theme_cowplot(font_size = 12) + #ylab("Location value") +
  geom_hline(yintercept = 0, size = 1.2, color = "black") +
  annotate(geom="text", x=1, y= -5, label="1 s ", size=fsize,  family=plot_font) +
  annotate(geom="text", x=1.9, y=- 4, label="2 s", size=fsize-2, family=plot_font, alpha = .5) +
  annotate(geom="text", x=2.9, y=-5, label="3 s", size=fsize-2,  family=plot_font, alpha = .5) +
  annotate(geom="text", x=3.85, y=-4, label="4 s    ", size=fsize-2,  family=plot_font, alpha = .5) +
  annotate(geom = "segment", x = 0, xend = 0, y = -4, yend = 15, size = 0.8, color="black", lineend="butt") +
  {if(time == "pre_update")
  annotate(geom="text", x=1, y= -20, label="You won\n110 points", size=fsize-2,  family=plot_font) +
  annotate(geom="point", x=1, y= 15, size=fsize + 2,  shape = 16, color = "darkgreen", fill = "darkgreen", alpha = 1)} +
  # annotate("segment", x = 0, y = -10, xend = 1.9, yend = -10,
  #          arrow = arrow(length = unit(0.5, "cm")), color = "black") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-3,4), "cm") 
  ) +
  ylim(c(-20,50)) + 
  xlim(c(0, 3.9)) +
  coord_polar(theta = "x")
dev.off()

##################
# entropy examples

gaussmat <- sapply(centers, function(x) {
  dvec <- dnorm(x=tvec, mean=x, sd=sig)
  dvec <- dvec/max(dvec) #renormalize to max=1
})

# pdf("gauss_basis_24basis.pdf", width=6, height=5)
# matplot(gaussmat, type="l", lty=1, lwd=5, col=colorRampPalette(c("blue", "red"))(13), ann=FALSE, xaxt='n', yaxt='n', bty='n')
# dev.off()

set.seed(1005)
weights <- runif(nbasis, min=0, max=8)
weights[1] <- 0
weights[24] <- 0
weights[23] <- 1


#initial values for basis elements
v <- t(sapply(1:nrow(gaussmat), function(r) {
  gaussmat[r,]*weights
}))

weights_df <- data.frame(time=centers/1000, weight=weights) %>% filter(time > 0 & time <= 4)
vm <- reshape2::melt(v, varnames=c("time", "basis")) %>% mutate(basis=factor(basis), time=time/1000) #seconds
vfunc <- apply(v, 1, sum) %>% as.data.frame() %>% setNames("value") %>% mutate(time=1:length(value)/1000)

plot_font <- "Helvetica"
setwd(plot_dir)
# pdf("hi_h.pdf", width=10, height=7)
pdf("hi_entropy_example_polar.pdf", height = 3, width = 3)
ggplot(vm, aes(x=time, y=value, group=basis, color = basis)) + scale_color_viridis(option = "magma", discrete = T) +
  geom_rect(data=weights_df, aes(xmin=time-.04, xmax=time+.04, ymin=0, ymax=weight, y=NULL, group=NULL, color=NULL),
            fill="grey70", show.legend = FALSE) + #weights
  geom_line(show.legend = FALSE, size=1.2) + #basis elements
  # geom_line(data=vfunc, aes(group=NULL), color="black", size=1.5) + #integrated value
  cowplot::theme_cowplot(font_size = 24) + ylab("Location value") +
  xlab("Location (seconds within the interval)") + 
  # annotate(geom="text", x=1.21, y=38.4, label="Reward", size=9, hjust=1, family=plot_font) +
  # annotate(geom="text", x=1.25, y=22, label="RPE+", size=9, hjust=1, family=plot_font) +
  annotate(geom="text", x=4, y=45, label="High\nentropy", size=12, hjust=1, family=plot_font, lineheight = 0.8) +
  # annotate(geom="point", x=1.4, y=38, size=9, color="darkblue") +
  # annotate(geom="segment", x=1.4, xend=1.4, y=8.5, yend=35.0, size=1.5, color="gray90", lineend="butt") +
  theme(text=element_text(family=plot_font), axis.title.x = element_text(margin=margin(t=15)), 
        axis.title.y = element_text(margin=margin(r=15))) + ylim(c(-10,12)) +
  theme(plot.title = element_text(face = "plain", size=32)) + coord_polar() + theme_minimal() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),panel.grid = element_blank())
# plot(g1)
dev.off()
# dev.off()

pdf("hi_entropy_example_linear.pdf", height = 2, width = 4)
ggplot(vm, aes(x=time, y=value, group=basis, color = basis)) + scale_color_viridis(option = "magma", discrete = T) +
  geom_rect(data=weights_df, aes(xmin=time-.04, xmax=time+.04, ymin=0, ymax=weight, y=NULL, group=NULL, color=NULL),
            fill="grey70", show.legend = FALSE) + #weights
  geom_line(show.legend = FALSE, size=1.2) + #basis elements
  geom_line(data=vfunc, aes(group=NULL), color="black", size=1.5) + #integrated value
  cowplot::theme_cowplot(font_size = 12) + ylab("Location value") +
  xlab("Location (seconds within the interval)") + 
  # annotate(geom="text", x=1.21, y=38.4, label="Reward", size=9, hjust=1, family=plot_font) +
  # annotate(geom="text", x=1.25, y=22, label="RPE+", size=9, hjust=1, family=plot_font) +
  annotate(geom="text", x=3, y=20, label="High entropy", size=8, hjust=1, family=plot_font, lineheight = 0.8) +
  # annotate(geom="point", x=1.4, y=38, size=9, color="darkblue") +
  # annotate(geom="segment", x=1.4, xend=1.4, y=8.5, yend=35.0, size=1.5, color="gray90", lineend="butt") +
  theme(text=element_text(family=plot_font), axis.title.x = element_text(margin=margin(t=15)), 
        axis.title.y = element_text(margin=margin(r=15))) + ylim(c(0,30)) +
  theme(plot.title = element_text(face = "plain", size=32)) + 
  theme_minimal() #+ 
# theme(axis.text = element_blank(),
#       axis.title = element_blank(),panel.grid = element_blank())
# plot(g1)
dev.off()
# dev.off()

# low-entropy example without compression


all <- list()
for (model in c("full", "selective")) {
  set.seed(1005)
  weights <- runif(nbasis, min=0, max=8)
  weights[1] <- 0
  weights[24] <- 0
  weights[23] <- 1
  if (model == "full") {
    weights[13] <- 10
    weights[14] <- 11
    weights[15] <- 20
    weights[16] <- 14
    weights[17] <- 11 } else if (model == "selective") {
      weights[1] <- 0
      weights[2] <- 0
      weights[3:13] <- runif(11, min = 1, max = 2.5)
      weights[14] <- 8
      weights[15] <- 20
      weights[16] <- 11
      weights[17:22] <- runif(6, min = 1, max = 2.5)
      weights[23] <- 1
    }
  v <- t(sapply(1:nrow(gaussmat), function(r) {
    gaussmat[r,]*weights
  }))
  weights_df <- data.frame(time=centers/1000, weight=weights) %>% filter(time > 0 & time <= 4)
  vm <- reshape2::melt(v, varnames=c("time", "basis")) %>% mutate(basis=factor(basis), time=time/1000) #seconds
  vfunc <- apply(v, 1, sum) %>% as.data.frame() %>% setNames("value") %>% mutate(time=1:length(value)/1000)
  all[[model]] <- list("weights_df" = weights_df, "vm" = vm, "vfunc" = vfunc)
}

# new version with single plot

colors = c("Information-\ncompressing RL" = "black", "Traditional RL" = "darkgrey", "Compression" = "darkgreen")

pdf(paste0("lo_entropy_example_both_models.pdf"), height = 2, width = 5)
ggplot(data = all$selective$vm, aes(x=time, y=value)) + 
  scale_color_manual(values = colors) +
  geom_line(data=all$selective$vfunc, aes(group=NULL, color="Information-\ncompressing RL"), size=1.2) + #integrated value
  geom_line(data=all$full$vfunc, aes(group=NULL , color="Traditional RL"), size=1.2) + #integrated value
  geom_rect(data=all$selective$weights_df, aes(xmin=time-.04, xmax=time+.04, ymin=0, ymax=weight, y=NULL, group=NULL, color=NULL),
            fill="grey70", show.legend = FALSE) + #weights
  labs(color = "") +
  new_scale_color() +
  geom_line(data = all$selective$vm, aes(group=basis, color = basis), alpha = 0.6, show.legend = FALSE, size=1.2) + #basis elements
  scale_color_viridis(option = "magma", discrete = T) +
  cowplot::theme_cowplot(font_size = 12) + ylab("Location value") +
  xlab("Location (seconds within the interval)") + 
  labs(color = "") +
  annotate(geom="text", x=4, y=26, label="Low\nentropy", size=5, hjust=1, family=plot_font, lineheight = 0.8) +
  annotate(geom="segment", x=0.5, xend=0.5, y=7, yend=4.3, size = 0.7, color="darkgreen",    arrow = arrow(length = unit(0.15,"cm"))) +
  annotate(geom="segment", x=2.05, xend=2.05, y=10, yend=4.6, size = 0.7, color="darkgreen",            arrow = arrow(length = unit(0.15,"cm"))) +
  annotate(geom="segment", x=3.05, xend=3.05, y=8.5, yend=5, size = 0.7, color="darkgreen",            arrow = arrow(length = unit(0.15,"cm"))) +
  annotate(geom="segment", x=3.4, xend=3.4, y=5.3, yend=4, size = 0.45, color="darkgreen",            arrow = arrow(length = unit(0.1,"cm"))) +
  annotate(geom="segment", x=3.9, xend=3.9, y=10, yend=4.1, size = 0.7, color="darkgreen", arrow = arrow(length = unit(0.15,"cm"))) +
  # annotate(geom="text", x=1.9, y=10, label="Compression", size=4.5, hjust=1, family=plot_font, lineheight = 0.8, color = "darkgreen", alpha = .8) +
  # annotate(geom="text", x=3.86, y=6.5, label="Compression", size=2.6, hjust=1, family=plot_font, lineheight = 0.8, color = "darkgreen", alpha = .7) +
  theme(text=element_text(family=plot_font), axis.title.x = element_text(margin=margin(t=15)), 
        axis.title.y = element_text(margin=margin(r=15))) + ylim(c(0,30)) + xlim(c(.1, 4)) +
  theme(plot.title = element_text(face = "plain", size=32)) + 
  theme_minimal() #+ 
# theme(axis.text = element_blank(),
#       axis.title = element_blank(),panel.grid = element_blank())
# plot(g1)
dev.off()




# # old version with separate plots
# pdf(paste0("lo_entropy_example_linear_", model, ".pdf"), height = 2, width = 4)
# ggplot(vm, aes(x=time, y=value, group=basis, color = basis)) + scale_color_viridis(option = "magma", discrete = T) +
#   geom_rect(data=weights_df, aes(xmin=time-.04, xmax=time+.04, ymin=0, ymax=weight, y=NULL, group=NULL, color=NULL),
#             fill="grey70", show.legend = FALSE) + #weights
#   geom_line(show.legend = FALSE, size=1.2) + #basis elements
#   geom_line(data=vfunc, aes(group=NULL), color="black", size=1.5) + #integrated value
#   cowplot::theme_cowplot(font_size = 12) + ylab("Location value") +
#   xlab("Location (seconds within the interval)") + 
#   # annotate(geom="text", x=1.21, y=38.4, label="Reward", size=9, hjust=1, family=plot_font) +
#   # annotate(geom="text", x=1.25, y=22, label="RPE+", size=9, hjust=1, family=plot_font) +
#   annotate(geom="text", x=4, y=45, label="High\nentropy", size=12, hjust=1, family=plot_font, lineheight = 0.8) +
#   # annotate(geom="point", x=1.4, y=38, size=9, color="darkblue") +
#   # annotate(geom="segment", x=1.4, xend=1.4, y=8.5, yend=35.0, size=1.5, color="gray90", lineend="butt") +
#   theme(text=element_text(family=plot_font), axis.title.x = element_text(margin=margin(t=15)), 
#         axis.title.y = element_text(margin=margin(r=15))) + ylim(c(0,30)) + xlim(c(.1, 4))
# theme(plot.title = element_text(face = "plain", size=32)) + 
#   theme_minimal() #+ 
# # theme(axis.text = element_blank(),
# #       axis.title = element_blank(),panel.grid = element_blank())
# # plot(g1)
# dev.off()


#first choice
weights[9] <- 22
weights[8] <- 11
weights[10] <- 12


#updated values for basis elements
v <- t(sapply(1:nrow(gaussmat), function(r) {
  gaussmat[r,]*weights
}))

weights_df <- data.frame(time=centers/1000, weight=weights) %>% filter(time > 0 & time <= 4)
vm <- reshape2::melt(v, varnames=c("time", "basis")) %>% mutate(basis=factor(basis), time=time/1000) #seconds
vfunc <- apply(v, 1, sum) %>% as.data.frame() %>% setNames("value") %>% mutate(time=1:length(value)/1000)

#dev.new()
#plot(vfunc, type="l")

#second choice
weights[16] <- 21
weights[17] <- 31
weights[18] <- 11

#updated values for basis elements
v_next <- t(sapply(1:nrow(gaussmat), function(r) {
  gaussmat[r,]*weights
}))

weights_df_next <- data.frame(time=centers/1000, weight=weights) %>% filter(time > 0 & time <= 4)
vm_next <- reshape2::melt(v, varnames=c("time", "basis")) %>% mutate(basis=factor(basis), time=time/1000) #seconds
vfunc_next <- apply(v_next, 1, sum) %>% as.data.frame() %>% setNames("value") %>% mutate(time=1:length(value)/1000)

pdf("g2.pdf", width=10, height=6)

g2 <- ggplot(vm, aes(x=time, y=value, group=basis)) + 
  geom_rect(data=weights_df, aes(xmin=time-.04, xmax=time+.04, ymin=0, ymax=weight, y=NULL, group=NULL, color=NULL), 
            fill="grey70", show.legend = FALSE) + #weights
  annotate(geom="segment", x=2.9, xend=2.9, y=5, yend=55, size=1.5, color="gray90", lineend="butt") +
  geom_line(show.legend = FALSE, size=1.2, color="red") + #basis elements
  geom_line(data=vfunc_next, aes(group=NULL), color="grey50", size=1.5, lty=6) + #integrated after update
  geom_line(data=vfunc, aes(group=NULL), color="black", size=1.5) + #integrated value
  cowplot::theme_cowplot(font_size = 24) + ylab("Expected Value (points)") +
  xlab("Time (seconds)") + annotate(geom="point", x=2.9, y=58, size=9, color="darkblue") +
  annotate(geom="text", x=2.7, y=58.5, label="Reward", size=9, hjust=1, family=plot_font) +
  #annotate(geom="text", x=3, y=30, label="PE+", size=9, hjust=1, family=plot_font) +
  annotate(geom="text", x=1.3, y=35, label=expression(RT[Vmax]), size=9, hjust=0.5, family=plot_font, parse=TRUE) +
  annotate(geom="segment", x=1.24, xend=2.75, y=41.5, yend=41.5, size=1.5, color="gray50", lineend="round", 
           arrow = arrow(length = unit(0.5,"cm"))) +
  annotate(geom="text", x=1.9, y=45.5, label=expression(Delta*RT[Vmax]), size=9, hjust=0.5, family=plot_font, color="gray50", parse=TRUE) +
  theme(text=element_text(family=plot_font), axis.title.x = element_text(margin=margin(t=15)), 
        axis.title.y = element_text(margin=margin(r=15))) + ylim(c(0,60)) +
  theme(plot.title = element_text(face = "plain", size=32))
plot(g2)
dev.off()

pdf("two_panel.pdf", width=12, height=5)
plot_grid(g1 + ggtitle("Early in learning"), g2 + ggtitle("Change in RT Vmax") + ylab(""), nrow=1)
dev.off()


#nth choice
weights <- runif(nbasis, min=0, max=4)
weights[1:3] <- 0
weights[22:24] <- 0
weights[13] <- 12
weights[14] <- 35
weights[15] <- 26
weights[16] <- 9
weights[17] <- 6

#updated values for basis elements
v <- t(sapply(1:nrow(gaussmat), function(r) {
  gaussmat[r,]*weights
}))

weights_df <- data.frame(time=centers/1000, weight=weights) %>% filter(time > 0 & time <= 4)
vm <- reshape2::melt(v, varnames=c("time", "basis")) %>% mutate(basis=factor(basis), time=time/1000) #seconds
vfunc <- apply(v, 1, sum) %>% as.data.frame() %>% setNames("value") %>% mutate(time=1:length(value)/1000)

# pdf("lo_h.pdf", width=10, height=6)

g3 <- ggplot(vm, aes(x=time, y=value, group=basis)) + 
  geom_rect(data=weights_df, aes(xmin=time-.04, xmax=time+.04, ymin=0, ymax=weight, y=NULL, group=NULL, color=NULL), 
            fill="grey70", show.legend = FALSE) + #weights
  geom_line(show.legend = FALSE, size=1.2, color="red") + #basis elements
  geom_line(data=vfunc, aes(group=NULL), color="black", size=1.5) + #integrated value
  cowplot::theme_cowplot(font_size = 24) + ylab("Location value") +
  xlab("Location (seconds within the interval)") + 
  #annotate(geom="point", x=3.10, y=58, size=9, color="darkblue") +
  #annotate(geom="text", x=2.95, y=58.5, label="Reward", size=9, hjust=1, family=plot_font) +
  #annotate(geom="text", x=3, y=30, label="PE+", size=9, hjust=1, family=plot_font) +
  annotate(geom="text", x=2.5, y=54, label=expression(RT[Vmax]), size=9, hjust=0.5, family=plot_font, parse=TRUE) +
  annotate(geom="text", x=0.1, y=45, label="Low\nentropy", size=12, hjust=0, family=plot_font, lineheight = 0.8) +
  #annotate(geom="segment", x=3.10, xend=3.10, y=5, yend=55, size=1.5, color="gray60", lineend="butt") +
  #annotate(geom="segment", x=1.24, xend=2.9, y=42, yend=42, size=1.5, color="gray60", lineend="round", 
  #         arrow = arrow(length = unit(0.5,"cm"))) +
  theme(text=element_text(family=plot_font), axis.title.x = element_text(margin=margin(t=15)), 
        axis.title.y = element_text(margin=margin(r=15))) + ylim(c(0,60)) + 
  theme(plot.title = element_text(face = "plain", size=32))

plot(g3)
# dev.off()

#pdf("three_panel.pdf", width=18, height=6)
gp <- plot_grid(g1 + ggtitle("Early in learning"), g2 + ggtitle("Change in RT Vmax") + ylab(""),
                g3 + ggtitle("Late in learning") + ylab(""), nrow=1)
# plot_grid(g1, g2 + ylab(""), g3 + ylab(""), label_size = 18, vjust=1, hjust=0.1,
#           nrow=1, labels=c("a) Early in learning", "b) Change in RT Vmax", "c) Late in learning"))
#dev.off()
ggsave("three_panel.pdf", gp, width=18, height=6, useDingbats=FALSE)

embed_fonts('three_panel.pdf', outfile='three_panel_embed.pdf')





#dev.new()
#plot(vfunc, type="l")

# pdf("g3.pdf", width=10, height=6)
# 
# g2 <- ggplot(vm, aes(x=time, y=value, group=basis)) + 
#   geom_rect(data=weights_df, aes(xmin=time-.04, xmax=time+.04, ymin=0, ymax=weight, y=NULL, group=NULL, color=NULL), 
#             fill="grey70", show.legend = FALSE) + #weights
#   geom_line(show.legend = FALSE, size=1.2, color="red") + #basis elements
#   geom_line(data=vfunc, aes(group=NULL), color="black", size=1.5) + #integrated value
#   cowplot::theme_cowplot(font_size = 24) + ylab("Expected Value (points)") +
#   xlab("Time (seconds)") + annotate(geom="point", x=3.1, y=58, size=9, color="darkblue") +
#   annotate(geom="text", x=2.90, y=58.5, label="Reward", size=9, hjust=1, family=plot_font) +
#   #annotate(geom="text", x=2.95, y=30, label="PE+", size=9, hjust=1, family=plot_font) +
#   annotate(geom="text", x=1.3, y=35, label="RT Vmax", size=9, hjust=0.5, family=plot_font) +
#   annotate(geom="segment", x=3.1, xend=3.1, y=5, yend=55, size=1.5, color="gray60", lineend="butt") +
#   theme(text=element_text(family=plot_font), axis.title.x = element_text(margin=margin(t=15)), 
#         axis.title.y = element_text(margin=margin(r=15))) + ylim(c(0,60))
# plot(g2)
# dev.off()

###OLD CODE


ntimesteps=13
maxt=4000
centers <- seq(0,4000, length=ntimesteps)
gaussmat <- sapply(centers, function(v) {
  dnorm(x=0:maxt, mean=v, sd=300)
})

#matplot(gaussmat, type="l", lty=1, lwd=5, col=colors(distinct=TRUE))
pdf("gauss_basis.pdf", width=6, height=5)
matplot(gaussmat, type="l", lty=1, lwd=5, col=colorRampPalette(c("blue", "red"))(13), ann=FALSE, xaxt='n', yaxt='n', bty='n')
dev.off()




##test out weights and centers to show example
weights<- 100*c(0.0776, 7.4801, 2.0792, 1.4008, 2.0000, 1.0000, 3.0000, 9.0000, 25.0000, 22.0024, 30.2286, 12.0326, 30.2332)
v <- 5* apply( sapply(1:nrow(gaussmat), function(r) {
  gaussmat[r,]*weights
}), 2, sum)



pdf("learned_ev.pdf", width=6, height=5)
df <- data.frame(time=seq(0,4,length=4001), ev=v)
ggplot(df, aes(x=time, y=ev)) + geom_line(size=1.5) + ylab("Expected value (learned)\n") + xlab("\nTime (s)") + theme_bw(base_size=24)
dev.off()


##uncertainty weights
uweights<- -10*c(0.3776, 1.4801, 0.9792, 1.4008, 2.0000, 1.0000, 3.0000, 5.0000, 7.0000, 8.0024, 10.2286, 3.0326, 1.2332)
u <- 5* apply( sapply(1:nrow(gaussmat), function(r) {
  gaussmat[r,]*uweights
}), 2, sum)


library(ggplot2)
pdf("learned_uncertainty.pdf", width=6, height=5)
df <- data.frame(time=seq(0,4,length=4001), u=u)
ggplot(df, aes(x=time, y=u)) + geom_line(size=1.5) + ylab("Uncertainty (experienced)\n") + xlab("\nTime (s)") + theme_bw(base_size=24)
dev.off()


pdf("weights.pdf", width=6, height=3.5)
par(mar=c(1, 5, 2, 1) + 0.1)

plot(c(0, 4000), c(0, max(weights) + 50), type = "n", xlab = "", ylab = "Basis weight (AU)", yaxs="i",
     cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, bty='n', xaxt='n')
rect(centers-50, 0, centers+50, weights, col="gray")
#axis(1, at=c(0,4), labels=c("0s", "4s"))#, pos=, lty=, col=, las=, tck=, ...)
dev.off()


##true underlying IEV contingency.
setwd("~/Data_Analysis/clock_analysis/td_model")
source("~/Data_Analysis/clock_analysis/td_model/getrew.R")
fm <- getMagFreq(0:4000, "IEV")
f <- fm[4002:8002]
m <- fm[1:4001]
ev <- f*m

library(ggplot2)
pdf("true_ev.pdf", width=6, height=5)
df <- data.frame(time=seq(0,4,length=4001), ev=ev)
ggplot(df, aes(x=time, y=ev)) + geom_line(size=1.5) + ylab("Expected value\n") + xlab("\nTime (s)") + theme_bw(base_size=24)
dev.off()

##real curve
pdf("iev_func.pdf", width=5, height=4)
ggplot(df, aes(x=time, y=ev)) + geom_line(size=1.5) + ylab("Expected value\n") + xlab("\nTime (ms)") + theme_bw(base_size=24)
dev.off()

df <- c()
#take CEVR out since CEV and CEVR are identical wrt EV (and we are not showing prob + freq)
for (cont in c("IEV", "DEV", "CEV", "CEVR")) { #, "CEVR"
  fm <- getMagFreq(0:4000, cont)
  #if (cont=="CEV") { cont="CEV/CEVR" } #for plot name
  df <- rbind(df, data.frame(contingency=cont, time=0:4000, mag=fm$Mag, freq=fm$Freq, ev=fm$Mag*fm$Freq))
}

#Figure: plot of EV in clock task
pdf("Clock contingencies.pdf", width=5, height=3.4)
ggplot(df, aes(x=time/1000, y=ev, color=contingency)) + geom_line(size=2) + ylab("Expected value (points)") + xlab("Time (seconds)") + scale_color_brewer("Contingency", palette="Dark2") +
  theme_bw(base_size=18) + theme(axis.title.x=element_text(margin = margin(t = 10)), axis.title.y=element_text(margin = margin(r = 10)), legend.margin = unit(0.15, "cm"), 
                                 plot.margin=margin(r=3, l=3, t=10, b=10))
dev.off()

#freq, prob, and ev
library(cowplot)
gcommon <- list(geom_line(size=2), xlab("Time (seconds)"), scale_color_brewer("Contingency", palette="Dark2"), theme_bw(base_size=18), 
                theme(axis.title.x=element_text(margin = margin(t = 10)), axis.title.y=element_text(margin = margin(r = 8)), 
                      legend.margin = unit(0.15, "cm"), 
                      plot.margin=margin(r=10, l=10, t=10, b=5)))


g1 <- ggplot(df, aes(x=time/1000, y=mag, color=contingency)) + gcommon + ylab("Reward magnitude (points)") + 
  theme(legend.position="none", plot.margin=margin(r=10, l=0, t=10, b=5))

g2 <- ggplot(df, aes(x=time/1000, y=freq, color=contingency)) + gcommon + ylab("Reward probability") + theme(legend.position="none")

df$ev[df$contingency=="CEVR"] <- df$ev[df$contingency=="CEVR"] + 0.5 #offset for plotting
df$ev[df$contingency=="CEV"] <- df$ev[df$contingency=="CEV"] - 0.5 #offset for plotting
g3 <- ggplot(df, aes(x=time/1000, y=ev, color=contingency)) + gcommon + ylab("Expected value (points)") + theme(legend.position="none")

pdf("Clock contingencies with freq mag.pdf", width=8.5, height=4)
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


setwd("~/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri")
#load(file="dataframe_for_entropy_analysis_Oct2016.RData")
#this contains data with 24 basis functions and post-Niv learning rule
#load(file="dataframe_for_entropy_analysis_Nov2016.RData")
load(file="dataframe_for_entropy_analysis_Mar2017.RData") #has the random priors entropy

bdf = bdf %>% rename(subject=rowID) %>% group_by(subject) %>% arrange(subject, run, trial) %>% mutate(totreward=sum(score), cumreward=cumsum(score)) %>% ungroup() %>%
  mutate(medreward=median(totreward), #between subjects
         msplit=factor(as.numeric(totreward > medreward), levels=c(0,1), labels=c("< Median", ">= Median")))

bdf <- bdf %>% mutate(msplit=recode(msplit, "< Median" = "Total~earnings<median", ">= Median"="Total~earnings>=median"))
bdf$rewFunc <- factor(bdf$rewFunc, levels=c("IEV", "DEV", "CEV", "CEVR")) #to match contingency plot

# Figure 1c
pdf("Fig_1c.pdf", width = 8.5, height = 3.75)
ggplot(bdf, aes(x=trial, y=rt/1000, color = rewFunc )) + stat_smooth(method="loess", size = 2) +
  scale_color_brewer("Contingency", palette="Dark2") +
  theme_bw(base_size=18) + facet_wrap(~msplit, labeller=label_parsed) + xlab("Trial") +
  ylab("Response time (seconds)") + scale_y_continuous(breaks=c(1.25, 1.5, 1.75, 2, 2.25)) +
  theme(axis.title.x=element_text(margin = margin(t = 12)), axis.title.y=element_text(margin = margin(r = 12)), 
        legend.margin = margin(t=0, r=2, b=0, l=5), plot.margin=margin(r=10, l=10, t=10, b=5)) + theme(legend.position="none")

dev.off()

# 1d
# pdf("Fig_1d.pdf", width = 10, height = 4)
# ggplot(subset(bdf), aes(x=trial, y=abstschange*100, color = rewFunc)) + stat_smooth(method="loess", size = 2) + theme_bw(base_size=25) + facet_wrap(~msplit) + ylab("RT swings, ms") + labs(colour = "Contingency") #facet_wrap(~msplit) #geom_jitter(alpha=0.2) +
# dev.off()

pdf("Fig_1d.pdf", width = 8.5, height = 3.75)
ggplot(bdf, aes(x=trial, y=abstschange/10, color = rewFunc )) + stat_smooth(method="loess", size = 2) +
  scale_color_brewer("Contingency", palette="Dark2") +
  theme_bw(base_size=18) + facet_wrap(~msplit, labeller=label_parsed) + xlab("Trial") +
  ylab("Change in RT (seconds)") +
  theme(axis.title.x=element_text(margin = margin(t = 12)), axis.title.y=element_text(margin = margin(r = 12)), 
        legend.margin = margin(t=0, r=2, b=0, l=5), plot.margin=margin(r=10, l=10, t=10, b=5)) + theme(legend.position="none")
dev.off()




#Supplementary figure: Sinusoidal contingency
ntimesteps=500

ev = 10*sin(2*pi*(1:ntimesteps)*1/ntimesteps) + 2.5*sin(2*pi*(1:ntimesteps)*2/ntimesteps) + 2.0*cos(2*pi*(1:ntimesteps)*4/ntimesteps)
ev = ev + abs(min(ev)) + 10;
prb = 25*cos(2*pi*(1:ntimesteps)*1/ntimesteps) + 10*cos(2*pi*(1:ntimesteps)*3/ntimesteps) + 6*sin(2*pi*(1:ntimesteps)*5/ntimesteps)
prb_max=0.7
prb_min=0.3
prb = (prb - min(prb))*(prb_max-prb_min)/(max(prb)-min(prb)) + prb_min

allshift = array(NA_real_, dim=c(ntimesteps, ntimesteps, 3))

for (i in 1:ntimesteps) {
  if (i > 1) {
    shift = c(i:ntimesteps, 1:(i-1))
  } else { shift <- 1:ntimesteps }
  evi = ev[shift]
  prbi = prb[shift]
  
  allshift[i,,1] = evi
  allshift[i,,2] = prbi
  allshift[i,,3] = evi/prbi
  
}

shift3 <- rbind(data.frame(time=1:ntimesteps/100, EV=allshift[,1,1], Probability=allshift[,1,2], Magnitude=allshift[,1,3], name="shift = 0"),
                data.frame(time=1:ntimesteps/100, EV=allshift[,100,1], Probability=allshift[,100,2], Magnitude=allshift[,100,3], name="shift = 100"),
                data.frame(time=1:ntimesteps/100, EV=allshift[,200,1], Probability=allshift[,200,2], Magnitude=allshift[,200,3], name="shift = 200"))

library(reshape2); library(ggplot2)
m3 <- melt(shift3, id.vars=c("time", "name"))

pdf("Sinusoid contingency.pdf", width=10, height=8)
ggplot(m3, aes(x=time, y=value)) + geom_line(size=1.5) + facet_grid(variable ~ name, scales="free_y") + theme_bw(base_size=24) + xlab("Time (seconds)") + ylab("") +
  theme(panel.margin = unit(20, "pt"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = "grey90", colour = "grey50", size = 0.2))
dev.off()

##BMC Figure
library(R.matlab)
setwd(file.path(getMainDir(), "temporal_instrumental_agent", "clock_task", "figures"))

scepticbmc <- readMat("finalicissimo_BMC_for_eLife_fig.mat")

#out.Ef contains estimated frequencies
#out.Vf contains the variance-covariance matrix of frequencies
#inside plotUncertainTimeSeries, which is called from VBA_groupBMC, it appears the SEs are derived by the sqrt of the diagonal of out.Vf

#somehow got mangled -- went into MATLAB and just saved these in a simpler .mat 
Ef <- scepticbmc$out[,,1]$Ef
scepticbmc$out[,,1]$Vf

#bmcef <- readMat("elife_bmc_frequencies.mat")
bmcef <- readMat("~/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri/figures/ploscompbio_bmc_frequencies.mat")

df <- data.frame(model=unlist(bmcef$modelnames), freq=bmcef$Ef, se=sqrt(diag(bmcef$Vf)))
#df$m_ordered <- ordered(df$model, levels=c("fixed", "fixed_uv", "fixed_decay", "kalman_softmax",
#        "kalman_uv_sum", "kalman_logistic", "kalman_processnoise", "kalman_sigmavolatility", "Qstep"),
#    labels=c("Fixed LR V", "Fixed LR U + V", "Fixed LR V Decay", "KF V", "KF U + V", 
#        "KF U -> V", "KF Process Noise", "KF Volatility", "TD"))

df$m_ordered <- ordered(df$model, levels=c("fixed", "fixed_uv", "fixed_decay", "kalman_softmax",
                                           "kalman_uv_sum", "Qstep"),
                        labels=c("Fixed LR V", "Fixed LR U + V", "Fixed LR V Sel. Maint.", "KF V", "KF U + V", "TD"))


#for ggplot with coord_flip, need to reverse
df$m_ordered <- factor(df$m_ordered, levels=rev(levels(df$m_ordered))) 

library(ggplot2)

#updated version with smaller model set for PLoS Comp Bio
pdf("SCEPTIC Main BMC v5 May2017.pdf", width=6, height=4)
ggplot(df, aes(x=m_ordered, y=freq, ymin=freq-se, ymax=freq+se)) + geom_pointrange(stat="identity", size=1.3, fatten=2.3) +
  annotate("text", x=4, y=0.65, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
  geom_label(mapping=aes(x=x,y=y,label=label, ymin=NULL, ymax=NULL), 
             data=data.frame(x=0.65, y=0.48, label=as.character(expression(paste("BOR < ",10^{-51})))),
             hjust=0, vjust=0, parse=TRUE, size=6, label.padding = unit(0.4, "lines")) +
  ylab("Estimated Model Frequency") + xlab("") + coord_flip() +
  theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)),
                                 panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(), plot.margin=margin(t=5, r=10, b=5, l=0)) +
  scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75), labels=c("0", ".25", ".5", ".75"))
dev.off()



pdf("SCEPTIC Main BMC.pdf", width=6, height=4)
ggplot(df, aes(x=m_ordered, y=freq, ymin=freq-se, ymax=freq+se)) + geom_pointrange(stat="identity") + geom_errorbar(width=0.5) +
  annotate("text", x=7, y=0.54, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
  annotate("text", x=1.2, y=0.30, label=as.character(expression(paste("BOR = ",8.03," x ",10^{-49}))), hjust=0, vjust=0, parse=TRUE, size=6) +
  ylab("Estimated Model Frequency") + xlab("") + coord_flip() +
  theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)))
dev.off()


library(ggplot2)
pdf("SCEPTIC Main BMC v2.pdf", width=6, height=4)
ggplot(df, aes(x=m_ordered, y=freq, ymin=freq-se, ymax=freq+se)) + geom_bar(stat="identity", fill="grey92", color="black") + geom_errorbar(width=0.5) +
  annotate("text", x=7, y=0.54, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
  annotate("text", x=1.2, y=0.30, label=as.character(expression(paste("BOR = ",8.03," x ",10^{-49}))), hjust=0, vjust=0, parse=TRUE, size=6) +
  ylab("Estimated Model Frequency") + xlab("") + coord_flip() +
  theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)))
dev.off()


library(ggplot2)
pdf("SCEPTIC Main BMC v3.pdf", width=6, height=4)
ggplot(df, aes(x=m_ordered, y=freq, ymin=freq-se, ymax=freq+se)) + geom_bar(stat="identity", fill="grey92", color="black") + geom_errorbar(width=0.5) +
  annotate("text", x=7, y=0.54, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
  geom_label(mapping=aes(x=x,y=y,label=label, ymin=NULL, ymax=NULL), 
             data=data.frame(x=0.8, y=0.47, label=as.character(expression(paste("BOR < ",10^{-49})))),
             hjust=0, vjust=0, parse=TRUE, size=6) +
  ylab("Estimated Model Frequency") + xlab("") + coord_flip() +
  theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)))
dev.off()

pdf("SCEPTIC Main BMC v4.pdf", width=6, height=4)
ggplot(df, aes(x=m_ordered, y=freq, ymin=freq-se, ymax=freq+se)) + geom_bar(stat="identity", fill="grey92", color="black") + geom_errorbar(width=0.5) +
  annotate("text", x=7, y=0.54, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
  geom_label(mapping=aes(x=x,y=y,label=label, ymin=NULL, ymax=NULL), 
             data=data.frame(x=0.8, y=0.47, label=as.character(expression(paste("BOR < ",10^{-49})))),
             hjust=0, vjust=0, parse=TRUE, size=6) +
  ylab("Estimated Model Frequency") + xlab("") + coord_flip() +
  theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)),
                                 panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank())
dev.off()

pdf("SCEPTIC Main BMC v5.pdf", width=6, height=4)
ggplot(df, aes(x=m_ordered, y=freq, ymin=freq-se, ymax=freq+se)) + geom_pointrange(stat="identity", size=1.3, fatten=2.5) +
  annotate("text", x=7, y=0.54, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
  geom_label(mapping=aes(x=x,y=y,label=label, ymin=NULL, ymax=NULL), 
             data=data.frame(x=0.8, y=0.47, label=as.character(expression(paste("BOR < ",10^{-43})))),
             hjust=0, vjust=0, parse=TRUE, size=6, label.padding = unit(0.4, "lines")) +
  ylab("Estimated Model Frequency") + xlab("") + coord_flip() +
  theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)),
                                 panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) +
  scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75), labels=c("0", ".25", ".5", ".75"))
dev.off()

#ar1 and schoenberg results
arfreqs <- readMat("ar_modelfreqs_Sep2016.mat")


#manual entry from Jon email 19Sep2016
mnames <- c("Fixed LR V",	"Fixed LR U + V", "Fixed LR V Sel. Maint.", "KF V", "KF Process Noise", "KF U + V", "KF Volatility")
df <- data.frame(model=ordered(mnames), freq=arfreqs$ar1Ef, se=sqrt(diag(arfreqs$ar1Vf)))

#for ggplot with coord_flip, need to reverse
df$m_ordered <- factor(df$model, levels=rev(levels(df$model))) 

pdf("Ar1 Main BMC v5.pdf", width=6, height=4.3)
ggplot(df, aes(x=m_ordered, y=freq, ymin=freq-se, ymax=freq+se)) + geom_pointrange(stat="identity", size=1.3, fatten=2.5) +
  annotate("text", x=5, y=0.45, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
  geom_label(mapping=aes(x=x,y=y,label=label, ymin=NULL, ymax=NULL), 
             data=data.frame(x=0.8, y=0.40, label=as.character(expression(paste("BOR < ",10^{-32})))),
             hjust=0, vjust=0, parse=TRUE, size=6, label.padding = unit(0.4, "lines")) +
  ylab("Estimated Model Frequency") + xlab("Includes AR(1) choice") + coord_flip() +
  theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)),
                                 panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) +
  scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75), labels=c("0", ".25", ".5", ".75"))
dev.off()

#manual entry from Jon email 19Sep2016
df <- data.frame(model=ordered(mnames), freq=arfreqs$schEf, se=sqrt(diag(arfreqs$schVf)))

#for ggplot with coord_flip, need to reverse
df$m_ordered <- factor(df$model, levels=rev(levels(df$model))) 

pdf("Scho Main BMC v5.pdf", width=6, height=4.3)
ggplot(df, aes(x=m_ordered, y=freq, ymin=freq-se, ymax=freq+se)) + geom_pointrange(stat="identity", size=1.3, fatten=2.5) +
  annotate("text", x=5, y=0.45, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
  geom_label(mapping=aes(x=x,y=y,label=label, ymin=NULL, ymax=NULL), 
             data=data.frame(x=0.8, y=0.43, label=as.character(expression(paste("BOR < ",10^{-37})))),
             hjust=0, vjust=0, parse=TRUE, size=6, label.padding = unit(0.4, "lines")) +
  ylab("Estimated Model Frequency") + xlab("Includes Schoenberg choice") + coord_flip() +
  theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)),
                                 panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) +
  scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75), labels=c("0", ".25", ".5", ".75"))
dev.off()


#frank TC (replicate v5 above)
frank_bmcef <- readMat("elife_bmc_franktc_frequencies.mat")

df <- data.frame(model=unlist(frank_bmcef$models), freq=frank_bmcef$Ef, se=sqrt(diag(frank_bmcef$Vf)))

df$modelmath <- ordered(df$model, levels=c("K", "K_Lambda", "K_Lambda_Nu", "K_Lambda_Nu_AlphaG",
                                           "K_Lambda_Nu_AlphaG_AlphaN", "K_Lambda_Nu_AlphaG_AlphaN_Rho", "K_Lambda_Nu_AlphaG_AlphaN_Rho_Epsilon"))

#for ggplot with coord_flip, need to reverse
df$modelmath <- factor(df$modelmath, levels=rev(levels(df$modelmath))) 


pdf("Frank TC BMC v5.pdf", width=6, height=4)
ggplot(df, aes(x=modelmath, y=freq, ymin=freq-se, ymax=freq+se)) + geom_pointrange(stat="identity", size=1.3, fatten=2.5) +
  annotate("text", x=5, y=0.63, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
  geom_label(mapping=aes(x=x,y=y,label=label, ymin=NULL, ymax=NULL), 
             data=data.frame(x=0.8, y=0.60, label=as.character(expression(paste("BOR < ",10^{-35})))),
             hjust=0, vjust=0, parse=TRUE, size=6, label.padding = unit(0.4, "lines")) +
  ylab("Estimated Model Frequency") + xlab("") + coord_flip() +
  theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)),
                                 axis.title.y=element_text(margin = margin(r = 15)),
                                 panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) +
  scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75), labels=c("0", ".25", ".5", ".75")) +
  scale_x_discrete("Parameter added to TC", labels=rev(expression(K, lambda, nu, alpha[G], alpha[N], rho, epsilon)))

#scale_x_discrete("test", labels=c(expression(alpha), expression(beta)))
dev.off()