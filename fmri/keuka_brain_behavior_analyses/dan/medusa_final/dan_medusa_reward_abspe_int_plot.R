##################
# little custom script for understanding whether PEs in the DAN are primarily signed or absolute
# first run dan_medusa_fmri_encode_meg_betas.R


# tdf <- trial_df %>% filter(abs_pe < 30)
# ggplot(trial_df, aes(abs_pe)) + geom_histogram() + facet_grid(~outcome)
if (filter_abspe) {fsuffix = "abspe_reward_interaction_filtered.pdf"} else {
  fsuffix = "abspe_reward_interaction_unfiltered.pdf"
}
fname <- file.path(out_dir, "rt_encode", fsuffix)

# df <- ddf$emmeans_list$emm_1 %>% setNames(make.names(names(.), unique = TRUE)) %>% 
#   select(-matches("*\\.[1-9]+$"))
# ggplot(df, aes_string( x="evt_time", y="estimate", ymin="estimate - std.error", shape = "as.factor(abs_pe)", lty = "as.factor(abs_pe)", ymax="estimate + std.error", color="vm_gradient")) + 
#   geom_line(size=1, position=position_dodge(width=0.4)) + 
#   geom_pointrange(position=position_dodge(width=0.4)) +
#   scale_color_brewer(palette="Dark2", labels=c("1" = "MT+, control", "2" = "Parieto-occipital", "3" = "Post. parietal", "4" = "Frontal")) +
#   geom_hline(yintercept = 0, size=1.5, alpha=0.6) +
#   geom_vline(xintercept = 0, size=1.5, alpha=0.6) +
#   # scale_size_manual(values=c(0.5, 0.8, 1.1, 1.4)) + 
#   theme_bw(base_size=15) + facet_wrap(reward ~ side)

pdf(fname, width = 6, height = 6)
ggplot(df, aes_string( x="evt_time", y="estimate", ymin="estimate - std.error", shape = "reward", lty = "reward", ymax="estimate + std.error", color="vm_gradient")) + 
  geom_line(size=1, position=position_dodge(width=0.4)) + 
  geom_pointrange(position=position_dodge(width=0.4)) +
  scale_color_brewer(palette="Dark2", labels=c("1" = "MT+, control", "2" = "Parieto-occipital", "3" = "Post. parietal", "4" = "Frontal")) +
  geom_hline(yintercept = 0, size=1.5, alpha=0.6) +
  geom_vline(xintercept = 0, size=1.5, alpha=0.6) +
  # scale_size_manual(values=c(0.5, 0.8, 1.1, 1.4)) + 
  theme_bw(base_size=15) + facet_wrap(as.factor(abs_pe) ~ side)
dev.off()
