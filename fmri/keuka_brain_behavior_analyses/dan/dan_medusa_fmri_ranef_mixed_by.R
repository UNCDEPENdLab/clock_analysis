library(tidyverse)
library(lme4)


medusa_dir = "~/Box/SCEPTIC_fMRI/dan_medusa/"
cache_dir = "~/Box/SCEPTIC_fMRI/dan_medusa/cache"
repo_directory <- "~/code/clock_analysis"
gc()

# load DAN RT-aligned data
# data loading options
reprocess = F # otherwise load data from cache
if (!reprocess) {
  wide_only = F
  tall_only = T# only load wide data (parcels and timepoints as variables)
}
replicate_compression = F
if(replicate_compression) {reprocess = T}

plots = F

# load MEDUSA deconvolved data
source(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/load_medusa_data_dan.R"))

# load behavioral data, get lags, merge with signal
load(file.path(repo_directory, '/fmri/keuka_brain_behavior_analyses/trial_df_and_vh_pe_clusters_u.Rdata'))
df <- df %>% select(id, run, run_trial, emotion, last_outcome, rt_next, pe_max, rt_vmax, score_csv,
                    v_max_wi, v_entropy_wi, v_entropy_b, v_entropy, v_max_b, u_chosen_quantile, u_chosen_quantile_lag, u_chosen_quantile_change, 
                    rt_vmax_lag_sc, rt_lag_sc,rt_lag2_sc, rt_csv_sc, trial_neg_inv_sc, Age, Female, kld3, kld4) %>% 
  group_by(id, run) %>% arrange(id, run, run_trial) %>% 
  mutate(rt_next = lead(rt_csv_sc),
         rt_change = rt_next - rt_csv_sc,
         rt_vmax_lead = lead(rt_vmax),
         rt_vmax_change_next = rt_vmax_lead - rt_vmax,
         v_entropy_wi_lead = lead(v_entropy_wi),
         v_entropy_wi_change = v_entropy_wi_lead-v_entropy_wi,
         v_entropy_wi_change_lag = lag(v_entropy_wi_change),
         u_chosen_quantile_next = lead(u_chosen_quantile),
         u_chosen_quantile_change_next = lead(u_chosen_quantile_change),
         kld3_lag = lag(kld3),
         kld3_lead = lead(kld3),
         outcome = case_when(
           score_csv>0 ~ 'Reward',
           score_csv==0 ~ "Omission"),
         abs_pe = abs(pe_max),
         abs_pe_lag = lag(abs_pe)
  ) %>% ungroup()
d <- merge(df, rt_comb, by = c("id", "run", "run_trial"))

# mixed_by call
splits = c("stream", "side", "evt_time")
encode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag)  + scale(rt_vmax_change) + 
                           v_entropy_wi + v_entropy_wi_change + kld3_lag  + v_max_wi  + scale(abs_pe) + outcome + (outcome + v_entropy_wi|id))
ddf <- mixed_by(d, outcomes = "decon_interp", rhs_model_formulae = encode_formula , split_on = splits,
                ncores = 20, refit_on_nonconvergence = 5,
                tidy_args = "ran_vals")
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_decode')
encode_results_fname = "rt_encode_output_streams_mixed_by_entropy_ranef.RDS"
saveRDS(file = encode_results_fname, ddf)

if (plots) {
  ## Check plots
  message("\nPlotting streams decoding")
  library(viridis)
  setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_decode')
  epoch_label = "Time relative to outcome, seconds"
  encode_results_fname = "rt_encode_output_streams_mixed_by_entropy_change_ranef.RDS"
  ddf <- as_tibble(ddf)
  ddf$t <- ddf$evt_time
  ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term, 
                         p_level_fdr = as.factor(case_when(
                           # p_fdr > .1 ~ '0',
                           # p_fdr < .1 & p_fdr > .05 ~ '1',
                           p_fdr > .05 ~ '1',
                           p_fdr < .05 & p_fdr > .01 ~ '2',
                           p_fdr < .01 & p_fdr > .001 ~ '3',
                           p_fdr <.001 ~ '4')))
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  ddf$`p, FDR-corrected` = ddf$p_level_fdr
  if ("visuomotor_grad" %in% splits) {ddf$visuomotor_grad <- factor(ddf$visuomotor_grad, labels=c("1" = "MT+, control", "2" = "Parieto-occipital", "3" = "Post. parietal", "4" = "Frontal"))}
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  for (fe in terms) {
    # fe <- terms[1] # test only
    edf <- ddf %>% filter(term == paste(fe) & t < 8) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    # plot stream gradients
    fname = paste("streams_", termstr, "_mixed_by.pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    print(ggplot(edf, aes(t, stream)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
            # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("") +
            scale_y_discrete(labels=c("visual-motion" = "MT+,\ncontrol", "ventro-dorsal" = "Ventro-dorsal\nstream",
                                      "oculomotor" = "Oculomotor\nstream", "dorso-dorsal" = "Dorso-dorsal\nstream")))
    dev.off()
    
    fname = paste("streams_line_", termstr, "_mixed_by.pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=stream, size=`p, FDR-corrected`)) + 
      geom_line(size = 1) + geom_point() +
      geom_errorbar() +
      # geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~visuomotor_grad) +
      geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
      scale_color_brewer(palette="Set1") + xlab(epoch_label) + 
      labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
    print(gg)
    dev.off()
    
  }
}
saveRDS(file = encode_results_fname, ddf)
