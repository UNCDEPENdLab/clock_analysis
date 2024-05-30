library(tidyverse)
trial_df <- trial_df %>% select(dataset, id, run, trial, rewFunc, rt_csv, magnitude, probability, ev, rt_vmax, score_csv)
save(trial_df, file = "~/OneDrive/collected_letters/papers/meg/sci_adv/zenodo_sci_adv_data/trial_data_compact.RData")
rt_comb <- rt_comb %>% select(id, run, run_trial, feedback_onset, rewFunc, atlas_value, label, decon_interp, side)
save(rt_comb, file = "~/OneDrive/collected_letters/papers/meg/sci_adv/zenodo_sci_adv_data/fig_3/rt_aligned_deconvolved_bold.RData")
library(data.table)
# fig. 4
df <- fread("/Users/alexdombrovski/Library/CloudStorage/OneDrive-Personal/collected_letters/papers/meg/sci_adv/zenodo_sci_adv_data/fig_4/Schaefer_DorsAttn_2.3mm_cope_l2.csv.gz") %>%
  select(id, l2_cope_name, l1_cope_name, mask_value, x, y, z, value) %>% filter(l2_cope_name == "overall" & l1_cope_name == "EV_entropy_change_feedback") %>%
  select(-l2_cope_name, -l1_cope_name)

fwrite(df, file = "/Users/alexdombrovski/Library/CloudStorage/OneDrive-Personal/collected_letters/papers/meg/sci_adv/zenodo_sci_adv_data/fig_4/entropy_change_betas.csv.gz")
# fig. 5
df <- readRDS("/Users/alexdombrovski/Library/CloudStorage/OneDrive-Personal/collected_letters/papers/meg/sci_adv/dryad_sci_adv_supplemental_data/fig_5/meg_ddf_wholebrain_entropy_change_ri.rds") %>% 
  filter(term == "entropy_change_t" & effect == "fixed") %>%
  select(Time, Freq, estimate, std.error, statistic, df, p.value, p_fdr)
fwrite(df, file = "/Users/alexdombrovski/Library/CloudStorage/OneDrive-Personal/collected_letters/papers/meg/sci_adv/dryad_sci_adv_supplemental_data/fig_5/meg_time_frequency_entropy_change_ri.rds")
