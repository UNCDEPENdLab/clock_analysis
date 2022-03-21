# merge in the DAN betas into the MEG trial_df

library(tidyverse)
library(psych)
library(corrplot)
library(data.table)
repo_directory <- "~/code/clock_analysis"
behavioral_data_file <- "~/code/clock_analysis/meg/MEG_n63_behavioral_data_preprocessed_trial_df.RDS"

trial_df <- readRDS(behavioral_data_file)

# read in the betas
labels <- as_tibble(read_table2("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/Schaefer2018_200Parcels_DAN_order_manual.txt", col_names = F)) %>% 
  mutate(side = case_when(
    str_detect(X2, "LH") ~ "L",
    str_detect(X2, "RH") ~ "R"
  ),
  label = paste0(X3, "_", side)) %>% select(c(X1, X3, side)) %>% rename(mask_value = X1, label_sym = X3)
betas <- read_csv("/Volumes/GoogleDrive/.shortcut-targets-by-id/1ukjK6kTlaR-LXIqX6nylYOPWu1j3XGyF/SCEPTIC_fMRI/final_betas/L1m-echange_abspe/Schaefer_DorsAttn_2.3mm_cope_l2.csv.gz")  %>%
  mutate(id  = as.integer(id)) %>% select(id, l1_model, l1_cope_name, l2_cope_name, mask_value, value) %>% filter(l1_cope_name == "EV_entropy_change_feedback", l2_cope_name == "overall")  %>%
  merge(labels) %>% as_tibble() %>% mutate(value = winsor(value, trim = .005))
wbetas <- betas %>% select(id, value, label_sym, side) %>% group_by(id, label_sym) %>% summarize(beta_bl = mean(value)) %>% 
  pivot_wider(names_from = label_sym, values_from = beta_bl) %>% ungroup() %>%
  rowwise() %>%
  mutate(ppc_ec_beta = mean(c(`2_ip_LIPd`,`2_ip_LIPv`,  `2_ip_VIP`, `3_sp_7AM`, `3_sp_7PC`)),
         mt_ec_beta = mean(c(`4_MT/V5_FST`, `4_MT/V5_MST`)),
         pfc_ec_beta = mean(c(`1_f_6a`, `1_f_FEF`, `1a_f_BA44`)))
just_betas <- wbetas %>% select(is.numeric) %>% select(-id)
cormat <- corr.test(just_betas, method = "pearson")
corrplot(cormat$r[1:13,1:13], order = "hclust", tl.pos = "n", cl.lim = c(0.5, 1), is.corr = T)

mec <- nfactors(cormat$r, n=5, rotate = "oblimin", diagonal = FALSE,fm = "pa", n.obs = 70, SMC = T)
ec.fa = psych::fa(just_betas, nfactors=3, rotate = "oblimin", fm = "pa")
# remove MT+
fp <- cormat$r[1:13,1:13]
mec_fp <- nfactors(fp, n=5, rotate = "oblimin", diagonal = FALSE,fm = "ml", n.obs = 70, SMC = T)
fp.fa <- psych::fa(r = fp, nfactors=2, rotate = "oblimin", fm = "ml")

# hbetas <- read_csv("/Volumes/GoogleDrive/.shortcut-targets-by-id/1ukjK6kTlaR-LXIqX6nylYOPWu1j3XGyF/SCEPTIC_fMRI/whole_brain/ebetas/L1m-entropy/Schaefer_DorsAttn_2.3mm_cope_l2.csv.gz")

# save parcel-wise and regional betas
saveRDS(wbetas, file = "~/code/clock_analysis/meg/data/fmri_ec_betas.RDS")

fmri_betas <- inner_join(trial_df, wbetas %>% select(c(ppc_ec_beta, pfc_ec_beta, mt_ec_beta, id)), by = "id")
# saveRDS(trial_df, file = behavioral_data_file)
