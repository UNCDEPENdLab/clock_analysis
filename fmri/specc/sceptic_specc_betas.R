#sceptic specc betas
library(readr)
library(dplyr)
library(ggplot2)
#group_dir <- "/Volumes/mnh5174_collab/Michael/specc_betas"
group_dir <- "~/Box/SCEPTIC_fMRI/specc_betas"
setwd(group_dir)

beta_files <- system(paste0("find ", group_dir, " -mindepth 2 -iname '*subj_betas.csv' -type f"), intern=TRUE)
beta_files <- grep("(clock_subj|feedback_subj|sceptic-clock-feedback-v_chosen-v_entropy-pe_max-d_auc-preconvolve_fse_groupfixed)", beta_files, value = TRUE, invert = TRUE)
# list.files(path = group_dir, pattern=".*betas.*.csv.*", full.names=TRUE, recursive=TRUE)

#model <- "Intercept-Age-BPD"
#model <- "Intercept-Age_c-BPD_c-BPD_Age"
#model <- "Intercept-Age_c-BPD_c-BPD_Age"
#"Intercept", 
# models <- c("Intercept-Age",
#             "Intercept-Age_c-BPD_c-BPD_Age",
#             "Intercept-Age-BPD",
#             "Intercept-I_Age_c-BPD_c-BPD_IAge")

#models <- "Intercept-I_Age_c-BPD_c-BPD_IAge"
models <- "Intercept-Age_c-BPD_c-BPD_Age"

for (model in models) {
  dir.create(file.path("figures", model), showWarnings = FALSE)
  
  for (b in beta_files) {
    betas <- read.csv(b) %>% 
      filter(model==!!model & l3_contrast != "Intercept" & !l2_contrast %in% c("scram", "fear", "happy")) %>%
      mutate(l3_contrast=recode(l3_contrast, BPD_c="BPD", Age_c="Age", I_Age_c="IAge"))
    
    metadata <- read.csv(sub("subj_betas", "cluster_metadata", b, fixed=TRUE)) %>% 
      filter(model==!!model & l3_contrast != "Intercept" & !l2_contrast %in% c("scram", "fear", "happy")) %>%
      mutate(l3_contrast=recode(l3_contrast, BPD_c="BPD", Age_c="Age", I_Age_c="IAge"))
    
    if (all(is.na(betas$Age))) { betas <- betas %>% select(-Age) %>% rename(Age="Age_c") }
    if (all(is.na(betas$BPD))) { betas <- betas %>% select(-BPD) %>% rename(BPD="BPD_c") }
    if (!is.null(betas$I_Age_c)) { betas <- betas %>% rename(IAge=I_Age_c) }
    
    bpd_effects <- metadata %>% filter(l3_contrast=="BPD")
    if (nrow(bpd_effects) > 0) { 
      pdf(paste0("figures/", model, "/", bpd_effects$l1_contrast[1L], "_BPD_effects.pdf"), width=6, height=6)
      for (cl in 1:nrow(bpd_effects)) {
        title <- with(bpd_effects[cl,], paste0(bpd_effects$l2_contrast[cl], "\n", sub("Focus point:  &nbsp;", "", label, fixed=TRUE), "\n", cluster_size, " vox, z > ", z_threshold, " x:", x, ", y:", y, ", z:", z))
        this_betaset <- bpd_effects %>% slice(cl) %>% left_join(betas) %>% #dplyr detects the merge correctly
          mutate(BPD=factor(if_else(BPD > 0, "BPD", "HC")))
        g <- ggplot(this_betaset, aes(x=BPD, y=cope_value)) + stat_summary(fun.data=mean_cl_boot, geom="pointrange") +
          ggtitle(title) + theme_bw(base_size=12)
        plot(g)
      }
      dev.off()
    }
    
    age_effects <- metadata %>% filter(l3_contrast=="Age")
    if (nrow(age_effects) > 0) { 
      pdf(paste0("figures/", model, "/", age_effects$l1_contrast[1L], "_age_effects.pdf"), width=6, height=6)
      for (cl in 1:nrow(age_effects)) {
        title <- with(age_effects[cl,], paste0(age_effects$l2_contrast[cl], "\n", sub("Focus point:  &nbsp;", "", label, fixed=TRUE), "\n", cluster_size, " vox, z > ", z_threshold, " x:", x, ", y:", y, ", z:", z))
        this_betaset <- age_effects %>% slice(cl) %>% left_join(betas)
        g <- ggplot(this_betaset, aes(x=Age, y=cope_value)) + geom_point() + stat_smooth(method=lm) +
          ggtitle(title) + theme_bw(base_size=12)
        plot(g)
      }
      dev.off()
    }
    
    int_effects <- metadata %>% filter(l3_contrast=="BPD_Age")
    if (nrow(int_effects) > 0) { 
      pdf(paste0("figures/", model, "/", int_effects$l1_contrast[1L], "_agexbpd_effects.pdf"), width=6, height=6)
      for (cl in 1:nrow(int_effects)) {
        title <- with(int_effects[cl,], paste0(int_effects$l2_contrast[cl], "\n", sub("Focus point:  &nbsp;", "", label, fixed=TRUE), "\n", cluster_size, " vox, z > ", z_threshold, " x:", x, ", y:", y, ", z:", z))
        this_betaset <- int_effects %>% slice(cl) %>% left_join(betas) %>%
          mutate(BPD=factor(if_else(BPD > 0, "BPD", "HC")))
        g <- ggplot(this_betaset, aes(x=Age, y=cope_value, color=BPD)) + geom_point() + stat_smooth(method=lm) +
          ggtitle(title) + theme_bw(base_size=12)
        plot(g)
      }
      dev.off()
    }
    
    #inverse age x BPD
    int_effects <- metadata %>% filter(l3_contrast=="BPD_IAge")
    if (nrow(int_effects) > 0) {
      pdf(paste0("figures/", model, "/", int_effects$l1_contrast[1L], "_iagexbpd_effects.pdf"), width=6, height=6)
      for (cl in 1:nrow(int_effects)) {
        title <- with(int_effects[cl,], paste0(int_effects$l2_contrast[cl], "\n", sub("Focus point:  &nbsp;", "", label, fixed=TRUE), "\n", cluster_size, " vox, z > ", z_threshold, " x:", x, ", y:", y, ", z:", z))
        this_betaset <- int_effects %>% slice(cl) %>% left_join(betas) %>%
          mutate(BPD=factor(if_else(BPD > 0, "BPD", "HC")))
        g <- ggplot(this_betaset, aes(x=100/IAge, y=cope_value, color=BPD)) + geom_point() + stat_smooth(method=lm) +
          ggtitle(title) + theme_bw(base_size=12)
        plot(g)
      }
      dev.off()
    }
    
  }
}
