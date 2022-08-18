# import explore subject-level data
library(tidyverse)
library(readr)
library(naniar) # missingness 

data_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/explore_medusa/data/"
setwd(data_dir)

sub_df <- read_csv("explore_clock_indv_diff_8_3.csv")
str(sub_df)


# compute groups:

sub_df <-  sub_df %>% mutate(
  Group = registration_group,
  Group = case_when(
    ideation_group == "low-ide" ~ "Depressed",
    ideation_group == "high-ide" ~ "Ideators",
    T ~ Group),
  Group = case_when(
    Group == "ATT" ~ "Attempters",
    Group == "HC" ~ "Controls",
    T ~ Group
  ),
  Group = factor(Group, levels = c("Controls","Depressed", "Ideators", "Attempters")),
  Group_d = factor(Group, levels = c("Depressed", "Controls", "Ideators", "Attempters")),
  Group_a = factor(Group, levels = c("Attempters", "Controls", "Depressed", "Ideators")),
  suds_dx = as.logical(suds_dx),
  anxiety_dx = as.logical(anxiety_dx),
  group_leth = case_when(
    max_lethality > 3  ~ "HL_Attempters",
    max_lethality < 4  ~ "LL_Attempters",
    T ~ as.character(Group)
  ),
  group_leth = factor(group_leth, levels = c("HL_Attempters", "Controls", "Depressed", "Ideators", "LL_Attempters")),
  id = as.integer(registration_redcapid),
  opioid = as.factor(ifelse(is.na(opioid), 0, opioid)),
  sedhyp = as.factor(ifelse(is.na(sedhyp), 0, sedhyp)),
  antipsychotic = as.factor(ifelse(is.na(antipsychotic), 0, antipsychotic))
)

# check for missingness
gg_miss_var(sub_df %>% select(-registration_lethality, -total_attempts, -max_lethality, -ideation_group, -ipipds_total), show_pct = TRUE,
            facet = registration_group)
gg_miss_var(sub_df %>% select(-registration_lethality, -total_attempts, -max_lethality, -ideation_group, -ipipds_total), show_pct = TRUE)


str(sub_df)

library(compareGroups)

t1 <- compareGroups(Group ~ ., data = sub_df[5:58])
createTable(t1)

t2 <- compareGroups(Group_d ~ ., data = sub_df[c(5:58)] %>% filter(Group!="Controls"))
createTable(t2)

t3 <- compareGroups(group_leth ~ ., data = sub_df[5:61] %>% filter(Group!="Controls"))
createTable(t3)

saveRDS(sub_df, file = "explore_n146.rds")
