# analyzes behavioral preliminary data from clock reversal task in explore
# first run run_fsl_pipeline_explore.R

df <- vba_output
subject_df$id <- as.character(subject_df$redcapid)
df <- inner_join(df, subject_df)
df <- df %>% mutate(case_when(trial <))
ggplot(df, (trial))