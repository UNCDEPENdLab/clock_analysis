library(cowplot)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

source("RewFunction_Reversed.R") #will create a cached RT - freq translation object called 'lookup' used in pmap_dfr

#just used for global diagnosis at the moment
df <- read_csv("~/Data_Analysis/clock_analysis/meg/data/mmclock_meg_decay_mfx_sceptic_global_statistics.csv")

ggplot(df, aes(x=beta, y=R2)) + geom_point() + stat_smooth()

#these are subjects who had all 0 RTs for unknown reasons
filter(df, R2 > .9) %>% pull(id)

#this is the meg file that generates (usually) slight differences between rt_csv and rt_backcalc
tdf <- read_csv("~/Data_Analysis/clock_analysis/meg/data/mmclock_meg_decay_mfx_trial_stats.csv.gz")

### PROBLEM: some rows are missing key information such as rewFunc (contingency)

miss_stuff <- tdf %>% filter(is.na(rewFunc)) %>% write_csv(path="meg_missing_info.csv")


%>%
  bind_cols(pmap_dfr(list(.$probability, .$rewFunc), RewFunction_Reversed, L=lookup)) #add on back-calculated RTs for all subjects

#for comparison, we don't see any discrepancies in MRI. back calculation is correct down to the millisecond
#tdf <- read_csv("~/Data_Analysis/clock_analysis/fmri/data/mmclock_fmri_decay_mfx_trial_stats.csv.gz")

zero_rts <- tdf %>% filter(rt_csv==0) %>% select(id, run, trial, rewFunc, emotion, rt_csv, score_csv,
                                                magnitude, probability, ev)

#cache for visual examination
write_csv(x=zero_rts, path="zero_rts.csv")

#look at subjects missing some, but not all, rts. These seem to merit further consideration, as this isn't 
#attributable to a bug that dropped all RTs.
zero_rts_some <- zero_rts %>%
  filter(!id %in% c("11229_20140212", "11280_20140926", "11321_20140904", "11322_20140918",
                    "11325_20140918", "11331_20141120", "11336_20141121", "11342_20150228", 
                    "11344_20141118", "11346_20141230")) %>%
  bind_cols(rt_backcalc=pmap_dfr(list(.$probability, .$rewFunc), RewFunction_Reversed, L=lookup))

tdf %>% filter(rt_csv > 4000) %>% nrow()

write_csv(x=zero_rts_some, path="zero_rts_not_all.csv")


#this subject has intact RTs,
#filter(id=="11343_20150110")
one_subj <- tdf %>% filter(id=="11262_20140312") %>% mutate(rewFunc=as.character(rewFunc)) %>%
  bind_cols(pmap_dfr(list(.$probability, .$rewFunc), RewFunction_Reversed, L=lookup))

one_subj <- tdf %>% filter(id==10637) %>% mutate(rewFunc=as.character(rewFunc)) %>%
  bind_cols(pmap_dfr(list(.$probability, .$rewFunc), RewFunction_Reversed, L=lookup))


tocorr <- one_subj %>% select(trial, rewFunc, rt_csv, rt_backcalc, rt_vba, probability, magnitude) %>% mutate(rt_vba=rt_vba*100) %>%
  mutate(rt_backcalc=if_else(rt_backcalc>=4000, NA_integer_, rt_backcalc))

xx <- split(tocorr, tocorr$rewFunc)

lapply(xx, function(df) { cor(df$rt_csv, df$rt_backcalc, use="pairwise.complete.obs")}) 

sum(is.na(tocorr$rt_backcalc))

tocorr %>% filter(is.na(rt_backcalc))

toplot <- tocorr %>% gather(key="rtsrc", value="rt", rt_csv, rt_backcalc) %>% select(-rt_vba) #, rt_vba

tocorr <- tocorr %>% mutate(rt_discrep=rt_csv - rt_backcalc) #%>% filter(abs(rt_csv -rt_backcalc) > 100)

filter(tocorr, abs(rt_discrep) > 100) %>% View()

ggplot(toplot, aes(x=trial, y=rt, color=rtsrc)) + geom_line() + xlim(c(1,150)) #facet_wrap(~rtsrc, ncol=1) +


### PROBLEM:   mismatch between rt_csv and rt_vba that appears to be related to trial numbering
### DIAGNOSIS: this occurs when a subject has a dropped run. The MATLAB export code numbers observed trials in
###            ascending order, which will not match the CSV containing correct trial numbering with gaps
### SOLUTION:  Add a ascend_trial column to the trial_stats data.frame and bind on that as the key before writing trial_stats.csv.gz files
#11313_20141104
#11320_20140908

#can poke around in a few weird subjects
questionable <- tdf %>% filter(id == "11313_20141104") %>% select(run, trial, rewFunc, emotion, rt_csv, rt_vba)

tdf <- tdf %>% mutate(rt_check=rt_vba*100) %>% filter(abs(rt_check - rt_csv) > 10) %>% droplevels() %>%
  select(id, run, trial, rewFunc, emotion, rt_csv, rt_vba)
xtabs(~id, tdf)

### PROBLEM:  there are small mismatches in MEG between the rt_csv column and the backcalculated RTs. These appear to occur
              in 
### 
#double check the probability -> RT lookup structure to make sure nothing is improper there
lookup_df <- as.data.frame(lookup) %>% gather(key="frqmag", value="value", -RT) %>%
  separate(frqmag, into=c("contingency", "frqmag"), sep="_") %>% spread(key=frqmag, value=value)

frqplot <- ggplot(lookup_df, aes(x=RT, y=frq, color=contingency)) + geom_line()
magplot <- ggplot(lookup_df, aes(x=RT, y=mag, color=contingency)) + geom_line()

#yep, looks right!
plot_grid(frqplot, magplot, align="h")
