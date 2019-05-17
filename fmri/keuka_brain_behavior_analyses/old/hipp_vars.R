library(vars)
load('~/Box Sync/SCEPTIC_fMRI/var/feedback_hipp_wide_ts.Rdata')
slices = 12
if (slices == 12) {df <- fb_wide; grp <- list(1:3,4:12)}
if (slices == 6) {df <- fb_wide6; grp <- list(1:2,3:6)}
rsort <- mixedsort(names(df[grep('_r', names(df))]))
lsort <- mixedsort(names(df[grep('_l', names(df))]))
plot.ts(df$hipp_1_l[1:10000])
plot.ts(df$hipp_1_r[1:10000])

# does not run with NAs
# adf1 <- summary(ur.df(df[,"hipp_1_l"], type = "trend", lags = 2))

VARselect(df[,5:28], lag.max = 5, type = "both")
