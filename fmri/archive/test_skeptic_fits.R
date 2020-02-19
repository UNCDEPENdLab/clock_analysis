#analyze SKEpTIC fits to behavior

library(R.matlab)
setwd(file.path(getMainDir(), "temporal_instrumental_agent", "clock_task"))

#behav <- readMat("behav.mat")
options(matlab="/Applications/MATLAB_R2015a.app/bin/matlab")

#Matlab$startServer()
#matlab <- Matlab()
#isOpen <- open(matlab)
#if (!isOpen) { throw("MATLAB server is not running: waited 30 seconds.") }
#print(matlab)
#
#evaluate(matlab, "cd '/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task'; \
#        load behav.mat; \
#				myvar=table2struct(behav{1,1}.data);
#")
#
#alpha <- getVariable(matlab, "myvar")[[1L]]
#dim(alpha) <- dim(alpha)[c(1,2)]
#
#df <- c()
#for (i in 1:dim(alpha)[1]) {
#  df <- cbind(df, data.frame(t(alpha[i,])))
#}
#
#library(plyr)
#result <- ldply(alpha, function(col) {
#      unlist(col)
#    })
#
#
#result <- adply(alpha, 2, function(col) {
#      data.frame(unlist(col))
#    })
#
#
#
#df <- do.call(data.frame,alpha)
#
#
#close(matlab)

#read csv files
csvfiles <- list.files("data", pattern=".*\\.csv", full.names=TRUE)
csvfiles <- csvfiles[!csvfiles=="data/11332.csv"] #only 300 trials screwing up split below
df <- c()
for (i in 1:length(csvfiles)) {
  thisguy <- read.csv(csvfiles[i], header=TRUE)
  thisid <- sub("^[^\\d]*(\\d+)[^\\d]*", "\\1", csvfiles[i], perl=TRUE)
  #if (nrow(thisguy) != 400) { browser()}
  thisguy$id <- factor(thisid)
  thisguy$full_rtpred_exploit <- thisguy$full_rtpred_exploit * 10 #convert to ms
  thisguy$full_rtpred_explore <- thisguy$full_rtpred_explore * 10 #convert to ms
  thisguy$limited_rtpred_exploit <- thisguy$limited_rtpred_exploit * 10 #convert to ms
  thisguy$limited_rtpred_explore <- thisguy$limited_rtpred_explore * 10 #convert to ms  
  df <- rbind(df, thisguy)
}

df_learnable <- gdata::drop.levels(subset(df, rewFunc %in% c("IEV", "DEV")))
#df_learnable <- df #above getting screwed up by the fact that CEV blocks fall in different places, so not fully crossed
byid <- split(df_learnable, list(df_learnable$id, df_learnable$run))

#split generates empty combinations due to CEV and CEVR falling at different run numbers depending on counterbalance order
invalid_combinations <- which(sapply(byid, function(x) { nrow(x) == 0}))
byid[invalid_combinations] <- NULL

options(error=recover)
allfits <- lapply(byid, function(iddf) {      
      m <- lm(rt ~ full_rtpred_explore + full_rtpred_exploit, iddf)
      sm <- summary(m)
      m_limited <- lm(rt ~ limited_rtpred_explore + limited_rtpred_exploit, iddf)
      sm_limited <- summary(m_limited)
      
      ret <- list(params_full=sm$coefficients, r2_full=sm$r.squared, f_full=sm$fstatistic,
          params_limited=sm_limited$coefficients, r2_limited=sm_limited$r.squared, f_limited=sm_limited$fstatistic)
      return(ret)
    })

r2_full <- sapply(allfits, "[[", "r2_full")
r2_limited <- sapply(allfits, "[[", "r2_limited")

summary(r2_full)
summary(r2_limited)

library(lattice)
histogram(r2_full)

b_exploit_full <- sapply(allfits, function(x) { x$params_full["full_rtpred_exploit", "Estimate"] })
b_exploit_limited <- sapply(allfits, function(x) { x$params_limited["limited_rtpred_exploit", "Estimate"] })

b_exploit_limited <- sapply(allfits, function(x) { nrow(x$params_limited) })
histogram(b_exploit_full)

high_corrs <- sapply(byid, function(x) { cor(x$limited_rtpred_explore, x$limited_rtpred_exploit)})
crazy_exploit <- high_corrs[which(high_corrs == 1)]

crazy_exploit <- which(b_exploit_full > 5)

pdf("figures/extreme_exploitfits.pdf", width=11, height=8)
for (i in 1:length(crazy_exploit)) {
  badguy <- byid[[ crazy_exploit[i] ]] #"11331.3"
  
  library(ggplot2)
  library(reshape2)
  m <- melt(badguy[,c("trial", "rt", "full_rtpred_explore", "full_rtpred_exploit")], id.vars="trial")
  
  scores <- badguy[,c("trial", "rt", "score")]
  scores <- plyr::rename(scores, c(rt="value"))
  
  g <- ggplot(m, aes(x=trial, y=value, color=variable)) + geom_line() + theme_bw(base_size=16) + 
      geom_text(data=scores, aes(x=trial, y=value + 100, label=score, color=NULL)) +
      ggtitle(paste(badguy$id[1], badguy$rewFunc[1]))
  
  plot(g)
}
dev.off()








b_explore_full <- sapply(allfits, function(x) { x$params_full["full_rtpred_explore", "Estimate"] })
b_explore_limited <- sapply(allfits, function(x) { x$params_limited["limited_rtpred_explore", "Estimate"] })


