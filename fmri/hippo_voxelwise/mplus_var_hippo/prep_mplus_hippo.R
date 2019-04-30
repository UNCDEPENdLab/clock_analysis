library(MplusAutomation)
library(dplyr)
load("/Users/mnh5174/Box/SCEPTIC_fMRI/var/feedback_hipp_wide_ts.Rdata")

setwd(file.path(getMainDir(), "clock_analysis", "fmri", "hippo_voxelwise", "mplus_var_hippo"))
str(fb_wide)
fb_wide <- fb_wide %>% mutate(block=factor(paste0(id, run))) #%>% dplyr::select(-evt_time)

#spread evt_time wide
fb_l <- fb_wide %>% dplyr::select(id, block, run, run_trial, evt_time, ends_with("_l"))

#shorten hippocampus labels to "lh" for left hippocampus (to avoid Mplus 8-character complaints)

names(fb_l) <- sub("hipp_(\\d+)_l", "lh\\1", names(fb_l), perl=TRUE)
fb_l_wide <- fb_l %>% dplyr::select(-block) %>% filter(evt_time >= 0) %>% 
  arrange(id, run_trial) %>%
  dplyr::group_by(id) %>% dplyr::mutate(trial=(run-1)*50+run_trial) %>% dplyr::select(-run_trial, -run) %>% ungroup() %>%
  reshape2::melt(id.vars=c("id", "trial", "evt_time")) %>%
  mutate(evt_time=factor(paste0("t", sprintf("%02d", evt_time)))) %>%
  reshape2::dcast(id + trial ~ variable + evt_time, value.var = "value")

#generate syntax of everything regressed on everything else... in super-wide format (time within trial as wide)
#hipp_slices <- 1:12
#hipp_slices <- 1:3
hipp_slices <- 1:3
#times <- 1:10
times <- 0:4
incl_contemp <- TRUE

m_string <- c()

#setup latent factors for each measurement -- these are used to capture within-trial temporal dependence through regressions
m_string <- c(m_string, "%WITHIN%", "! indicator factor definitions")
for (i in hipp_slices) {
  for (j in times) {
    m_string <- c(m_string,
      paste0("fh", hipp_slices[i], "_t", sprintf("%02d", j), " BY lh", hipp_slices[i], "_t", sprintf("%02d", j), "@1;"),
      paste0("lh", hipp_slices[i], "_t", sprintf("%02d", j), "@0; !no resid var"),
      ""
    )
  }
  m_string <- c(m_string, "") #empty line
}

#Setup equally weighted factors to capture trial-level dependence
m_string <- c(m_string, "! between-trial correlation using unit-weighted factor with lag")
for (i in hipp_slices) {
  # m_string <- c(m_string, paste0("fh", hipp_slices[i], "trial BY fh", hipp_slices[i], "_t", sprintf("%02d", min(times)),
  #                                "-fh", hipp_slices[i], "_t", sprintf("%02d", max(times)), "@1 (&1);")
  
  #the observed variables are indicators of the trial, not, the latent indicators (parcellation)
  m_string <- c(m_string, paste0("fh", hipp_slices[i], "trial BY lh", hipp_slices[i], "_t", sprintf("%02d", min(times)),
                                 "-lh", hipp_slices[i], "_t", sprintf("%02d", max(times)), "@1 (&1);")
  )
}
m_string <- c(m_string, "") #empty line

#setup between-trial AR
m_string <- c(m_string, "! setup trial AR")
for (i in hipp_slices) {
  m_string <- c(m_string, paste0("fh", hipp_slices[i], "trial ON fh", hipp_slices[i], "trial&1;"))
}
m_string <- c(m_string, "") #empty line

for (i in 1:length(hipp_slices)) {
  for (j in times) {

    if (j == min(times)) { next }
    
    #setup regressions of every measure on the preceding measures
    m_string <- c(m_string, paste("! slice", hipp_slices[i], "time", sprintf("%02d", j), "regressed on time", sprintf("%02d", j-1)))
    for (k in hipp_slices) {
      #lag 1
      m_string <- c(
        m_string, paste0("fh", hipp_slices[i], "_t", sprintf("%02d", j), " ON ", paste0("fh", k, "_t", sprintf("%02d", j-1), " (a", hipp_slices[i], k, ")", collapse=" "), ";")
      )
    }
    m_string <- c(m_string, "")
    
    #do not allow within-trial and between-trial factors to correlate (undermines variance partition)
    if (j==max(times)) {
      m_string <- c(m_string, paste0("fh", hipp_slices[i], "trial WITH fh", hipp_slices, "_t", sprintf("%02d", j), "@0; !no corr between within trial and bw trial factors"), "")
    }
    
    #lag 1
    # m_string <- c(
    #   m_string, paste0("lh", hipp_slices[i], "_t", sprintf("%02d", j), " ON ", paste0("lh", hipp_slices, "_t", sprintf("%02d", j-1), " (a", hipp_slices[i], ")\n", collapse=" "), ";")
    # )
    
  }   
}
m_string <- c(m_string, "")

#contemporaneous correlations (structural VAR)
if (incl_contemp) {
  pairs <- gtools::combinations(length(hipp_slices), 2, hipp_slices)
  for (j in times) {
    m_string <- c(m_string, paste("! time", sprintf("%02d", j), "contemp correlations"))
    for (p in 1:nrow(pairs)) {
      m_string <- c(
        #m_string, paste0("fh", pairs[p,1], "_t", sprintf("%02d", j), " WITH ", paste0("fh", pairs[p,2], "_t", sprintf("%02d", j), " (c", pairs[p,1], pairs[p,2], ")", collapse=" "), ";")
        m_string, paste0("fh", pairs[p,1], "_t", sprintf("%02d", j), " WITH ", paste0("fh", pairs[p,2], "_t", sprintf("%02d", j), ";"))
      )
    }
    m_string <- c(m_string, "")
  }
  # m_string <- c(m_string, paste("! slice", hipp_slices[i], "time", sprintf("%02d", j), "contemp regressions on", sprintf("%02d", j)))
  # for (k in hipp_slices) {
  #   if (k == i) { next } #don't regress on oneself
  #   m_string <- c(
  #     m_string, paste0("fh", hipp_slices[i], "_t", sprintf("%02d", j), " ON ", paste0("fh", k, "_t", sprintf("%02d", j), " (c", hipp_slices[i], k, ")", collapse=" "), ";")
  #   )
  # }
  m_string <- c(m_string, "")
  
  # m_string <- c("! contemporaneous",
  #   m_string, paste0("fh", hipp_slices[i], "_t", sprintf("%02d", j), " ON ", paste0("fh", hipp_slices[-i], "_t", sprintf("%02d", j), collapse=" "), ";"),
  #   ""
  # )
}

#not really sure we need a between-subs model, but doesn't hurt...
m_string <- c(m_string, "%BETWEEN%", "! BW subjects correlations")
for (i in hipp_slices) {
  m_string <- c(m_string, paste0("bh", hipp_slices[i], " BY lh", hipp_slices[i], "_t", sprintf("%02d", min(times)),
                                 "-lh", hipp_slices[i], "_t", sprintf("%02d", max(times)), "@1;")
  )
}
m_string <- c(m_string, "") #empty line

m_string <- c(m_string, "! zero means")
for (i in hipp_slices) {
  m_string <- c(
    m_string, 
    #paste0("[lh", hipp_slices[i], "_t", min(times), "-lh", hipp_slices[i], "_t", max(times), "] (ic", hipp_slices[i], ");") #equal intercepts (blows up)
    paste0("[lh", hipp_slices[i], "_t", sprintf("%02d", min(times)), "-lh", hipp_slices[i], "_t", sprintf("%02d", max(times)), "@0];")
  )
}
m_string <- c(m_string, "") #empty line


cat(m_string, file="test_syntax.txt", sep="\n")

#scale up variables to make parameter estimates easier to see
#prepareMplusData(fb_l_wide %>% mutate_at(vars(starts_with("lhipp")), list(~.*10)), "fb_l_superwide.dat")
prepareMplusData(fb_l_wide, "fb_l_superwide.dat")

######
#plot demo output
prototype_out <- readModels("mega_wide.out")
pars <- prototype_out$parameters$unstandardized

library(stringr)
lag1 <- pars %>% filter(str_detect(paramHeader, "^FH\\d+_T01.ON") & str_detect(param, "^FH\\d+_T00"))
#matrix of lagged relationships
mat <- matrix(NA, nrow=3, ncol=3)
for (i in 1:nrow(mat)) {
  for (j in 1:ncol(mat)) {
    #column is 'from', row is 'to'
    mat[i,j] <- lag1 %>% filter(str_detect(paramHeader, paste0("^FH", i, "_T02.ON")) & str_detect(param, paste0("^FH", j, "_T01"))) %>% pull(est)
  }
}

library(qgraph)
lag1_graph <- qgraph(mat, layout="spring", edge.labels=TRUE, 
                            nonsig = "show", vsize=5, esize=2, asize=2, edge.label.cex=0.5, filetype="jpeg", filename="lag1_test",
                            fade=FALSE, directed=TRUE, mar=c(8,8,8,8))

#mlvar help!
# doesn't want to actually write mplus syntax
# Model <- mlVARsim(nPerson = 50, nNode = 3, nTime = 50, lag=1)
# fit1 <- mlVAR(Model$Data, vars = Model$vars, idvar = Model$idvar, lags = 1, 
#               temporal = "orthogonal", estimator="Mplus", verbose=TRUE)
# 
# fit1 <- mlVAR:::Mplus_mlVAR(Model$Data, vars = Model$vars, idvar = Model$idvar, lags = 1, 
#                estimator="Mplus", verbose=TRUE)
