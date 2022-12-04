library(pracma)
library(tidyverse)
library(entropy)
library(parallel)

source("~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R")
trial_df <- get_trial_data(repo_directory = "~/Data_Analysis/clock_analysis", dataset="mmclock_fmri", groupfixed = TRUE) %>%
  select(id, run, trial, rt_csv, rew_om) %>%
  rename(reward=rew_om) %>%
  mutate(rt_csv = rt_csv*10) # deciseconds for alignment with basis

#calculate entropy of the basis weights (as in the Cognition paper)
basis <- read.csv("~/Data_Analysis/clock_analysis/fmri/data/mmclock_fmri_sceptic_decay_fits/mmclock_fmri_decay_factorize_selective_psequate_fixedparams_ffx_sceptic_basis.csv") %>%
  as.matrix()

# port matlab
gaussmf <- function(x=NULL, sigma, mu) {
  exp(-(x - mu)^2/(2*sigma^2))
}

# args to setup_rbf
nbasis = 24
trial_length = 4000
ntimesteps = 40
#if contains(so.model, 'psequate'), so.max_prop_spread = -1; end %equate SD to underlying basis
tvec <- 1:ntimesteps
sig_spread <- 1.3945 # ps equate
refspread <- 3.4954

get_past_rts <- function(rt_vec, rewards, rewom=NULL, n=5) {
  if (!is.null(rewom)) rt_vec <- rt_vec[rewards==rewom] # subset to relevant rewards or omissions
  nrt <- length(rt_vec)
  rt_vec[max(1, nrt-n+1):nrt] # don't allow negative index
}

# basis is basis x timesteps (24 x 40)
basis_mult <- function(times, tvec = 1:40, basis, sig_spread = 1.3945) {
  tmult <- sapply(times, function(tt) {
    elig <- gaussmf(tvec, sigma = sig_spread, mu = tt)
    elig <- elig/sum(elig) # normalize by AUC
    e <- rowSums(repmat(elig, nbasis, 1)*basis)
  })
  
  # sum each eligibility function for individual times to get combined representation
  return(rowSums(tmult))
}

# wmat is trials x basis elements and reflects trialwise weights
calc_entropy <- function(wmat) {
  apply(wmat, 1, function(basis_weights) {
    w_norm <- basis_weights/sum(basis_weights)
    nz <- w_norm[w_norm > 0]
    entropy <- -sum(nz * log10(nz))
    return(entropy)
  })
}

# split dataset by id (and eventually run for real data)
gdf <- trial_df %>%
  group_by(id, run)

gkeys <- group_keys(gdf)
glist <- group_split(gdf, .keep = FALSE)

elist <- mclapply(glist, function(ss) {
  
  ntrials <- diff(range(ss$trial)) + 1
  choice_v <- matrix(NA, nrow=ntrials, ncol=nbasis)
  reward_v <- matrix(NA, nrow=ntrials, ncol=nbasis)
  omission_v <- matrix(NA, nrow=ntrials, ncol=nbasis)
  rewom_entropy <- rep(NA, ntrials)
  
  for (t in 2:ntrials) {
    # only feed responses and outcomes up to t-1 since we want entropy representation on trial t prior to making the choice
    past_choices <- get_past_rts(ss$rt_csv[1:t-1], ss$reward[1:t-1], rewom=NULL, n=4)
    choice_v[t,] <- basis_mult(past_choices, basis = basis)
    
    past_rewards <- get_past_rts(ss$rt_csv[1:t-1], ss$reward[1:t-1], rewom=1, n=4)
    reward_v[t,] <- basis_mult(past_rewards, basis = basis)
    
    past_omissions <- get_past_rts(ss$rt_csv[1:t-1], ss$reward[1:t-1], rewom=0, n=4)
    omission_v[t,] <- basis_mult(past_omissions, basis = basis)
    
    rewom_entropy[t] <- entropy.empirical(ss$reward[max(1, t-4):(t-1)])
  }
  
  choice_entropy <- calc_entropy(choice_v)
  reward_entropy <- calc_entropy(reward_v)
  omission_entropy <- calc_entropy(omission_v)
  
  # return(fmri.pipeline:::named_list(choice_v, reward_v, omission_v, 
  #                                   choice_entropy, reward_entropy, omission_entropy, rewom_entropy))
  
  return(data.frame(trial=ss$trial, choice_entropy, reward_entropy, omission_entropy, rewom_entropy))
}, mc.cores = 6)

# add list column
gkeys$edata <- elist

edata_expand <- gkeys %>% unnest(cols = edata)
saveRDS(edata_expand, file="~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/wm_entropy/mmclock_wm_entropy.rds")


# matplot(t(choice_v), type = "l")

# choice_entropy <- calc_entropy(choice_v)
# reward_entropy <- calc_entropy(reward_v)
# omission_entropy <- calc_entropy(omission_v)

# summary(lm(choice_entropy~reward_entropy + omission_entropy))

# add binary reward/omission entropy to set of scalars that are retained -- max at 50/50, min at 100/0

# plot(choice_v[2,], type="l") # trial 2 -- just one trial in the history
# plot(choice_v[3,], type="l") # trial 3
# plot(choice_v[5,], type="l") # trial 5



# test dataset
# rt_df <- tibble::tribble(
#   ~rt_csv, ~trial, ~reward, ~magnitude,
#   1.4, 1, 1, 100,
#   1.8, 2, 0, 0,
#   2.5, 3, 1, 50,
#   3.1, 4, 1, 20,
#   1.2, 5, 0, 0,
#   1.1, 6, 0, 0,
#   3.9, 7, 1, 18,
#   2.0, 8, 1, 55,
#   2.2, 9, 1, 60,
#   2.4, 10, 0, 0,
#   2.7, 11, 1, 90
# )
# rt_df$id <- 1
# rt_df$rt_csv <- rt_df$rt_csv*10 # go to deciseconds for consistency with standard sceptic



# spm hrf
#y = spm_hrf(0.2);
#x = 0:0.2:32
# dd <- R.matlab::readMat("/Users/hallquist/Documents/MATLAB/spm12/spm0p2.mat")
# 
# y <- spm_hrf(0.2)$hrf
# x <- seq(0, length.out=length(y), by=0.2)
# 
# plot(x,y, type="l", col="blue")
# lines(dd$x,dd$y, type="l", col="orange")
# summary(y - dd$y)
# 
# 

# x <- seq(1, 100, 0.2)
# plot(x, gaussmf(x, 4, 50))
