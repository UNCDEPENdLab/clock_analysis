library(data.table)
library(tidyverse)
library(lme4)
analysis_dir <- "~/code/clock_analysis/fmri/explore/deepmreye/"
setwd(analysis_dir)
df <- read_csv("2023-05-17-DeepMReye-N=50.csv")[-1] %>% filter(online == T) %>% 
  select(id, trial, run_trial, run, rewFunc, contains("theta"), rt_csv, rt_vmax, TR_num, dtheta, dr, c(124:135)) %>%
  mutate(tr_num_f = as.factor(TR_num)) %>% group_by(id, trial) %>% 
  mutate(trs_until_rt = max(TR_num) - TR_num,
         trs_until_rt_f = as.factor(trs_until_rt))

# something wrong with TR_num, filter big values for now

ggplot(df %>% filter(theta_rt > 1 & theta_rt < 5.5), aes(theta_rt, theta)) + geom_smooth(method = "gam")

ggplot(df %>% filter(theta_rt_vmax > 1 & theta_rt_vmax < 5.5), aes(theta_rt_vmax, theta)) + geom_smooth(method = "gam")

ggplot(df, aes(TR_num, theta, color = rewFunc)) +
  # geom_smooth(method = "loess") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 8, bs = "cs")) +
  coord_polar(theta = "y", clip = "off") +
  scale_y_continuous(limits = c(0, 6.33)) + facet_wrap(~run) +
  theme_bw()

ggplot(df, aes(-trs_until_rt, theta, color = rewFunc)) +
  # geom_smooth(method = "loess") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 8, bs = "cs")) +
  coord_polar(theta = "y", clip = "off") +
  # scale_y_continuous(limits = c(0, 6.33)) + #facet_wrap(~run) +
  theme_bw()



ggplot(df, aes(run_trial, theta, color = rewFunc)) +
  # geom_smooth(method = "loess") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3, bs = "cs")) +
  coord_polar(theta = "y", clip = "off") +
  scale_y_continuous(limits = c(0, 6.33)) + facet_wrap(~run) +
  theme_bw()


m1 <- lmer(theta ~ theta_rt*TR_num + theta_rt_vmax*TR_num + TR_num + (1|id), df)
summary(m1)
car::Anova(m1)
m1a <- lmer(theta ~ theta_rt*trs_until_rt + theta_rt_vmax*trs_until_rt + TR_num + (TR_num|id), df)
summary(m1a)
car::Anova(m1a)
anova(m1, m1a)

m2 <- lmer(theta ~ theta_rt*tr_num_f + rewFunc*tr_num_f + (1|id), df)
summary(m2)
car::Anova(m2)
m2a <- lmer(theta ~ theta_rt*trs_until_rt_f + rewFunc*trs_until_rt_f + (1|id), df)
summary(m2a)
car::Anova(m2)

anova(m2, m2a)

library(circular)
set.seed(1234)
x <- cbind(rnorm(10), rep(1, 10))
y <- circular(2*atan(c(x%*%c(5,1))))+rvonmises(10, mu=circular(0), kappa=100)
lm.circular(y=y, x=x, init=c(5,1), type='c-l', verbose=TRUE)

cm1 <- lm.circular(circular(df$theta), df$TR_num, init = 1, type = "c-l", verbose = T)
cm2 <- lm.circular(circular(df$theta), circular(df$theta_rt),  type = "c-c")

x = cbind(rep(1, length(df$TR_num)), df$TR_num, df$theta_rt, df$theta_rt_vmax, scale(df$TR_num)*df$theta_rt)
cm3 <- lm.circular(circular(df$theta), x, init = c(rep(1, dim(x)[2])), type = "c-l", verbose = T)

# try to learn it in brms

library(rethinking)
data(Howell1)
d <- Howell1
rm(Howell1)
detach(package:rethinking, unload = T)
library(brms)
d %>% select(height) %>% glimpse()
n <- 200
d2 <- 
  d %>%
  filter(age >= 18)
d_grid <-
  # we'll accomplish with `tidyr::crossing()` what McElreath did with base R `expand.grid()`
  crossing(mu    = seq(from = 140, to = 160, length.out = n),
           sigma = seq(from = 4, to = 9, length.out = n))

grid_function <- function(mu, sigma) {
  
  dnorm(d2$height, mean = mu, sd = sigma, log = T) %>% 
    sum()
  
}
d_grid <-
  d_grid %>% 
  mutate(log_likelihood = map2(mu, sigma, grid_function)) %>%
  unnest(log_likelihood) %>% 
  mutate(prior_mu    = dnorm(mu, mean = 178, sd = 20, log = T),
         prior_sigma = dunif(sigma, min = 0, max = 50, log = T)) %>% 
  mutate(product = log_likelihood + prior_mu + prior_sigma) %>% 
  mutate(probability = exp(product - max(product)))

head(d_grid)
library(brms)
b4.5 <- 
  brm(data = d2, 
      family = gaussian,
      height ~ 1,
      prior = c(prior(normal(178, 10), class = Intercept),
                prior(uniform(0, 50), class = sigma, ub = 50)),
      iter = 2000, warmup = 1000, chains = 4, cores = 10,
      seed = 4,
      file = "fits/b04.05")
print(b4.5)
