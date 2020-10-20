library(lme4)
data <- read.csv("/Users/mnh5174/Downloads/Data_for_testing.csv")
data$Pow_dB = 10*log10(data$Pow + 101)
data$Trial_z <- as.vector(scale(data$Trial))
data$Trial_rel <- data$Trial %% 63 #63 trials per run, right? This will be slightly off with truncations...

#For python comparison. Doesn't match. Not totally sure wht
mm_compare <- lmer(Pow ~ Faces  + Age + Faces*Age + Trial + Rewarded + Rewarded*Faces + (1 + Trial | Subject/Run), data=data, REML=TRUE)
summary(mm_compare)

#problems:
#1. Need to log transform power to have outcome be approximately normal

hist(data$Pow)
hist(data$Pow_dB)

#2. Convergence is hindered or will fail when variables are on very different scales.
#   This is the source of problems with trial

var(data$Age)
var(data$Trial)

#basic model with trial fails because of scaling
mm <- lmer(Pow_dB ~ Faces  + Age + Faces*Age + Trial + Rewarded + Rewarded*Faces + (1 + Trial | Subject/Run), data=data)

#3. Random effect specification
#There is essentially no variation in Trial slopes across runs compared to overall slope effects,
#  and it leads to excessive correlations for within-run random effects.
#Altogether, the trial slope is not moderated by run. This seems overparameterized conceptually.

#using relative trial is conceptually appealing, but doesn't solve the random effects problem
mm2 <- lmer(Pow_dB ~ Faces  + Age + Faces*Age + Trial_rel + Rewarded + Rewarded*Faces + (1 + Trial_rel | Subject/Run), data=data)

#ultimately, it appears we don't need a nested slope for trial
#here, we get a) random slope of trial per subject, b) random intercept per subject,
#             c) random intercept per run within subject (i.e., run-specific intercepts)
#This specification also omit correlations among random effects, which seems fine.
mm_mod <- lmer(Pow_dB ~ Faces  + Age + Faces*Age + Trial + Rewarded + Rewarded*Faces + 
                 (0 + Trial | Subject) + (1 | Subject/Run), data=data, REML=FALSE)

#z-scoring trial will fix convergence and speed up estimation massively (scaling problems)
mm_mod2 <- lmer(Pow_dB ~ Faces  + Age + Faces*Age + Trial_z + Rewarded + Rewarded*Faces + 
                 (0 + Trial_z | Subject) + (1 | Subject/Run), data=data, REML=FALSE)

#what about removing overall RE intercept and just having intercepts within runs
mm_mod3 <- lmer(Pow_dB ~ Faces  + Age + Faces*Age + Trial_z + Rewarded + Rewarded*Faces + 
                 (0 + Trial_z | Subject) + (1 | Subject:Run), data=data, REML=FALSE)

#in this case, we get a small boost from overall intercept.
anova(mm_mod2, mm_mod3)

#here's the RE structure of the winning model. Makes sense.
str(ranef(mm_mod2))

#winning model
summary(mm_mod2)