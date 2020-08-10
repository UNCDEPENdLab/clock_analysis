library(emmeans)
library(multcomp)
dmat <- data.frame(emotion=factor(
  c("Scram", "Fear", "Scram", "Fear", 
    "Scram", "Happy", "Scram", "Happy")), run=0:7, dummy=rnorm(8) #run=1:8, dummy=rnorm(8)
)
dmat$emotion <- relevel(dmat$emotion, ref="Scram")
mm <- lm(dummy ~ emotion + run, dmat)
summary(mm)
model.matrix(mm)

#emotion condition averages
xx <- lsmeans(mm, "emotion")
xx@linfct #emotion effects with run at the mean

cm <- contrast(xx, "pairwise")
cm@linfct #matches my setup in FSL
cm

#trying to figure out overall average
mm_norun <- lm(dummy ~ emotion, dmat)
xx <- lsmeans(mm_norun, "emotion")
(cm <- contrast(xx, "eff"))
cm@linfct


#figure out run effect
(xx <- lsmeans(mm, "run"))
xx@linfct #contrast for effect of run at the average emotion

summary(mm)

summary(lm(dummy~run, dmat))
cm_manual <- rbind(c(1, 0.25, 0.25, 1))
summary(glht(mm, cm_manual))


dmat$res <- resid(lm(dummy ~ emotion, dmat))
summary(lm(res ~ run, dmat))



#the 0.33 coding is not making sense to me from the perspective of weighted averaging
#take a crack at some ridiculous mean differences in contrast...

dmat <- data.frame(emotion=factor(
  c("Scram", "Fear", "Scram", "Fear", 
    "Scram", "Happy", "Scram", "Happy")), run=0:7,
  dv=c(10, 2, 11, 3, 9, 20, 10, 22)
)

dmat$emotion <- relevel(dmat$emotion, ref="Scram")
mm <- lm(dv ~ emotion, dmat)
summary(mm)
model.matrix(mm)


(xx <- lsmeans(mm, "emotion"))
xx@linfct #run at the mean

(cm <- contrast(xx, "pairwise"))
cm@linfct #matches my setup in FSL

(cm <- contrast(xx, "eff"))
cm@linfct

cm_manual <- rbind(c(1, 0.25, 0.25),
                   c(1, 0.33, 0.33))
library(multcomp)
summary(glht(mm, cm_manual))

tapply(dmat$dv, dmat$emotion, mean)
mean(dmat$dv)

#As anticipated, the 1, 0.25, 0.25 weighted coding does yield the
#grand mean in the sample: 10.875, whereas the 0.33 coding yields 11.155,
#overweighting the fear and happy cells.


