setwd("/Users/michael/CogEmoFaceReward/analysis")
library(plyr)
library(ggplot2)
library(reshape2)
library(gdata)

source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))

#questionnaire data
behav <- read.xls("clock_questionnaires_n36.xlsx", sheet="Sheet1", skip=1)
behav$LunaID <- factor(behav$LunaID)


png("AgeHist.png", width=700, height=500, res=300)
ggplot(behav, aes(x=AgeAtVisit)) + geom_histogram(binwidth=2) + xlab("Age at Visit") + ylab("Count") + theme_bw(base_size=12)
dev.off()

##compare fits across subjects and models
fitFiles <- list.files(path="../fit_behavior", pattern="SubjsSummary.*", full.names = TRUE)
allM <- list()
for (f in fitFiles) {
    model <- sub(pattern="^.*SubjsSummary_(.*)\\.txt", replacement="\\1", x=f, perl=TRUE)
    p <- read.table(f, header=TRUE, sep="\t")
    p$nparams <- length(which(names(p) %notin% c("Subject", "Session", "ignore", "SSE", "model")))
    ##p <- subset(p, select=c("Subject", "SSE", "nparams"))
    if (model == "noemo_scram" || model=="noemosticky_scram") {
      p$ntrials <- 168 #fixed: 4 runs of 42 trials
    } else { 
      p$ntrials <- 504 #fixed: 12 runs of 42 trials 
    }
    p$Subject <- factor(as.integer(p$Subject))
    p$model <- factor(model)
    allM[[f]] <- p
}

allM <- do.call(rbind.fill, allM)

allM$AIC <- with(allM, ntrials*(log(2*pi*(SSE/ntrials))+1) + 2*nparams)

nosticky <- gdata::drop.levels(subset(allM, model %in% levels(allM$model)[!grepl("sticky", levels(allM$model), fixed=TRUE)]))
noscram <- gdata::drop.levels(subset(allM, model %in% levels(allM$model)[!grepl("scram", levels(allM$model), fixed=TRUE)]))

##identify subjects who had a positive epsilon in at least one non-sticky model
posEpsSubjects <- unique(unlist(subset(nosticky, (explore > 0 | (explore_scram > 0 | explore_fear > 0 | explore_happy > 0)), select=Subject)))

posEps <- gdata::drop.levels(subset(nosticky, Subject %in% posEpsSubjects))


##for compatibility with spm_BMS.m (Bayesian Model Selection),
##need a subjects x models AIC matrix
#AICmat <- do.call(cbind, lapply(split(allM, allM$model), "[[", "AIC"))
AICmat <- do.call(cbind, lapply(split(noscram, noscram$model), "[[", "AIC"))

AICmat_posEps <- do.call(cbind, lapply(split(posEps, posEps$model), "[[", "AIC"))

##per discussion with Michael Frank (and review of Stephan 2009),
##the use of AIC as an approximation of log-evidence is the LL - nparams
##Thus, for use with spm_BMS, we need to negate the values of AIC computed
AICmat <- -1*AICmat
AICmat_posEps <- -1*AICmat_posEps

library(R.matlab)
#writeMat(con="AICmatrix_n36.mat", AICmat=AICmat, mnames=levels(allM$model))
writeMat(con="AICmatrix_n36.mat", AICmat=AICmat, mnames=levels(noscram$model))
#writeMat(con="AICmatrix_n36.mat", AICmat=AICmat_posEps, mnames=levels(posEps$model))

##run SPM BMS in MATLAB
#system("matlab -nodisplay < computeBMSprobs.m")
system("/Applications/MATLAB_R2013b.app/bin/matlab -nodisplay < computeBMSprobs.m")

BMSresults <- readMat(con="AICresults_n36.mat")

##form into a data.frame
#BMSresults <- with(BMSresults, data.frame(model=unlist(mnames), alpha=as.vector(alpha), expr=as.vector(expr), xp=as.vector(xp), mAIC=apply(AICmat, 2, mean)))
BMSresults <- with(BMSresults, data.frame(model=unlist(mnames), alpha=as.vector(alpha), expr=as.vector(expr), xp=as.vector(xp), mAIC=apply(AICmat, 2, mean)))

library(xtable)
print(xtable(BMSresults[order(BMSresults$mAIC),]), include.rownames=FALSE)

allM <- ddply(allM, .(Subject), function(subdf) {
    minAIC <- min(subdf$AIC)
    minSSE <- min(subdf$SSE)
    subdf$AICdiff <- subdf$AIC - minAIC
    subdf$SSEdiff <- subdf$SSE - minSSE
    subdf
})

png("Subject_SSE_byModel.png", width=8, height=6, units="in", res=300)
ggplot(allM, aes(x=Subject, y=SSE, color=model)) + geom_jitter(size=2, position = position_jitter(width = .3, height=0)) + xlab("Subject") + ylab("SSE") +
    theme(axis.text.x=element_text(angle=90)) + scale_color_brewer(palette="Dark2") + coord_flip() + theme_bw(base_size=14)
dev.off()

png("Avg_SSE_byModel.png", width=8, height=6, units="in", res=300)
ggplot(allM, aes(x=model, y=SSE)) + geom_boxplot() + xlab("Model") + ylab("SSE") +
    theme(axis.text.x=element_text(angle=90))
dev.off()

png("Avg_AIC_byModel.png", width=8, height=6, units="in", res=300)
ggplot(allM, aes(x=model, y=AIC)) + geom_boxplot() + xlab("Model") + ylab("AIC") +
    theme(axis.text.x=element_text(angle=90))
dev.off()

sort(tapply(allM$AIC, allM$model, mean))
tapply(allM$AIC, allM$model, sd)

library(nlme)
modelDiffs <- lme(AIC ~ model, random = ~1 | Subject, data=allM)
#library(lme4)
#modelDiffs <- lmer(AIC ~ model + (1 | Subject), data=allM)
anova(modelDiffs)
summary(modelDiffs)
library(multcomp)
summary(glht(modelDiffs, linfct=mcp(model="Tukey")))

png("Subject_AICdiff_byModel.png", width=8, height=6, units="in", res=300)
ggplot(allM, aes(x=Subject, y=AICdiff, color=model)) + geom_jitter(size=2, position = position_jitter(width = .3, height=0)) + xlab("Subject") + ylab("AIC difference (from best)") +
    theme(axis.text.x=element_text(angle=90)) + scale_color_brewer(palette="Dark2") + coord_flip() + theme_bw(base_size=14)
dev.off()

png("Subject_SSEdiff_byModel.png", width=8, height=6, units="in", res=300)
ggplot(allM, aes(x=Subject, y=AICdiff, color=model)) + geom_jitter(size=2, position = position_jitter(width = .3, height=0)) + xlab("Subject") + ylab("SSE difference (from best)") +
    theme(axis.text.x=element_text(angle=90)) + scale_color_brewer(palette="Dark2") + coord_flip() + theme_bw(base_size=14)
dev.off()



##look at how allowing epsilon to be negative improves fits
noemo <- gdata::drop.levels(subset(allM, model %in% c("noemo", "noemosticky")))
noemo <- ddply(noemo, .(Subject), function(subdf) {
    negeps <- subset(subdf, model=="noemosticky")
    poseps <- subset(subdf, model=="noemo")

    data.frame(Subject=negeps$Subject[1L],
               epsChange=(poseps$explore - negeps$explore),
               AICchange=(poseps$AIC - negeps$AIC),
               stickyEps=negeps$explore,
               constrainedEps=poseps$explore)
})

cor.test(~ AICchange + epsChange, noemo)
plot(~ AICchange + epsChange, noemo)


#learningParams <- read.table("../fit_behavior/SubjsSummary_emoexplore.txt", header=TRUE)

learning_emoexplore <- read.table("../fit_behavior_matlab/SubjsSummary_emoexploresticky.txt", header=TRUE)
learning_noemo <- read.table("../fit_behavior_matlab/SubjsSummary_noemosticky.txt", header=TRUE)
#learning_noemo <- read.table("../fit_behavior_matlab/SubjsSummary_noemo.txt", header=TRUE)
#learning_noemo <- read.table("../fit_behavior/SubjsSummary_noemo_scram.txt", header=TRUE)

learning_emoexplore <- rename(learning_emoexplore, c(Subject="LunaID"))
learning_emoexplore <- merge(learning_emoexplore, behav[,c("LunaID", "AgeAtVisit", 
            "UPPS_Urg", "UPPS_PosUrg", "UPPS_SS", "UPPS_Prem", "UPPS_Pers", 
            "NN", "NSR", "EPA", "ESA", "OAI", "OII", "OU", "ANO", "APO", "CO", "CGS", "CD",
            "RIST.INDEZ", "SSS_Total", "TAS", "ES", "DIS", "BS", 
            "ADI_Total", "DERS_Total", "DERS_Clarity", "DERS_Strategies", "DERS_Awareness", "DERS_Impulse", "DERS_NonAccept", 
            "STAI_Score", "ADI_Emotional", "ADI_Behavioral", "ADI_Cognitive", "RPI", "RSE")], by="LunaID")

learning_noemo <- rename(learning_noemo, c(Subject="LunaID"))
learning_noemo <- merge(learning_noemo, behav[,c("LunaID", "AgeAtVisit", 
            "UPPS_Urg", "UPPS_PosUrg", "UPPS_SS", "UPPS_Prem", "UPPS_Pers",
            "NN", "NSR", "EPA", "ESA", "OAI", "OII", "OU", "ANO", "APO", "CO", "CGS", "CD",
            "RIST.INDEZ", "SSS_Total", "TAS", "ES", "DIS", "BS", 
            "ADI_Total", "DERS_Total", "DERS_Clarity", "DERS_Strategies", "DERS_Awareness", "DERS_Impulse", "DERS_NonAccept", 
            "STAI_Score", "ADI_Emotional", "ADI_Behavioral", "ADI_Cognitive", "RPI", "RSE")], by="LunaID")

learning_noemo$explorePos <- sapply(learning_noemo$explore, function(x) { ifelse(x > 0, x, NA) })
learning_noemo$exploreGt0 <- sapply(learning_noemo$explore, function(x) { ifelse(x > 0, 1, 0) })

cor.test(~ exploreGt0 + AgeAtVisit, learning_noemo)
cor.test(~ explore + AgeAtVisit, learning_noemo)
cor.test(~ explorePos + AgeAtVisit, learning_noemo)
cor.test(~ explore + DERS_Total, learning_noemo)
cor.test(~ alphaG + SSS_Total, learning_noemo)
cor.test(~ alphaN + SSS_Total, learning_noemo)

#go and no go learning rates
cor.test(~ alphaG + AgeAtVisit, learning_noemo)
cor.test(~ alphaN + AgeAtVisit, learning_noemo)



ggplot(learning_noemo, aes(x=AgeAtVisit, y=alphaG)) + geom_point() + stat_smooth()
ggplot(learning_noemo, aes(x=AgeAtVisit, y=alphaN)) + geom_point() + stat_smooth()


ggplot(learning_noemo, aes(x=AgeAtVisit, y=exploreGt0)) + geom_point()

corstarsl(learning_noemo, omit=c("LunaID", "Session", "ignore"))
corstarsl(learning_emoexplore, omit=c("LunaID", "Session", "ignore"))
corstarsl(learning_emoexplore, omit=c("LunaID", "Session", "ignore"))

params <- c("lambda", "alphaN", "alphaG", "explore_scram", "explore_fear", "explore_happy", "K", "sticky_decay", "rho", "SSE")
#params <- c("lambda", "alphaN", "alphaG", "explore", "K", "sticky_decay", "rho", "SSE")
selfreports <- c("UPPS_Urg", "UPPS_PosUrg", "UPPS_SS", "UPPS_Prem", "UPPS_Pers", 
"NN", "NSR", "EPA", "ESA", "OAI", "OII", "OU", "ANO", "APO", "CO", "CGS", "CD",
"RIST.INDEZ", "SSS_Total", "TAS", "ES", "DIS", "BS",
"DERS_Total", "DERS_Clarity", "DERS_Strategies", "DERS_Awareness", "DERS_Impulse", "DERS_NonAccept",
"ADI_Total", "DERS_Total", "STAI_Score", "ADI_Emotional", "ADI_Behavioral", "ADI_Cognitive",
"STAI_Score", "RSE", "RSI")

sigrs <- corwithtarget(learning_emoexplore, pmin=.05, omit=c("LunaID", "Session", "ignore", "explore_HappyMScramble", "explore_FearMScramble"), target=params)
#sigrs <- corwithtarget(learning_noemo, pmin=.05, omit=c("LunaID", "Session", "ignore", "explore_HappyMScramble", "explore_FearMScramble"), target=params)

#sigrs <- corwithtarget(learning_noemo, pmin=.05, target="AgeAtVisit", with=c(params, selfreports))

#follow-up to see whether there are linear or quadratic interactions with age
pdf("learning_by_age_interactions.pdf", width=11, height=7)
for (p in 1:length(sigrs)) {
  if (is.null(sigrs[[p]])) { next }
  p_name <- names(sigrs)[p]
  
  #get variables with which this is correlated
  sreports <- dimnames(sigrs[[p]])[[2]]
  
  for(s in 1:length(sreports)) {
    s_name <- sreports[s]
#    df_test <- data.frame(
#        as.vector(scale(learning_emoexplore[[ p_name ]], scale=FALSE)),
#        learning_emoexplore[[ s_name ]],
#        as.vector(scale(learning_emoexplore$AgeAtVisit, scale=FALSE))
#    )
    
    df_test <- data.frame(
        learning_emoexplore[[ p_name ]],
        learning_emoexplore[[ s_name ]],
        learning_emoexplore$AgeAtVisit
    )
#    df_test <- data.frame(
#        as.vector(scale(learning_noemo[[ p_name ]], scale=FALSE)),
#        learning_noemo[[ s_name ]],
#        as.vector(scale(learning_noemo$AgeAtVisit, scale=FALSE))
#    )    
    
    names(df_test) <- c(p_name, s_name, "AgeAtVisit")
    #form <- as.formula(paste(s_name, "~", p_name, "*AgeAtVisit + ", p_name, "*I(AgeAtVisit^2)"))
    form <- as.formula(paste(s_name, "~", p_name, "*AgeAtVisit"))
    model <- lm(form, data=df_test)
    pvals <- summary(model)$coefficients[,c("Pr(>|t|)")][-1L] #-1L to drop p-value for intercept
    #if (any(pvals < .05)) {
    #if (pvals[paste0(p, ":I(AgeAtVisit^2)")] < .05) {
    if (pvals[paste0(p_name, ":AgeAtVisit")] < .05) {
      
      cat("------\n\nParameter:", p_name, ", Self-report:", s_name, "\n")
      print(summary(model))
      cat("\n-----\n\n")
      
      p <- plot(effect(term=paste0(p_name, ":AgeAtVisit"),mod=model,default.levels=5, x.var="AgeAtVisit"),multiline=TRUE)
      print(p)
    }
  }
}
dev.off()

learning_emoexplore$alphaG.c <- learning_emoexplore$alphaG - mean(learning_emoexplore$alphaG, na.rm=TRUE)
learning_emoexplore$AgeAtVisit.c <- learning_emoexplore$AgeAtVisit - mean(learning_emoexplore$AgeAtVisit, na.rm=TRUE)
#adiModel <- lm(ADI_Total ~ alphaG.c*AgeAtVisit.c, learning_emoexplore)
adiModel <- lm(ADI_Total ~ alphaG*AgeAtVisit, learning_emoexplore)

gr <- expand.grid(
    AgeAtVisit = seq(14, 31, by=0.5),
    alphaG = c(-0.1544, 0, .1544))
    #ADI_Total=0)

gr$pred <- predict(adiModel, gr)
gr$alphaG <- ordered(gr$alphaG)
ggplot(gr, aes(x=AgeAtVisit, y=pred, color=alphaG)) + geom_line()

upps <- lm(UPPS_PosUrg ~ alphaG*AgeAtVisit, learning_emoexplore)

sigrs <- corwithtarget(learning_emoexplore, omit=c("LunaID", "Session", "ignore", "explore_HappyMScramble", "explore_FearMScramble"), target="AgeAtVisit")
sigrs <- corwithtarget(learning_noemo, omit=c("LunaID", "Session", "ignore", "explore_HappyMScramble", "explore_FearMScramble"), target="AgeAtVisit")
    #withvars=c("UPPS_SS", "UPPS_PosUrg", "NN"))

print(sigrs)

pdf("param_age_interactions.pdf", width=9, height=9)
for (p in params) {
  #g <- ggplot(learning_emoexplore, aes_string(x="AgeAtVisit", y=p)) + geom_point() + stat_smooth(width=1.5, method="loess") + ggtitle(p)
  g <- ggplot(learning_noemo, aes_string(x="AgeAtVisit", y=p)) + geom_point() + stat_smooth(width=1.5, method="loess") + ggtitle(p)
  print(g)
}
dev.off()


#just look at at age interactions with self-reports
pdf("sig_AgeSelfReport_Assoc.pdf", width=9, height=9)
for (s in selfreports) {
  df_test <- data.frame(
      learning_emoexplore[[s]],
      as.vector(scale(learning_emoexplore$AgeAtVisit, scale=FALSE)) #center age to reduce linear-quadratic collinearity
  )
  names(df_test) <- c(s, "AgeAtVisit")
  
  form <- as.formula(paste(s, "~ AgeAtVisit + I(AgeAtVisit^2)"))
  model <- lm(form, data=df_test)
  
  pvals <- summary(model)$coefficients[,c("Pr(>|t|)")][-1L] #-1L to drop p-value for intercept
  if (any(pvals < .05)) {
    
    predData <- list()
    
    predData[["AgeAtVisit"]] <- seq(min(df_test$AgeAtVisit, na.rm=TRUE), max(df_test$AgeAtVisit, na.rm=TRUE), length=cont.points)
    
    predData[[s]] <- 0 #dv
    
    pred_df <- do.call(expand.grid, predData)
    
    pred_df$`I(AgeAtVisit^2)` <- pred_df$AgeAtVisit^2
    
    pred_df[[s]] <- predict(model, pred_df)
    pred_df[["AgeAtVisit"]] <- pred_df[["AgeAtVisit"]] + mean(learning_emoexplore[["AgeAtVisit"]]) #add mean back in for plotting 
    
    g<- ggplot(pred_df, aes_string(x="AgeAtVisit", y=s)) + geom_line() + ggtitle(s) + geom_point(data=learning_emoexplore)
    print(g)
    
    cat("Self-report:", s, "\n")
    print(summary(model)) 
  }
}
dev.off()

#conclusion: age reductions for 
# UPPS_PosUrg (quadratic down)
# UPPS_SS (quadratic down)
# UPPS_Prem (linear down)
# NSR (Neuroticism self-reproach) quad down
# CO: orderliness linear up
# ADI_Total: linear down
# STAI_SCORE: quad down


#check age x process interactions for each variable
n.divide <- 3
cont.points <- 20
library(ggplot2)
pdf("sigplots.pdf", width=9, height=9)
for (p in params) {
  for (s in selfreports) {
    df_test <- data.frame(
        as.vector(scale(learning_emoexplore[[p]], scale=FALSE)),
        learning_emoexplore[[s]],
        as.vector(scale(learning_emoexplore$AgeAtVisit, scale=FALSE))
    )
    names(df_test) <- c(p, s, "AgeAtVisit")
    
    #form <- as.formula(paste(s, "~", p, "*AgeAtVisit + ", p, "*I(AgeAtVisit^2)"))
    form <- as.formula(paste(s, "~", p, "*AgeAtVisit"))
    model <- lm(form, data=df_test)
    pvals <- summary(model)$coefficients[,c("Pr(>|t|)")][-1L] #-1L to drop p-value for intercept
    #if (any(pvals < .05)) {
    #if (pvals[paste0(p, ":I(AgeAtVisit^2)")] < .05) {
    if (pvals[paste0(p, ":AgeAtVisit")] < .05) {
      #generate plot
      #p_sd <- sd(learning_emoexplore[[p]], na.rm=TRUE)
      #p_m <- mean(learning_emoexplore[[p]], na.rm=TRUE)
      
      p_sd <- sd(df_test[[p]], na.rm=TRUE)
      p_m <- mean(df_test[[p]], na.rm=TRUE) #should be zero, but whatever
      
      predData <- list()
      
      predData[[p]] <- if (n.divide==3) { c(p_m-p_sd, p_m, p_m+p_sd)
          } else { c(p_m-p_sd*2, p_m-p_sd, p_m, p_m+p_sd, p_m+p_sd*2) }     
      
      predData[["AgeAtVisit"]] <- seq(min(df_test$AgeAtVisit, na.rm=TRUE), max(df_test$AgeAtVisit, na.rm=TRUE), length=cont.points)
      predData[[s]] <- 0 #dv
      pred_df <- do.call(expand.grid, predData)
      #pred_df$`I(AgeAtVisit^2)` <- pred_df$AgeAtVisit^2
      
      pred_df[[s]] <- predict(model, pred_df)
      
      pred_df[[p]] <- factor(plyr::round_any(pred_df[[p]], .001)) #since we discretized parameter into SD ranges, need to store as factor for visual display
      pred_df[["AgeAtVisit"]] <- pred_df[["AgeAtVisit"]] + mean(learning_emoexplore[["AgeAtVisit"]]) #add mean back in for plotting
      
      print(ggplot(pred_df, aes_string(x="AgeAtVisit", y=s, color=p)) + geom_line() + ggtitle(paste(s, p)))
      
      cat("Parameter:", p, ", Self-report:", s, "\n")
      print(summary(model)) 
    }
  }
}
dev.off()

plot(effect(term="lambda:AgeAtVisit",mod=model,default.levels=10),multiline=TRUE)

#explore by emotion
explore.melt <- melt(learning_emoexplore[,c("LunaID", "AgeAtVisit", "UPPS_SS", "UPPS_Urg", "UPPS_PosUrg", 
            "explore_scram", "explore_fear", "explore_happy")], id.vars=c("LunaID", "AgeAtVisit", "UPPS_SS", "UPPS_Urg", "UPPS_PosUrg"))

learning_emoexplore$explore_HappyMScramble <- with(learning_emoexplore, explore_happy - explore_scram)
learning_emoexplore$explore_FearMScramble <- with(learning_emoexplore, explore_fear - explore_scram)


hms <- ggplot(learning_emoexplore, aes(x=factor(1:36), y=sort(explore_HappyMScramble))) + geom_bar(stat="identity") +
        coord_flip() + theme_bw(base_size=16) + xlab("Subject") + ylab("Explore (Happy - Scrambled)") +
        theme(axis.text.y=element_blank())
print(hms)

fms <- ggplot(learning_emoexplore, aes(x=factor(1:36), y=sort(explore_FearMScramble))) + geom_bar(stat="identity") +
        coord_flip() + theme_bw(base_size=16) + xlab("Subject") + ylab("Explore (Fear - Scrambled)") + 
        theme(axis.text.y=element_blank())
print(fms)



library(gridExtra)
pdf("Within-subjects EmoExplore Diffs from Scrambled.pdf", width=11, height=6)
grid.arrange(hms, fms, ncol=2)
dev.off()

#non-parametric repeated measures test for explore parameter
friedman.test(value ~ variable | LunaID, data = explore.melt)
library(ez)

library(lme4)

summary(glmer(value > 0 ~ variable | LunaID, data = explore.melt, family=binomial))

tapply(explore.melt$value, explore.melt$variable, mean)
ezANOVA(explore.melt, dv="value", wid="LunaID", within="variable")
summary(lm(value ~ variable, data = explore.melt))

summary(lm(value ~ variable*AgeAtVisit*UPPS_Urg, data = explore.melt))
anova(lm(value ~ variable*AgeAtVisit*UPPS_Urg, data = explore.melt))

library(nlme)
anova(lme(value ~ variable, random=~1 | LunaID, data=explore.melt))
anova(lme(value ~ variable*AgeAtVisit, random=~1 | LunaID, data=explore.melt))

summary(lme(value ~ variable*AgeAtVisit*UPPS_Urg, random=~1 | LunaID, data=explore.melt))
anova(lme(value ~ variable*scale(AgeAtVisit)*scale(UPPS_Urg), random=~1 | LunaID, data=explore.melt))
anova(lme(value ~ variable*AgeAtVisit*UPPS_PosUrg, random=~1 | LunaID, data=explore.melt))
anova(lme(value ~ variable*AgeAtVisit*UPPS_SS, random=~1 | LunaID, data=explore.melt))

#ageUrg <- lmer(value ~ variable*AgeAtVisit*UPPS_Urg + (1 | LunaID), explore.melt)
ageUrg <- lmer(value ~ variable*UPPS_Urg + (1 | LunaID), explore.melt)
summary(ageUrg)
anova(ageUrg)

cm <- lmerCellMeans(ageUrg, divide="AgeAtVisit")
cm <- lmerCellMeans(ageUrg)

ggplot(cm, aes(x=UPPS_Urg, y=value, shape=AgeAtVisit, color=variable)) + scale_color_brewer("Emotion") + ylab("Explore") +
    xlab("Negative Urgency") + scale_shape("Age At Visit") + geom_line() + geom_point()

pdf("UPPS Urgency by Emotion Exploration.pdf", width=9, height=6)
ggplot(cm, aes(x=UPPS_Urg, y=value, color=variable, ymin=value-se, ymax=value+se)) + scale_color_brewer("Emotion", palette="Dark2") + ylab("Explore") +
    xlab("Negative Urgency") + geom_line() + geom_point() + geom_errorbar() + theme_bw(base_size=16) +
    geom_hline(yintercept=0)
dev.off()


png("Explore_by_emotion.png", width=12, height=10, units="in", res=300)
ggplot(explore.melt, aes(x=value)) + geom_histogram(binwidth=2000) + 
        facet_wrap(~variable) + ggtitle("Exploration parameter by emotion")
dev.off()

#plot scatter plot of exploration by condition and age
ggplot(explore.melt, aes(x=AgeAtVisit, y=value)) + geom_point() + stat_smooth(se=FALSE, method="lm") + stat_smooth(se=FALSE) +
    facet_wrap(~variable) + ggtitle("Exploration parameter by emotion and age")

ggplot(explore.melt, aes(x=UPPS_SS, y=value)) + geom_point() + stat_smooth(se=FALSE, method="lm") + facet_wrap(~variable) + ggtitle("Exploration parameter by emotion and sensation seeking")

ggplot(explore.melt, aes(x=UPPS_Urg, y=value)) + geom_point() + stat_smooth(se=FALSE, method="lm") + facet_wrap(~variable) + ggtitle("Exploration parameter by emotion and negative urgency")
ggplot(explore.melt, aes(x=UPPS_PosUrg, y=value)) + geom_point() + stat_smooth(se=FALSE, method="lm") + facet_wrap(~variable) + ggtitle("Exploration parameter by emotion and negative urgency")

ggplot(explore.melt, aes(x=AgeAtVisit, y=value)) + geom_point() + stat_smooth(se=FALSE, method="lm") + facet_wrap(~variable) + ggtitle("Exploration parameter by emotion and sensation seeking")

ggplot(behav, aes(x=AgeAtVisit, y=alphaG)) + geom_point() + stat_smooth(se=FALSE, method="lm") + ggtitle("Go param by age")
ggplot(behav, aes(x=AgeAtVisit, y=alphaN)) + geom_point() + stat_smooth(se=FALSE, method="lm") + ggtitle("NoGo param by age")
ggplot(behav, aes(x=AgeAtVisit, y=rho)) + geom_point() + stat_smooth(se=FALSE, method="lm") + ggtitle("rho param by age")

##More simple checks and analyses of RTs
##Build big data frame of subjects, trials, RTs, and conditions

##sadly, readMat is very slow!
##tcMats <- list.files(path="../subjects", pattern=".*_tc\\.mat", full.names=TRUE)
## for (f in tcMats) {
##     sdata <- readMat(f)
##     sdf <- data.frame(subject=as.vector(sdata$subject[[1L]]),
##                       sex=as.vector(sdata$subject[[2L]]),
##                       condition=unlist(lapply(sdata$order, "[[", 1), use.names=FALSE),
                      

##                       )
##     browser()
## }

tcFiles <- list.files(path="../subjects", pattern=".*_tc\\.txt", full.names=TRUE)

allData <- list()
for (f in tcFiles) {
    subject <- sub("^.*/(\\d+)_tc\\.txt$", "\\1", f, perl=TRUE)
    sdata <- read.table(f, header=TRUE, comment.char="#")
    sdata$Null <- NULL #delete dummy column
    sdata$Subject <- factor(subject)
    sdata <- ddply(sdata, .(Emotion, Func), function(subdf) {
        bestRew <- -1
        bestEV <- -1
        bestRewRT <- subdf[1, "RT"]
        bestEVRT <- subdf[1, "RT"]
        subdf$bestRewRT <- NA_real_
        subdf$bestEVRT <- NA_real_
        for (i in 1:nrow(subdf)) {
            if (subdf[i,"ScoreInc"] >= bestRew) {
                bestRewRT <- subdf[i,"bestRewRT"] <- subdf[i,"RT"]
                bestRew <- subdf[i,"ScoreInc"]
            } else {
                subdf[i,"bestRewRT"] <- bestRewRT
            }
            if (subdf[i,"ScoreInc"] > 0 && subdf[i,"EV"] >= bestEV) {
                bestEVRT <- subdf[i,"bestEVRT"] <- subdf[i,"RT"]
                bestEV <- subdf[i,"EV"]
            } else {
                subdf[i,"bestEVRT"] <- bestEVRT
            }
        }
        subdf
    })
    allData[[f]] <- sdata
}

##rearrange column headers for readability
allData <- do.call(rbind, allData)
row.names(allData) <- NULL
allData <- allData[,c("Subject", "Run", "Block", "Trial", "Func", "Emotion", "Mag", "Freq", "ScoreInc", "EV", "RT", "bestRewRT", "bestEVRT", "Image")]
allData <- plyr::rename(allData, c(Subject="LunaID"))

##verify that each subject completed 42 trials for each Func x Emotion condition
with(allData, table(LunaID, Func, Emotion))

##compute trial within a block
allData$TrialRel <- unlist(lapply(split(allData, f=list(allData$LunaID, allData$Func, allData$Emotion)), function(l) { return(1:nrow(l)) } ))
##allData$TrialRel2 <- 1:42 ##shouldn't this be identical and easier? :) Just use recycling

allData <- merge(allData, behav[,c("LunaID", "AgeAtVisit")], by="LunaID")
allData$Half <- factor(sapply(allData$TrialRel, function(x) { ifelse(x > 21, "H2", "H1")}))



##compute median reaction time across subjects
allSubjAgg <- ddply(allData, .(LunaID, Func), function(subdf) {
    subdf <- subset(subdf, RT > 80)
    med <- median(subdf$RT, na.rm=TRUE)
    return(c(med=med))
})

with(allSubjAgg, tapply(med, Func, median))
with(allSubjAgg, tapply(med, Func, max))
max(allSubjAgg$med)


pdf("AllSubjRTs.pdf", width=11, height=8)
for (s in split(allData, allData$LunaID)) {
    g <- ggplot(s, aes(x=TrialRel, y=RT)) + geom_line() + facet_grid(Emotion ~ Func) + ggtitle(s$LunaID[1L])
    print(g)
}
dev.off()

pdf("AllSubjDensity.pdf", width=11, height=8)
for (s in split(allData, allData$LunaID)) {
  g <- ggplot(s, aes(x=RT)) + geom_density() + facet_grid(Emotion ~ Func) + ggtitle(s$LunaID[1L])
  print(g)
}
dev.off()


pdf("AllSubjRTs_withMax.pdf", width=11, height=8)
for (s in split(allData, allData$LunaID)) {
    sm <- reshape2::melt(s[,c("TrialRel", "RT", "Emotion", "Func", "bestRewRT", "bestEVRT")], id.vars=c("TrialRel", "Emotion", "Func"))
    g <- ggplot(sm, aes(x=TrialRel, y=value, color=variable)) + geom_line() + facet_grid(Emotion ~ Func) + ggtitle(s$LunaID[1L]) + ylab("RT") + xlab("Trial")
    print(g)
}
dev.off()



pdf("AllSubjRTHist.pdf", width=11, height=8)
ggplot(allData, aes(x=RT)) + geom_histogram(binwidth=250) + facet_grid(Emotion ~ Func) + ggtitle("All Subjects")
for (s in split(allData, allData$LunaID)) {
  g <- ggplot(s, aes(x=RT)) + geom_histogram(binwidth=250) + facet_grid(Emotion ~ Func) + ggtitle(s$LunaID[1L])
  print(g)
}
dev.off()

#Look at expected value as a function of emotion and reward contingency.
#This is pure expected value (i.e., for a given RT, just frequency * magnitude) and is not based on subjects' reward history.
#Consequently, this is our best estimate of "good performance"
evMixedAll <- lmer(EV ~ Func*Emotion + (1|LunaID), allData)
evMixedH1 <- lmer(EV ~ Func*Emotion + (1|LunaID), subset(allData, TrialRel <= 21))
evMixedH2 <- lmer(EV ~ Func*Emotion + (1|LunaID), subset(allData, TrialRel > 21))
cmAll <- lmerCellMeans(evMixedAll)
cmH1 <- lmerCellMeans(evMixedH1)
cmH2 <- lmerCellMeans(evMixedH2)

pdf("Average EV.pdf", width=9, height=7)
ggplot(cmAll, aes(x=Func, y=EV, color=Emotion, group=Emotion, ymin=plo, ymax=phi)) + geom_point(size=5) + geom_errorbar(width=0.5) + ggtitle("All Trials")
ggplot(cmH1, aes(x=Func, y=EV, color=Emotion, group=Emotion, ymin=plo, ymax=phi)) + geom_point(size=5) + geom_errorbar(width=0.5) + ggtitle("First Half")
ggplot(cmH2, aes(x=Func, y=EV, color=Emotion, group=Emotion, ymin=plo, ymax=phi)) + geom_point(size=5) + geom_errorbar(width=0.5) + ggtitle("Second Half")
dev.off()

#look at how EV is modulated by Age
evMixedAllAge <- lmer(EV ~ Func*Emotion*AgeAtVisit + (1|LunaID), allData)
evMixedAllAgeHalf <- lmer(EV ~ Func*Emotion*AgeAtVisit*Half + (1|LunaID), allData)
evMixedH1Age <- lmer(EV ~ Func*Emotion*AgeAtVisit + (1|LunaID), subset(allData, TrialRel <= 21))
evMixedH2Age <- lmer(EV ~ Func*Emotion*AgeAtVisit + (1|LunaID), subset(allData, TrialRel > 21))
cmAllAge <- lmerCellMeans(evMixedAllAge, n.cont=10)
cmAllAgeHalf <- lmerCellMeans(evMixedAllAgeHalf, n.cont=10)
cmH1Age <- lmerCellMeans(evMixedH1Age, n.cont=10)
cmH2Age <- lmerCellMeans(evMixedH2Age, n.cont=10)


pdf("Average EV with Age.pdf", width=9, height=7)
ggplot(cmAllAgeHalf, aes(x=AgeAtVisit, y=EV, color=Emotion, ymin=plo, ymax=phi)) + facet_grid(Func~Half, scales="free_y") + geom_point(size=5) + geom_errorbar(width=0.5) + geom_line(size=1) + ggtitle("All Trials")
ggplot(cmAllAge, aes(x=AgeAtVisit, y=EV, color=Func, ymin=plo, ymax=phi)) + facet_wrap(~Emotion) + geom_point(size=5) + geom_errorbar(width=0.5) + geom_line(size=1) + ggtitle("All Trials")
ggplot(cmH1Age, aes(x=AgeAtVisit, y=EV, color=Func, ymin=plo, ymax=phi)) + facet_wrap(~Emotion) + geom_point(size=5) + geom_errorbar(width=0.5) + geom_line(size=1) + ggtitle("First Half")
ggplot(cmH2Age, aes(x=AgeAtVisit, y=EV, color=Func, ymin=plo, ymax=phi)) + facet_wrap(~Emotion) + geom_point(size=5) + geom_errorbar(width=0.5) + geom_line(size=1) + ggtitle("Second Half")
dev.off()

allData <- ddply(allData, .(LunaID, Block), function(subdf) {
      #for complete safety, re-sort by trial so that smooth is proper
      subdf <- subdf[order(subdf$TrialRel),]
      subdf$RTSpline <- smooth.spline(x=subdf$TrialRel, y=subdf$RT)$y
      subdf$RTLoess0p2 <- lowess(x=subdf$TrialRel, y=subdf$RT, f = 0.2)$y
      subdf$EVSpline <- smooth.spline(x=subdf$TrialRel, y=subdf$EV)$y
      subdf$EVLoess0p2 <- lowess(x=subdf$TrialRel, y=subdf$EV, f = 0.2)$y
      return(subdf)
    })

#look at trial-by-trial changes in smoothed EV as a function of emotion, age, and contingency
#treat trial as factor for a moment just to get the cell means by trial (not assuming a smooth interaction
#with age)
allData$TrialRelFac <- factor(allData$TrialRel) 
evMixedTrial <- lmer(EVSpline ~ Func*Emotion*TrialRelFac*AgeAtVisit + (1|LunaID), allData)
cmMixedTrial <- lmerCellMeans(evMixedTrial, n.cont=10, divide="AgeAtVisit")

pdf("Estimated EV over Trials by Age, Function, and Emotion.pdf", width=14, height=11)
ggplot(cmMixedTrial, aes(x=TrialRelFac, y=EVSpline, color=AgeAtVisit, ymin=plo, ymax=phi)) +
    geom_point(size=4) +
    geom_errorbar(width=0.25) +
    facet_grid(Func~Emotion)
dev.off()
    
#similar idea for RT     
rtMixedTrial <- lmer(RTSpline ~ Func*Emotion*TrialRelFac*AgeAtVisit + (1|LunaID), allData)
cmMixedTrial <- lmerCellMeans(rtMixedTrial, n.cont=10, divide="AgeAtVisit")

pdf("Estimated RT over Trials by Age, Function, and Emotion.pdf", width=14, height=11)
ggplot(cmMixedTrial, aes(x=TrialRelFac, y=RTSpline, color=AgeAtVisit, ymin=plo, ymax=phi)) +
    geom_point(size=4) +
    geom_errorbar(width=0.25) +
    facet_grid(Func~Emotion) +
    ggtitle("Average RT over trials by Age, Function and Emotion")
dev.off()


##look at first-half versus last-half RTs for each condition per LunaID
RTagg <- ddply(allData, .(LunaID, Func, Emotion), function(subdf) {
      firstHalf <- subset(subdf, TrialRel <= 0.5*floor(nrow(subdf)))
      lastHalf <- subset(subdf, TrialRel > 0.5*floor(nrow(subdf)))
      h1 <- data.frame(mRT=mean(firstHalf$RT, na.rm=TRUE),
          mRTSpline=mean(firstHalf$RTSpline, na.rm=TRUE),
          mRTLoess0p2=mean(firstHalf$RTLoess0p2, na.rm=TRUE), half=factor("1"))
      
      h2 <- data.frame(mRT=mean(lastHalf$RT, na.rm=TRUE),
          mRTSpline=mean(lastHalf$RTSpline, na.rm=TRUE),
          mRTLoess0p2=mean(lastHalf$RTLoess0p2, na.rm=TRUE), half=factor("2"))              

      agg <- rbind(h1, h2)
      return(agg)  ##subsetting factors automatically added to data.frame
    })

for (s in split(RTagg, list(RTagg$Emotion, RTagg$Func))) {
    cat("Func: ", as.character(s$Func[1L]), ", Emo: ", as.character(s$Emotion[1L]), "\n")
    print(t.test(mRT ~ half, s, paired=TRUE))
}

for (s in split(RTagg, list(RTagg$Emotion, RTagg$Func))) {
  cat("Func: ", as.character(s$Func[1L]), ", Emo: ", as.character(s$Emotion[1L]), "\n")
  print(t.test(mRTSpline ~ half, s, paired=TRUE))
}

for (s in split(RTagg, list(RTagg$Emotion, RTagg$Func))) {
  cat("Func: ", as.character(s$Func[1L]), ", Emo: ", as.character(s$Emotion[1L]), "\n")
  print(t.test(mRTSpline ~ half, s, paired=TRUE))
}

aggMelt <- melt(RTagg, id.vars=c("LunaID", "Func", "Emotion", "half"))
png("Split block RT averages.png", width=6, height=6, units="in", res=300)
ggplot(aggMelt, aes(x=half, y=value)) + geom_boxplot() + facet_grid(variable ~ Func*Emotion) + ylab("Average RT") + xlab("1st half or 2nd half of block")
dev.off()

library(tables)

tabular(Emotion*Func*half ~ mRT*(mean+sd+min+max), data=RTagg)
tabular(Emotion*Func*half ~ mRTSpline*(mean+sd+min+max), data=RTagg)
tabular(Emotion*Func*half ~ mRTLoess0p2*(mean+sd+min+max), data=RTagg)

library(gplots)

library(lme4)

RTagg <- plyr::rename(RTagg, c(LunaID="LunaID"))
RTagg <- merge(RTagg, behav[,c("LunaID", "AgeAtVisit", "ADI_Emotional", "ADI_Behavioral", "ADI_Cognitive", 
            "UPPS_Urg", "UPPS_Prem", "UPPS_Pers", "UPPS_SS", "UPPS_PosUrg", "UPPS_Total", "RIST.INDEZ")], by="LunaID")

summary(lmer(mRTSpline ~ Emotion*Func*half*ADI_Emotional + (1 | LunaID), RTagg))

rtModel <- lmer(RTSpline ~ Emotion*Func*TrialRel + (1 + TrialRel | LunaID), allData)
anova(rtModel)
rtPred <- lmerCellMeans(rtModel)

ggplot(rtPred, aes(x=TrialRel, y=RTSpline, color=Func)) + facet_wrap(~Emotion) + geom_point() + geom_line()

png("Descriptives of RTs by Half.png", width=8, height=4, units="in", res=300)
textplot(tabular(Emotion*Func ~ (RTh1 + RTh2)*(mean+sd+min+max), data=RTagg), show.rownames = FALSE, show.colnames = FALSE)
dev.off()

#look at CEV vs. CEVR probability-magnitude tradeoff
cevH1 <- gdata::drop.levels(subset(RTagg, Func=="CEV" & half=="1"))
cevrH1 <- gdata::drop.levels(subset(RTagg, Func=="CEVR" & half=="1"))

pmBias <- data.frame(cevH1[,c("LunaID", "Emotion", "AgeAtVisit")], pmTradeoff=cevrH1$mRT - cevH1$mRT)
tapply(pmBias$pmTradeoff, pmBias$Emotion, mean)

summary(aov(lm(pmTradeoff ~ Emotion, pmBias)))
summary(lm(pmTradeoff ~ Emotion, pmBias))

summary(lm(pmTradeoff ~ Emotion*AgeAtVisit, pmBias))
library(car)
Anova(lm(pmTradeoff ~ Emotion*AgeAtVisit, pmBias))

