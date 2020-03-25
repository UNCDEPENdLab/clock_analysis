
##look at rpe- fear > happy etc. in BPD versus control. Lots of action in the amygdala

setwd("/Users/michael/ics/SPECC/specc_3dmema/all_results")

library(oro.nifti)
cmask <- readAFNI("fear_gt_happy_bpd_control_clusters_p005+tlrc.HEAD")

maskvals <- sort(unique(as.vector(cmask@.Data)))
maskvals <- maskvals[maskvals != 0] #zero not of interest

glmfiles <- list.files(path="../glm_stat_files", full.names=TRUE, pattern=".*\\.HEAD")

maskinterest <- 10 #amygdala is cluster number 10

#38 voxels
mi <- which(cmask==maskinterest, arr.ind=TRUE)

briksinterest <- c("m_rpe_neg_scram#0_Coef", "m_rpe_neg_happy#0_Coef", "m_rpe_neg_fear#0_Coef")

statmat <- matrix(NA_real_, nrow=length(glmfiles), ncol=length(briksinterest))

for (i in 1:length(glmfiles)) {
    glm_stats <- readAFNI(glmfiles[i])
    briknames <- strsplit(glm_stats@BRICK_LABS, "~")[[1]]
    
    for (j in 1:length(briksinterest)) {
        dim4 <- which(briknames == briksinterest[j])
        this_mi <- cbind(mi[,1:3], dim4) #cbind the 3d mask indices with this sub-brik
        statmat[i,j] <- mean(glm_stats[this_mi])
    }
}

save(statmat, file="rpe_neg_conditionmeans.RData")
load(file=file.path(getMainDir(), "clock_analysis", "fmri", "rpe_neg_conditionmeans.RData"))
subid <- sub("^.*/glm_stat_files/([^_]+)_.*$", "\\1", glmfiles, perl=TRUE)
df <- data.frame(subid, statmat)
names(df) <- c("subid", "Scrambled", "Happy", "Fear")
df$bpd <- factor(as.numeric(grepl("^\\d{3}[a-z]{2}", df$subid, perl=TRUE)), levels=c(0,1), labels=c("Control", "BPD"))

library(reshape2)
library(ggplot2)
library(plyr)
m <- melt(df, id.vars=c("subid", "bpd"))

ggplot(m, aes(x=variable, y=value, color=bpd)) + stat_summary(fun.data= "mean_cl_boot", geom = "bar", width = 0.3, position="dodge")

library(ez)
a <- ezANOVA(data=m, dv=value, within=variable, wid=subid, between=bpd)
a
#library(multcomp)
#glht(a, linfct = mcp(variable="Tukey"))

aggm <- ddply(m, .(variable, bpd), function(subdf) {
    return(data.frame(m=mean(subdf$value), se=plotrix::std.error(subdf$value)))
})

pdf("rpe_neg_emotion.pdf", width=8, height=5)
negplot <- ggplot(aggm, aes(x=variable, y=m, fill=bpd)) + geom_bar(position=position_dodge(width=0.8), stat="identity", width=0.8) + geom_errorbar(aes(ymin=m-se, ymax=m+se), position=position_dodge(width=0.8), width=0.3) +
    ylab("R Amy. Neg. RPE (AU)\n") + xlab("\nEmotion") + scale_fill_brewer("Group", palette="Set1") + theme_bw(base_size=24)
plot(negplot)
dev.off()

t.test(value ~ bpd, subset(m, variable=="Scrambled"))
t.test(value ~ bpd, subset(m, variable=="Happy"))
t.test(value ~ bpd, subset(m, variable=="Fear"))



##compare against RPE +
briksinterest <- c("m_rpe_pos_scram#0_Coef", "m_rpe_pos_happy#0_Coef", "m_rpe_pos_fear#0_Coef")

statmat_pos <- matrix(NA_real_, nrow=length(glmfiles), ncol=length(briksinterest))

for (i in 1:length(glmfiles)) {
    glm_stats <- readAFNI(glmfiles[i])
    briknames <- strsplit(glm_stats@BRICK_LABS, "~")[[1]]
    
    for (j in 1:length(briksinterest)) {
        dim4 <- which(briknames == briksinterest[j])
        this_mi <- cbind(mi[,1:3], dim4) #cbind the 3d mask indices with this sub-brik
        statmat_pos[i,j] <- mean(glm_stats[this_mi])
    }
}

save(statmat_pos, file="rpe_pos_conditionmeans.RData")

load(file=file.path(getMainDir(), "clock_analysis", "fmri", "rpe_pos_conditionmeans.RData"))

subid <- sub("^.*/glm_stat_files/([^_]+)_.*$", "\\1", glmfiles, perl=TRUE)
df_pos <- data.frame(subid, statmat_pos)
names(df_pos) <- c("subid", "Scrambled", "Happy", "Fear")
df_pos$bpd <- factor(as.numeric(grepl("^\\d{3}[a-z]{2}", df$subid, perl=TRUE)), levels=c(0,1), labels=c("Control", "BPD"))


library(reshape2)
library(ggplot2)
library(plyr)
m_pos <- melt(df_pos, id.vars=c("subid", "bpd"))

a_pos <- ezANOVA(data=m_pos, dv=value, within=variable, wid=subid, between=bpd)
a_pos

aggm_pos <- ddply(m_pos, .(variable, bpd), function(subdf) {
    return(data.frame(m=mean(subdf$value), se=plotrix::std.error(subdf$value)))
})

pdf("rpe_pos_emotion.pdf", width=8, height=5)
ggplot(aggm_pos, aes(x=variable, y=m, fill=bpd)) + geom_bar(position=position_dodge(width=0.8), stat="identity", width=0.8) + geom_errorbar(aes(ymin=m-se, ymax=m+se), position=position_dodge(width=0.8), width=0.3) +
    ylab("R Amy. Pos. RPE (AU)\n") + xlab("\nEmotion") + scale_fill_brewer("Group", palette="Set1") + theme_bw(base_size=24) + coord_cartesian(ylim = c(-0.002318872,  0.002807108)) #identical to neg emo
dev.off()


##So, this is specific to RPE- in amygdala


##
##similar effects in amygdala for relative uncertainty
##clust 13 is left amygdala
##clust 3 is R hippocampus

cmask_unc <- readAFNI("rel_uncertainty_fear_gt_happy_bpd_control_clusters_p005+tlrc.HEAD")

maskvals <- sort(unique(as.vector(cmask_unc@.Data)))
maskvals <- maskvals[maskvals != 0] #zero not of interest

glmfiles <- list.files(path="../glm_stat_files", full.names=TRUE, pattern=".*\\.HEAD")

maskinterest <- 3 #r hippocampus

#38 voxels
mi <- which(cmask_unc==maskinterest, arr.ind=TRUE)

briksinterest <- c("m_rel_uncertainty_scram#0_Coef", "m_rel_uncertainty_happy#0_Coef", "m_rel_uncertainty_fear#0_Coef")

statmat_unc <- matrix(NA_real_, nrow=length(glmfiles), ncol=length(briksinterest))

for (i in 1:length(glmfiles)) {
    glm_stats <- readAFNI(glmfiles[i])
    briknames <- strsplit(glm_stats@BRICK_LABS, "~")[[1]]
    
    for (j in 1:length(briksinterest)) {
        dim4 <- which(briknames == briksinterest[j])
        this_mi <- cbind(mi[,1:3], dim4) #cbind the 3d mask indices with this sub-brik
        statmat_unc[i,j] <- mean(glm_stats[this_mi])
    }
}


save(statmat_unc, file="rel_unc_conditionmeans.RData")


subid <- sub("^.*/glm_stat_files/([^_]+)_.*$", "\\1", glmfiles, perl=TRUE)
df_unc <- data.frame(subid, statmat_unc)
names(df_unc) <- c("subid", "Scrambled", "Happy", "Fear")
df_unc$bpd <- factor(as.numeric(grepl("^\\d{3}[a-z]{2}", df_unc$subid, perl=TRUE)), levels=c(0,1), labels=c("Control", "BPD"))


library(reshape2)
library(ggplot2)
library(plyr)
m_unc <- melt(df_unc, id.vars=c("subid", "bpd"))

a_unc <- ezANOVA(data=m_unc, dv=value, within=variable, wid=subid, between=bpd)
a_unc

aggm_unc <- ddply(m_unc, .(variable, bpd), function(subdf) {
    return(data.frame(m=mean(subdf$value), se=plotrix::std.error(subdf$value)))
})

pdf("rel_unc_emotion_rhipp.pdf", width=8, height=5)
ggplot(aggm_unc, aes(x=variable, y=m, fill=bpd)) + geom_bar(position=position_dodge(width=0.8), stat="identity", width=0.8) + geom_errorbar(aes(ymin=m-se, ymax=m+se), position=position_dodge(width=0.8), width=0.3) +
    ylab("R Hipp. Rel. Unc. (AU)\n") + xlab("\nEmotion") + scale_fill_brewer("Group", palette="Set1") + theme_bw(base_size=24)
dev.off()


##cluster 13 is l amygdala

maskinterest <- 13

mi <- which(cmask_unc==maskinterest, arr.ind=TRUE)

briksinterest <- c("m_rel_uncertainty_scram#0_Coef", "m_rel_uncertainty_happy#0_Coef", "m_rel_uncertainty_fear#0_Coef")

statmat_unc_lamy <- matrix(NA_real_, nrow=length(glmfiles), ncol=length(briksinterest))

for (i in 1:length(glmfiles)) {
    glm_stats <- readAFNI(glmfiles[i])
    briknames <- strsplit(glm_stats@BRICK_LABS, "~")[[1]]
    
    for (j in 1:length(briksinterest)) {
        dim4 <- which(briknames == briksinterest[j])
        this_mi <- cbind(mi[,1:3], dim4) #cbind the 3d mask indices with this sub-brik
        statmat_unc_lamy[i,j] <- mean(glm_stats[this_mi])
    }
}


save(statmat_unc_lamy, file="rel_unc_conditionmeans_lamy.RData")

subid <- sub("^.*/glm_stat_files/([^_]+)_.*$", "\\1", glmfiles, perl=TRUE)
df_unc_lamy <- data.frame(subid, statmat_unc_lamy)
names(df_unc_lamy) <- c("subid", "Scrambled", "Happy", "Fear")
df_unc_lamy$bpd <- factor(as.numeric(grepl("^\\d{3}[a-z]{2}", df_unc_lamy$subid, perl=TRUE)), levels=c(0,1), labels=c("Control", "BPD"))


library(reshape2)
library(ggplot2)
library(plyr)
m_unc_lamy <- melt(df_unc_lamy, id.vars=c("subid", "bpd"))

a_unc_lamy <- ezANOVA(data=m_unc_lamy, dv=value, within=variable, wid=subid, between=bpd)
a_unc_lamy




aggm_unc_lamy <- ddply(m_unc_lamy, .(variable, bpd), function(subdf) {
    return(data.frame(m=mean(subdf$value), se=plotrix::std.error(subdf$value)))
})

pdf("rel_unc_emotion_lamy.pdf", width=8, height=5)
ggplot(aggm_unc_lamy, aes(x=variable, y=m, fill=bpd)) + geom_bar(position=position_dodge(width=0.8), stat="identity", width=0.8) + geom_errorbar(aes(ymin=m-se, ymax=m+se), position=position_dodge(width=0.8), width=0.3) +
    ylab("L Amy. Rel. Unc. (AU)\n") + xlab("\nEmotion") + scale_fill_brewer("Group", palette="Set1") + theme_bw(base_size=24)
dev.off()


