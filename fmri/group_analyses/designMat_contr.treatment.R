library(lsmeans)
nSubs <- 37
group <- factor(c(rep("BPD", 17), rep("HC",20))) #17 BPD, 21 Control
subject <- factor(paste0("sub", sprintf("%02d", 1:nSubs)))
run <- 1:4
so <- c("self", "other")
ft <- c("first", "third")
design <- expand.grid(subject=subject, run=run, so=so, ft=ft)

design$cope <- apply(design[,c("so", "ft")], 1, function(r) {
  if ( r["so"] == "self" && r["ft"] == "first") {
    cope <- 1
  } else if ( r["so"] == "self" && r["ft"] == "third") {
    cope <- 3
  } else if ( r["so"] == "other" && r["ft"] == "first") {
    cope <- 2
  } else if ( r["so"] == "other" && r["ft"] == "third") {
    cope <- 4
  } else { stop("Unmatched") }
  
  return(cope)
})

design$group <- group #add group status factor

#dummy DV to get model.matrix
design$dummy <- rnorm(nrow(design), 0, 1)

##lookup file name
#setwd("~/Desktop")
runDirs <- scan("SPBP_RunDirs.txt", what="char", sep="\n")
group <- sub("^.*/data/\\d+_[^_]+_(\\w+).*$", "\\1", runDirs, perl=TRUE)
runDirs <- runDirs[order(group)] #BPD first, control second
group <- group[order(group)] #BPD first, control second
#fMRIID <- sub("^.*/data/(\\d+).*$", "\\1", runDirs, perl=TRUE)
#groupNum <- sapply(group, function(x) { ifelse(x=="HC", -1, 1)})
subject <- factor(paste0("sub", sprintf("%02d", 1:length(runDirs))))

fileDF <- data.frame(runDirs, subject) #group=groupNum, fMRIID


#for testing the merge
#copeRun <- expand.grid(subject=1:37, run=1:4, cope=1:4)
#bringTogether <- merge(copeRun, fileDF, by="subject")

#the real thing
#N.B.: The file paths must be sorted to be BPD first for the ascending subject ids to match as expected
#design$subject <- factor(design$subject)
#fileDF$subject <- factor(fileDF$subject)
designCombine <- merge(design, fileDF, by="subject")
designCombine <- designCombine[order(designCombine$subject, designCombine$run, designCombine$cope),]
designCombine$fullPath <- with(designCombine, paste0(runDirs, "/spbp_r", run, ".feat/stats/cope", cope, ".nii.gz"))
write.table(designCombine$fullPath, file="copeList.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#include interations at possible expense of ME interpretability?
designmat.int <- lm(dummy ~ group*so*ft + subject, data=designCombine)

mm <- model.matrix(designmat.int)
mm <- mm[,-which(dimnames(mm)[[2]] == "subjectsub37")]

#correlation among effects
#cmat <- cor(mm)
options(contrasts=c("contr.helmert", "contr.poly"))
options(contrasts=c("contr.treatment", "contr.poly"))

cor(mm[,c("group1", "so1", "ft1", "group1:so1", "group1:ft1", "so1:ft1", "group1:so1:ft1")])

cor(mm[,c("groupHC", "soother", "ftthird", "groupHC:soother", "groupHC:ftthird", "soother:ftthird", "groupHC:soother:ftthird")])

write.table(mm, file="fslDesign_RcontrTx_int.mat", row.names=FALSE, col.names=FALSE, quote=FALSE)
dimnames(mm)[2]

#adjusted group means
l <- lsmeans(designmat.int, ~ group)
#manually remove the subject37 dummy code, which is linearly dependent
#should just remove it without tweaking contrast code for subjects, since it should remain 1/37.
linfct <- l@linfct
linfct <- linfct[,-1*which(colnames(linfct)=="subjectsub37")]
write.table(file="groupME_adjMeans_RcontrTx_int.txt", l@linfct, row.names=FALSE, col.names=FALSE, quote=FALSE)

l <- lsmeans(designmat.int, ~ so)
linfct <- l@linfct
linfct <- linfct[,-1*which(colnames(linfct)=="subjectsub37")]
write.table(file="soME_adjMeans_RcontrTx_int.txt", linfct, row.names=FALSE, col.names=FALSE, quote=FALSE)

l <- lsmeans(designmat.int, ~ ft)
linfct <- l@linfct
linfct <- linfct[,-1*which(colnames(linfct)=="subjectsub37")]
write.table(file="ftME_adjMeans_RcontrTx_int.txt", linfct, row.names=FALSE, col.names=FALSE, quote=FALSE)

#lsmeans(designmat.int, ~ group | so*ft )
#with(designCombine, tapply(dummy, group, mean))

#additive model only
#generate a big ME of group
#for testing that contrast code is basically correct
designCombine$dummy[which(designCombine$group=="BPD")] <- designCombine$dummy[which(designCombine$group=="BPD")] + rnorm(nrow(subset(designCombine, group=="BPD")), 10, 1)

#the contr.treatment setup here generates a rank deficient design
#this occurs because the contrast for HC vs. BPD must have a reference condition wrt the subject effect.
#if we drop the subject37 dummy, I believe we implicitly set subject37 as the reference for the group ME?
#anyhow, this is totally sensible.
#need to drop last subject dummy code from model matrix and figure out the lsmeans from that.
designmat.add <- lm(dummy ~ group + so + ft + subject, data=designCombine, singular.ok=TRUE)
mm <- model.matrix(designmat.add)
mm <- mm[,-which(dimnames(mm)[[2]] == "subjectsub37")]


#last column is linearly dependent, and is dropped by lm
#same results produced when manually dropping column 40
#test <- lm(designCombine$dummy~-1+mm[,-40])
#summary(test)

#buttttt, we still have to get the lsmeans call right. :)

write.table(mm, file="fslDesign_RcontrTx_add.mat", row.names=FALSE, col.names=FALSE, quote=FALSE)
dimnames(mm)[2]

#check for which column is most squirrelly
#rankifremoved <- sapply(1:ncol(mm), function (x) qr(mm[,-x])$rank)
#which(rankifremoved == max(rankifremoved))

#adjusted group means
l <- lsmeans(designmat.add, ~ group)

#tweak
## ltweak <- l@linfct
## ltweak <- ltweak[,-40]

## designmat.addTweak <- lm(designCombine$dummy~-1+mm[,-40])
## mult <- ltweak %*%  matrix(data=coef(designmat.addTweak), ncol=1)
## g1 <-  sum(coef(designmat.addTweak)*ltweak[1,])
## g2 <-  sum(coef(designmat.addTweak)*ltweak[2,])
## sum(mult[1,])
## sum(mult[2,])

## with(designCombine, tapply(dummy, group, mean))

linfct <- l@linfct
linfct <- linfct[,-1*which(colnames(linfct)=="subjectsub37")]
write.table(file="groupME_adjMeans_RcontrTx_add.txt", linfct, row.names=FALSE, col.names=FALSE, quote=FALSE)

l <- lsmeans(designmat.add, ~ so)
linfct <- l@linfct
linfct <- linfct[,-1*which(colnames(linfct)=="subjectsub37")]
write.table(file="soME_adjMeans_RcontrTx_add.txt", linfct, row.names=FALSE, col.names=FALSE, quote=FALSE)

l <- lsmeans(designmat.add, ~ ft)
linfct <- l@linfct
linfct <- linfct[,-1*which(colnames(linfct)=="subjectsub37")]
write.table(file="ftME_adjMeans_RcontrTx_add.txt", linfct, row.names=FALSE, col.names=FALSE, quote=FALSE)
