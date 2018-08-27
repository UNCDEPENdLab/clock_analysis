#1) load spatial maps RData object
#2) rebuild into 4d cube (where fourth dimension is run/subject)
#3) clusterize each effect of interest in stats outputs using 3dclust and generating mask
library(ggplot2)
library(abind)
library(oro.nifti)
library(reshape2)
library(robust)
library(car)
library(dplyr)

fslgroupdir <- "/gpfs/group/mnh5174/default/MMClock/fsl_sceptic_group_jun2017"

#cope 6 is relative uncertainty in LVL1 runs
l2dirs <- data.frame(dir=dirname(readLines(file.path(fslgroupdir, "sceptic_vchosen_ventropy_dauc_pemax_preconvolve_cope1_inputs")))) #use dirname to pull off cope
if (as.character(l2dirs[nrow(l2dirs), "dir"]) == "") { l2dirs <- l2dirs[-nrow(l2dirs),,drop=FALSE] } #there is an empty trailing row sometimes
l2dirs$ID <- as.numeric(sub(".*/MR_Proc/(\\d{5})_.*", "\\1", l2dirs$dir, perl=TRUE))

#master cov list
#n73_covs <- read.table("/Volumes/Serena/MMClock/fsl_group/groupcov_ageexplore_n73.txt", header=TRUE)
##n73_covs <- read.table("/Users/michael/TresorSync/fmri/fsl_group/groupcov_ageexplore_n73.txt", header=TRUE)

subinfo <- read.table("/gpfs/group/mnh5174/default/clock_analysis/fmri/subinfo_db", header=TRUE)
subinfo <- subinfo %>% rename(ID=lunaid)

l2dirs <- merge(l2dirs, subinfo, by="ID", all.x=TRUE)
l2dirs$age.c <- l2dirs$age - mean(l2dirs$age, na.rm=TRUE)
l2dirs$female.c <- l2dirs$female - mean(l2dirs$female, na.rm=TRUE)
l2dirs$age_female <- l2dirs$age.c * l2dirs$female.c #numerical interaction
l2dirs$female_fac <- factor(l2dirs$female, levels=c(0,1), labels=c("Male", "Female"))

#list of LVL2 copes
#cope1: m_scram
#cope2: m_fear
#cope3: m_happy
#cope4: m_overall
#cope5: fear > scram
#cope6: happy > scram
#cope7: fear > happy
#cope8: m_run (linear run effect)

source(file.path(getMainDir(), "clock_analysis", "fmri", "glm_helper_functions.R"))
afnidir <- "/opt/aci/sw/afni/17.0.02/bin"
setwd(file.path(getMainDir(), "clock_analysis/fmri/group_analyses/sceptic_cluster_plots_fsl_Jun2017"))

#for explorer x age interaction group map, generate clusters using 3dclust, then extract average betas within clusters to generate plots
#use 3dclustsim 2-sided threshold with voxelwise p = .005 and NN=1
# 3dClustSim -mask /Volumes/Serena/MMClock/fsl_group/clock_onset.gfeat/mask.nii.gz -fwhmxyz 5.857861 6.010121 5.735949 -pthr .02 .01 .005 .001
# 2-sided thresholding
# Grid: 84x100x84 2.30x2.30x2.30 mm^3 (141450 voxels in mask)
#
# CLUSTER SIZE THRESHOLD(pthr,alpha) in Voxels
# -NN 1  | alpha = Prob(Cluster >= given size)
#  pthr  | .10000 .05000 .02000 .01000
# ------ | ------ ------ ------ ------
#0.020000    75.2   83.1   94.7  102.6
#0.010000    49.9   55.5   63.0   69.3
#0.005000    35.2   39.5   45.7   50.0
#0.001000    17.6   20.2   23.7   25.7

#so, 40 voxel cluster minimum

l1copes <- c("clock_onset", "feedback_onset", "vchosen", "ventropy", "dauc", "pemax")
l2copes <- c("m_scram", "m_fear", "m_happy", "m_overall", "fear_gt_scram", "happy_gt_scram", "fear_gt_happy", "m_run")
l3copes <- c("intercept", "female", "age", "female_x_age")

#create a multi-dimensional list to store all graphs
#then we can extract specific graphs of interest

graphlist <- vector("list", length(l1copes)*length(l2copes)*length(l3copes))
dim(graphlist) <- c(length(l1copes), length(l2copes), length(l3copes))
dimnames(graphlist) <- list(l1cope=l1copes, l2cope=l2copes, l3cope=l3copes)

zthresh <- 2.807 #2-tailed p=.005 for z stat
clustsize <- 34 #based on 3dClustSim using ACFs for first-level FEAT runs

for (l1 in 1:length(l1copes)) {
  for (l2 in 1:length(l2copes)) {
    #load relevant single-subject statistics from LVL2 directories
    #note that these do not vary by level 3 covariate (because each represents the activation for a given event (e.g., clock_onset) and multi-run contrast (e.g., fear > scram)
    copefiles <- file.path(l2dirs$dir, paste0("cope", l1, ".feat"), "stats", paste0("cope", l2, ".nii.gz"))
    
    suppressMessages(library(Rniftilib))
    imgdims <- Rniftilib::nifti.image.read(copefiles[1], read_data=0)$dim
    detach("package:Rniftilib", unload=TRUE) #necessary to avoid dim() conflict with oro.nifti
    
    #generate concatenated cope file image
    copeconcat <- array(0, dim=c(imgdims, length(copefiles)))
    for (i in 1:length(copefiles)) {
      copeconcat[,,,i] <- readNIfTI(copefiles[i], reorient=FALSE)@.Data
    }
    
    for (l3 in 1:length(l3copes)) {
      #groupmap <- file.path(fslgroupdir, paste0(l1copes[l1], ".gfeat"), paste0("cope", l2, ".feat"), "stats", paste0("zstat", l3, ".nii.gz")) #this line works if .gfeat directories are named by the contrast
      groupmap <- file.path(fslgroupdir, paste0("cope", l1, ".gfeat"), paste0("cope", l2, ".feat"), "stats", paste0("zstat", l3, ".nii.gz")) #this is for numeric naming of .gfeat dirs
      
      #gdat <- readNIfTI(groupmap, reorient=FALSE)
      #generate cluster mask
      runAFNICommand(paste0("3dclust -overwrite -1Dformat -nosum -1dindex 0 -1tindex 0",
                            " -1thresh ", zthresh, " -dxyz=1 -savemask tmpclust 1.01 ", clustsize, " ", groupmap), 
                     afnidir=afnidir, stdout="tmpclust.1D")
      
      #get coordinates and names of regions
      lookup <- runAFNICommand(paste0("whereami -coord_file tmpclust.1D'[1,2,3]' -space MNI -lpi -atlas CA_ML_18_MNIA"),
                               afnidir=afnidir, stderr="/dev/null", intern=TRUE) #/Volumes/Serena/bars_ica/output/einfomax_60_30_prenorm_28Jan2014/ic11behavLMER+tlrc.HEAD
      exitstatus <- attr(lookup, "status")  
      if (!is.null(exitstatus) && exitstatus != 0) next #whereami failed, which occurs when there are no clusters. Skip to next tbrik
      
      #get voxel sizes of clusters
      vsizes <- read.table("tmpclust.1D")$V1
      
      atlaslines <- grep("Atlas CA_ML_18_MNIA: Macro Labels (N27)", lookup, fixed=TRUE)
      bestguess <- sub("(^\\s*|\\s*$)", "", lookup[atlaslines+1], perl=TRUE) #first match after atlas for each cluster
      
      coordlines <- grep("Focus point (LPI)=", lookup, fixed=TRUE)
      coords <- lookup[coordlines+2] #first line after header is TLRC, second is MNI
      #coords <- sub("<a href=.*$", "", coords, perl=TRUE)
      coords <- sub("^\\s*(-?\\d+\\s*mm.*\\{MNI\\})\\s*<a href=.*$", "\\1", coords, perl=TRUE)
      
      #plot title for each of k ROIs
      plottitles <- paste(bestguess, paste(coords, paste(vsizes, "vox"), sep="; "), sep="\n")
      
      roimask <- readAFNI("tmpclust+tlrc.HEAD", vol=1)
      #afni masks tend to read in as 4D matrix with singleton 4th dimension. Fix this
      if (length(dim(roimask)) == 4L) {
        roimask@.Data <- roimask[,,,,drop=T]    
      }
      
      maskvals <- sort(unique(as.vector(roimask)))
      maskvals <- maskvals[!maskvals == 0]
      
      #generate a list of roi averages across subjects
      roimats <- lapply(maskvals, function(v) {
        mi <- which(roimask==v, arr.ind=TRUE)
        nsubj <- length(copefiles)
        nvox <- nrow(mi)
        mi4d <- cbind(pracma::repmat(mi, nsubj, 1), rep(1:nsubj, each=nvox))
        
        mat <- matrix(copeconcat[mi4d], nrow=nvox, ncol=nsubj) #need to manually reshape into matrix from vector
        
        #for each subject, compute huber m-estimator of location/center Winsorizing at 2SD across voxels (similar to voxel mean)
        #clusavg <- apply(mat, 2, function(x) { MASS::huber(x, k=2)$mu })
        clusavg <- apply(mat, 2, mean)
      })      
      
      curoutdir <- file.path(getwd(), l1copes[l1], l2copes[l2])
      dir.create(curoutdir, showWarnings=FALSE, recursive=TRUE)
      pdf(file.path(curoutdir, paste0(l1copes[l1], "_", l2copes[l2], "_", l3copes[l3], ".pdf")), width=11, height=8)
      clustgraphs <- vector("list", length(bestguess))
      for (k in 1:length(bestguess)) {
        df <- data.frame(l2dirs, vox=roimats[[k]])
        
        if (l3copes[l3] == "age") {            
          robcorr <- covRob(na.omit(df[,c("age", "vox")]), estim="mcd", corr=TRUE)$cov[1,2]
          corr <- cor.test(df$age, df$vox)
          plotnote <- paste0("r(", corr$parameter, ") = ", round(corr$estimate, 2), ", p = ", round(corr$p.value, 3), ", rob r = ", round(robcorr, 2))
          
          g <- ggplot(df, aes(x=age, y=vox)) + geom_point() + stat_smooth(method="loess", color="red", size=2, se=FALSE) + stat_smooth(method="lm") +
            ggtitle(plottitles[k]) + annotate("text", x = min(df$age), y = max(df$vox), label = plotnote, hjust=0) + theme_bw(base_size=18) +
            ylab(paste(l1copes[l1], l2copes[l2])) + xlab("Age (years)")
          plot(g)
        } else if (l3copes[l3] == "female") {
          m <- t.test(vox ~ female_fac, df)
          plotnote <- paste0("t(", round(m$parameter, 2), ") = ", round(m$statistic, 2), ", p = ", round(m$p.value, 4))
          g <- ggplot(df, aes(x=female_fac, y=vox)) + geom_boxplot() + 
            ggtitle(plottitles[k]) + annotate("text", x = 1, y = max(df$vox)+0.1*max(df$vox), label = plotnote, hjust=0) + theme_bw(base_size=18) +
            ylab(paste(l1copes[l1], l2copes[l2])) + xlab("Sex")
          plot(g)
        } else if (l3copes[l3] == "intercept") {
          m <- t.test(df$vox, mu=0) #test against 0
          plotnote <- paste0("t(", round(m$parameter, 2), ") = ", round(m$statistic, 2), ", p = ", round(m$p.value, 4))
          g <- ggplot(df, aes(x=factor(0), y=vox)) + geom_boxplot() + 
            ggtitle(plottitles[k]) + annotate("text", x = 0.5, y = max(df$vox)+0.1*max(df$vox), label = plotnote, hjust=0) +
            theme_bw(base_size=18) + ylab(paste(l1copes[l1], l2copes[l2])) +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank())
          plot(g)          
        } else if (l3copes[l3] == "female_x_age") {
          m <- summary(lm(vox ~ age*female, df))
          eff_b <- m$coefficients["age:female","Estimate"]
          eff_t <- m$coefficients["age:female","t value"]
          eff_p <- m$coefficients["age:female","Pr(>|t|)"]
          plotnote <- paste0("b = ", round(eff_b, 3), ", t = ", round(eff_t, 3), ", p = ", round(eff_p, 4))
          
          g <- ggplot(df, aes(x=age, y=vox, color=female_fac)) + geom_point() + stat_smooth(method="lm") +
            ggtitle(plottitles[k]) + annotate("text", x = min(df$age), y = max(df$vox), label = plotnote, hjust=0) +
            theme_bw(base_size=18) +
            ylab(paste(l1copes[l1], l2copes[l2])) + xlab("Age (years)")
          plot(g)
        }
        
        clustgraphs[[k]] <- g
      }
      dev.off()
      
      graphlist[[l1,l2,l3]] <- clustgraphs
      
    }
  }
}

save(graphlist, file="cluster_ggplot_objs_Jun2017.RData")


load(file="cluster_ggplot_objs_Jun2017.RData")
dimnames(graphlist)

#L and R IFG fear > scram age effect. L IFG is cluster 3
df <- graphlist[["ventropy", "fear_gt_scram", "age"]][[3]]$data

corr <- cor.test(df$age, df$vox)
plotnote <- paste0("r(", corr$parameter, ") = ", round(corr$estimate, 2), ", p = ", format(corr$p.value, digits=3))

pdf("ventropy_fear_gt_scram_age_increases_lifg.pdf", width=8, height=6)
g <- ggplot(df, aes(x=age, y=vox)) + geom_point(size=3) + stat_smooth(method="lm", se=FALSE, size=3, color="darkblue") +
  annotate("text", x = min(df$age), y = max(df$vox), label = plotnote, hjust=0, size=8) + #ggtitle("Age-related increases in L IPL activation to pos. RPEs") + 
  theme_bw(base_size=24) + ylab("Value entropy activity in L IFG (AU)\n") + xlab("Age (years)")
plot(g)          
dev.off()


#MPFC age increase for happy > scram. Cluster 1 is mPFC/ACC
df <- graphlist[["ventropy", "happy_gt_scram", "age"]][[1]]$data

corr <- cor.test(df$age, df$vox)
plotnote <- paste0("r(", corr$parameter, ") = ", round(corr$estimate, 2), ", p = ", format(corr$p.value, digits=3))

pdf("ventropy_happy_gt_scram_age_increases_mpfc.pdf", width=8, height=6)
g <- ggplot(df, aes(x=age, y=vox)) + geom_point(size=3) + stat_smooth(method="lm", se=FALSE, size=3, color="darkblue") +
  annotate("text", x = min(df$age), y = max(df$vox), label = plotnote, hjust=0, size=8) + #ggtitle("Age-related increases in L IPL activation to pos. RPEs") + 
  theme_bw(base_size=24) + ylab("Value entropy activity in mPFC (AU)\n") + xlab("Age (years)")
plot(g)          
dev.off()


#age-related increases in SMA for PEs
df <- graphlist[["pemax", "m_overall", "age"]][[3]]$data

corr <- cor.test(df$age, df$vox)
plotnote <- paste0("r(", corr$parameter, ") = ", round(corr$estimate, 2), ", p = ", round(corr$p.value, digits=3))

pdf("rpe_pos_age_increases_sma.pdf", width=8, height=6)
g <- ggplot(df, aes(x=age, y=vox)) + geom_point(size=3) + stat_smooth(method="lm", se=FALSE, size=3, color="darkblue") +
  annotate("text", x = min(df$age), y = max(df$vox) - .01, label = plotnote, hjust=0, size=8) + #ggtitle("Age-related increases in L IPL activation to pos. RPEs") + 
  theme_bw(base_size=24) + ylab("RPE activity in SMA (AU)\n") + xlab("Age (years)")
plot(g)          
dev.off()

#
##age-related increase for clock onset in L IFG
#df <- graphlist[["clock_onset", "happy_gt_scram", "age"]][[6]]$data
#
#corr <- cor.test(df$age, df$vox)
#plotnote <- paste0("r(", corr$parameter, ") = ", round(corr$estimate, 2), ", p = ", round(corr$p.value, 4))
#
#pdf("clock_age_increases_lifg.pdf", width=8, height=6)
#g <- ggplot(df, aes(x=age, y=vox)) + geom_point() + stat_smooth(method="lm", se=FALSE, size=3) +
#    annotate("text", x = min(df$age), y = max(df$vox), label = plotnote, hjust=0) + #ggtitle("Age-related increases in L IPL activation to pos. RPEs") + 
#    theme_bw(base_size=24) + ylab("Activation to clock in L IFG (AU)\n") + xlab("Age (years)")
#plot(g)          
#dev.off()
#
#
##age-related decrease for feedback onset in R IFG
#df <- graphlist[["feedback_onset", "happy_gt_scram", "age"]][[4]]$data
#
#corr <- cor.test(df$age, df$vox)
#plotnote <- paste0("r(", corr$parameter, ") = ", round(corr$estimate, 2), ", p = ", round(corr$p.value, 4))
#
#pdf("feedback_age_decreases_rifg.pdf", width=8, height=6)
#g <- ggplot(df, aes(x=age, y=vox)) + geom_point() + stat_smooth(method="lm", se=FALSE, size=3) +
#    annotate("text", x = min(df$age), y = max(df$vox), label = plotnote, hjust=0) + #ggtitle("Age-related increases in L IPL activation to pos. RPEs") + 
#    theme_bw(base_size=24) + ylab("Activation to feedback in R IFG (AU)\n") + xlab("Age (years)")
#plot(g)          
#dev.off()
#
#
##age x explorer interaction in L precentral/premotor cortex
#df <- graphlist[["rel_unc", "m_overall", "age_explorer"]][[1]]$data
#
##generate age-explorer graph
#m <- summary(lm(vox ~ female.c + age*explorer, df))
#int_b <- m$coefficients["age:explorerExplorer","Estimate"]
#int_t <- m$coefficients["age:explorerExplorer","t value"]
#int_p <- m$coefficients["age:explorerExplorer","Pr(>|t|)"]
#plotnote <- paste0("b = ", round(int_b, 3), ", t = ", round(int_t, 3), ", p = ", round(int_p, 5))
#
#pdf("rel_unc_age_explorer_lpremotor.pdf", width=8, height=6)
#g <- ggplot(df, aes(x=age, y=vox, color=explorer)) + geom_point() + stat_smooth(method="lm", se=FALSE, size=3) +
#    annotate("text", x = min(df$age), y = 10, label = plotnote, hjust=0) + ylim(-10, 12) +
#    theme_bw(base_size=23) + ylab("Rel. unc. activation in L premotor (AU)\n") + xlab("Age (years)") +
#    scale_color_brewer("", palette="Set2")
#plot(g)
#dev.off()
#
