####
##SCEPTIC L3 ANALYSIS
##models <- c("sceptic_dauc_preconvolve", "sceptic_pemax_preconvolve", "sceptic_vmax_preconvolve", "sceptic_vchosen_preconvolve", "sceptic_ventropy_preconvolve")
library(dplyr)
options(dplyr.width=200)
#preproc_dirname <- "mni_5mm_aroma"
preproc_dirname <- "mni_5mm_3ddespike"
fsl_dirname <- "FEAT_LVL2_runtrend.gfeat"
##models <- c("sceptic_vchosen_ventropy_decay_matlab_dauc_pemax_preconvolve")
models <- c("sceptic_vchosen_ventropy_dauc_pemax_preconvolve")

#dataset <- "MMY3"
dataset <- "SPECC"

if (dataset=="MMY3") {
  subinfo <- read.table("/gpfs/group/mnh5174/default/clock_analysis/fmri/subinfo_db", header=TRUE)
  subinfo <- subinfo %>% rename(ID=lunaid, Age=age, Female=female, ScanDate=scandate)
  subinfo$mr_dir <- with(subinfo, paste0("/gpfs/group/mnh5174/default/MMClock/MR_Proc/", ID, "_", format((as.Date(ScanDate, format="%Y-%m-%d")), "%Y%m%d"))) #convert to Date, then reformat YYYYMMDD
  outdir <- "/gpfs/group/mnh5174/default/MMClock/fsl_sceptic_group_jun2017"
} else if (dataset=="SPECC") {
  subinfo <- read.csv("/gpfs/group/mnh5174/default/SPECC/SPECC_Participant_Info.csv", header=TRUE)
  subinfo <- subinfo %>% rowwise() %>% mutate(mr_dir=ifelse(LunaMRI==1,
    paste0("/gpfs/group/mnh5174/default/MMClock/MR_Proc/", Luna_ID, "_", format((as.Date(ScanDate, format="%Y-%m-%d")), "%Y%m%d")), #convert to Date, then reformat YYYYMMDD
    paste0("/gpfs/group/mnh5174/default/SPECC/MR_Proc/", tolower(SPECC_ID), "_", tolower(format((as.Date(ScanDate, format="%Y-%m-%d")), "%d%b%Y")))))
  subinfo <- subinfo %>% rename(ID=SPECC_ID, Age=AgeAtScan)
  outdir <- "/gpfs/group/mnh5174/default/SPECC/fsl_sceptic_group"  
}


#verify that mr_dir is present as expected
subinfo$dirfound <- file.exists(subinfo$mr_dir)
subset(subinfo, dirfound==FALSE)

dir.create(outdir, showWarnings=FALSE)
setwd(outdir)

#cope structure for preconvolve models
#1 = clock_onset
#2 = feedback_onset
#3 = regressor of interest
            
ncopes <- 6 #clock, feedback, sceptic vchosen, ventropy, dauc, pemax

fsl_dirname <- "FEAT_LVL2_runtrend.gfeat"
##models <- c("sceptic_vchosen_ventropy_decay_matlab_dauc_pemax_preconvolve")
models <- c("sceptic_vchosen_ventropy_dauc_pemax_preconvolve")

##rework using subinfo structure as the authoritative guide (rather than repeated searches)
copedf <- c()
for (m in models) {
  for (s in 1:nrow(subinfo)) {
    for (cope in 1:ncopes) {
      expectdir <- file.path(subinfo[s,"mr_dir"], preproc_dirname, m, fsl_dirname, paste0("cope", cope, ".feat"))
      if (dir.exists(expectdir)) {
        copedf <- rbind(copedf, data.frame(ID=subinfo[s,"ID"], model=m, cope=cope, fsldir=expectdir))
      } else {
        message("could not find expected directory: ", expectdir)
      }
    }
  }
}


mdf <- merge(subinfo, copedf, by="ID", all.y=TRUE)
badids <- c(11335, #low IQ, ADHD Hx, loss of consciousness
            11332, #should be excluded, but scan was terminated early due to repeated movement
            11282, #RTs at the floor for essentially all runs. Not appropriate
            11246, #huge movement and RTs at floor
            10637, #large and many movements in later runs (need to revisit to confirm)
            10662  #I think there are reconstruction problems here -- need to revisit
            ) 

mdf <- mdf %>% filter(!ID %in% badids)
mdf <- arrange(mdf, ID, model, cope)

##fsl constructs models by cope
bycope <- split(mdf, mdf$cope)

##create input files and design matrices for L3 models by cope
##separate models for inverse age from age
if (dataset=="MMY3") {
  lapply(bycope, function(cdf) {
    outname <- paste0(cdf$model[1], "_cope", cdf$cope[1])
    cdf$female.c <- cdf$Female - mean(cdf$Female)
    cdf$age.c <- cdf$Age - mean(cdf$Age)
    cdf$iage <- 1/cdf$Age * -100 #scale up age to make regressors not tiny and multiply by -1 so that larger values indicate developmental increase
    cdf$iage.c <- cdf$iage - mean(cdf$iage)
    cdf$femage <- cdf$age.c * cdf$female.c
    cdf$ifemage <- cdf$iage.c * cdf$female.c
    cat(as.character(cdf$fsldir), quote="", sep="\n", file=file.path(outdir, paste0(outname, "_inputs")))
    write.table(cbind(1, cdf[,c("female.c", "age.c", "femage")]), file=file.path(outdir, paste0(outname, "_design")), sep="\t", col.names=FALSE, row.names=FALSE)
    write.table(cbind(1, cdf[,c("female.c", "iage.c", "ifemage")]), file=file.path(outdir, paste0(outname, "_design_iage")), sep="\t", col.names=FALSE, row.names=FALSE)
  })
 
} else if (dataset=="SPECC") {

  lapply(bycope, function(cdf) {
    outname <- paste0(cdf$model[1], "_cope", cdf$cope[1])
    cdf$female.c <- cdf$Female - mean(cdf$Female)
    cdf$bpd.c <- cdf$BPD - mean(cdf$BPD)
    cdf$age.c <- cdf$Age - mean(cdf$Age)
    cdf$iage <- 1/cdf$Age * 100 #scale up age to make regressors not tiny
    cdf$iage.c <- cdf$iage - mean(cdf$iage)
    cdf$bpdage <- cdf$age.c * cdf$bpd.c
    cdf$bpdiage <- cdf$age.c * cdf$bpd.c
    cdf$femage <- cdf$age.c * cdf$female.c
    cdf$ifemage <- cdf$iage.c * cdf$female.c
    cat(as.character(cdf$fsldir), quote="", sep="\n", file=file.path(outdir, paste0(outname, "_inputs")))
    write.table(cbind(1, cdf[,c("bpd.c", "female.c", "age.c", "bpdage", "femage")]), file=file.path(outdir, paste0(outname, "_design")), sep="\t", col.names=FALSE, row.names=FALSE)
    write.table(cbind(1, cdf[,c("bpd.c", "female.c", "iage.c", "bpdiage", "ifemage")]), file=file.path(outdir, paste0(outname, "_design_iage")), sep="\t", col.names=FALSE, row.names=FALSE)
  })

}
####

#Older Approach that searches file system and aligns with subject info
procdir <- "/gpfs/group/mnh5174/default/MMClock/MR_Proc"
for (m in models) {
  for (cope in 1:ncopes) {
    outname <- paste0(m, "_cope", cope)
    subjdirs <- system(paste0("find ", procdir, " -mindepth 1 -maxdepth 1 -type d | sort -n"), intern=TRUE)
    copedirs <- c()
    for (s in subjdirs) {
      if (file.exists(checkpath <- file.path(s, preproc_dirname, m, fsl_dirname, paste0("cope", cope, ".feat")))) {
        copedirs <- c(copedirs, checkpath)
      }
    }
    
    ##copedirs <- system(paste0("find ", procdir, " -iname 'cope", cope, ".feat' -ipath '*", m, "/FEAT_LVL2.gfeat*' -type d | sort -n"), intern=TRUE)
    ##copedirs <- read.table(outname)$V1
    copedf <- data.frame(fsldir=copedirs, lunaid=as.numeric(sub("^.*/MR_Proc/(\\d{5})_\\d+/.*$", "\\1", copedirs, perl=TRUE)))
    mdf <- merge(subinfo, copedf, by="lunaid", all.y=TRUE) #should probably do a setdiff to look for discrepancies
    mdf <- filter(mdf, !lunaid %in% badids) %>% arrange(lunaid)
    ##dplyr::anti_join(subinfo, copedf, by="lunaid")
    ##dplyr::anti_join(copedf, subinfo, by="lunaid")
    mdf$female.c <- mdf$female - mean(mdf$female)
    mdf$age.c <- mdf$age - mean(mdf$age)
    mdf$agefem <- mdf$age.c * mdf$female.c
    cat(as.character(mdf$fsldir), quote="", sep="\n", file=file.path(outdir, paste0(outname, "_inputs")))
    write.table(cbind(1, mdf[,c("female.c", "age.c", "agefem")]), file=file.path(outdir, paste0(outname, "_design")), sep="\t", col.names=FALSE, row.names=FALSE)
  }
}

####
##BEGIN L3

#need to pull in age here
n73_covs <- read.table("/Volumes/Serena/MMClock/fsl_group/groupcov_ageexplore_n73.txt", header=TRUE)
outdir <- "/Volumes/Serena/MMClock/fsl_group"
dir.create(outdir, showWarnings=FALSE)

#generate cope directories
#system("find /Volumes/Serena/MMClock/MR_Proc -iname 'cope1.feat' -ipath '*fsl_tc_nomeanunc/FEAT_LVL2.gfeat*' -type d > cope1dirs")

###setup L3 models for FSL no mean uncertainty
#1 = clock_onset
#2 = feedback_onset
#3 = ev
#4 = rpe_pos
#5 = rpe_neg
#6 = rel_unc
for (cope in 1:6) {
  copedirs <- read.table(paste0("/Volumes/Serena/MMClock/MR_Proc/cope", cope, "dirs"))$V1
  copedf <- data.frame(fsldir=copedirs, lunaid=as.numeric(sub("^.*/MR_Proc/(\\d{5})_\\d+/.*$", "\\1", copedirs, perl=TRUE)))
  m <- merge(n73_covs, copedf, by="lunaid", all.y=TRUE)
  
  #mean center covs
  m$female.c <- m$female - mean(m$female)
  m$age.c <- m$age - mean(m$age)
  m$alpha_diff.c <- m$alpha_diff - mean(m$alpha_diff)
  m$exp_age.c <- m$explorer*m$age.c 
  
  #write out vector of runs and covariates to model
  cat(as.character(copedf$fsldir), quote="", sep="\n", file=file.path(outdir, paste0("cope", cope, "inputs")))
  write.table(cbind(1, m[,c("female.c", "age.c", "explorer", "exp_age.c", "alpha_diff.c")]), file=file.path(outdir, paste0("cope", cope, "covs_fem_age_exp_expage_alphadiff")), sep="\t", col.names=FALSE, row.names=FALSE)
}


##need to check mask coverage of FSL level 1 files
#these are generated (reasonably) by a -Tmin over the first level runs.
library(oro.nifti)

mask_files <- system("find /Volumes/Serena/MMClock/MR_Proc -iname 'mask.nii.gz' -type f -ipath '*fsl_tc_nomeanunc/FEAT_LVL1_run[0-9].feat*' -type f", intern=TRUE)
#mask_files <- system("find /Volumes/Serena/MMClock/MR_Proc -iname 'subject_mask.nii.gz' -type l -ipath '*mni_5mm_wavelet/clock[0-9]*'", intern=TRUE)

#pulling in reg_standard mask, too. filter out
mask_files <- mask_files[!grepl("reg_standard", mask_files)]
#group_mask <- "/Users/michael/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_mask_2.3mm.nii"

#actually, first-level runs are masked by a -thrP 10 application to the 2.3mm brain. This is slightly dilated wrt the mask file above
#recreate here
runFSLCommand("fslmaths /Users/michael/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_brain_2.3mm -thrP 10 -bin /Users/michael/templateMask_2.3 -odt char", fsldir="/usr/local/ni_tools/fsl")
group_mask <- "/Users/michael/templateMask_2.3.nii.gz"

gdat <- readNIfTI(group_mask, reorient=FALSE)
subid <- as.numeric(sub(".*/MR_Proc/(\\d{5})_\\d+/.*", "\\1", mask_files, perl=TRUE))
runnum <- as.numeric(sub(".*/MR_Proc/\\d{5}_\\d+/mni_5mm_wavelet/fsl_tc_nomeanunc/FEAT_LVL1_run([0-9])\\.feat/.*", "\\1", mask_files, perl=TRUE))
#runnum <- as.numeric(sub(".*/MR_Proc/\\d{5}_\\d+/mni_5mm_wavelet/clock([0-9])/.*", "\\1", mask_files, perl=TRUE))

#matrix of missing voxels. voxmiss is how many voxels are missing in the subject image relative to mask. Overflow is for (unlikely) possibility of more than MNI mask
voxCheck <- matrix(NA_real_, nrow=length(mask_files), ncol=4, dimnames=list(NULL, c("id", "run", "voxmiss", "overflow")))
for (m in 1:length(mask_files)) {  
  mdat <- readNIfTI(mask_files[m], reorient=FALSE)
  mdiff <- gdat - mdat
  missVox <- sum(mdiff == 1)
  weirdVox <- sum(mdiff == -1)
  voxCheck[m,] <- c(subid[m], runnum[m], missVox, weirdVox)
}

#fortunately, overflow is all 0 -- so single subject mask was consistent in preprocessing
#look at runs > 75%ile of missingness
cutoff <- quantile(voxCheck[,"voxmiss"], 0.75)
highmiss <- voxCheck[voxCheck[,"voxmiss"] > cutoff,]
highmiss[order(highmiss[,"voxmiss"], decreasing=TRUE),]
table(highmiss[,"id"])

voxCheckNew <- voxCheck
load("MissingVoxWRTTemplate_BeforeBBRFix.RData")
save(voxCheck, file="MissingVoxWRTTemplate_AfterBBRFix.RData")

diff <- cbind(voxCheck[,1:2], voxCheck[,3] - voxCheckNew[,3])

diff <- merge(as.data.frame(voxCheck), as.data.frame(voxCheckNew), by=c("id", "run"), all=TRUE)

head(diff)

diff$afterMbefore <- diff$voxmiss.y - diff$voxmiss.x

#worse
diff[diff$afterMbefore > 100,]

#better
diff[diff$afterMbefore < -100,]

#save(voxCheck, file="MissingVoxWRTTemplate_BeforeBBRFix.RData")






#soooo.... fsl is very slow and I don't think it handles the multiple sessions per subject very well (since we have dummy regressors for subject)
#what about using weighted lmer approach ala David Paulsen?

#loop over copes and varcopes and load into 4-d array (x y z subject)
library(oro.nifti)
library(doSNOW)
library(foreach)
library(iterators)

load(file="Feat_runinfo.RData")

#about 1.29 GB in RAM per array
exampCope <- readNIfTI(file.path(featL1Df$featRun[1], "stats", "cope1.nii.gz"), reorient=FALSE)@.Data
mask <- readNIfTI("/Users/michael/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_mask_2.3mm.nii", reorient=FALSE)@.Data

#generate a group cope and varcope image for all subjects
for (copenum in 1:5) {
  cope_array <- array(NA_real_, c(dim(exampCope),nrow(featL1Df)))
  varcope_array <- array(NA_real_, c(dim(exampCope),nrow(featL1Df)))
  
  for (i in 1:nrow(featL1Df)) {
    cope <- readNIfTI(file.path(featL1Df$featRun[i], "stats", paste0("cope", copenum, ".nii.gz")), reorient=FALSE)@.Data
    varcope <- readNIfTI(file.path(featL1Df$featRun[i], "stats", paste0("varcope", copenum, ".nii.gz")), reorient=FALSE)@.Data
    
    #apply mask
    cope[which(mask==0)] <- NA_real_
    varcope[which(mask==0)] <- NA_real_
    
    cope_array[,,,i] <- cope
    varcope_array[,,,i] <- varcope
    
  }
  
  save(cope_array, varcope_array, file=paste0("allsubj_cope", copenum, ".RData"))
}


#now have 4d files for copes and varcopes in same order as featL1Df
#loop over voxels that are within the mask

#only consider unmasked voxels
vindices <- which(mask==1, arr.ind=TRUE)

copenum <- 1
load(paste0("allsubj_cope", copenum, ".RData"))

#each 2-D array is ~283MB
cope_voxels <- apply(cope_array, 4, '[', vindices) #obtain 154985 x 240 matrix (voxels x sessions)
varcope_voxels <- apply(varcope_array, 4, '[', vindices) #obtain 154985 x 240 matrix (voxels x sessions)

rm(cope_array)
rm(varcope_array)

#maskData <- maskData[1:1000,] #for testing

#setDefaultClusterOptions(master="localhost", port=10290) #move away from 10187 to avoid collisions
#clusterobj <- makeSOCKcluster(8)
#registerDoSNOW(clusterobj)

#as a memory-saving hack, cbind cope and varcope matrices together so that we can iterate by row (and avoid exporting full data to each worker)
vanalyze <- cbind(cope_voxels, varcope_voxels)
rownames(vanalyze) <- 1:nrow(vanalyze) #allows for lookup of mask indices inside dopar loop

#parallel loop over mask indices
fourDMat <- foreach(vox=iter(vanalyze, by="row"), .inorder=TRUE, .packages=c("lme4", "pbkrtest"), .combine=rbind, 
        .noexport=c("cope_array", "varcope_array", "cope_voxes", "varcope_voxels")) %do% {
      
      #each v is 1 x nsessions*2 matrix
      #first half of columns are betas, second half are variances
      vdf <- cbind(featL1Df[,c("subid", "emotion", "contingency", "runnums")], data.frame(b=vox[1, 1:(ncol(vox)/2)], v=vox[1,(ncol(vox)/2 + 1):ncol(vox)]))
      
      #for some edge voxels, some subjects to not have any values. In this case, make a note of the voxel coordinates and analyze those with data present... 
      if (any(vdf$v==0.0)) {
        cat("Vmiss: ", paste(vindices[rownames(vox),], collapse=" "), "indices: ", paste(which(vdf$v==0.0), collapse=","), "\n", file="misslog", append=TRUE)
        vdf <- subset(vdf, v > 0)
      }
      
      browser()
      require(pbkrtest)
      
      #normalize reciprocal of variances to have sd = 1.0 (and mean of ~0.5-0.7)
      #vdf$weight <- scale(1/vdf$v, center=FALSE)
      
      #ldiv <- lmer(b ~ 1 + emotion + (1 | subid), data=vdf, weights=1/(v/1000000))
      l <- lmer(b ~ 1 + emotion + (1 | subid), data=vdf, weights=weight)
      #l <- lmer(b ~ 1 + emotion + (1 | subid), data=vdf, weights=1/(v/10000))
      #l <- lmer(b ~ 1 + emotion + (1 | subid), data=vdf)
      
      #Kenward-Roger approximation to df
      #note: this only produces one value, not a per-parameter estimate
      #df.KR <- get_ddf_Lb(l, fixef(l))
      # get p-values from the t-distribution using the t-values and approximated
      # degrees of freedom
      #coefs$p.KR <- 2 * (1 - pt(abs(coefs$t.value), df.KR))
      #coefs
      
      ltest <- summary(lmerTest::lmer(b ~ 1 + emotion + (1 | subid), data=vdf))$coefficients[,"df"]
      l2 <- lmerTest::lmer(b ~ 1 + emotion + (1 | subid), data=vdf, weights=1/v)
      #mod <- lmer(b ~ 1 + emotion + (1 | subid), data=vdf, weights=1/v)
      rl <- rlmer(b ~ 1 + emotion + (1 | subid), data=vdf)
      wt_huber <- getME(rl, "w_e")
      
      #nlme needs ~v, not 1/v, to get variance weighting correct
      test <- nlme::lme(b ~ 1 + emotion, random=~1|subid, data=vdf, weights=~v)
      ranef(test)
      
      #ltest <- lmer(b ~ 1 + emotion + (1 | subid), data=vdf, weights=wt_huber)
      
      af <- afex::mixed(b ~ 1 + emotion + (1 | subid), data=vdf, method="KR", test.intercept=TRUE, per.parameter=".")
      
      af <- mixed(b ~ 1 + emotion + (1 | subid), data=vdf, method="KR", test.intercept=TRUE, per.parameter=".")
      
      library(varComp)
      
      test2 <- varComp::varComp(b ~ 1 + emotion, vdf, ~subid)
      
      vout <- ltest$coefficients["(Intercept)", c("Estimate", "t value", "df", "Pr(>|t|)")]
      return(vout)      
      #afex::mixed breaks down for intercept only model
      #voxDf$dummy <- 0
      #ltest <- mixed(v ~ 1 + dummy + (1 | LunaID), voxDf, method="KR", test.intercept=TRUE, per.parameter=".")      
    }


colnames(fourDMat) <- c("b", "t", "df", "p")

#allocate memory
#nifti.image.alloc.data(out_hdr)

#recycle each row of the mask indices 4 times
#this is ugly
#maskIndicesMod <- cbind(mi[rep(1:nrow(mi), each=ncol(fourDMat)),], dim4=1:ncol(fourDMat))
#fourDUnmat <- gdata::unmatrix(fourDMat, byrow=TRUE) #convert into a vector by row (b, t, df, p) to match indices above

#insert results into 4d array
#system.time(out_hdr[maskIndicesMod] <- fourDUnmat) #as.matrix(fourDMat))

#aha! Use abind to build 4-d stats mat
#lapply over columns of matrix (take advantage of data.frame inheriting list)
out_hdr@.Data <- abind(lapply(data.frame(fourDMat), function(col) {
          m <- array(0, dim(mask)[1:3]) #empty 3d matrix matching dims of images
          m[mi] <- col
          return(m)
        }), along=4)



#from here: http://mindingthebrain.blogspot.com/2014/02/three-ways-to-get-parameter-specific-p.html
require(pbkrtest)
# get the KR-approximated degrees of freedom
df.KR <- get_ddf_Lb(m.sem, fixef(m.sem))
# get p-values from the t-distribution using the t-values and approximated
# degrees of freedom
coefs$p.KR <- 2 * (1 - pt(abs(coefs$t.value), df.KR))
coefs
