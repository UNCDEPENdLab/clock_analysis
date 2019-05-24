# This script sets up the .fsf files to run a group analysis
## FSL Feat Level 3 analysis -- that is, mixed effects combination of subjects

#used for testing: group fixed entropy
#Sys.setenv(fsl_pipeline_file="/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline/configuration_files/MMClock_aroma_preconvolve_fse_groupfixed.RData")
#Sys.setenv(run_model_index=2)

#load the master configuration file
to_run <- Sys.getenv("fsl_pipeline_file")

run_model_index <- as.numeric(Sys.getenv("run_model_index")) #which variant to execute
if (nchar(to_run) == 0L) { stop("Cannot locate environment variable fsl_pipeline_file") }
if (!file.exists(to_run)) { stop("Cannot locate configuration file", to_run) }
if (is.na(run_model_index)) { stop("Couldn't identify usable run_model_index variable.") }

load(to_run)

library(tidyverse)
library(dependlab)

#verify that mr_dir is present as expected
subinfo <- fsl_model_arguments$subject_covariates
id_col <- fsl_model_arguments$id_col
feat_run_outdir <- fsl_model_arguments$outdir[run_model_index] #the name of the subfolder for the current run-level model
feat_lvl3_outdir <- file.path(fsl_model_arguments$group_output_dir, feat_run_outdir) #output directory for this run-level model
n_l1_copes <- fsl_model_arguments$n_l1_copes[run_model_index] #number of l1 copes determines number of FEAT LVL3 analyses to run (1 per LVL1 cope)
l1_cope_names <- fsl_model_arguments$l1_cope_names[[run_model_index]] #names of l1 copes (used for folder naming)

subinfo$dir_found <- file.exists(subinfo$mr_dir)

rerun <- FALSE #TODO: move to fsl_model_arguments list

cat("The following subjects were in the covariate file, but not the processed MRI data\n")
print(subset(subinfo, dir_found==FALSE))

dir.create(feat_lvl3_outdir, showWarnings=FALSE, recursive=TRUE)
setwd(feat_lvl3_outdir)

#cope structure for preconvolve models
#1 = clock_onset
#2 = feedback_onset
#3 = regressor of interest (in single-param models)

feat_lvl2_dirname <- "FEAT_LVL2.gfeat" #should populate this to the structure at some point
models <- fsl_model_arguments$group_model_variants #different covariate models for the current run-level model (run_model_index)

##rework using subinfo structure as the authoritative guide (rather than repeated searches)
copedf <- c()
nodatadf<-c()
for (s in 1:nrow(subinfo)) {
  for (cope in 1:n_l1_copes) {
    expectdir <- file.path(subinfo[s,"mr_dir"], fsl_model_arguments$expectdir, feat_run_outdir, feat_lvl2_dirname, paste0("cope", cope, ".feat"))
    if (dir.exists(expectdir)) {
      #print("yes")
      copedf <- rbind(copedf, data.frame(id=subinfo[s,id_col], model=feat_run_outdir, cope=cope, fsldir=expectdir,stringsAsFactors = F))
    } else {
      message("could not find expected directory: ", expectdir)
      nodatadf <- rbind(nodatadf, data.frame(id=subinfo[s,id_col], model=feat_run_outdir, cope=cope, fsldir=expectdir,stringsAsFactors = F))
    }
  }
}

names(copedf)[1] <- id_col #for matching

mdf <- merge(subinfo, copedf, by=id_col, all.y=TRUE)

#remove bad ids
mdf <- mdf %>% filter(!id %in% fsl_model_arguments$badids)
mdf <- arrange(mdf, id, model, cope) #should really use !!id_col here?

##fsl constructs models by cope
bycope <- lapply(split(mdf, mdf$cope), droplevels)

#loop over group-level models, setup the design matrix and spawn a FSL Level 3 job
l3template <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "fsf_templates", "feat_lvl3_explore_template.fsf"))

#loop over copes and group models, setting up .fsf files for each combination
rundf<-do.call(rbind,lapply(1:length(bycope),function(cope){
  if (is.null(bycope[[cope]]$Intercept)) { bycope[[cope]]$Intercept <- 1 } #add the column of ones

  copename <- l1_cope_names[cope]

  #cope-level subfolder
  #currently organized by run_model_name/l1_cope/l3_model (promotes comparisons of alternative l3 models of a given run-level effect)
  #could reorganize as run_model_name/l3_model/l1_cope (promotes comparisons of maps within an l3 model)
  model_output_dir <- file.path(feat_lvl3_outdir, copename) # paste0("cope", cope))
  dir.create(model_output_dir, showWarnings=FALSE)

  do.call(rbind,lapply(models,function(this_model) {

    model_df <- bycope[[cope]]
    fsf_syntax <- l3template #copy shared ingredients

    classesx<-sapply(this_model,function(m){class(model_df[[m]])})
    #We have to reorder the groups....it turns out
    if ("factor" %in% classesx){this_model<-names(c(which(classesx!="factor"),which(classesx=="factor")))}
    if (!"Intercept" %in% this_model ) { this_model <- c("Intercept", this_model) } #at present, force an intercept column

    if (fsl_model_arguments$center_l3_predictors) {
      for (p in this_model) {
        if (p != "Intercept" && is.numeric(model_df[[p]])) {
          model_df[[p]] <- model_df[[p]] - mean(model_df[[p]], na.rm=TRUE)
        }
      }
    }

    model_df$dummy_ <- rnorm(nrow(model_df))
    mform <- as.formula(paste("dummy_ ~ -1 + ", paste(this_model, collapse=" + ")))


    if ("factor" %in% classesx){
      model_df_og<-model_df
      this_model<-this_model[this_model!="Intercept"]
      factorvariname<-names(which(classesx=="factor"))
      nonfactorvariname<-names(which(classesx!="factor"))
      nonfactorvariname<-c(nonfactorvariname,"Intercept")
      if(length(factorvariname)>1){stop("Currently more than one group variable is not supported!")}
      glvels_og<-levels(model_df_og[[names(which(classesx=="factor"))]])
      pairsdf<-expand.grid(glvels_og,glvels_og,stringsAsFactors = F)
      pairsdf<-pairsdf[pairsdf$Var1!=pairsdf$Var2,]
      indx <- !duplicated(t(apply(pairsdf, 1, sort))) # finds non - duplicates in sorted rows
      pairsdf<-rev(pairsdf[indx, ])
      do.call(rbind,lapply(1:nrow(pairsdf),function(srx){
        groupdf<-pairsdf[srx,]

        model_df<-droplevels(model_df_og[which(model_df_og[[factorvariname]] %in% unlist(groupdf,use.names = F)),])
        glvels<-levels(model_df[[factorvariname]])
        fit_lm <- lm(mform, model_df)
        dmat <- model.matrix(fit_lm) #eventually allow interactions and so on??
        model_df$dummy_ <- NULL #clean up
        dmat<-dmat[,!colnames(dmat) %in% "Intercept"]

          sr<-groupdf
          ngrp<-length(sr)
          cmat<-matrix(nrow = 1,ncol = ngrp,data = 0)
          cmat[,match(sr[2],glvels)]<-1
          cmat[,match(sr[1],glvels)]<- -1
          cmat_name<-paste(rev(sr),collapse = ">")
        #cmat<-do.call(rbind,lapply(cmat_ls,function(r){r$cmat}))



        #add the other subjectwise co-variates
        horicmat<-matrix(0,nrow = length(nonfactorvariname),ncol = ncol(cmat))
        cmat<-rbind(cmat,horicmat)
        rownames(cmat)<-c(cmat_name,nonfactorvariname)
        vericmat<-matrix(0,ncol = length(nonfactorvariname[nonfactorvariname!="Intercept"]),nrow = nrow(cmat))
        cmat<-cbind(cmat,vericmat)
        colnames(cmat)<-c(paste0(factorvariname,glvels),nonfactorvariname[nonfactorvariname!="Intercept"])


        for (rnamex in nonfactorvariname){
          xa<-match(rnamex,rownames(cmat))
          if (rnamex=="Intercept") {
            refdfx<-as.data.frame(table(model_df[[names(which(classesx=="factor"))]]))
            for(rg in 1L:nrow(refdfx)){
              cmat[xa,which(grepl(paste0(refdfx$Var1[rg],"$"),colnames(cmat)))]<-round(refdfx$Freq[rg] / sum(refdfx$Freq),2)
            }
          } else {
            gra<-strsplit(rnamex,split = " > ")[[1]]
            cmat[xa,which(grepl(paste0(gra[1],"$"),colnames(cmat)))]<-1
            cmat[xa,which(grepl(paste0(gra[2],"$"),colnames(cmat)))]<-(-1)
          }
        }
        cmat<-cmat[,match(colnames(dmat),colnames(cmat))]

        fsf_syntax <- c(fsf_syntax, generate_fsf_contrast_syntax(cmat))
        fsf_syntax <- c(fsf_syntax, generate_fsf_ev_syntax(inputs=model_df$fsldir, dmat=dmat))
        fsf_syntax <- gsub(".OUTPUTDIR.", file.path(model_output_dir, paste0(copename, "-", paste(c(this_model,paste(unlist(groupdf,use.names = F),collapse = "and")), collapse="-"))), fsf_syntax, fixed=TRUE)

        #write the FSF to file
        out_fsf <- file.path(model_output_dir, paste0(copename, "-", paste(c(this_model,paste(unlist(groupdf,use.names = F),collapse = "and")), collapse="-"), ".fsf"))

        writeLines(fsf_syntax, con=out_fsf)
        write.table(table(model_df[[factorvariname]]),file=sub(".fsf", "_subj_cout.txt", out_fsf, fixed=TRUE),row.names=FALSE)
        dmat_file <- sub(".fsf", "_design.txt", out_fsf, fixed=TRUE)
        #write the design matrix to file for matching with extracted betas later
        model_df$feat_input_id <- 1:nrow(model_df) #for matching with extracted betas
        model_df <- model_df %>% select(-dir_found, -mr_dir) %>% select(id, feat_input_id, model, cope, fsldir, everything())
        write.table(model_df, file=dmat_file, row.names=FALSE)
        if (!file.exists(sub(".fsf", ".gfeat", out_fsf, fixed=TRUE)) || rerun) {
          runthis<-T
        } else {runthis<-F}
        data.frame(path=out_fsf,ifrun=runthis,cope=copename,stringsAsFactors = F)
      })
      )
   #This chunk of codes allows interaction between continous and factor regressors
      # dmat<-cbind(dmat[,which(!colnames(dmat) %in% names(which(classesx!="factor")))] ,do.call(cbind,lapply(glvels,function(r){
      #   dmat[,which(colnames(dmat) == names(which(classesx!="factor")))] * dmat[,grep(r,colnames(dmat))]
      # }))
      # )
      #colnames(dmat)[colnames(dmat)==""]<-NA
      #colnames(dmat)[is.na(colnames(dmat))]<-paste(glvels,names(which(classesx!="factor")),sep = "")

      #Dealing with factors:



      #This one allows group with their own variance intercept
      #fsf_syntax <- c(fsf_syntax, generate_fsf_ev_syntax(inputs=model_df$fsldir, dmat=dmat,group_membership = as.numeric(model_df[[names(which(classesx=="factor"))]])))
    } else {
      fit_lm <- lm(mform, model_df)
      dmat <- model.matrix(fit_lm) #eventually allow interactions and so on??
      model_df$dummy_ <- NULL #clean up
      cmat <- diag(length(this_model))
      rownames(cmat) <- this_model #just name the contrasts after the EVs themselves
      fsf_syntax <- gsub(".OUTPUTDIR.", file.path(model_output_dir, paste0(copename, "-", paste(this_model, collapse="-"))), fsf_syntax, fixed=TRUE)
      fsf_syntax <- c(fsf_syntax, generate_fsf_contrast_syntax(cmat))
      fsf_syntax <- c(fsf_syntax, generate_fsf_ev_syntax(inputs=model_df$fsldir, dmat=dmat))
      #write the FSF to file
      out_fsf <- file.path(model_output_dir, paste0(copename, "-", paste(this_model, collapse="-"), ".fsf"))

      writeLines(fsf_syntax, con=out_fsf)
      write.table(table(model_df[[factorvariname]]),file=sub(".fsf", "_subj_cout.txt", out_fsf, fixed=TRUE),row.names=FALSE)
      dmat_file <- sub(".fsf", "_design.txt", out_fsf, fixed=TRUE)
      #write the design matrix to file for matching with extracted betas later
      model_df$feat_input_id <- 1:nrow(model_df) #for matching with extracted betas
      model_df <- model_df %>% select(-dir_found, -mr_dir) %>% select(id, feat_input_id, model, cope, fsldir, everything())
      write.table(model_df, file=dmat_file, row.names=FALSE)


      if (!file.exists(sub(".fsf", ".gfeat", out_fsf, fixed=TRUE)) || rerun) {
        runthis<-T
      } else {runthis<-F}
      data.frame(path=out_fsf,ifrun=runthis,cope=copename,stringsAsFactors = F)
    }





  }))

}))

rundf[rundf$ifrun,]
print(rundf)

rundf<-rundf[which(rundf$cope==unique(rundf$cope)[3]),]
require(parallel)
cla<- makeForkCluster(fsl_model_arguments$ncpus,verbose=TRUE)
cat("Start Now")
runfeat <- function(fsf) {
  runname <- basename(fsf)
  runFSLCommand(args = paste("feat", fsf),stdout=file.path(dirname(fsf), paste0("feat_stdout_", runname)), stderr=file.path(dirname(fsf), paste0("feat_stderr_", runname)))
}
NX<-parallel::clusterApply(cla, rundf$path[which(rundf$ifrun)], runfeat)
parallel::stopCluster(cla)
