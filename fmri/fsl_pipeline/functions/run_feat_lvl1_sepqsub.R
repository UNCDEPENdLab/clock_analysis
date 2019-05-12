#this version of the FSL LVL1 feat estimation creates multiple qsub scripts in a temporary directory
##where each script has a number of feat calls that are forked, then the script waits for completion
##This circumvents the ICS limit on multiple nodes in R using a SNOW cluster.
##The primary calculation is how many files there are to run relative to processors and files per processor (chunking)

## to_run <- Sys.getenv("fsl_pipeline_file")

## run_model_index <- as.numeric(Sys.getenv("run_model_index")) #which variant to execute
## if (nchar(to_run) == 0L) { stop("Cannot locate environment variable fsl_pipeline_file") }
## if (!file.exists(to_run)) { stop("Cannot locate configuration file", to_run) }
## if (is.na(run_model_index)) { stop("Couldn't identify usable run_model_index variable.") }

## load(to_run)
## wait_for <- Sys.getenv("WAIT_FOR")

## #call function below
## run_feat_lvl1_sepqsub(feat_model_arguments, run_model_index, rerun=FALSE, wait_for=wait_for)

run_feat_lvl1_sepqsub <- function(fsl_model_arguments, run_model_index, rerun=FALSE, wait_for="") {

  #This function is now called within run_fsl_pipeline, rather than being run in its own qsub
  #to_run <- Sys.getenv("fsl_pipeline_file")
  #run_model_index <- as.numeric(Sys.getenv("run_model_index")) #which variant to execute
  #if (nchar(to_run) == 0L) { stop("Cannot locate environment variable fsl_pipeline_file") }
  #if (!file.exists(to_run)) { stop("Cannot locate configuration file: ", to_run) }
  #if (is.na(run_model_index)) { stop("Couldn't identify usable run_model_index variable.") }

  #load(to_run)

  subject_covariates <- fsl_model_arguments$subject_covariates
  expectdir <- fsl_model_arguments$expectdir
  model_match <- fsl_model_arguments$outdir[run_model_index]
  workdir <- fsl_model_arguments$workdir[run_model_index]

  cpusperjob <- 4 #number of cpus per qsub
  runsperproc <- 2 #number of feat calls per processor

  #look in the subfolder for each subject for fsf files
  fsfFiles <- do.call(c, lapply(subject_covariates$mr_dir, function(s) {
    system(paste0("find ", file.path(s, expectdir), " -mindepth 2 -iname \"FEAT_LVL1_*.fsf\" -ipath \"*/", model_match, "/*\" -type f"), intern=TRUE)
  }))

  #figure out which fsf files have already been run
  dirExpect <- gsub("\\.fsf$", ".feat", fsfFiles, perl=TRUE)
  torun <- c()

  for (f in 1:length(fsfFiles)) {
    if (file.exists(dirExpect[f])) {
      if (rerun) {
        cmd <- paste0("rm -rf \"", dirExpect[f], "\"")
        cat("Removing old directory: ", cmd, "\n")
        system(cmd)
        torun <- c(torun, fsfFiles[f]) #add to queue
      } else {
        cat("Skipping existing directory: ", dirExpect[f], "\n")
      }
    } else {
      torun <- c(torun, fsfFiles[f])
    }
  }

  if (length(torun) == 0) {
    cat("No LVL1 .fsf files to execute.\n\n")
    return(NULL)
  }

  cat("About to run the following fsf files in parallel:\n\n")
  cat(torun, sep="\n")

  if (!file.exists(workdir)) { dir.create(workdir, recursive=TRUE) }

  njobs <- ceiling(length(torun)/(cpusperjob*runsperproc))

  #use length.out on rep to ensure that the vectors align even if chunks are uneven wrt files to run
  df <- data.frame(fsf=torun, job=rep(1:njobs, each=cpusperjob*runsperproc, length.out=length(torun)), stringsAsFactors=FALSE)
  df <- df[order(df$job),]

  #We run this locally;
  require(parallel)
  cla<-parallel::makeCluster((cpusperjob*runsperproc),type="FORK",outfile="lvl1_seqsub_out.txt")
  runfeat <- function(fsf) {
    runname <- basename(fsf)
    runFSLCommand(paste("feat", fsf),stdout=file.path(dirname(fsf), paste0("feat_stdout_", runname)), stderr=file.path(dirname(fsf), paste0("feat_stderr_", runname)),fsldir = "/usr/local/ni_tools/fsl")
  }
  NX<-clusterApply(cl_fork, allFeatFiles, runfeat)
  stopCluster(cla)

}
