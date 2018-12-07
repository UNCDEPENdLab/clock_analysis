library(parallel)
source("glm_helper_functions.R")

target=Sys.getenv("TARGET")
if (target=="") { stop("Must pass in TARGET as environment variable") }
stopifnot(file.exists(target))

rerun <- FALSE
ncpus <- 40
setwd(target)
##stop("exiting in dir: ", target) #leftover for testing
fsfFiles <- system(paste("find", getwd(), "-mindepth 3 -iname \"FEAT_LVL1_*.fsf\" -ipath \"*sceptic_vchosen*\" -type f"), intern=TRUE)
#    list.files(path=getwd(), pattern="FEAT_LVL1_.*\\.fsf", recursive=TRUE, full.names=TRUE)

#figure out which fsf files have already been run
dirExpect <- gsub("\\.fsf$", ".feat", fsfFiles, perl=TRUE)
torun <- list()

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

cat("About to run the following fsf files in parallel:\n\n")
cat(unlist(torun), sep="\n")

cl_fork <- makeForkCluster(nnodes=ncpus)

runfeat <- function(fsf) {
    runname <- basename(fsf)
    runFSLCommand(paste("feat", fsf), stdout=file.path(dirname(fsf), paste0("feat_stdout_", runname)),
                  stderr=file.path(dirname(fsf), paste0("feat_stderr_", runname)), fsldir="/gpfs/group/mnh5174/default/sw/fsl")
}

clusterApply(cl_fork, torun, runfeat)
stopCluster(cl_fork)

system(paste0("bash /gpfs/group/mnh5174/default/clock_analysis/fmri/gen_feat_reg_dir.bash ", target))
