library(parallel)

scripts_dir <- "/gpfs/group/mnh5174/default/clock_analysis/fmri"

source(file.path(scripts_dir, "glm_helper_functions.R"))


##this version of the FSL LVL1 feat estimation creates multiple qsub scripts in a temporary directory
##where each script has a number of feat calls that are forked, then the script waits for completion
##This circumvents the ICS limit on multiple nodes in R using a SNOW cluster.
##The primary calculation is how many files there are to run relative to processors and files per processor (chunking)

target=Sys.getenv("TARGET")
model_match=Sys.getenv("MODEL_MATCH")
wait_for=Sys.getenv("WAIT_FOR")
if (target=="") { stop("Must pass in TARGET as environment variable") }
if (model_match=="") { model_match <- "sceptic-v_chosen-v_entropy-d_auc" }
stopifnot(file.exists(target))

rerun <- FALSE
cpusperjob <- 8 #number of cpus per qsub
runsperproc <- 3 #number of feat calls per processor
setwd(target)

fsfFiles <- system(paste0("find ", getwd(), " -mindepth 3 -iname \"FEAT_LVL1_*.fsf\" -ipath \"*/", model_match, "/*\" -type f"), intern=TRUE)
# list.files(path=getwd(), pattern="FEAT_LVL1_.*\\.fsf", recursive=TRUE, full.names=TRUE)

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

preamble <- c(
  "#PBS -A mnh5174_a_g_hc_default",
  paste0("#PBS -l nodes=1:ppn=", cpusperjob, ":himem"),
  ifelse(wait_for != "", paste0("#PBS -W depend=afterok:", wait_for), ""), #allow job dependency on upstream setup
  "#PBS -l walltime=4:00:00",
  "#PBS -j oe",
  "#PBS -M michael.hallquist@psu.edu",
  "#PBS -m abe",
  "#PBS -W group_list=mnh5174_collab",
  "",
  "",
  "source /gpfs/group/mnh5174/default/lab_resources/ni_path.bash",
  "module load \"fsl/5.0.11\" >/dev/null 2>&1",
  "",
  "cd $PBS_O_WORKDIR"
)

workdir <- "/gpfs/scratch/mnh5174/feat_qsub_temp"
if (!file.exists(workdir)) { dir.create(workdir) }

njobs <- ceiling(length(torun)/(cpusperjob*runsperproc))

#use length.out on rep to ensure that the vectors align even if chunks are uneven wrt files to run
df <- data.frame(fsf=unlist(torun), job=rep(1:njobs, each=cpusperjob*runsperproc, length.out=length(torun)), stringsAsFactors=FALSE)
df <- df[order(df$job),]
for (j in 1:njobs) {
    outfile <- paste0(workdir, "/qsub_featsep_", j, "_", basename(tempfile()), ".pbs")
    cat(preamble, file=outfile, sep="\n")
    thisrun <- with(df, fsf[job==j])
    cat(paste("feat", thisrun, "&"), file=outfile, sep="\n", append=TRUE)
    cat("wait\n\n", file=outfile, append=TRUE)
    cat(paste0("bash /gpfs/group/mnh5174/default/clock_analysis/fmri/gen_feat_reg_dir.bash ", unique(dirname(thisrun))), sep="\n", file=outfile, append=TRUE)
    system(paste0("qsub ", outfile))
}
