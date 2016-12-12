#!/usr/bin/env Rscript

#read in command line arguments.
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0L) {
    message("subset_afni_runtiming.R expects -timing_dir and -good_runs <1 .. n> as inputs\n")
    quit(save="no", 1, FALSE)
}

argpos <- 1
goodruns <- c()
timing_dir <- c()
while (argpos <= length(args)) {
    if (args[argpos] == "-good_runs") {
        argpos <- argpos + 1
        while (argpos <= length(args) && substr(args[argpos], 1, 1) != "-") {
            goodruns <- c(goodruns, as.numeric(args[argpos]))
            argpos <- argpos + 1
        }
    } else if (args[argpos] == "-timing_dir") {
        timing_dir <- args[argpos + 1]
        argpos <- argpos + 2
        stopifnot(file.exists(timing_dir))
    } else { stop("Unrecognized option: ", args[argpos]) }
}

regs <- c("clock", "feedback", "ev", "rpe_pos", "rpe_neg", "mean_uncertainty", "rel_uncertainty")
emos <- c("", "_happy", "_scram", "_fear")

##determine number of volumes in each run
setwd(timing_dir)
runvols <- sapply(1:8, function(r) {
    nvols <- R.utils::countLines(paste0("run", r, "_clock.1D"))
})

lastpos <- cumsum(runvols)
tokeep <- do.call(c, lapply(goodruns, function(x) {
    if (x == 1) startpos <- 1
    else startpos <- lastpos[x - 1] + 1
    stoppos <- lastpos[x]
    return(startpos:stoppos)
}))

for (r in regs) {
    for (e in emos) {
        txt <- readLines(paste0(r, e, "_concat.1D"))
        writeLines(text=txt[tokeep], con=paste0(r, e, "_concat_runs", paste(goodruns, collapse=""), ".1D"))
    }
}

##motion pcs
mot <- readLines("../motion_pcs.txt")
writeLines(text=mot[tokeep], con=paste0("../motion_pcs_runs", paste(goodruns, collapse=""), ".txt"))
##censor file
cens <- readLines("../censor_intersection_concat.1D")
writeLines(text=cens[tokeep], con=paste0("../censor_intersection_concat_runs", paste(goodruns, collapse=""), ".1D"))
