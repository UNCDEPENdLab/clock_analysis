#!/usr/bin/Rscript

#read in command line arguments.
#current format:
#arg 1: directory to process (default current directory)
#arg 2: number of parallel jobs (default 8)
#arg 3: folder containing MRRC reconstructed MB data (rsync from meson)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0L) {
    goto <- args[1L]
    if (! file.exists(goto)) { stop("Cannot find directory: ", goto) }
    setwd(goto)
}

basedir <- getwd() #root directory for processing

if (length(args) > 1L) {
    njobs <- as.numeric(args[2L])
} else {
    njobs <- 8
}

if (length(args) > 2L) {
    MB_src <- args[3L] #folder containing MRRC reconstructed data
} else {
    MB_src <- normalizePath(Sys.glob("../WPC-*_MB")) #assume that MB data are up one directory in folder called WPC-XXXX_MB
}

#pull in cfg environment variables from bash script
mprage_dirpattern=Sys.getenv("mprage_dirpattern")
preprocessed_dirname=Sys.getenv("preprocessed_dirname")
paradigm_name=Sys.getenv("paradigm_name")
n_expected_funcruns=Sys.getenv("n_expected_funcruns")

##handle all mprage directories
##overload built-in list.dirs function to support pattern match
list.dirs <- function(...) {
    args <- as.list(match.call())[-1L] #first argument is call itself

    if (! "recursive" %in% names(args)) { args$recursive <- TRUE } #default to recursive
    if (! ("full.names" %in% names(args))) { args$full.names <- TRUE } #default to full names
    if (! "path" %in% names(args)) { args$path <- getwd() #default to current directory
                                 } else { args$path <- eval(args$path) }
    args$include.dirs <- TRUE

    flist <- do.call(list.files, args)

    oldwd <- getwd()
    if (args$full.names == FALSE) {
        #cat("path: ", args$path, "\n")
        setwd(args$path)
    }
    ##ensure that we only have directories (no files)
    ##use unlist to remove any NULLs from elements that are not directories
    dlist <- unlist(sapply(flist, function(x) { if (file.info(x)$isdir) { x } else { NULL } }, USE.NAMES = FALSE))
    setwd(oldwd)
    return(dlist) #will be null if no matches
}

#find original mprage directories to rename
mprage_dirs <- list.dirs(pattern=mprage_dirpattern)

if (!is.null(mprage_dirs)) {
    cat("Renaming original mprage directories to \"mprage\"\n")
    for (d in mprage_dirs) {
        mdir <- file.path(dirname(d), "mprage")
        file.rename(d, mdir) #rename to mprage
    }
}

#find all renamed mprage directories for processing
#use beginning and end of line markers to force exact match
#use getwd to force absolute path since we setwd below
mprage_dirs <- list.dirs(pattern="^mprage$", path=getwd()) 

library(foreach)
library(doMC)

registerDoMC(njobs) #setup number of jobs to fork

#for (d in mprage_dirs) {
f <- foreach(d=mprage_dirs, .inorder=FALSE) %dopar% {
    setwd(d)
    #call preprocessmprage
    if (file.exists(".mprage_complete")) {
        return("complete") #skip completed mprage directories
    } else {
        if (file.exists("mprage.nii.gz")) {
            args <- "-delete_dicom archive -template_brain MNI_2mm -nifti mprage.nii.gz"
        } else {
            args <- "-delete_dicom archive -template_brain MNI_2mm"
        }
        
        ret_code <- system2("preprocessMprage", args, stderr="preprocessMprage_stderr", stdout="preprocessMprage_stdout")
        if (ret_code != 0) { stop("preprocessMprage failed.") }

        #echo current date/time to .mprage_complete to denote completed preprocessing
        sink(".mprage_complete")
        cat(as.character(Sys.time()))
        sink()

        if (file.exists("need_analyze")) { unlink("need_analyze") } #remove dummy file
        if (file.exists("analyze")) { unlink("analyze") } #remove dummy file

        if (file.exists("mprage_bet.nii.gz")) {
            file.symlink("mprage_bet.nii.gz", "mprage_brain.nii.gz") #symlink to _brain for compatibility with FEAT/FSL
        }
    }
    return(d)
}

#get list of subject directories in root directory
subj_dirs <- list.dirs(path=basedir, recursive=FALSE)

#make run processing parallel, not subject processing
#f <- foreach(d=subj_dirs, .inorder = FALSE) %dopar% {
all_funcrun_dirs <- list()
for (d in subj_dirs) {
    cat("Processing subject: ", d, "\n")
    setwd(d)

    ##create paradigm_run1-paradigm_run8 folder structure and copy raw data
    if (!file.exists(preprocessed_dirname)) { #create preprocessed folder if absent
        dir.create(file.path(d, preprocessed_dirname), showWarnings = FALSE)
    } else {
        ##preprocessed folder exists, check for .preprocessfunctional_complete files
        extant_funcrundirs <- list.dirs(path=file.path(d, preprocessed_dirname), pattern=paste0(paradigm_name,"[0-9]+"), full.names=TRUE, recursive=FALSE)
        if (length(extant_funcrundirs) > 0L &&
            length(extant_funcrundirs) >= n_expected_funcruns &&
            all(sapply(extant_funcrundirs, function(x) { file.exists(file.path(x, ".preprocessfunctional_complete")) }))) {
            cat("   preprocessing already complete for all functional run directories\n\n")
            next
        }
    }

    #identify original reconstructed flies for this subject
    subid <- basename(d)
    mbraw_dirs <- list.dirs(path=MB_src, recursive = FALSE, full.names=FALSE) #all original recon directories, leave off full names for grep

    #approximate grep is leading to problems with near matches!!
    #example: 11263_20140307; WPC5640_11253_20140308
    #srcmatch <- agrep(subid, mbraw_dirs, max.distance = 0.1, ignore.case = TRUE)[1L] #approximate id match in MRRC directory
    srcmatch <- grep(subid, mbraw_dirs, ignore.case = TRUE)[1L] #id match in MRRC directory
    
    if (is.na(srcmatch)) {
        warning("Unable to identify reconstructed images for id: ", subid, " in MB source directory: ", MB_src)
        next #skip this subject
    }

    srcdir <- file.path(MB_src, mbraw_dirs[srcmatch])
    cat("Matched with MB src directory: ", srcdir, "\n")
    mbfiles <- list.files(path=srcdir, pattern=".*ep2d_MB_E?clock.*_MB.hdr$", full.names = TRUE) #images to copy

    ##figure out run numbers based on file names
    ##there is some variability in how files are named.
    ## v1: ep2d_MB_clock1_MB.hdr
    ## v2: ep2d_MB_clock1_8_MB.hdr (ambiguous!)
    ## v3: ep2d_MB_clock_1_MB.hdr
    ## occasionally "Eclock"?

    runnums <- sub("^.*ep2d_MB_E?clock(\\d?)_?(\\d?)_?(_FID)*.*_MB.hdr$",
                   "\\1 \\2", mbfiles, perl=TRUE, ignore.case = TRUE)

    run_split <- strsplit(runnums, "\\s+", perl=TRUE)
    run_lens <- sapply(run_split, length)

    if (any(run_lens > 1L)) {
        #at least one file name contains two potential run numbers
        #if any file has just one run number, duplicate it for comparison
        run_split <- lapply(run_split, function(x) { if(length(x) == 1L) { c(x,x) } else { x } } )

        #determine which potential run number contains unique information
        R1 <- unique(sapply(run_split, "[[", 1))
        R2 <- unique(sapply(run_split, "[[", 2))

        if (length(unique(R1)) > length(unique(R2))) {
            runnums <- R1
        } else {
            runnums <- R2
        }            
    }
            
    if (length(runnums) > length(unique(runnums))) {
        print(mbfiles)
        stop("Duplicate run numbers detected.")
    }

    runnums <- as.numeric(runnums)
    if (any(is.na(runnums))) { stop ("Unable to determine run numbers:", runnums) }

    cat("Detected run numbers, MB Files:\n")
    print(cbind(runnum=runnums, mbfile=mbfiles))
   
    #loop over files and setup run directories in preprocessed_dirname
    for (m in 1:length(mbfiles)) {
        #only copy data if folder does not exist
        if (!file.exists(file.path(d, preprocessed_dirname, paste0(paradigm_name, runnums[m])))) {
            dir.create(file.path(d, preprocessed_dirname, paste0(paradigm_name, runnums[m])))
            
            ##use 3dcopy to copy dataset as .nii.gz
            system(paste0("3dcopy \"", mbfiles[m], "\" \"", file.path(d, preprocessed_dirname, paste0(paradigm_name, runnums[m]), paste0(paradigm_name, runnums[m])), ".nii.gz\""))
        }
    }

    #now that preprocessed_dirname files are copied, preprocess all
    setwd(preprocessed_dirname)

    all_funcrun_dirs[[d]] <- list.dirs(pattern=paste0(paradigm_name, ".*"), path=getwd(), recursive = FALSE)
}

all_funcrun_dirs <- unname(unlist(all_funcrun_dirs)) #generate vector of all functional runs to process

#loop over directories to process
##for (cd in clockdirs) {
f <- foreach(cd=all_funcrun_dirs, .inorder=FALSE) %dopar% {
    setwd(cd)
    
    ##determine phase versus magnitude directories for fieldmap
    ##in runs so far, magnitude comes first. preprocessFunctional should handle properly if we screw this up...
    fmdirs <- sort(normalizePath(Sys.glob("../../gre_field_mapping*")))
    magdir <- file.path(fmdirs[1], "MR*")
    phasedir <- file.path(fmdirs[2], "MR*")

    cfile <- Sys.glob(paste0(paradigm_name, "*.nii.gz"))
    ##run preprocessFunctional
    args <- paste0("-4d ", cfile, " -tr 1.0 -mprage_bet ../../mprage/mprage_bet.nii.gz -warpcoef ../../mprage/mprage_warpcoef.nii.gz -threshold 98_2 ",
                   "-hp_filter 100 -rescaling_method 10000_globalmedian -template_brain MNI_2.3mm -func_struc_dof bbr -warp_interpolation spline ",
                   "-constrain_to_template y -wavelet_despike -4d_slice_motion -custom_slice_times /Volumes/Serena/SPECC/MR_Raw/speccMBTimings.1D ",
                   "-fm_phase \"", phasedir, "\" -fm_magnitude \"", magdir, "\" -fm_cfg clock ", #quote phasedir and magdir to prevent wildcard expansion
                   "-mc_movie -motion_censor fd=0.9,dvars=20")
    
    ret_code <- system2("preprocessFunctional", args, stderr="preprocessFunctional_stderr", stdout="preprocessFunctional_stdout")
    if (ret_code != 0) { stop("preprocessFunctional failed.") }
}
