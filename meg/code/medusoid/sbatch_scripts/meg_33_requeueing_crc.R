# R script for handling requeuing of 0.33s epoch

epochs <- c("RT", "clock")

# FINISHING UP CLOCK ENTROPY CHANGE
# epochs <- c( "clock")

# effectively 3rd run of 1-sensor RT prediciton
# for fourth run, consider increasing number of cores/job to 4 or 8
# epochs <- c("RT")
# regressors_of_interest <- c("entropy")
regressors_of_interest <- c("entropy")
basedir <- "/bgfs/adombrovski/tfr_rds1"
sbatch_dir <- "~/code/clock_analysis/meg/code/medusoid/sbatch_scripts"
setwd(basedir)
test <- F
step_up <- tibble::tribble(
  ~gb, ~time,
  15, "7:00:00",
  15, "1-00:00:00",
  30, "1-00:00:00",
  30, "1-00:00:00",
  40, "1-00:00:00",
    40, "1-00:00:00",
      40, "1-00:00:00"

)

for (ee in epochs) {
  epochdir <- file.path(basedir, ee)
  flist <- list.files(pattern = ".*freq_t.*", path = epochdir)
  fnum <- grep("0.33", flist)
  flist33 <- flist[fnum]
  #file.create(paste0(epochdir,"/","tempfile.csv"))
  #this_f <- "tempfile.csv"
  #writeLines(as.character(1), this_f)
  for (rr in regressors_of_interest) {
    # out_exists <- list.files(pattern = paste0("meg_mixed_by_tf_ddf_wholebrain_", rr, "_rs.*"), path = epochdir)
    # entropy change:
    out_expect <- file.path(epochdir, paste0("meg_mixed_by_tf_ddf_wholebrain_", rr, "_rs_single_", ee, fnum))
    # RT:
  #  out_expect <- file.path(epochdir, paste0("meg_tf_rdf_wholebrain_", rr, "_rs_single_sensor33_", ee, fnum))
    compute_expect <- file.path(epochdir, paste0(".", rr, "_it", fnum, "_compute33"))
    out_exists <- file.exists(out_expect)
    if (any(out_exists)) {
      #message("The following files have been analyzed: ")
      #print(flist[out_exists])
      message("Number of files analyzed:\n ")
      print(length(flist33[out_exists]))
      message("Number remaining:\n ")
      print(length(flist33[!out_exists]))
    }
    
    to_run <- flist33[!out_exists]
    out_run <- out_expect[!out_exists]
    it_run <- fnum[!out_exists]
    compute_run <- compute_expect[!out_exists]

    if (length(to_run) > 0) {
      for (ff in seq_along(to_run)) {
        #look at prior compute file
        this_f <- compute_run[ff]
        if (file.exists(compute_run[ff])) {
          curLevel <- as.integer(readLines(this_f, n = 1))
          level <- curLevel + 1 
        } else {
          level <- 1
        }

        if (level > nrow(step_up)) {
          warning("Maximum compute used for: ", to_run[ff])
          next
        }
        if (!test) {
          system(
          paste0(
            "cd ~/code/clock_analysis/meg/code/medusoid/sbatch_scripts; ",
            "sbatch -t ", step_up$time[level], " --mem=", step_up$gb[level], "g",
            " --export=epoch=", ee, ",sourcefilestart=", it_run[ff], ",regressor=", rr,
             " sbatch_meg_mixed_wholebrain.bash"
          )
        )
        #write compute level to temporary file
        setwd(epochdir)
        writeLines(as.character(level), this_f)
        }
        
         cat(
           paste0(
             "sbatch -t ", step_up$time[level], " --mem=", step_up$gb[level], "g",
             " --export=epoch=", ee, ",sourcefilestart=", it_run[ff], ",regressor=", rr,
             " sbatch_meg_mixed_wholebrain.bash\n"
           )
         )

        
      }
      
    } else {
      message("All files have been analyzed for this iteration")
      next
    }

  }
  #file.remove(paste0(epochdir,"/","tempfile.csv"))
}
