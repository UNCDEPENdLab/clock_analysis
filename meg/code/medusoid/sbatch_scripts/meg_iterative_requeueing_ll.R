# R script for handling requeuing

epochs <- c("RT")

# FINISHING UP CLOCK ENTROPY CHANGE
# epochs <- c( "clock")

# effectively 3rd run of 1-sensor RT prediciton
# for fourth run, consider increasing number of cores/job to 4 or 8
# epochs <- c("RT")
# regressors_of_interest <- c("entropy_rs", "entropy_ri")
# regressors_of_interest <- c("entropy_change_full_ri", "entropy_change_pos_ri", "entropy_change_neg_ri")
regressors_of_interest <- c("abspe_by_rew") #negative is still running, only 77 done
basedir <- "/proj/mnhallqlab/projects/Clock_MEG/atfr_rds"
sbatch_dir <- "/nas/longleaf/home/dnpl/code/clock_analysis/meg/code/medusoid/sbatch_scripts"
setwd(basedir)
test <- T
# try random-intercept models with lower RAM
# modified to run subsets based on CRC version
start_at = 0
end_at = 170
step_up <- tibble::tribble(
  ~gb, ~time,
  30, "4-00:00:00",
  35, "4-00:00:00",
  40, "4-00:00:00",
  70, "4-00:00:00"
)

for (ee in epochs) {
  epochdir <- file.path(basedir, ee)
  flist <- list.files(pattern = "^freq_t", path = epochdir)
 fnum <- seq_along(flist)
 run_already <- fnum <= start_at | fnum >= end_at
  #file.create(paste0(epochdir,"/","tempfile.csv"))
  #this_f <- "tempfile.csv"
  #writeLines(as.character(1), this_f)
  for (rr in regressors_of_interest) {
    # out_exists <- list.files(pattern = paste0("meg_mixed_by_tf_ddf_wholebrain_", rr, "_rs.*"), path = epochdir)
    #out_expect <- file.path(epochdir, paste0("meg_mixed_by_tf_ddf_wholebrain_entropy_change_rs_single_", rr, "_rs_single_", ee, fnum))
    # entropy change:
  #  out_expect <- file.path(epochdir, paste0("meg_mixed_by_tf_ddf_wholebrain_", rr, "_rs_single", ee, fnum))
    # RT:
    if (rr=="rt") {
    out_expect <- file.path(epochdir, paste0("meg_tf_rdf_wholebrain_", rr, "_rs_single_sensor_", ee, fnum))
    } else {
      if (stringr::str_detect(rr, "_ri")){
    out_expect <- file.path(epochdir, paste0("meg_mixed_by_tf_ddf_wholebrain_", rr, "_single_", ee, fnum))  
      } else {
     out_expect <- file.path(epochdir, paste0("meg_mixed_by_tf_ddf_wholebrain_", rr, "_rs_single_", ee, fnum))  
      }
    }
    compute_expect <- file.path(epochdir, paste0(".", rr, "_it", fnum, "_compute"))
    out_exists <- file.exists(out_expect)
    if (any(out_exists)) {
      #message("The following files have been analyzed: ")
      #print(flist[out_exists])
      message("Number of files analyzed: ")
      print(length(flist[out_exists]))
      message("Number remaining: ")
      print(length(flist[!out_exists]))
    }
    to_run <- flist[!out_exists & !run_already]
    it_run <- fnum[!out_exists  & !run_already]
    compute_run <- compute_expect[!out_exists & !run_already]
    # out_run <- out_expect[!out_exists]
    
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
            "cd /nas/longleaf/home/dnpl/code/clock_analysis/meg/code/medusoid/sbatch_scripts; ",
            "sbatch -t ", step_up$time[level], " --mem=", step_up$gb[level], "g",
            " --export=epoch=", ee, ",sourcefilestart=", it_run[ff], ",regressor=", rr,
            " sbatch_meg_mixed_wholebrain_ll.bash"
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
             " sbatch_meg_mixed_wholebrain_ll.bash\n"
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
