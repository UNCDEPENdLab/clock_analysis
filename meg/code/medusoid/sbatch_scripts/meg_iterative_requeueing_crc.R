# R script for handling requeuing
epochs <- c("RT")
# epochs <- c("RT")

# epochs <- c("RT", "clock")
regressors_of_interest <- c("abspe_by_rew") # some duplication for v_max, of which 92 are still running
# regressors_of_interest <- c("entropy_change_ri", "entropy_change",  "abspe_by_rew")
 # HAVEN'T yet done a second round on abspe*reward
# regressors_of_interest <- c("entropy_change")
# regressors_of_interest <- c("entropy_change_fmr2")

basedir <- "/bgfs/adombrovski/tfr_rds1"
sbatch_dir <- "~/code/clock_analysis/meg/code/medusoid/sbatch_scripts"
setwd(basedir)
test <- F
silent <- T
# running in smaller batches since CRC queues up ≤1000 at a time
start_at = 0
end_at = 685
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
    if (stringr::str_detect(rr, "_ri")) {
    out_exists <- file.exists(pattern = paste0("meg_mixed_by_tf_ddf_wholebrain_", rr), path = epochdir)
    out_expect <- file.path(epochdir, paste0("meg_mixed_by_tf_ddf_wholebrain_", rr, "_single_", ee, fnum))} else {
    out_exists <- file.exists(pattern = paste0("meg_mixed_by_tf_ddf_wholebrain_", rr, "_rs.*"), path = epochdir)
    out_expect <- file.path(epochdir, paste0("meg_mixed_by_tf_ddf_wholebrain_", rr, "_rs_single_", ee, fnum))

    }

    compute_expect <- file.path(epochdir, paste0(".", rr, "_it", fnum, "_compute"))
    out_exists <- file.exists(out_expect)
    #"|",paste0("meg_mixed_by_tf_ddf_combined_", rr, "_rs_", ee, "__", gsub("_freq_t", "", basename(flist[fnum]))))
    compute_expect <- file.path(epochdir, paste0(".", rr, "_it", fnum, "_compute"))
    if (any(out_exists)) {
      #message("The following files have been analyzed: ")
      #print(flist[out_exists])
      message("Total number of files to re-run: ")
      print(sum(!out_exists))
    }
  # for i in meg_mixed_by_tf_ddf_combined*; do [[ -e ${i/__/_freq_t_all_} ]] || mv "$i" "${i/__/_freq_t_all_}"; done 
    to_run <- flist[!out_exists & !run_already]
    #out_run <- out_expect1[!out_exists]
    it_run <- fnum[!out_exists  & !run_already]
    compute_run <- compute_expect[!out_exists & !run_already]

    if (length(to_run) > 0) {
     message("Number of files to re-run in this batch: ")
      print(length(to_run))
    
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
        
    #  njobs <- as.integer(system("squeue --m | wc -l", intern=T)) - 1
    #    while (njobs < 1000)
    #    {
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
        writeLines(as.character(level), this_f)
        }
        if (!silent) {
        cat(
           paste0(
             "sbatch -t ", step_up$time[level], " --mem=", step_up$gb[level], "g",
            " --export=epoch=", ee, ",sourcefilestart=", it_run[ff], ",regressor=", rr,
             " sbatch_meg_mixed_wholebrain.bash\n",
             to_run[ff], "\n"
           )
         )
        }
        
      }
      
    } else {
      message("All files have been analyzed for this iteration")
      next
    }

  }
  #file.remove(paste0(epochdir,"/","tempfile.csv"))
}
