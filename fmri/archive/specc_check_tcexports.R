#data check for SPECC fMRI
library(gdata)
library(dplyr)
basedir <- "/Users/michael/Dropbox/Hallquist_K01/Data/fMRI"
pinfo <- read.xls("/Users/michael/Box_Sync/DEPENd/Projects/SPECC/ID Management/SPECC_Participant_Info.xlsx")
pinfo <- filter(pinfo, HasClock==1)

for (s in 1:nrow(pinfo)) {
  thisid <- as.character(pinfo[s, "SPECC_ID"])
  dir <- system(paste0("find ", basedir, " -iname '*", thisid, "*' -type d"), intern=TRUE)
  if (length(dir)==0) { 
    message("Unable to locate fMRI dir for subject: ", thisid)  
    next
  } else if (length(dir) > 1) {
    message("Multiple matches found for subject: ", thisid)
    print(dir)
    next
  }
  
  file <- system(paste0("find ", dir, " -iname '*fMRIEmoClock*tcExport.csv' -type f"), intern=TRUE)
  
  if (length(file)==0) { 
    message("Unable to locate tcExport.csv file in dir: ", dir)  
    next
  } else if (length(dir) > 1) {
    message("Multiple CSV files found in dir: ", dir)
    print(dir)
    next
  }
}