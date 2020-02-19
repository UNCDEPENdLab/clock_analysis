mrfile <- "/Volumes/Serena/MMClock/MR_Proc/cope1dirs"
subfile <- "subinfo_db"

subinfo <- read.table(subfile, header=TRUE)
mrinfo <- read.table(mrfile)$V1
mrorder <- 1:length(mrinfo) #check fidelity of row order to ensure match in FEAT
mrdf <- data.frame(n=mrorder, file=mrinfo, lunaid=sub("^.*/MR_Proc/(\\d+)_\\d+/.*$", "\\1", mrinfo, perl=TRUE))

submerge <- merge(mrdf, subinfo, by="lunaid", all.x=TRUE)

submerge$age.c <- submerge$age - mean(submerge$age)
submerge$female.c <- submerge$female - mean(submerge$female)
    
##FSL design with age and sex centered
design <- cbind(1, submerge$age.c, submerge$female.c)
write.table(design, file="fsl_agec_femalec_design.txt", row.names = FALSE, col.names=FALSE)
#write.table(submerge$file, file="fsl_cope.txt", row.names = FALSE, col.names=FALSE)
