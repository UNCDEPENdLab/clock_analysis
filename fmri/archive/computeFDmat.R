#compute a new fd censoring mat at 0.9mm
fdTxt <- list.files(path="/Volumes/Serena/MMClock/MR_Raw", pattern="fd.txt", full.names=TRUE, recursive = TRUE)

thresh = 0.9 #Siegel et al. 2014 HBM

for (f in fdTxt) {
    fd <- read.table(f)$V1

    censormat=do.call(cbind, sapply(1:length(fd), function(x) { if (abs(fd[x]) > thresh) {
        v <- rep(0, length(fd)); v[x] <- 1; v } else { NULL } }))

    if (!is.null(censormat)) {
        write.table(censormat, file=file.path(dirname(f), paste0("fd_", thresh, ".mat")), col.names=FALSE, row.names=FALSE, sep="\t")
        afnicensor <- 1 - rowSums(censormat)
    } else {
        afnicensor <- rep(1, length(fd)) #no bad volumes
    }
    write.table(afnicensor, file=file.path(dirname(f), paste0("fd_", thresh, "_censor.1D")), col.names=FALSE, row.names=FALSE, sep="\t")

    ##leave out reference rms change for now because there is always a big spike around the reference (middle volume)....
    ##stopifnot(file.exists(filepath(dirname(f), "refrms.txt")))
    ##refrms <- read.table("refrms.txt")$V1
    ##high_cutoff <- quantile(refrms, 0.75) + 10.0*IQR(refrms) #rms is always positive, so only have upper cutoff
    ##which(refrms > cutoff)

    #power et al. 2014 recommends DVARS of 20 as threshold for pre-BP data, assuming mode 1000 normalization
    
}
