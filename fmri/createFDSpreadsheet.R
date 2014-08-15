fdfiles <- list.files(path="/Volumes/Serena/MMClock/MR_Raw", pattern="fd\\.txt", full.names=TRUE, recursive=TRUE)

subids <- unique(gsub("^.*/MR_Raw/(\\d+)_\\d+/.*$", "\\1", fdfiles, perl=TRUE))
nsubj <- length(subids)
df <- data.frame(LunaID=as.numeric(subids),
                 clock1nspikes_fd0.9=rep(NA_integer_, nsubj),
                 clock1include=rep(0L, nsubj),
                 clock1notes=rep("", nsubj),
                 clock1spikevols=rep(NA_character_, nsubj),
                 clock2nspikes_fd0.9=rep(NA_integer_, nsubj),
                 clock2include=rep(0L, nsubj),
                 clock2notes=rep("", nsubj),
                 clock2spikevols=rep(NA_character_, nsubj),
                 clock3nspikes_fd0.9=rep(NA_integer_, nsubj),
                 clock3include=rep(0L, nsubj),
                 clock3notes=rep("", nsubj),
                 clock3spikevols=rep(NA_character_, nsubj),
                 clock4nspikes_fd0.9=rep(NA_integer_, nsubj),
                 clock4include=rep(0L, nsubj),
                 clock4notes=rep("", nsubj),
                 clock4spikevols=rep(NA_character_, nsubj),
                 clock5nspikes_fd0.9=rep(NA_integer_, nsubj),
                 clock5include=rep(0L, nsubj),
                 clock5notes=rep("", nsubj),
                 clock5spikevols=rep(NA_character_, nsubj),
                 clock6nspikes_fd0.9=rep(NA_integer_, nsubj),
                 clock6include=rep(0L, nsubj),
                 clock6notes=rep("", nsubj),
                 clock6spikevols=rep(NA_character_, nsubj),
                 clock7nspikes_fd0.9=rep(NA_integer_, nsubj),
                 clock7include=rep(0L, nsubj),
                 clock7notes=rep("", nsubj),
                 clock7spikevols=rep(NA_character_, nsubj),
                 clock8nspikes_fd0.9=rep(NA_integer_, nsubj),
                 clock8include=rep(0L, nsubj),
                 clock8notes=rep("", nsubj),
                 clock8spikevols=rep(NA_character_, nsubj),
                 stringsAsFactors = FALSE
                 )
                            
                 
for (f in fdfiles) {
    subid <- as.numeric(gsub("^.*/MR_Raw/(\\d+)_\\d+/.*$", "\\1", f, perl=TRUE))
    runid <- gsub("^.*/MR_Raw/\\d+_\\d+/MBclock_recon/(clock[0-9])/.*$", "\\1", f, perl=TRUE)
    fd <- read.table(f, header=FALSE)$V1
    fdspikes <- which(fd > 0.9)
    nspikes <- length(fdspikes)
    df[which(df$LunaID == subid), paste0(runid, "nspikes_fd0.9")] <- nspikes
    df[which(df$LunaID == subid), paste0(runid, "include")] <- as.integer(nspikes < floor(length(fd)/10)) #fewer than 10% spikes
    df[which(df$LunaID == subid), paste0(runid, "spikevols")] <- paste(fdspikes, collapse="; ")
}

library(xlsx)
if(!file.exists("MMY3_motion_checks.xlsx")) { write.xlsx(x=df, file="MMY3_motion_checks.xlsx", sheetName="fdchecks", row.names=FALSE) }

