###Cleaning pipeline:

rootdir="/Volumes/bek/explore/MR_Proc"
paradname="clockRev_proc"
paradname2="clockrev"
procpatt="nfswudktm_clockrev[0-9]_7.nii.gz"
allprocID<-list.dirs(rootdir,recursive = F,full.names = F)
behavID<-unique(trial_df$id)
allIDs<-unique(c(allprocID,behavID))
alldirsdf<-do.call(rbind,lapply(allIDs,function(d){
  hasbehav<-d %in% behavID

  procexists<-dir.exists(file.path(rootdir,d))
  if(procexists){
    dx<-list.files(path = file.path(file.path(rootdir,d),paradname),pattern = paradname2,full.names = T)

    if(length(dx)>0){clockproc<-T} else {clockproc<-F;dx<-NA}
  } else {clockproc<-F;dx<-NA}



  rx<-data.frame(paths=dx,stringsAsFactors = F)
  rx$ID<-d
  rx$behavexists<-hasbehav
  rx$procexists<-procexists
  rx$clockproc<-clockproc
  return(rx)
})
)


alldirsdf$preprocfinished<-sapply(alldirsdf$paths,function(gra){
  length(list.files(path = gra,pattern = procpatt,recursive = F,all.files = T))>0
})

alldirsdf$LVL1_COMPLETE<-sapply(alldirsdf$paths,function(gra){
  runnum<-sub('.*(?=.$)', '', gra, perl=T)
  length(list.files(file.path(dirname(gra),"sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed"),
             pattern = paste0("^FEAT_LVL1_run",runnum,".feat"),recursive = F,include.dirs = T))>0
 })

alldirsdf$LVL2_COMPLETE<-sapply(alldirsdf$paths,function(gra){
  #runnum<-sub('.*(?=.$)', '', gra, perl=T)
  length(list.files(file.path(dirname(gra),"sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed"),
                    pattern ="FEAT_LVL2.gfeat",recursive = F,include.dirs = T))>0
})
#Get these folks:
donelv1<-alldirsdf[which(alldirsdf$behavexists & alldirsdf$preprocfinished & !alldirsdf$LVL1_COMPLETE),]
donelv2<-alldirsdf[which(alldirsdf$behavexists & alldirsdf$preprocfinished & alldirsdf$LVL1_COMPLETE & alldirsdf$LVL2_COMPLETE),]





moinfo<-do.call(rbind,lapply(1:nrow(alldirsdf),function(i){

  fd<-read.table( file.path(alldirsdf$paths[i], "motion_info", "fd.txt") )$V1
  propSpikes_0p9 <- sum(as.integer(fd > 0.9))/length(fd)
  spikeExclude <- if (propSpikes_0p9 > .10) 1 else 0
  maxFD <- max(fd)
  meanFD <- mean(fd)
  maxMotExclude <- if (maxFD > 5) 1 else 0
  anyExclude <- if (as.logical(spikeExclude) | as.logical(maxMotExclude)) 1 else 0
  ExcludeReason<-if (as.logical(spikeExclude) && as.logical(maxMotExclude)) "Both" else if (as.logical(spikeExclude)) "SpikeExclude" else if (as.logical(maxMotExclude)) "maxMotExclude" else "None"
  data.frame(propSpikes_0p9, spikeExclude, meanFD, maxFD, maxMotExclude, anyExclude, ExcludeReason,stringsAsFactors=FALSE)
}))

alldirsdf_a<-cbind(alldirsdf,moinfo)



torundf<-alldirsdf[which(alldirsdf$preprocver=="2018/7/30"),]
torundf$donzo<-FALSE

#We will let throndike do the first 60 and mine do the rest
fslpipe::fsl_2_sys_env()
cl<-parallel::makeCluster(12,type = "FORK")
NX<-parallel::parSapply(cl,61:nrow(torundf),function(i){
  torundf$paths[i]->ipath
  message(ipath)
  cmda<-gsub("REPLACETHISPATH",ipath,
            "cd REPLACETHISPATH
            echo REPLACETHISPATH
fslmaths swudktm_REPLACETHISNUM_7 -mul wudktm_REPLACETHISNUM_extents_mask swudktm_REPLACETHISNUM_7 -odt float
fslmaths 'swudktm_REPLACETHISNUM_7' -Tmean tempMean
fslmaths 'swudktm_REPLACETHISNUM_7' -bptf 70.7714 -1 -add tempMean 'fswudktm_REPLACETHISNUM_7'
date > .temporal_filtering_complete
fslmaths 'fswudktm_REPLACETHISNUM_7' -Tmean 'fswudktm_mean_float' -odt float
fslmaths 'fswudktm_REPLACETHISNUM_7' -mul 100 -div fswudktm_mean_float 'nfswudktm_REPLACETHISNUM_7' -odt float
date > .rescaling_complete"
  )
  cmd<-gsub("REPLACETHISNUM",basename(ipath),cmda)
  #if(!file.exists(file.path(ipath,"tempMean.nii.gz"))){
  system(cmd,intern = F)
  #}
  #Cleanup clock;
  cpaths<-c(list.files(path = dirname(ipath),pattern = "^FEAT_LVL1_run[0-9].feat",recursive = T,include.dirs = T,full.names = T),
            list.files(path = dirname(ipath),pattern = "^FEAT_LVL2_runtrend.gfeat$",recursive = T,include.dirs = T,full.names = T)  )
  lapply(cpaths,unlink,recursive=T)
  torundf$donzo[i]<-TRUE
})
parallel::stopCluster(cl)

###Do aroma:

alldirsdf$ica_aroma<-sapply(alldirsdf$paths,function(gra){
  dir.exists(file.path(gra,".ica_aroma_complete"))
})
alldirsdf$wavelet_despike<-sapply(alldirsdf$paths,function(gra){
  file.exists(file.path(gra,".despike_complete"))
})
alldirsdf$preproccompleted<-sapply(alldirsdf$paths,function(gra){
  file.exists(file.path(gra,"nfaswuktm_mean_func_6.nii.gz"))
})
alldirsdf$motioninfo<-sapply(alldirsdf$paths,function(gra){
  file.exists(file.path(gra,"motion_info", "fd.txt"))
})

alldirsdf$motioninfo<-sapply(alldirsdf$paths,function(gra){
  file.exists(file.path(gra,"motion_info", "fd.txt"))
})
percent_now<-length(which(alldirsdf$motioninfo))/nrow(alldirsdf)

while(percent_now != 1.0) {
  alldirsdf$motioninfo<-sapply(alldirsdf$paths,function(gra){
    file.exists(file.path(gra,"motion_info", "fd.txt"))
  })
  percent_now<-length(which(alldirsdf$motioninfo))/nrow(alldirsdf)
  print(percent_now)
  Sys.sleep(120)
}



lapply(alldirsdf$paths,function(gra){
  a1<-list.files(gra,"dktm_perserveoutput",include.dirs = T,full.names = T)
  if(length(a1)>0){
    a2<-list.files(a1,full.names = T)
    a3<-unique(dirname(dirname(a2)))
    a4<-list.files(a1,full.names = F)
    file.rename(from = a2,to = file.path(a3,a4))
  }
})


alldirsdf$LVL1_runnum<-sapply(alldirsdf$paths,function(gra){

  length(list.files(path = file.path(dirname(gra),"sceptic-clock-feedback-v_chosen-preconvolve_fse_groupfixed"),
                    pattern = paste0("^FEAT_LVL1_run",last(last(strsplit(gra,split = ""))),".feat"),recursive = F,include.dirs = T,full.names = F))
})





sp_alldirsdf<-split(alldirsdf_a,alldirsdf$ID)

allout<-do.call(rbind,lapply(sp_alldirsdf,function(dfx){
  LVL1_runnum=length(which(as.logical(dfx$LVL1_runnum)))
  lvl2runs=length(which(!as.logical(dfx$anyExclude)))
  anylvl2=any(!as.logical(dfx$anyExclude))
  RanLvL1=if(LVL1_runnum==2) T else F
  BothLvl2 = if (lvl2runs == 2) T else F
  Group = explore_demo$registration_group[which(explore_demo$registration_redcapid %in% unique(dfx$ID))]
  if(length(Group)<1){Group<-NA}
  data.frame(ID=unique(dfx$ID),anylvl2,
       RanLvL1,BothLvl2,Group,stringsAsFactors=FALSE)
}))












