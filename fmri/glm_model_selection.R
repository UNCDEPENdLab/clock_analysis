##concatenate design matrices across runs and provide unique baseline per run
##also apply .01 Hz high-pass filter to be consistent with fMRI preprocessing
##depends on rs-fcMRI_Functions.R
concatDesign <- function(d, hpass=.01) {
    require(plyr)
    
    d_allruns <- do.call(rbind.fill, lapply(1:length(d$design.convolve), function(r) {
        thisrun <- d$design.convolve[[r]]
        basecols <- grepl("base", names(thisrun))
        ##note that this will rename repeated names into something like ev and ev.1, which is good
        if (!is.null(hpass)) {
            thisrun <- as.data.frame(cbind(lapply(thisrun[,!basecols], function(col) { lmBandpass(col, dt=1.0, hpass, 5) }), thisrun[,basecols])) #.01 Hz high-pass filter
        }
        names(thisrun) <- gsub("base", paste0("run", r, "base"), names(thisrun))
        thisrun
    }))

    d_allruns[which(is.na(d_allruns), arr.ind=TRUE)] <- 0
    d_allruns <- as.matrix(d_allruns) #needed for lm

    d_allruns

}

visualizeDesignMatrix <- function(d, outfile=NULL, runboundaries=NULL, events=NULL, includeBaseline=TRUE) {
  require(ggplot2)
  require(reshape2)

  if (!includeBaseline) {
      d <- d[,!grepl("run[0-9]+base", colnames(d))]
  }

  print(round(cor(d), 3))
  d <- as.data.frame(d)
  d$volume <- 1:nrow(d)
  d.m <- melt(d, id.vars="volume")
  g <- ggplot(d.m, aes(x=volume, y=value)) + geom_line(size=1.2) + theme_bw(base_size=15) + facet_grid(variable ~ ., scales="free_y")
  
  colors <- c("black", "blue", "red", "orange") #just a hack for color scheme right now
  
  if (!is.null(runboundaries)) {
    g <- g + geom_vline(xintercept=runboundaries, color=colors[1L])
  }
  
  #browser()
  
  if (!is.null(events)) {
    for (i in 1:length(events)) {
      g <- g + geom_vline(xintercept=events[[i]], color=colors[i+1])
    }
  }
  
  if (!is.null(outfile)) {
      ggsave(filename=outfile, plot=g, width=21, height=9)
  }
  return(invisible(g))
}

incr_fit <- function(designmat, sparsemat, mask3d, modelname=NULL, outputDf=FALSE, outputPvals=FALSE, outputBetas=TRUE, njobs=20, models="auto", add_derivs=FALSE) {
    require(oro.nifti)
    require(foreach)
    require(iterators)
    require(doSNOW)
    require(abind)
    require(compiler)
    
    on.exit(stopCluster(clusterobj))
    
    setDefaultClusterOptions(master="localhost", port=10290) #move away from 10187 to avoid collisions
    clusterobj <- makeSOCKcluster(njobs)
    registerDoSNOW(clusterobj)
    
    #expect baseline columns to be named runXbaseY
    baselineCols <- grepl("run[0-9]+base", colnames(designmat))
    dmat_baseline <- designmat[,baselineCols]
    dmat_interest <- designmat[,!baselineCols]
    comb3d <- function(...) { abind(..., along=3) }

    if (is.null(modelname)) {
      modelname <- paste(colnames(dmat_interest), collapse="_")
    }

    visualizeDesignMatrix(dmat_interest, outfile=paste0(modelname, "_design.pdf"), includeBaseline=FALSE)
    
    rownames(sparsemat) <- 1:nrow(sparsemat) #set rownames to track progress inside foreach
    ##pb <- txtProgressBar(min = 1, max = nrow(sparsemat), style = 3)

    ##core subfunction that fits incremental regression models based on a single dependent variable and a design matrix of substantive interest
    incr_regress <- function(y, dmat_interest, outputDf=FALSE) {
        if (models=="auto") {
            models <- as.list(colnames(dmat_interest))
            message("Incremental fits determined by adding columns sequentially to design matrix: ", paste(colnames(dmat_interest), collapse=", "))
            dmat_incrlist <- lapply(1:ncol(dmat_interest), function(col) {
                dmat_interest[,1:col, drop=FALSE]
            })
        } else {
            ##alternative specification is to add multiple columns to design at a time 
            dmat_incrlist <- lapply(1:length(models), function(m) {
                ##each vector within the models list specifies additional columns to add
                ##example: list(c("clock", "feedback"), c("rpe_pos", "rpe_neg"))
                dmat[,models[[m]], drop=FALSE]
            })
        }

        if (add_derivs) {
            message("Adding temporal derivatives of each substantive regressor orthogonalized against design matrix")
            message("These are added in the same step as the regressor itself")

            ##generate derivative matrix for each column
            dmat_derivatives <- do.call(cbind, lapply(1:ncol(dmat_interest), function(col) {
                dx <- c(0, diff(dmat_interest[,col]))

                ##orthogonalize wrt design
                return(residuals(lm(dx ~ dmat_interest)))                
            }))

            colnames(dmat_derivatives) <- paste0("d", colnames(dmat_interest))

            #add relevant derivatives in each model step
            dmat_incrlist <- lapply(dmat_incrlist, function(m) {
                cbind(m, dmat_derivatives[, paste0("d", colnames(m)), drop=FALSE])
            })

            ##update expected columns for each model to include derivatives
            models <- lapply(models, function(m) c(m, paste0("d", m)) )
        }

        names(dmat_incrlist) <- as.character(1:length(dmat_incrlist))

        lmlist <- lapply(0:length(dmat_incrlist), function(n) {
            if (n == 0) {
                lm(y ~ 1) #intercept only null model  
            } else {
                dat <- data.frame(y, dmat_incrlist[[n]])
                fmodel <- as.formula(paste("y ~ 1 + ", paste(dimnames(dmat_incrlist[[n]])[[2]], collapse=" + ")))
                #lm(y ~ 1 + dmat_incrlist[[n]])
                lm(fmodel, dat)
            }
        })

        #retain incremental R^2
        rsq <- diff(sapply(lmlist, function(x) { summary(x)$r.squared} ))
        
        an <- do.call(anova, lmlist)
        
        ##remove null model from list since this is just for comparison to additional models
        an <- an[2:nrow(an),]
        lmlist <- lmlist[2:length(lmlist)]

        ##add full R2, F, and p values for full model
        fullM <- summary(lmlist[[length(lmlist)]])
        
        ## basic incremental statistics: rsq, F
        stats <- cbind(rsq=c(rsq, fullM$r.squared), F=c(an$F, unname(fullM$fstatistic[1])))
        
        #get the betas for the added effects
        if (outputBetas) { betas <- c(lapply(1:length(lmlist), function(m) { coefficients(lmlist[[m]])[ models[[m]] ] }), NULL) } #add null to list at the end because fullm has no betas of interest
        if (outputPvals) { stats <- cbind(stats, invp=c(1-an$`Pr(>F)`, do.call(pf, as.list(unname(fullM$fstatistic))))) }
        if (outputDf) { stats <- cbind(stats, numDf=c(an$Df, unname(fullM$fstatistic[2])), denDf=c(an$Res.Df, unname(fullM$fstatistic[3]))) }

        ##rownames(stats) <- c(sapply(models, function(m) { paste0("add_", m, collapse="_") }), "fullm")  ##too long: vars added should be clear from betas
        rownames(stats) <- c(paste0("m", 1:length(models)), "fullm")
        
        ##new approach: unmatrix the stats outputs here, such that we return a concatenated vector of stats for each model
        ##this allows for a variable number of betas per model to be returned (otherwise, we would have a non-rectangular matrix)
        stat_comb <- lapply(1:nrow(stats), function(r) {
            if (r < nrow(stats)) c(betas[[r]], stats[r,]) else stats[r,]
        })
        names(stat_comb) <- rownames(stats)
        stat_vec <- unlist(stat_comb)

        return(stat_vec)
    }

    #incr_regress <- cmpfun(incr_regress) #byte compile central function for speed
    
    #figure out degrees of freedom for incremental models, since these will not change per voxel
    #note: we cannot get away with simply adding a vector of 1s and explicitly modeling intercept because R changes
    #its calculation of R^2 and F in this case in ways that make things non-comparable.
    #see here: http://stats.stackexchange.com/questions/26176/removal-of-statistically-significant-intercept-term-boosts-r2-in-linear-model

    ##for adding df stats to AFNI header (so that p-values are computed within AFNI), make into a models x (numDf, denDf) matrix
    modeldf <- incr_regress(sparsemat[1,], dmat_interest, outputDf=TRUE)
    modeldf <- matrix(modeldf[grepl("(numDf|denDf)", names(modeldf))], ncol=2, byrow = TRUE, dimnames=list(NULL, c("numDf", "denDf")))

    cat("", file="log.txt") #clear
    #need to explicitly output settings to be used within incr_regress since it operates as a subfunction and foreach misses the variables
    res <- foreach(vox=iter(sparsemat, by="row"), .combine=rbind, .multicombine=TRUE,
                   .export=c("models", "outputDf", "outputPvals", "outputBetas", "add_derivs"), .noexport=c("sparsemat", "concatMR"),
                   .packages="forecast", .inorder=TRUE) %dopar% {
        ##cat(rownames(vox), "\n")
        ##setTxtProgressBar(pb, as.integer(rownames(vox)))
        ##flush.console()
        vnum <- as.integer(rownames(vox))
        if (vnum %% 100 == 0) cat("Voxel ", vnum, "\n", file="log.txt", append=TRUE)

        ##pull out baseline model
        vox <- vox[1,] #convert to vector
        mbase <- lm(vox ~ 0 + dmat_baseline) #intercept is included for each run. need to exclude typical intercept to avoid rank degeneracy
        voxnobase <- residuals(mbase)
        
        ##fit initial lm to residuals after baseline adjustment
        m <- lm(voxnobase ~ dmat_interest)
        
        ##fit arma(2,2) to residuals. looks like most voxels have approximately this structure
        a <- Arima(residuals(m), c(2, 0, 2))
        
        ##prewhiten voxel and design matrix using ARIMA coefficients from residuals, then re-run lm
        dmat_whiten <- apply(dmat_interest, 2, function(col) {
            as.vector(Arima(col, model=a)$residuals)
        })
        
        voxwhite <- Arima(voxnobase, model=a)$residuals

        stats <- incr_regress(voxwhite, dmat_whiten, outputDf)

        return(stats)
    }

    ##close(pb)

    ##new approach returns a voxel x stat 2d matrix
    ##permute to voxel x model x stat
    ##res <- aperm(res, c(3,1,2))
    ##dimnames(res)[[2]] <- c("fullm", paste0("add_", colnames(dmat_interest)))
    ##rescollapse <- t(apply(res, 1, function(x){ gdata::unmatrix(x, byrow=TRUE) }))

    ##generate empty cube to fill
    dat <- array(0, c(dim(runmask), dim(res)[2]))
    
    ##fill matrix
    ind4d <- cbind(pracma::repmat(mask3d, ncol(res), 1), rep(1:ncol(res), each=nrow(res)))
    
    dat[ind4d] <- res
        
    #use a 3dREMLfit output as a starting point for an AFNI header
    afniout <- readAFNI("glm_hrf_clockRT_feedback0_rpepos0_rpeneg0_evtmaxRPE_stats+tlrc.HEAD")
    
    ##hard to create an empty one apparently
#    new("afni", dat, dim(dat), c("x", "y", "z", "stat"),
#        SCENE_DATA=c(2L, 11L, 1L),
#        TYPESTRING="3DIM_HEAD_FUNC",
#        DATASET_RANK=c(3L, ncol(res))
#    )
    
    #documentation: http://afni.nimh.nih.gov/pub/dist/src/README.attributes
    afniout@.Data <- dat
    afniout@ORIGIN <- c(94.5, 132, -76.5)
    afniout@SCENE_DATA=c(2L, 11L, 1L) #2=tlrc view, 11=anat_buck_type, 1=3dim_head_func typestring
    afniout@ORIENT_SPECIFIC <- c(1L,2L,4L) #LPI
    afniout@DELTA <- c(-3, -3, 3) #voxel size in xyz
    afniout@DATASET_RANK[1:2] <- c(3L, ncol(res))
    afniout@DATASET_DIMENSIONS[1:3] <- dim(dat)[1:3] #x y z dimensions
    afniout@BRICK_TYPES <- rep(3L, ncol(res)) #float
    afniout@BRICK_FLOAT_FACS <- rep(0, ncol(res)) #no scaling of brik value
    afniout@BRICK_LABS <- paste(colnames(res), collapse="~") #labels in AFNI viewer
    afniout@BRICK_STATS <- sapply(as.vector(apply(dat, 4, range)), function(x) { if(is.na(x)) 0 else x }) #range of values in each sub brick
    
    #add partial F-statistic information for use in AFNI
    #http://afni.nimh.nih.gov/afni/doc/source/parser__int_8c.html
    #FUNC_FT_TYPE = 4 is type for F-statistic
    fstats <- grep(".F$", colnames(res))
    afniout@BRICK_STATAUX <- as.vector(sapply(1:length(fstats), function(x) {
              c(fstats[x]-1, #brick number (0-based)
                  4, #FUNC_FT_TYPE
                  2, #number of additional parameters (numDF and denDF)
                  modeldf[x,"numDf"], modeldf[x,"denDf"]
              )
            }))
    
    
    afniout@IDCODE_STRING <- modelname
    
    writeAFNI(afniout, paste0(modelname, "_rsq+tlrc"), verbose=TRUE)

    return(res)
}

##  lmrefit <- lm(voxwhite ~ dmat)
##  test <- as.data.frame(cbind(dmat, voxwhite))
##  system.time(cand <- bestglm(test, TopModels=20))
##  lm1 <- lm(voxwhite ~ dmat[,1])
##  lm2 <- lm(voxwhite ~ dmat[,c(1,2)])
        
