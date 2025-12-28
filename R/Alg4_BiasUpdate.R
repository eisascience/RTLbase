

#Mahyari, Eisa: Implementation of Lee et al.'s TL-FC into R
#March, 2017



##### Bias Compensation (Algorithm 4) ##########
###
### Input: hyperplane (w,b), target task data T
##  1: Compute: z_i = <w,x_t_i> +b, for all i
##  2: Build a Grid: s_i = sort (z_i)
##  3: for j = 1 to N_t
##  4: cj = sum(II{abs(z_i - s_j)/||w|| < 1}
##  5: end for
##  6: h = kernelbandwidth ({(s_j,c_j)}_j)
##  7: Smooth: p(z) = sum(c_j k_h (z, s_j))_j
##  8: z∗ = gradientdescent(p(z), 0)
##  9: b_new = b - z∗
##  Output: b_new or f_b(x)= <w,x_i> + b_new
###
############################################################
### Lee, Gyemin, L Stoolman, and C Scott.
### "Transfer Learning for Auto-Gating of Flow Cytometry Data.”
### JMLR(workshop), 2012, 155–66.
############################################################
###
############################################################
###
############################################################
#########
######
###

#' Update classification bias on target tasks
#'
#' @description
#' Implements Algorithm 4 to refine decision boundaries on target tasks by
#' estimating density-informed bias adjustments. Supports optional data
#' manipulations and RCS (rare-cell subset) focusing.
#'
#' @param task_list List of target task feature matrices.
#' @param alg1_result Output list from [alg1_baselineClass()].
#' @param alg2_result Output list from [alg2_rob_meanNCov()].
#' @param alg3_result Vector of bias updates from [alg3_shiftComp()].
#' @param goodColumns Optional column indices to retain.
#' @param save2file Logical; persist diagnostic plots to disk.
#' @param save_dir Optional directory where diagnostic plots are written when
#'   `save2file` is `TRUE`. Defaults to a temporary directory.
#' @param Marg Margin parameter controlling neighborhood width.
#' @param alg4MinFx Mode selection for minima detection (`"gd"`, `"mean"`, `"win"`).
#' @param ADM Logical; apply additional data manipulations.
#' @param useMedian Logical; smooth minima search using medians.
#' @param ZnormMappingBL Logical; use z-normalized mapping for bias updates.
#' @param datatyp Character flag indicating data type (e.g., `"FC"`).
#' @param RCSmodeBL Logical; enable rare-cell subset adjustments.
#' @param RCSfreqSet Numeric vector specifying frequency thresholds for RCS mode.
#' @param CoreClassifier Core classifier identifier (e.g., `"LinSVM"`).
#' @param verbose Logical; print progress updates.
#' @param use_parallel Logical; parallelize per-task updates when `TRUE`.
#' @param parallel_cores Optional integer overriding detected core count.
#' @param wide_data_threshold Integer; coerce wide feature sets via `data.table`
#'   to reduce copies during mapping.
#'
#' @return A list with updated predictions, density summaries, and adjusted
#'   intercepts for each target task.
#' @export
alg4_BiasUpdate <- function(task_list, alg1_result, alg2_result,
                            alg3_result, goodColumns,
                            save2file = FALSE, Marg, alg4MinFx, ADM=F,
                            useMedian = T, ZnormMappingBL=F, datatyp="FC",
                            RCSmodeBL = F, RCSfreqSet = c(0,0),
                            CoreClassifier = "LinSVM", verbose = FALSE, save_dir = NULL,
                            use_parallel = FALSE, parallel_cores = NULL, wide_data_threshold = 200){



  # #Generic inputs for Dev needs; see kfold code as well
  # task_list = lapply(TestXY.ls, function(x){x$X.test});
  # alg1_result = alg1_res;
  # alg2_result = alg2_res;
  # alg3_result = alg3_res;
  # goodColumns = "";
  # Marg = 1; #[.1; 1]
  # alg4MinFx = "gd"; #mean; gd; win
  # save2file = F; RCSmodeBL = F;CoreClassifier = "LinSVM"; ADM=F; ZnormMappingBL = F; datatyp="FC"





  verbose <- isTRUE(verbose)
  save2file <- isTRUE(save2file)
  if (save2file && is.null(save_dir)) {
    save_dir <- file.path(tempdir(), "RTLbase_alg4")
  }
  if (save2file) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  }

  palette <- ColorTheme()
  col_vector <- palette$col_vector

  n_testSets = length(task_list)

  alg3_intercepts <- if (inherits(alg3_result, "rtl_alg3_result")) alg3_result$intercepts else alg3_result
  if (is.null(alg3_intercepts)) {
    stop("alg3_result must supply intercept updates.", call. = FALSE)
  }
  if (length(alg3_intercepts) < n_testSets) {
    stop("alg3_result must contain an intercept for each task.", call. = FALSE)
  }


  if(verbose) {
    message("starting Alg 4, Bias Update Base Version. No. of test sets:")
    message(n_testSets)
  }


  task_hat_outputs <- map_with_backend(seq_len(n_testSets), use_parallel = use_parallel, parallel_cores = parallel_cores, FUN = function(m){

    TASK <- coerce_feature_frame(task_list[[m]], goodColumns, wide_data_threshold = wide_data_threshold)
    if(anyNA(TASK)) stop("Task data contains missing values; please clean inputs before bias update.")

    if(ADM) TASK <- coerce_feature_frame(AllDataManipulations(TASK,
                                                              Scale=ScalePerCh,
                                                              Center=MaxPerCh,
                                                              X_cols2Keep = goodColumns,
                                                              globalASINhTrans = ifelse(datatyp=="FC", T, F),
                                                              globalRange1010Scale = F,
                                                              globalScaleByVal = F)[,goodColumns],
                                        goodColumns, wide_data_threshold = wide_data_threshold)

    if(ncol(TASK)==1){
      TASK <- as.data.frame(cbind(TASK, rep(-1, length(TASK))))
      colnames(TASK) <- c("X", "interc_multp")
    } else {
      TASK$interc_multp <- rep(-1, nrow(TASK))
    }




    if(CoreClassifier == "LinSVM"){
      task_hat_norm_upd <- as.matrix(TASK) %*% c(-alg2_result$U_robust[-length(alg2_result$U_robust)],
                                                 alg3_intercepts[m])
    }


    tempYhat.orgi <- task_hat_norm_upd

    if(RCSmodeBL) {
      xyDF <- as.data.frame(cbind(density(task_hat_norm_upd, n = nrow(task_hat_norm_upd))$x, density(task_hat_norm_upd, n = nrow(task_hat_norm_upd))$y))
      colnames(xyDF) <- c("x", "y")
      rownames(xyDF) <- rownames(task_hat_norm_upd)[order(task_hat_norm_upd)]
      DescTools::AUC(x=xyDF$x, y=xyDF$y)

      aprioriMedianFreq <- max(RCSfreqSet)

      aprioriMedianFreq <- ifelse(aprioriMedianFreq < 0.05, 0.05, aprioriMedianFreq)


      thresholds <- seq(min(xyDF$x), max(xyDF$x), length.out = 100)

      AUCthr <- as.numeric(lapply(thresholds, function(xcut){
        tempCutSubDF <- subset(xyDF, x > xcut)
        DescTools::AUC(x=tempCutSubDF$x, y=tempCutSubDF$y)
      }))

      xyDF.sub <- subset(xyDF, x >= thresholds[max(which(AUCthr >= aprioriMedianFreq))] )
      xyDF.sub$x.true <- task_hat_norm_upd[rownames(xyDF.sub),]

      if(length(xyDF.sub$y) > 3){
        if (save2file){
          fileID = file.path(save_dir, paste0("alg4_RCS_DensFocPT_", m, ".png"))
          png(fileID, width = 1024*2, height = 768*2, units = "px")
          par(mfrow=c(2,2))
          plot(xyDF, xlab="yhat", ylab="density", main=paste("PT thresholded set in red\nMarginal Freq from training: %", round(aprioriMedianFreq*100,3), sep=""), pch=19, col="grey")
          lines(xyDF, col="navy")
          points(xyDF.sub, pch=20, col="red")
          abline(v=thresholds[max(which(AUCthr >= aprioriMedianFreq))], lw=2, lty=3)
          plot(thresholds, AUCthr, pch=20, col="skyblue", main="RCS thresholding", ylab="Marginal Freq.", ylim=range(0,1), xlab="yhat")
          abline(v=thresholds[max(which(AUCthr >= aprioriMedianFreq))], lw=2, lty=3)
          abline(h=aprioriMedianFreq, lw=2, lty=3)
          plot(xyDF.sub[,1:2], main="focused region", col="red", lwd=2, pch=20)
          plot(density(xyDF.sub$x.true), main="density-curve of focused region", col="red", lwd=2)


          dev.off()
        }
        task_hat_norm_upd <- xyDF.sub$x.true
        remove(tempYhat.orgi)

      } else{
        task_hat_norm_upd <- tempYhat.orgi
        remove(tempYhat.orgi)}



    }


    z.norm <- (task_hat_norm_upd-mean(task_hat_norm_upd))/sd(task_hat_norm_upd) ## standardized data
    list(task_hat_norm_upd = task_hat_norm_upd, task_hat_Z = z.norm)
  })

  task_hat_norm_upd <- lapply(task_hat_outputs, function(x) x$task_hat_norm_upd)
  task_hat_Z <- lapply(task_hat_outputs, function(x) x$task_hat_Z)

  bias_outputs <- map_with_backend(seq_len(n_testSets), use_parallel = use_parallel, parallel_cores = parallel_cores, FUN = function(tn){
    if(verbose) message(paste("Starting Bias Update on test set #", tn, sep=""))

    if(!ZnormMappingBL) z_i <- task_hat_norm_upd[[tn]]
    if(ZnormMappingBL) z_i <- task_hat_Z[[tn]]

    s_j <- sort(z_i) #s_i ...grid of biases
    names(s_j) <- 1:length(s_j)

    w_t_euc_mag <- alg2_result$w_euc_mag #euclidean norm of base

    distance_grid <- abs(outer(z_i, s_j, "-"))/w_t_euc_mag
    c_j <- colSums(distance_grid < Marg)

    names(c_j) <- 1:length(c_j)
    if(verbose) message("counting within margin complete")

    h_bandwidth_val <- round(KBand_fx(s_k=s_j, c_k=c_j), 9)

    OneD_OptiBW <- density(s_j, n=length(s_j))$bw
    MeanBW = mean(c( h_bandwidth_val, OneD_OptiBW))


    Gaus_Ker_Smooth_sj_cj = ksmooth(s_j,c_j,kernel="normal", h_bandwidth_val,
                                    n.points = max(c(1000L, length(s_j))))



    Gaus_Ker_Smooth_sj_cj.optbw = ksmooth(s_j,c_j,kernel="normal", OneD_OptiBW,
                                          n.points = max(c(1000L, length(s_j))))
    Gaus_Ker_Smooth_sj_cj.meanbw = ksmooth(s_j,c_j,kernel="normal", MeanBW,
                                           n.points = max(c(1000L, length(s_j))))





    Gaus_Ker_Smooth_sj_cj.DF         <- as.data.frame(Gaus_Ker_Smooth_sj_cj)
    Gaus_Ker_Smooth_sj_cj.optbw.DF   <- as.data.frame(Gaus_Ker_Smooth_sj_cj.optbw)
    Gaus_Ker_Smooth_sj_cj.meanbw.DF  <- as.data.frame(Gaus_Ker_Smooth_sj_cj.meanbw)


    Gaus_Ker_Smooth_sj_cj.DF        <- Gaus_Ker_Smooth_sj_cj.DF[complete.cases(Gaus_Ker_Smooth_sj_cj.DF), ]
    Gaus_Ker_Smooth_sj_cj.optbw.DF  <- Gaus_Ker_Smooth_sj_cj.optbw.DF[complete.cases(Gaus_Ker_Smooth_sj_cj.optbw.DF), ]
    Gaus_Ker_Smooth_sj_cj.meanbw.DF <- Gaus_Ker_Smooth_sj_cj.meanbw.DF[complete.cases(Gaus_Ker_Smooth_sj_cj.meanbw.DF), ]


    Gaus_Ker_Smooth_sj_cj.DF         <- Gaus_Ker_Smooth_sj_cj.DF[ min( which ( Gaus_Ker_Smooth_sj_cj.DF$y != 0 )):max( which( Gaus_Ker_Smooth_sj_cj.DF$y != 0 )), ]
    Gaus_Ker_Smooth_sj_cj.optbw.DF   <- Gaus_Ker_Smooth_sj_cj.optbw.DF[ min( which ( Gaus_Ker_Smooth_sj_cj.optbw.DF$y != 0 )):max( which( Gaus_Ker_Smooth_sj_cj.optbw.DF$y != 0 )), ]
    Gaus_Ker_Smooth_sj_cj.meanbw.DF  <- Gaus_Ker_Smooth_sj_cj.meanbw.DF[ min( which ( Gaus_Ker_Smooth_sj_cj.meanbw.DF$y != 0 )):max( which( Gaus_Ker_Smooth_sj_cj.meanbw.DF$y != 0 )), ]

    if(datatyp == "FC"){
      MinDensThreshold = round(max(Gaus_Ker_Smooth_sj_cj.DF$y)*ifelse(max(Gaus_Ker_Smooth_sj_cj.DF$y)>100, 0.1, 0.01))


      Gaus_Ker_Smooth_sj_cj.DF         <- subset(Gaus_Ker_Smooth_sj_cj.DF, y > MinDensThreshold)
      Gaus_Ker_Smooth_sj_cj.optbw.DF   <- subset(Gaus_Ker_Smooth_sj_cj.optbw.DF, y > MinDensThreshold)
      Gaus_Ker_Smooth_sj_cj.meanbw.DF  <- subset(Gaus_Ker_Smooth_sj_cj.meanbw.DF, y > MinDensThreshold)


    }
    if(datatyp == "scRNASeqLogNormCount"){

      if(RCSmodeBL) {
        MinDensThreshold = 1

      } else {
        MinDensThreshold = round( max(Gaus_Ker_Smooth_sj_cj.DF$y)*ifelse(max(Gaus_Ker_Smooth_sj_cj.DF$y)>1000, 0.1, 0.05))

      }

      if(nrow(subset(Gaus_Ker_Smooth_sj_cj.DF, y > MinDensThreshold))>2) Gaus_Ker_Smooth_sj_cj.DF         <- subset(Gaus_Ker_Smooth_sj_cj.DF, y > MinDensThreshold)
      if(nrow(subset(Gaus_Ker_Smooth_sj_cj.optbw.DF, y > MinDensThreshold))>2) Gaus_Ker_Smooth_sj_cj.optbw.DF   <- subset(Gaus_Ker_Smooth_sj_cj.optbw.DF, y > MinDensThreshold)
      if(nrow(subset(Gaus_Ker_Smooth_sj_cj.meanbw.DF, y > MinDensThreshold))>2) Gaus_Ker_Smooth_sj_cj.meanbw.DF  <- subset(Gaus_Ker_Smooth_sj_cj.meanbw.DF, y > MinDensThreshold)


    }




    Gaus_Ker_Smooth_sj_cj.DF        <- try(as.data.frame(approx(Gaus_Ker_Smooth_sj_cj.DF, n=nrow(Gaus_Ker_Smooth_sj_cj.DF))))
    Gaus_Ker_Smooth_sj_cj.optbw.DF  <- try(as.data.frame(approx(Gaus_Ker_Smooth_sj_cj.optbw.DF, n=nrow(Gaus_Ker_Smooth_sj_cj.optbw.DF))))
    Gaus_Ker_Smooth_sj_cj.meanbw.DF <- try(as.data.frame(approx(Gaus_Ker_Smooth_sj_cj.meanbw.DF, n=nrow(Gaus_Ker_Smooth_sj_cj.meanbw.DF))))




    if(class(Gaus_Ker_Smooth_sj_cj.DF)=="try-error"){
      if(class(Gaus_Ker_Smooth_sj_cj.optbw.DF)=="try-error"){
        if(class(Gaus_Ker_Smooth_sj_cj.meanbw.DF)=="try-error"){
          Gaus_Ker_Smooth_sj_cj.DF         <- as.data.frame(Gaus_Ker_Smooth_sj_cj)
        } else{
          Gaus_Ker_Smooth_sj_cj.DF         <- Gaus_Ker_Smooth_sj_cj.meanbw.DF
        }
      } else {

        Gaus_Ker_Smooth_sj_cj.DF         <- Gaus_Ker_Smooth_sj_cj.optbw.DF

      }

    }


    SMmin.try   <- try(SmartMinimaAk(Gaus_Ker_Smooth_sj_cj.DF, print2screen = T, learnRate=ifelse(datatyp=="FC", 0.0000001 ,0.00001), mode=ifelse(RCSmodeBL, "MedianAllMin", ifelse(datatyp=="FC","MedianAll", "MedianAll")) ), silent = T)


    if(class(SMmin.try)=="try-error"){
      SMmin.try   <- try(SmartMinimaAk(Gaus_Ker_Smooth_sj_cj.DF, print2screen = F, learnRate=ifelse(datatyp=="FC", 0.0000001 ,0.0001), mode=ifelse(RCSmodeBL, "LeftOfPeakNPavg", ifelse(datatyp=="FC","MedianAll", "MedianAll")) ), silent = T)
      if(class(SMmin.try)=="try-error"){
        SMmin.try   <- try(SmartMinimaAk(Gaus_Ker_Smooth_sj_cj.DF, print2screen = F, learnRate=ifelse(datatyp=="FC", 0.0000001 ,0.0001), mode=ifelse(RCSmodeBL, "NP", ifelse(datatyp=="FC","NP", "NP")) ), silent = T)
      }
    }

    if(class(SMmin.try)=="try-error"){

      SMmin.meanbw.try  <- try(SmartMinimaAk(Gaus_Ker_Smooth_sj_cj.meanbw.DF, print2screen = T, learnRate=ifelse(datatyp=="FC", 0.0000001 ,0.00001), mode=ifelse(RCSmodeBL, "LeftOfPeakNPavg", ifelse(datatyp=="FC","NPZAvg", "NP")) ), silent = T)

      if(class(SMmin.meanbw.try)=="try-error"){

        SMmin.optbw.try   <- try(SmartMinimaAk(Gaus_Ker_Smooth_sj_cj.optbw.DF, print2screen = F, learnRate=ifelse(datatyp=="FC", 0.0000001 ,0.00001), mode=ifelse(RCSmodeBL, "LeftOfPeak", ifelse(datatyp=="FC","NPZAvg", "NPZAvg")) ), silent = T)


        if(class(SMmin.optbw.try)=="try-error"){
          SMmin.heur <- try(Gaus_Ker_Smooth_sj_cj.DF[find_peaks(-1*Gaus_Ker_Smooth_sj_cj.DF$y),])

          if(class(SMmin.heur)=="try-error"){
            SMmin.heur <- Gaus_Ker_Smooth_sj_cj.DF[which.min(abs(Gaus_Ker_Smooth_sj_cj.DF$x - median(Gaus_Ker_Smooth_sj_cj.DF$x))),]
          } else {
            if(nrow(SMmin.heur)==0) {
              SMmin.heur <- Gaus_Ker_Smooth_sj_cj.DF[which.min(abs(Gaus_Ker_Smooth_sj_cj.DF$x - median(Gaus_Ker_Smooth_sj_cj.DF$x))),]
            }
          }

          SMmin.ls <- list(PeakA = Gaus_Ker_Smooth_sj_cj.DF[which.max(Gaus_Ker_Smooth_sj_cj.DF$y),],
                           PeakB = Gaus_Ker_Smooth_sj_cj.DF[which.max(Gaus_Ker_Smooth_sj_cj.DF$y),],
                           Minima = SMmin.heur)

        } else{
          SMmin.ls <- SMmin.optbw.try
        }
      } else {
        SMmin.ls <- SMmin.meanbw.try
      }
    } else {
      SMmin.ls <- SMmin.try
    }


    if(save2file){
      fileID = file.path(save_dir, paste0("alg4_KeyOptima_", names(task_list)[tn], "_", round(h_bandwidth_val,4), ".png"))
      png(fileID, width = 1024, height = 768, units = "px")

      plot(as.data.frame(cbind(s_j, c_j)), typ="l", lwd=2, col="grey")
      points((Gaus_Ker_Smooth_sj_cj.DF), col="red", pch=19)
      abline(h=MinDensThreshold, lty=3, col="navy")

      lines(Gaus_Ker_Smooth_sj_cj.DF, type="l", col="dodgerblue", lwd=2,
            main=paste("Key Local Optima on test sample: ",names(task_list)[tn] , sep=""),
            xlim=range(min(Gaus_Ker_Smooth_sj_cj$x),max(Gaus_Ker_Smooth_sj_cj$x)));
      abline(v=-alg3_intercepts[tn], col="dodgerblue", lwd=2)
      abline(v=-alg2_result$U_robust["b.int"], col="forestgreen", lty=2, lwd=2);
      abline(v=SMmin.ls$PeakA, col="orange", lwd=2);
      abline(v=SMmin.ls$PeakB, col="orange", lwd=2);
      abline(v=SMmin.ls$Minima$x, col="purple", lwd=2, lty=3)


      legend("topleft",
             c("alg2_robustMean", "alg3_shifted", "Peak Positions", "Found Local minima", "Counts above 1% of total"),
             lwd=c(2,2,2,2,2),
             col=c("forestgreen", "dodgerblue","orange","purple", "red"),
             lty=c(2,2,2,3,2))

      dev.off()
    }


    s_j_KeyOptima <- c(LocMin = SMmin.ls$Minima, peakYB = SMmin.ls$PeakB, peakYA = SMmin.ls$PeakA )
    names(s_j_KeyOptima) <- c("LocMin","peakYB","peakYA")

    as.data.frame(cbind(alg3_b_upd=alg3_intercepts[tn], alg4_b_upd = (s_j_KeyOptima[c("LocMin")])))
  })

  SVM_b_all.df <- rbindlist(bias_outputs)


  #correction based on position
  SVM_b_all.df$alg4_b_upd <- as.numeric(SVM_b_all.df$alg3_b_upd) + as.numeric(SVM_b_all.df$alg4_b_upd)

  SVM_b_all.df <- as.data.frame( cbind(rep(alg2_result$U_robust[length(alg2_result$U_robust)],
                                           nrow(SVM_b_all.df)), SVM_b_all.df))
  colnames(SVM_b_all.df) <- c("b_alg2_norm", "b_alg3_norm", "b_alg_norm")
  rownames(SVM_b_all.df) <- paste("Task", 1:nrow(SVM_b_all.df), sep="")

  class(SVM_b_all.df) <- c("rtl_alg4_result", class(SVM_b_all.df))
  attr(SVM_b_all.df, "metadata") <- list(margin = Marg, core_classifier = CoreClassifier)

  return(SVM_b_all.df)


}
