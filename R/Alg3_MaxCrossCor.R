
#Mahyari, Eisa: Implementation of Lee et al.'s TL-FC into R
#March, 2017

#computing the kernel density estimate (KDE) of the projections onto w_0
#align the target data to each source data using maximum cross-correlation
# in a leave-one-out setting

##### Shift Compensation (Algorithm 3) ##########
###
### Input: hyperplane (w, b),
###    source task data {D_m} for m=1 to M, and
###    target task data Tau
###
###    z_t_i <- <w,x_m_i> + b, for all i
###
###    for m = 1 to M do
###        z_m_i <- <w,x_m_i> + b, for all i
###        e_m <- argmax_z KDE(z,z_t_i) *max_cross-cor* KDE(z, z_m_i)
###    end for
###
###    b <- b - median(e_m)
###
### Output: b
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

ccfmaxV3 <- function(a, b, e=0, maxLag, useAbsCor = T)
{
  # a = c(rnorm(5000,1,1), rnorm(5000,6,1))
  # b = c(rnorm(5000,2,1), rnorm(5000,8,1))
  # plot(density(a))
  # lines(density(b))
  # maxLag = 10


  #if(useAbsCor) print("AbsCor")
  #d <- ccf(a, b, plot = FALSE, lag.max = length(a)/2)
  d <- ccf(a, b, plot = FALSE, lag.max = maxLag)
  cor = round(d$acf[,,1], 3)
  if(useAbsCor) cor = abs(cor)
  lag = d$lag[,,1]
  res = data.frame(cor, lag)
  res_max = res[which.max(res$cor),]



  return(res_max)
}

#' Apply shift compensation using maximum cross-correlation
#'
#' @description
#' Aligns target task projections to source datasets by estimating robust lag
#' shifts via maximum cross-correlation between kernel density estimates. The
#' updated intercepts are returned for downstream bias correction.
#'
#' @param source_list List of source feature matrices used to train the
#'   classifiers.
#' @param task_list List of target task feature matrices to be aligned.
#' @param alg1_result Output list from [alg1_baselineClass()].
#' @param alg2_result Output list from [alg2_rob_meanNCov()].
#' @param verbose Logical; print progress updates.
#' @param save2file Logical; save diagnostic plots to disk.
#' @param save_dir Optional directory where diagnostic plots are written when
#'   `save2file` is `TRUE`. Defaults to a temporary directory.
#' @param maximumLag Maximum lag to consider in cross-correlation.
#' @param ImpFeats Optional vector of feature indices to retain.
#' @param ADM Logical; apply additional data manipulations.
#' @param datatyp Character flag indicating data type (e.g., `"FC"`).
#' @param useAbsCor Logical; use absolute correlations when selecting lags.
#' @param medianMediansBL Logical; use the median of medians strategy for shift
#'   estimation.
#' @param CoreClassifier Core classifier identifier (e.g., `"LinSVM"`).
#' @param use_parallel Logical; parallelize per-task shift estimation when `TRUE`.
#' @param parallel_cores Optional integer overriding the detected core count.
#' @param wide_data_threshold Integer; when data are wide, coerce with
#'   `data.table` to reduce copies before projection.
#'
#' @return A structured list containing updated intercepts per task and
#'   cross-correlation diagnostics.
#' @export
alg3_shiftComp <- function(source_list, task_list, alg1_result, alg2_result,
                           save2file = FALSE, maximumLag = 0, ImpFeats = NA,
                           ADM=F, datatyp="FC", useAbsCor = T, medianMediansBL = F,
                           CoreClassifier, use_parallel = FALSE, parallel_cores = NULL,
                           wide_data_threshold = 200, verbose = FALSE, save_dir = NULL,
                           print2screen = NULL){

  # task_list      = TestXls.t;
  # source_list    =  TrainXls.t;
  # alg2_result    =  alg2_res;
  # alg1_result    = alg1_res
  # print2screen   =  T;
  # save2file      =  F;
  # maximumLag     =  0;
  # ImpFeats = ImpFeats

  verbose <- isTRUE(verbose) || isTRUE(print2screen)
  n_testSets = length(task_list)
  M = length(source_list)
  palette <- ColorTheme()
  col_vector <- palette$col_vector

  if (length(source_list) != length(alg1_result$baselineSVM) && !is.null(alg1_result$baselineSVM)) {
    warning("Source list length does not match alg1_result; proceeding with provided data.", call. = FALSE)
  }

  if (n_testSets == 0 || M == 0) {
    stop("source_list and task_list must contain at least one entry.", call. = FALSE)
  }

  if (save2file && is.null(save_dir)) {
    save_dir <- file.path(tempdir(), "RTLbase_alg3")
  }
  if (save2file) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  }
  log_dir <- if (is.null(save_dir)) tempdir() else save_dir

  if(verbose) message("starting Alg 3, shift compensation")

  #alg1task_hat <- vector(mode="list")
  #alg1source_hat <- vector(mode="list")

  w_vec.alg2 <- alg2_result$U_robust[1:(length(alg2_result$U_robust)-1)]
  bint.alg2  <- alg2_result$U_robust[length(alg2_result$U_robust)]



  ###########################CLASSIFIER Selection
  if(CoreClassifier=="LinSVM") hyp.alg2  <- -c(w_vec.alg2, b.int=bint.alg2)
  ###########################CLASSIFIER Selection



  source_hat <- map_with_backend(seq_len(M), use_parallel = use_parallel, parallel_cores = parallel_cores, FUN = function(m){
    if(verbose) message(paste("Inference of class by the datasets's baseline hyperplane #", m,sep=""))

    SOURCE <- coerce_feature_frame(source_list[[m]], ImpFeats, wide_data_threshold = wide_data_threshold)

    if(ADM) {
      SOURCE_proc <- AllDataManipulations(SOURCE,
                                          Scale=ScalePerCh,
                                          Center=MaxPerCh,
                                          globalASINhTrans = ifelse(datatyp=="FC", T, F),
                                          globalRange1010Scale = F,
                                          globalScaleByVal = F)
      if(!is.na(ImpFeats[1])) SOURCE_proc <- SOURCE_proc[,ImpFeats]
      SOURCE <- coerce_feature_frame(SOURCE_proc,
                                     ImpFeats, wide_data_threshold = wide_data_threshold)
    }

    if(ncol(SOURCE)==1){
      SOURCE <- as.data.frame(cbind(SOURCE, rep(ifelse(CoreClassifier=="LinSVM", -1, 1), length(SOURCE))))
      colnames(SOURCE) <- c("X", "interc_multp")
    } else {
      SOURCE$interc_multp <- rep(ifelse(CoreClassifier=="LinSVM", -1, 1), nrow(SOURCE))
    }
    as.data.frame(as.matrix(SOURCE) %*% hyp.alg2)
  })
  names(source_hat) <- names(source_list)

  task_hat <- map_with_backend(seq_len(n_testSets), use_parallel = use_parallel, parallel_cores = parallel_cores, FUN = function(j){
    if(verbose) message(paste("Mapping to Y_hat for Task: ", j, sep=""))

    TASK <- coerce_feature_frame(task_list[[j]], ImpFeats, wide_data_threshold = wide_data_threshold)

    if(ADM) {
      TASK_proc <- AllDataManipulations(TASK,
                                        Scale=ScalePerCh,
                                        Center=MaxPerCh,
                                        globalASINhTrans = ifelse(datatyp=="FC", T, F),
                                        globalRange1010Scale = F,
                                        globalScaleByVal = F)
      if(!is.na(ImpFeats[1])) TASK_proc <- TASK_proc[,ImpFeats]
      TASK <- coerce_feature_frame(TASK_proc,
                                   ImpFeats, wide_data_threshold = wide_data_threshold)
    }

    if(ncol(TASK)==1){
      TASK <- as.data.frame(cbind(TASK, rep(ifelse(CoreClassifier=="LinSVM", -1, 1), length(TASK))))
      colnames(TASK) <- c("X", "interc_multp")
    } else {
      TASK$interc_multp <- rep(ifelse(CoreClassifier=="LinSVM", -1, 1), nrow(TASK))
    }

    as.data.frame(as.matrix(TASK) %*% hyp.alg2)
  })
  names(task_hat) <- names(task_list)

  if(verbose) message("step3 reached, task and source are mapped, starting Max-Cross Corr")

  per_task_results <- map_with_backend(seq_len(n_testSets), use_parallel = use_parallel, parallel_cores = parallel_cores, FUN = function(t){
    if(verbose) message(paste("starting cross-cor for Task: ", t, sep=""))
    maxLagLocal <- maximumLag
    if (maxLagLocal == 0){
      if(length(source_hat[[1]][,1])<50){
        maxLagLocal   =  length(source_hat[[1]][,1])
      } else{
        maxLagLocal   =  length(source_hat[[1]][,1])/10
      }

    }

    e <- vector()
    cor_e <- vector()
    MS <- vector()

    for (i in 1:M) {

      DensA <- density(task_hat[[t]][,1], n=length(task_hat[[t]][,1]))
      bwsig = DensA$bw
      if(verbose) message(paste("with source: ", i, sep=""))

      CorThresh = 0.9
      rep4Roba = 3

      RobustCrossCorLagALL <- lapply(1:rep4Roba, function(id){
        set.seed(abs(round(rnorm(id)*30000))[id])

        DensB1 <- density(source_hat[[i]][,1][sample(1:length(source_hat[[i]][,1]), min(c(length(source_hat[[i]][,1]), length(task_hat[[t]][,1]))))], bw=bwsig, n=length(task_hat[[t]][,1]))

        e_resamp <- ccfmaxV3(DensA$y,
                             DensB1$y,
                             maxLag = maxLagLocal, useAbsCor = useAbsCor);
        if(e_resamp$cor >= CorThresh) return(e_resamp)

      })

      RobustCrossCorLag <- as.numeric(as.character(lapply(RobustCrossCorLagALL, function(x){
        ifelse(is.null(x$lag), NA, x$lag)
      })))

      RobustCrossCorLag <- round(mean(RobustCrossCorLag[!is.na(RobustCrossCorLag)]))

      while(is.na(RobustCrossCorLag)){
        if(verbose) message(paste("robust cross-cor failed w. threshold: ", CorThresh, sep=""))
        CorThresh <- CorThresh - 0.01
        rep4Robb = 10
        RobustCrossCorLagALL <- lapply(1:rep4Robb, function(id){
          set.seed(abs(round(rnorm(1)*30000)))
          DensB1 <- density(source_hat[[i]][,1][sample(1:length(source_hat[[i]][,1]), min(c(length(source_hat[[i]][,1]), length(task_hat[[t]][,1]))))], bw=bwsig, n=length(task_hat[[t]][,1]))
          e_resamp <- ccfmaxV3(DensA$y,
                               DensB1$y,
                               maxLag = maxLagLocal, useAbsCor = useAbsCor);
          if(e_resamp$cor > CorThresh) return(e_resamp)

        });
        if(verbose) message(CorThresh)

        RobustCrossCorLag.t <- as.numeric(as.character(lapply(RobustCrossCorLagALL, function(x){
          ifelse(is.null(x$lag), NA,x$lag)
        })))

        RobustCrossCorLag <- Mode(RobustCrossCorLag.t[!is.na(RobustCrossCorLag.t)])

      }
      if(verbose) message(paste("robust cross-cor found @ threshold: ", CorThresh, sep=""))



      RobustCrossCor <- as.numeric(as.character(lapply(RobustCrossCorLagALL, function(x){
        x$cor
      })))
      RobustCrossCor <- Mode(RobustCrossCor[!is.na(RobustCrossCor)])

      DensB1 <- density(source_hat[[i]][,1], n=length(source_hat[[i]][,1]))


      if(RobustCrossCorLag <  0) deltaLag = (DensB1$x[abs(RobustCrossCorLag)] - DensA$x[1])
      if(RobustCrossCorLag >  0) deltaLag = -(DensB1$x[abs(RobustCrossCorLag)] - DensA$x[1])
      if(RobustCrossCorLag == 0) deltaLag = 0

      MedianShift <- median(task_hat[[t]][,1]) - median(source_hat[[i]][,1])


      if (save2file) {

        fileID <- file.path(save_dir, paste0("alg3Shifts_", "t_", t, "_s_", i, ".png"))

        png(fileID, width = 1024, height = 768, units = "px")

        par(mfrow=c(2,3))
        boxplot(list(task=task_hat[[t]][,1], source=source_hat[[i]][,1]),
                main=paste("task: ", t, " | source: ", i, "\n", names(task_hat)[t], " | ",
                           names(source_hat)[i], sep=""),
                ylab = "Prediction (y_hat)")

        ccfig <- ccf(DensA$y,
                     DensB1$y,
                     plot = F, lag.max = maxLagLocal)
        corTemp = round(ccfig$acf[,,1], 3)
        if(useAbsCor) corTemp = abs(corTemp)
        lagTemp = ccfig$lag[,,1]

        plot(data.frame(lagTemp, corTemp), main=paste("Maximum cross-correlation\n@ robust-lag:", RobustCrossCorLag, " | t_",t,"_s_",i,sep=""), xlab="Lag", ylab="abs(cor)")
        abline(v=RobustCrossCorLag, lwd=2, lty=2, col="red")

        plot(DensA, type='l',lwd=2, col="dodgerblue", main="pre-shift")
        lines(DensB1, lwd=2, col='firebrick')


        D2<- as.data.frame(cbind(x=DensA$x - deltaLag,
                                 y=DensA$y))

        plot(D2, lwd=2, type='l', col="dodgerblue",
             main=paste("post-shift w/ delta : ",round(deltaLag,3),
                        "\nMax(abs(Cor)) >= thr",
                        " and lag : ", round(RobustCrossCorLag,3),sep=""))
        lines(DensB1, lwd=2, col="firebrick")


        D2<- as.data.frame(cbind(x=DensA$x - MedianShift, y=DensA$y))

        plot(D2, type='l', lwd=2, col="dodgerblue",
             main=paste("post-shift w/ delta : ",round(MedianShift,3),
                        "\nMedian-based"))
        lines(DensB1, lwd=2, col="firebrick")


        dev.off()

      }

      e[i] <- deltaLag
      cor_e[i] <- RobustCrossCor
      MS[i] <- MedianShift
    }
    list(e = e, cor_e = cor_e, MS = MS)
  })

  e_med_mat <- lapply(per_task_results, function(res) res$e)
  cor_e_med_mat <- lapply(per_task_results, function(res) res$cor_e)
  MS_med_mat <- lapply(per_task_results, function(res) res$MS)

  #names(e_med_mat) <- names(task_hat)
  #names(cor_e_med_mat) <- names(task_hat)
  #names(MS_med_mat) <- names(task_hat)

  medianShifts <- as.numeric(lapply(e_med_mat, function(x){
    median(x, na.rm = T)}))
  medianMedianShifts <- as.numeric(lapply(MS_med_mat, median))
  e_med_matDF <- as.data.frame(e_med_mat)
  #rownames(e_med_matDF) <- names(source_hat)
  cor_e_med_matDF <- as.data.frame(cor_e_med_mat)
  #rownames(cor_e_med_matDF) <- names(source_hat)
  MS_med_matDF <- as.data.frame(MS_med_mat)
  #rownames(MS_med_matDF) <- names(source_hat)

  remove(t)
  if (save2file) {


    fileID <- file.path(save_dir, "alg3ShiftsCombo_.png")

    png(fileID, width = 1024, height = 768, units = "px")
    par(mfrow=c(1,3))
    i=1
    D2<- as.data.frame(cbind(x=density(task_hat[[i]][,1], bw=bwsig)$x, y=density(source_hat[[i]][,1], bw=.1)$y))
    plot(D2, col=col_vector[i], type='l', main="Task_hats pre-shift", lwd=2)
    if(length(task_hat)>1){
      for (i in 2:length(medianShifts)){
        D2<- as.data.frame(cbind(x=density(task_hat[[i]][,1], bw=bwsig)$x, y=density(task_hat[[i]][,1], bw=.1)$y))
        lines(D2, col=col_vector[i], type='l', main="Task_hats\nmedian shifts", lwd=2)
      }
    }


    i=1
    D2<- as.data.frame(cbind(x=density(task_hat[[i]][,1], bw=bwsig)$x-medianShifts[i], y=density(source_hat[[i]][,1], bw=.1)$y))

    plot(D2, col=col_vector[i], type='l', main="Task_hats\nmedian of lag-shifts", lwd=2)
    if(length(task_hat)>1){
      for (i in 2:length(medianShifts)){
        D2<- as.data.frame(cbind(x=density(task_hat[[i]][,1], bw=.5)$x-medianShifts[i], y=density(task_hat[[i]][,1], bw=.1)$y))
        lines(D2, col=col_vector[i], type='l', main="Task_hats\nmedian of lag-shifts", lwd=2)
      }
    }


    i=1
    D2<- as.data.frame(cbind(x=density(task_hat[[i]][,1], bw=.5)$x-medianMedianShifts[i], y=density(source_hat[[i]][,1], bw=.1)$y))
    plot(D2, col=col_vector[i], type='l', main="Task_hats\nmedian of Median-shifts", lwd=2, xlim=range(-10:10))
    if(length(task_hat)>1){
      for (i in 2:length(medianMedianShifts)){
        D2<- as.data.frame(cbind(x=density(task_hat[[i]][,1], bw=.5)$x-medianMedianShifts[i], y=density(task_hat[[i]][,1], bw=.1)$y))
        lines(D2, col=col_vector[i], type='l', main="Task_hats\nmedian of Median-shifts", lwd=2)
      }
    }

    dev.off()

    e_med_matDF$median <- apply(e_med_matDF, 1, median)
    print(xtable(t(cor_e_med_matDF), caption='Bias Shift Correlations in Test (rows) Vs. Training (columns) ',
                 label='tab:cor_shifts', digits=3),
          file=paste(log_dir,
                     "/cor_shifts.txt", sep=""))
    print(xtable(t(e_med_matDF), caption='Bias Shift Values in Test (rows) Vs. Training (columns) ',
                 label='tab:cor_shifts', digits=3),
          file=paste(log_dir,
                     "/val_shifts.txt", sep=""))

    print(xtable(t(MS_med_matDF), caption='Bias Shift Values in Test (rows) Vs. Training (columns) ',
                 label='tab:medMed_shifts', digits=3),
          file=paste(log_dir,
                     "/val_Medshifts.txt", sep=""))




  }

  if(CoreClassifier=="LinSVM") {
    #b_updated <- alg2_result$U_robust_norm[length(alg2_result$U_robust_norm)] + -1*abs(as.numeric(lapply(e_med_mat, median)))
    b_updated <- hyp.alg2[length(alg2_result$U_robust)] - (as.numeric(lapply(e_med_mat, median)))
    names(b_updated) <- paste("Task", 1:length(b_updated), "_MedianOfTasknSourceLags", sep="")
    if(medianMediansBL){
      b_updated <- alg2_result$U_robust[length(alg2_result$U_robust)] - medianMedianShifts
    }
    ##as.numeric(lapply(e_med_mat, median)))
  }


  if(verbose) message("Shift Compesation complete for all tasks")
  if(verbose) message("the normed intercepts (b_t) have been updated by median lag in max cross-correlation")

  if(verbose) message("the shifts for each task:")
  if(verbose) message(as.numeric(lapply(e_med_mat, median)))
  if(verbose) message("Final intercepts for all tasks:")
  if(verbose) message(b_updated)

  result <- list(
    intercepts = b_updated,
    lag_matrix = as.data.frame(e_med_mat),
    correlation_matrix = as.data.frame(cor_e_med_mat),
    median_shifts = medianShifts,
    median_of_medians = medianMedianShifts,
    median_shift_matrix = as.data.frame(MS_med_mat),
    metadata = list(
      core_classifier = CoreClassifier,
      use_abs_cor = useAbsCor,
      median_of_medians = medianMediansBL,
      maximum_lag = maximumLag
    )
  )
  class(result) <- "rtl_alg3_result"
  return(result)

}
