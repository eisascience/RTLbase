




#' Train baseline SVM classifiers for multiple source datasets
#'
#' @description
#' Fits a linear SVM for each source dataset and returns the learned hyperplane
#' parameters along with cross-validation diagnostics. This is the entry point
#' for initializing the RTL workflow before transfer learning adjustments are
#' applied.
#'
#' @param TrainXls List of training feature matrices (one per source dataset).
#' @param TrainYls List of training labels aligned to `TrainXls`.
#' @param TestXls List of test feature matrices (one per target task).
#' @param TestYls List of test labels aligned to `TestXls`.
#' @param K_forCrossV Number of folds for internal SVM cross-validation.
#' @param svmGamma Radial basis gamma parameter passed to `e1071::svm`.
#' @param svmCost Cost parameter passed to `e1071::svm`.
#' @param verbose Logical; print progress to the console.
#' @param X_cols2Keep Column indices to retain from the feature matrices.
#' @param transX Logical; apply additional transformations via `RTL::AllDataManipulations`.
#' @param sampleRed Optional integer sample size for down-sampling during training.
#' @param doParalellSVM Logical; use a parallelized SVM fit when `TRUE`.
#' @param datatyp Character flag indicating data type (e.g., `"FC"` for flow cytometry).
#' @param use_parallel Logical; parallelize per-dataset training when `TRUE`.
#' @param parallel_cores Optional integer overriding the detected core count.
#' @param wide_data_threshold Integer; when the number of columns exceeds this
#'   threshold, inputs are coerced via `data.table` to minimize copies.
#'
#' @return A list containing the baseline hyperplanes, per-dataset results, and
#'   metadata (number of datasets and dimensions).
#' @export
alg1_baselineClass <- function(TrainXls,
                               TrainYls,
                               TestXls,
                               TestYls,
                               K_forCrossV,
                               svmGamma,
                               svmCost,
                               X_cols2Keep,
                               transX=F, sampleRed=F, doParalellSVM=F, datatyp="FC",
                               use_parallel = FALSE, parallel_cores = NULL,
                               wide_data_threshold = 200,
                               verbose = FALSE,
                               prnt2scr = NULL){


  #initializations
  errorFlag = F
  verbose <- isTRUE(verbose) || isTRUE(prnt2scr)


  #Source task data Dm for m = 1, ...., M datasets

  if (is.null(TrainXls) || is.null(TrainYls)) {
    stop("Training feature and label lists must be provided.", call. = FALSE)
  }

  M = length(TrainXls); if(verbose) message(paste("# of datasets = ", M, sep=""))

  if (length(TrainYls) != M) {
    stop("Training feature and label list lengths must match.", call. = FALSE)
  }



  #Assuming that all df will have the same # of cols;
  if(class(TrainXls[[1]]) == "numeric"){
    nDims = 1
    if(verbose) message("datasets have 1 dim/feat each")

  } else {
    nDims = ncol(as.data.frame((TrainXls[[1]]))[,X_cols2Keep])
    if(verbose) message(paste("datasets have ", nDims, " dims/feats each", sep=""))
  }


  #this will hold the hyperplane parameters for each SVM
  baselineSVM <- matrix(0, nrow=M, ncol=nDims+1)

  #empty vectors for storage of results
  kfoldSVM_ls    <- vector(mode="list")
  results.all    <- vector(mode="list")

  #-------------------end of inputs and initializations


  if(verbose) message(paste("Starting training with a SVM classifier || cost: ",
                            svmCost, ", gamma: ", svmGamma,
                            ", Kernel: linear" , ", cross: ", K_forCrossV, sep=""))

  train_indices <- seq_len(M)
  dataset_results <- map_with_backend(train_indices, use_parallel = use_parallel, parallel_cores = parallel_cores, FUN = function(m){
    if(verbose) message(m)

    Dm.train  <- coerce_feature_frame(TrainXls[[m]], X_cols2Keep, wide_data_threshold = wide_data_threshold)

    if(transX){
      Dm.train <- as.data.frame(RTL::AllDataManipulations(Dm.train,
                                                          Scale=MaxPerCh,
                                                          Center=ScalePerCh,
                                                          X_cols2Keep = X_cols2Keep,
                                                          globalASINhTrans = ifelse(datatyp=="FC", T, F),
                                                          globalRange1010Scale = F,
                                                          globalScaleByVal = F))[,X_cols2Keep]
      Dm.train <- coerce_feature_frame(Dm.train, X_cols2Keep, wide_data_threshold = wide_data_threshold)
    }

    if(!is.factor(TrainYls[[m]])) {
      TrainYls[[m]] <- factor(TrainYls[[m]])
    }

    label_levels <- levels(TrainYls[[m]])
    if("TRUE" %in% label_levels) {
      TrainYls[[m]] <- factor(ifelse(TrainYls[[m]]==T,1,-1))
      label_levels <- levels(TrainYls[[m]])
    }

    if(any(label_levels %in% c("Pos", "Neg"))) {
      POS <- c("Neg", "Pos")
    } else if(any(label_levels %in% c("1", "-1", 1, -1))) {
      POS <- c(-1, 1)
    } else {
      POS <- label_levels
    }

    if(ncol(Dm.train)==1){
      Dm.train <- as.data.frame(cbind(Dm.train, factor(TrainYls[[m]], levels = POS)))
      colnames(Dm.train) <- c("X", "labelSVM")
    } else {
      Dm.train$labelSVM        <- factor(TrainYls[[m]], levels = POS)
    }

    if(nDims == 1) colnames(Dm.train) <- c(paste(rep("X", nDims), 1:nDims, sep=""), "labelSVM")

    if(!(sampleRed==F)){
      if(nrow(Dm.train)>sampleRed){
        Dm.train <- Dm.train[sample(1:nrow(Dm.train), sampleRed),]
      }
    }

    if(doParalellSVM==T){
      Dm.train <<- Dm.train
      Dm.train_model <- parallelSVM(labelSVM ~. , data=Dm.train,
                                    cost = svmCost , gamma = svmGamma,
                                    type="C-classification",
                                    kernel = 'linear',
                                    cross = K_forCrossV,
                                    scale = F, numberCores = 4)
    } else {

      Dm.train_model <- svm(labelSVM ~. , data=Dm.train,
                            cost=svmCost, gamma=svmGamma,
                            type="C-classification",
                            kernel = 'linear',
                            cross = K_forCrossV,
                            scale = F)
    }


    Dm.train.pred  <- predict(Dm.train_model, Dm.train)

    if("Pos" %in% levels(Dm.train.pred)) POS = "Pos"
    if("1" %in% levels(Dm.train.pred)) POS = "1"


    train_stats <- list(train=conf.mat.stats(conf.mat.pred=Dm.train.pred, conf.mat.truth=Dm.train$labelSVM, POS),
                        DmTrainPred = Dm.train.pred, train_model=Dm.train_model)


    if(doParalellSVM){
      Wm <- drop(t(Dm.train_model[[1]]$coefs) %*% Dm.train_model[[1]]$SV)
      bm <- drop(Dm.train_model[[1]]$rho)
    } else {
      Wm <- drop(t(Dm.train_model$coefs) %*% Dm.train_model$SV)
      bm <- drop(Dm.train_model$rho)

    }

    list(hyperplane = as.vector(cbind(t(Wm), bm)), stats = train_stats)
  })

  baselineSVM <- do.call(rbind, lapply(dataset_results, function(res) res$hyperplane))
  results.all <- lapply(dataset_results, function(res) res$stats)
  colnames(baselineSVM) <- c(paste(rep("Wm", nDims), 1:nDims, sep=""), "b.int")


  rownames(baselineSVM) <- names(TrainXls)
  names(results.all) <- names(TrainXls)

  structure(list(baselineSVM = baselineSVM,
                 results.all = results.all,
                 M = M,
                 nDims = nDims),
            class = "rtl_alg1_result")

}





