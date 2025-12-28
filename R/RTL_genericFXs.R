#' Build a qualitative color palette for RTL visuals
#'
#' @description
#' Generates qualitative and sequential color palettes used across RTL plotting
#' routines.
#'
#' @return A list containing `col_vector` (qualitative palette) and
#'   `scaleyellowred` (sequential palette).
#' @export
ColorTheme <- function(){
  #scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)
  scaleyellowred <- colorRampPalette(c("dodgerblue", "lightyellow", "red"), space = "rgb")(30)

  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- col_vector[-4]
  return(list(col_vector=col_vector, scaleyellowred=scaleyellowred))
}

conf.mat.stats <- function(conf.mat.pred, conf.mat.truth, POS="1"){

  #conf.mat.pred   =  factor(sign(y_hat[,x_m]), levels=c(-1,1))
  #conf.mat.truth  =  TASK$TrueClass

  temp.confMatdat <- confusionMatrix(data=conf.mat.pred, reference=conf.mat.truth, positive=POS, mode="everything")
  return(list(overall= temp.confMatdat$overall, table=temp.confMatdat$table, byClass=temp.confMatdat$byClass))

}


#' Compute the statistical mode of a vector
#'
#' @param x Vector of values.
#'
#' @return The most frequent value in `x`.
#' @export
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Convert logical flags to indicator values
#'
#' @param x Logical value to convert.
#'
#' @return `1` when `x` is `TRUE`, `0` when `FALSE`; otherwise prints a warning.
#' @export
IndicatorFX <- function(x){
  if(x==T) return(1)
  if(x==F) return(0)
  if(!(x %in% c(T, F))) {
    print("Not T/F")}
}


is.odd <- function(x){
  x%%2 == 1
}
is.even <- function(x){
  x%%2 == 0
}

#' Apply a sequential or parallel map with consistent ordering
#'
#' @param indices Index vector to iterate over.
#' @param FUN Function to apply.
#' @param use_parallel Logical; when `TRUE` attempt to use `parallel::mclapply`.
#' @param parallel_cores Optional integer to override the detected core count.
#'
#' @return A list matching the length and order of `indices`.
map_with_backend <- function(indices, FUN, use_parallel = FALSE, parallel_cores = NULL){
  if (!length(indices)) {
    return(list())
  }
  if (!use_parallel || length(indices) == 1) {
    return(lapply(indices, FUN))
  }

  cores <- if (is.null(parallel_cores)) {
    max(1L, parallel::detectCores() - 1L)
  } else {
    max(1L, parallel_cores)
  }

  if (.Platform$OS.type == "windows") {
    warning("Parallel backend is not available on Windows; falling back to sequential execution.", call. = FALSE)
    return(lapply(indices, FUN))
  }

  parallel::mclapply(indices, FUN, mc.cores = cores, mc.preschedule = TRUE)
}

#' Coerce feature inputs while minimizing copies for wide data
#'
#' @param x Input data frame or matrix.
#' @param cols Optional column indices to keep.
#' @param wide_data_threshold Number of columns above which `data.table` is used.
#'
#' @return A data.frame or data.table with selected columns.
coerce_feature_frame <- function(x, cols, wide_data_threshold = 200){
  if(!is.na(cols[1]) && !(length(cols) == 1 && identical(cols, ""))){
    x <- x[, cols, drop = FALSE]
  }
  if (ncol(x) >= wide_data_threshold) {
    return(data.table::as.data.table(x))
  }
  as.data.frame(x)
}

#' Split a labelled data frame into cross-validated train/test lists
#'
#' @description
#' Generates stratified K-fold splits for labelled feature data and produces
#' paired train/test lists compatible with RTL algorithms.
#'
#' @param X Data frame or matrix of features.
#' @param Y Vector or factor of class labels aligned to `X`.
#'
#' @return A list with `TrainXY.ls` and `TestXY.ls` elements containing per-fold
#'   feature/label pairs.
#' @export
DF2TrainTestls <- function(X, Y){
  #X = BCC$PBMC
  #Y = BCC$PBMC.Class
  DF2KfoldBYClass <- function(tSDF, tSYV, Ksplits=10){
    #use factor to split first to make sure both cases are in each set

    # tSDF     = BCC$PBMC
    # tSYV     = BCC$PBMC.Class
    # Ksplits  = 10

    if(!is.factor(tSYV)) tSYV <- factor(tSYV)

    ClassA.DF <- subset(cbind(tSDF, tSYV), tSYV == levels(tSYV)[1])
    ClassB.DF <- subset(cbind(tSDF, tSYV), tSYV == levels(tSYV)[2])
    err = F
    if(nrow(ClassA.DF) < Ksplits) {
      print("cant split, too few A cases")
      err = T
    }
    if(nrow(ClassB.DF) < Ksplits) {
      print("cant split, too few B cases")
      err = T
    }
    if(!err){

      tempSampAIDs <- 1:nrow(ClassA.DF)
      tempSampBIDs <- 1:nrow(ClassB.DF)
      options(warn=-1)
      tempSampAIDs.sp <- split(sample(tempSampAIDs), 1:Ksplits)
      tempSampBIDs.sp <- split(sample(tempSampBIDs), 1:Ksplits)
      options(warn=0)
      Y.split.ls <- lapply(1:length(tempSampBIDs.sp), function(n){
        #n = 1
        tempYV <- c(ClassA.DF[tempSampAIDs.sp[[n]], "tSYV"],
                    ClassB.DF[tempSampBIDs.sp[[n]],"tSYV"])
        factor(tempYV, levels = c(1,2), labels = levels(tSYV))
      })

      X.split.ls <- lapply(1:length(tempSampBIDs.sp), function(n){
        #n = 1

        rbind(ClassA.DF[tempSampAIDs.sp[[n]],
                        setdiff(colnames(ClassA.DF), "tSYV")],
              ClassB.DF[tempSampBIDs.sp[[n]],
                        setdiff(colnames(ClassA.DF), "tSYV")])[,]
      })

      return(lapply(1:length(X.split.ls), function(leng){
        #leng = 1
        list(x=X.split.ls[[leng]],y=Y.split.ls[[leng]], KselectRows=c(tempSampAIDs.sp[[leng]], tempSampBIDs.sp[[leng]]))

      }))
    }

  }

  TestXY.ls <- DF2KfoldBYClass(X, Y)
  names(TestXY.ls) <- paste("TestSubSec", 1:length(TestXY.ls), sep="")

  Nminus1Train <- function(testXYlist){
    #testXYlist = TestXY.ls

    lapply(names(testXYlist), function(Nm){
      #Nm=names(testXYlist)[1]
      list(x = as.data.frame(rbindlist(lapply(names(testXYlist)[which(!(names(testXYlist) %in% Nm))], function(xNm){
        #xNm=names(testXYlist)[2]
        as.data.frame(testXYlist[[xNm]]$x)
      }))),
      y = unlist(lapply(names(testXYlist)[which(!(names(testXYlist) %in% Nm))], function(xNm){
        #xNm=names(testXYlist)[2]
        (testXYlist[[xNm]]$y)
      })))

    })

  }

  TrainXY.ls <- Nminus1Train(TestXY.ls)
  names(TrainXY.ls) <- paste("TrainSubSec", 1:length(TrainXY.ls), sep="")


  return(list(TrainXY.ls = TrainXY.ls, TestXY.ls = TestXY.ls))
}
