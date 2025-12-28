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

#' Normalize and validate RNG seeds
#'
#' @param seed Optional numeric seed.
#'
#' @return A single integer seed or `NULL` when not supplied.
#' @export
normalize_seed <- function(seed){
  if (is.null(seed)) return(NULL)

  if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed)) {
    stop("seed must be a single finite numeric value.", call. = FALSE)
  }

  as.integer(seed)
}

#' Evaluate an expression while restoring the incoming RNG state
#'
#' @param seed Optional numeric seed used for reproducibility.
#' @param expr Expression to evaluate.
#'
#' @return The result of `expr`.
#' @export
scoped_seed <- function(seed, expr){
  seed <- normalize_seed(seed)

  if (is.null(seed)) {
    return(force(expr))
  }

  seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (seed_existed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    if (seed_existed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  }, add = FALSE)

  set.seed(seed)
  force(expr)
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


#' Validate cell-level metadata alongside an assay matrix
#'
#' @param feature_matrix Matrix or data.frame with rows representing cells and
#'   columns representing features.
#' @param cell_metadata Data.frame containing per-cell metadata.
#' @param cell_id_column Name of the column storing unique cell identifiers.
#' @param class_column Optional name of the column storing class labels.
#' @param required_metadata Optional character vector of columns that must be
#'   present in `cell_metadata`. Missing columns are created with `NA`.
#'
#' @return A list with harmonized `features` (data.frame) and `cell_metadata`.
#' @export
validate_assay_inputs <- function(feature_matrix,
                                  cell_metadata,
                                  cell_id_column = "cell_id",
                                  class_column = NULL,
                                  required_metadata = NULL){
  if (is.null(cell_metadata)) {
    stop("cell_metadata must be provided for validation.", call. = FALSE)
  }

  features <- as.data.frame(feature_matrix)
  metadata <- as.data.frame(cell_metadata)

  if (!nrow(features)) {
    stop("No cells provided for validation.", call. = FALSE)
  }

  if (!(cell_id_column %in% colnames(metadata))) {
    metadata[[cell_id_column]] <- paste0("cell_", seq_len(nrow(metadata)))
  }

  if (nrow(features) != nrow(metadata)) {
    stop("feature_matrix and cell_metadata must have the same number of rows.", call. = FALSE)
  }

  if (anyDuplicated(metadata[[cell_id_column]]) > 0) {
    stop("cell_metadata must contain unique cell identifiers.", call. = FALSE)
  }

  if (!is.null(class_column) && !(class_column %in% colnames(metadata))) {
    metadata[[class_column]] <- NA
  }

  if (!is.null(required_metadata)) {
    missing_cols <- setdiff(required_metadata, colnames(metadata))
    for (missing_col in missing_cols) {
      metadata[[missing_col]] <- NA
    }
  }

  rownames(features) <- metadata[[cell_id_column]]

  bad_cells <- !is.finite(rowSums(features))
  if (any(bad_cells)) {
    warning("Removing cells with non-finite values after validation.", call. = FALSE)
    features <- features[!bad_cells, , drop = FALSE]
    metadata <- metadata[!bad_cells, , drop = FALSE]
  }

  list(
    features = features,
    cell_metadata = metadata
  )
}


read_tabular_if_path <- function(x){
  if (is.character(x) && length(x) == 1 && file.exists(x)) {
    return(data.table::fread(x))
  }
  as.data.frame(x)
}


coerce_numeric_columns <- function(df, columns){
  df[, columns] <- lapply(columns, function(col){
    values <- df[[col]]
    if (is.numeric(values)) {
      return(values)
    }
    suppressWarnings(as.numeric(values))
  })
  df
}


#' Load CyTOF-style assay tables with harmonized metadata
#'
#' @param exprs A data.frame/matrix of expression values or a path to a delimited file.
#' @param cell_metadata Optional data.frame of per-cell metadata. If `NULL`,
#'   non-marker columns in `exprs` are treated as metadata.
#' @param marker_columns Optional character vector specifying which columns are
#'   assay markers. Defaults to all numeric columns after removing metadata.
#' @param arcsinh_cofactor Cofactor used for arcsinh transformation. Set to
#'   `NULL` to skip transformation.
#' @param class_column Name of the metadata column containing class labels.
#' @param cell_id_column Name of the metadata column containing cell identifiers.
#' @param qc_max_na Maximum allowed fraction of missing values per cell before
#'   filtering.
#'
#' @return A list containing transformed `exprs`, `cell_metadata`,
#'   `feature_metadata`, and an `assay_type` flag.
#' @export
load_cytof_assay <- function(exprs,
                             cell_metadata = NULL,
                             marker_columns = NULL,
                             arcsinh_cofactor = 5,
                             class_column = "label",
                             cell_id_column = "cell_id",
                             qc_max_na = 0.05){
  exprs_df <- read_tabular_if_path(exprs)

  if (is.null(marker_columns)) {
    marker_columns <- colnames(exprs_df)[vapply(exprs_df, is.numeric, logical(1))]
  }

  marker_columns <- setdiff(marker_columns, c(class_column, cell_id_column))
  marker_columns <- marker_columns[marker_columns %in% colnames(exprs_df)]

  if (!length(marker_columns)) {
    stop("No marker columns detected for CyTOF assay.", call. = FALSE)
  }

  exprs_df <- coerce_numeric_columns(exprs_df, marker_columns)

  markers <- exprs_df[, marker_columns, drop = FALSE]

  if (!is.null(qc_max_na) && qc_max_na < 1) {
    keep_cells <- rowMeans(is.na(markers)) <= qc_max_na
    if (any(!keep_cells)) {
      warning("Removing cells exceeding missing value threshold.", call. = FALSE)
      markers <- markers[keep_cells, , drop = FALSE]
      if (!is.null(cell_metadata)) {
        cell_metadata <- cell_metadata[keep_cells, , drop = FALSE]
      } else {
        exprs_df <- exprs_df[keep_cells, , drop = FALSE]
      }
    }
  }

  if (!is.null(arcsinh_cofactor)) {
    markers <- as.data.frame(asinh(as.matrix(markers) / arcsinh_cofactor))
  }

  if (is.null(cell_metadata)) {
    non_marker_cols <- setdiff(colnames(exprs_df), marker_columns)
    cell_metadata <- exprs_df[, non_marker_cols, drop = FALSE]
  }

  validated <- validate_assay_inputs(
    feature_matrix = markers,
    cell_metadata = cell_metadata,
    cell_id_column = cell_id_column,
    class_column = class_column,
    required_metadata = c(class_column)
  )

  feature_metadata <- data.frame(
    marker = marker_columns,
    channel = marker_columns,
    measurement_type = "protein",
    stringsAsFactors = FALSE
  )

  list(
    exprs = validated$features,
    cell_metadata = validated$cell_metadata,
    feature_metadata = feature_metadata,
    assay_type = "cytof",
    transform = list(method = "arcsinh", cofactor = arcsinh_cofactor)
  )
}


#' Load scRNA-seq count matrices with minimal QC
#'
#' @param counts Matrix/data.frame of raw counts (genes as rows, cells as
#'   columns) or a path to a delimited file.
#' @param cell_metadata Optional data.frame with per-cell annotations. If `NULL`,
#'   a metadata frame with only `cell_id_column` is created.
#' @param feature_metadata Optional data.frame describing genes.
#' @param normalization One of `"log1p"`, `"cpm"`, or `"none"` indicating the
#'   normalization applied to counts.
#' @param cell_id_column Name of the metadata column containing cell identifiers.
#' @param gene_id_column Column name (if present) that stores gene identifiers
#'   when reading from tabular files.
#' @param min_cells Minimum number of cells a gene must appear in to be retained.
#' @param min_features Minimum number of detected genes required for each cell.
#' @param class_column Optional metadata column storing class labels.
#'
#' @return A list containing `counts` (genes x cells), `normalized_counts`
#'   (cells x genes), `cell_metadata`, `feature_metadata`, and `assay_type`.
#' @export
load_scrnaseq_counts <- function(counts,
                                 cell_metadata = NULL,
                                 feature_metadata = NULL,
                                 normalization = c("log1p", "cpm", "none"),
                                 cell_id_column = "cell_id",
                                 gene_id_column = "gene_id",
                                 min_cells = 5,
                                 min_features = 200,
                                 class_column = "label"){
  normalization <- match.arg(normalization)
  counts_df <- read_tabular_if_path(counts)

  if (!is.null(gene_id_column) && gene_id_column %in% colnames(counts_df)) {
    gene_ids <- counts_df[[gene_id_column]]
    counts_df[[gene_id_column]] <- NULL
    rownames(counts_df) <- gene_ids
  }

  counts_mat <- as.matrix(counts_df)

  if (is.null(rownames(counts_mat))) {
    rownames(counts_mat) <- paste0("gene_", seq_len(nrow(counts_mat)))
  }

  if (is.null(colnames(counts_mat))) {
    colnames(counts_mat) <- paste0("cell_", seq_len(ncol(counts_mat)))
  }

  genes_keep <- rowSums(counts_mat > 0) >= min_cells
  cells_keep <- colSums(counts_mat > 0) >= min_features

  counts_mat <- counts_mat[genes_keep, cells_keep, drop = FALSE]

  norm_counts <- switch(
    normalization,
    log1p = log1p(counts_mat),
    cpm = {
      scaling <- colSums(counts_mat)
      scaling[scaling == 0] <- 1
      log1p(t(t(counts_mat) / scaling * 1e6))
    },
    none = counts_mat
  )

  if (is.null(cell_metadata)) {
    cell_metadata <- data.frame(
      cell_id = colnames(counts_mat),
      stringsAsFactors = FALSE
    )
  }

  validated <- validate_assay_inputs(
    feature_matrix = t(norm_counts),
    cell_metadata = cell_metadata,
    cell_id_column = cell_id_column,
    class_column = class_column,
    required_metadata = c(class_column)
  )

  if (is.null(feature_metadata)) {
    feature_metadata <- data.frame(
      gene_id = rownames(counts_mat),
      detected_cells = rowSums(counts_mat > 0),
      stringsAsFactors = FALSE
    )
  }

  list(
    counts = counts_mat,
    normalized_counts = validated$features,
    cell_metadata = validated$cell_metadata,
    feature_metadata = feature_metadata,
    assay_type = "scrna",
    normalization = normalization
  )
}
