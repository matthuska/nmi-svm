#' Fit a spectrum kernel SVM using Shogun
#'
#' description
#'
#' @param x A character vector of DNA sequences
#' @param y A factor with two levels, "NMI" and "Background"
#' @param param SVM parameters (C, order)
#' @return A fit model to be used with predict_spectrum_svm
#' @export
#' @examples
#' param <- list(C=1.0, order=4, alphabet="SNP")
#' fit <- train_spectrum_svm(x, y, 750, param)
fit_spectrum_svm <- function(x, y, param) {
  if (class(y) != "factor" || length(levels(y)) != 2) {
    cat("y must be a factor with two levels, aborting\n")
  }

  if (class(x[,1]) != "character") {
    cat("x must be a character vector, aborting\n")
  }

  lev <- levels(y)
  y <- as.double(ifelse(y == lev[[1]], 1L, -1L))
  labels <- BinaryLabels()
  labels$set_labels(y)

  if (param$alphabet == "SNP")
    x[,1] <- chartr("N", "0", x[,1])

  charfeat <- StringCharFeatures(param$alphabet)
  charfeat$set_features(x)
  charfeat$io$set_loglevel("MSG_WARN")

  alpha <- charfeat$get_alphabet()
  if (alpha$get_num_bits() * param$k <= 16) {
    # Can use "Word"
    feats <- StringWordFeatures(alpha)
    preproc <- SortWordString()
  }  else {
    # Need to use Ulong
    feats <- StringUlongFeatures(alpha)
    preproc <- SortUlongString()
  }

  gap <- 0L
  reverse <- FALSE
  feats$obtain_from_char(charfeat, param$k - 1, param$k, gap, reverse)
  preproc$init(feats)
  feats$add_preprocessor(preproc)
  feats$apply_preprocessor()

  if (alpha$get_num_bits() * param$k <= 16) {
    kernel <- CommWordStringKernel(feats, feats, param$usesign)
  } else {
    kernel <- CommUlongStringKernel(feats, feats, param$usesign)
  }

  svm <- SVMLight(param$C, kernel, labels)

  svm$parallel$set_num_threads(getOption("cores"))
  svm$set_epsilon(param$eps)
  retval <- svm$train()

  stopifnot(retval)

  return(svm)
}
