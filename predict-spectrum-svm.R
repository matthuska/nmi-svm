#' Predict regions based on a trained Spectrum Kernel SVM
#'
#' description
#'
#' @param fit A model trained using train_spectrum_svm
#' @param x A character vector of DNA sequences
#' @param window The window width for binning (in bp)
#' @param param A list of SVM parameters (C, order)
#' @return A list with scores and probabilites for each window
#' @export
#' @examples
#' param <- list(C=1.0, order=4, alphabet="SNP")
#' fit <- train_spectrum_svm(x, y, 750, param)
#' scores <- predict_spectrum_svm(fit, newx, 750, param)
predict_spectrum_svm = function(fit, x, window, param) {
  require(shogun)

  x <- toupper(x)
  if (param$alphabet == "SNP")
    x <- chartr("N", "0", x)

  fit$parallel$set_num_threads(getOption("cores"))

  gap <- 0L
  reverse <- FALSE
  order <- param$k

  charfeat <- StringCharFeatures(param$alphabet)
  charfeat$set_features(x)
  charfeat$io$set_loglevel("MSG_WARN")
  charfeat$obtain_by_sliding_window(window, window)

  alpha <- charfeat$get_alphabet()

  if (alpha$get_num_bits() * order <= 16) {
    feats <- StringWordFeatures(alpha)
    preproc <- SortWordString()
  } else {
    feats <- StringUlongFeatures(alpha)
    preproc <- SortUlongString()
  }

  feats$obtain_from_char(charfeat, order - 1, order, gap, reverse)
  preproc$init(feats)
  feats$add_preprocessor(preproc)
  feats$apply_preprocessor()

  labels_out <- fit$apply(feats)
  scores <- labels_out$get_values()
  labels_out$scores_to_probabilities()
  prob <- labels_out$get_values()

  return(list(scores = scores, prob = prob))
}
