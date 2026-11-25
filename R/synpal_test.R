#' Test one variant with SynPAL
#'
#' @param y Phenotype; `NA` for the samples in which it was not measured.
#' @param s Synthetic phenotype of every sample.
#' @param g Genotype of the variant (allele count or dosage).
#' @param x Covariate matrix with one row per sample and no intercept column.
#' @param model `"linear"` for a continuous phenotype, `"logistic"` for a binary phenotype coded 0/1
#'   (`s` is a probability between 0 and 1) or `"poisson"` for a count phenotype (`s` is non-negative).
#' @param w Calibration weight; computed from `y`, `s` and `x` when `NULL`.
#' @return A named vector: `BETA` (effect per allele), `SE` and `P`.
#' @export
synpal_test <- function(y, s, g, x, model = c("linear", "logistic", "poisson"), w = NULL) {
  fun <- engine(match.arg(model))
  r <- ifelse(is.na(y), 0, 1)
  x <- scale(as.matrix(x))
  if (is.null(w)) w <- fun$weight(r, y, s, x)

  g[is.na(g)] <- mean(g, na.rm = TRUE)
  g <- scale(g)
  res <- fun$test(y, s, x, r, g, w)
  sd_g <- attr(g, "scaled:scale")

  c(BETA = res[[1]] / sd_g, SE = res[[2]] / sd_g, P = 2 * stats::pnorm(-abs(res[[1]] / res[[2]])))
}

engine <- function(model) {
  switch(model,
         linear = list(weight = weight_linear, test = test_linear),
         logistic = list(weight = weight_logistic, test = test_logistic),
         poisson = list(weight = weight_poisson, test = test_poisson))
}
