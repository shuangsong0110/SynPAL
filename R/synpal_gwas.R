#' Genome-wide association analysis with SynPAL
#'
#' @param bfile Prefix of the PLINK `.bed`, `.bim` and `.fam` files.
#' @param pheno File with the columns `IID`, `pheno` (the observed phenotype; `NA` where it was not
#'   measured) and `s_pheno` (the synthetic phenotype of every sample).
#' @param covar File with the column `IID` followed by the covariates.
#' @param model `"linear"` for a continuous phenotype, `"logistic"` for a binary phenotype coded 0/1
#'   (`s_pheno` is a probability between 0 and 1) or `"poisson"` for a count phenotype (`s_pheno` is
#'   non-negative).
#' @param screen_p Variants whose P-value in the labeled samples is at most `screen_p` are tested with
#'   SynPAL; `1` tests every variant.
#' @param ncores Number of cores.
#' @param out File to which the results are written.
#' @return A data frame with one row per variant: `CHR`, `SNP`, `BP`, `A1` (the effect allele), `BETA`
#'   (effect per copy of `A1`), `SE`, `P` and `P_labeled` (the P-value in the labeled samples).
#' @examples
#' d <- system.file("extdata", package = "SynPAL")
#' res <- synpal_gwas(file.path(d, "demo"), file.path(d, "demo_pheno_continuous.txt"),
#'                    file.path(d, "demo_covar.txt"), model = "linear")
#' head(res[order(res$P), ])
#' @export
synpal_gwas <- function(bfile, pheno, covar, model = c("linear", "logistic", "poisson"),
                        screen_p = 5e-4, ncores = 1, out = NULL) {
  model <- match.arg(model)
  fun <- engine(model)

  fam <- data.table::fread(paste0(bfile, ".fam"), header = FALSE, data.table = FALSE)
  bim <- data.table::fread(paste0(bfile, ".bim"), header = FALSE, data.table = FALSE,
                           colClasses = list(character = c(1, 2, 5, 6)))
  bed <- BEDMatrix::BEDMatrix(paste0(bfile, ".bed"), n = nrow(fam), p = nrow(bim))
  pheno <- data.table::fread(pheno, data.table = FALSE)
  pheno <- pheno[!is.na(pheno$s_pheno), ]
  covar <- data.table::fread(covar, data.table = FALSE)

  id <- intersect(intersect(as.character(fam$V2), as.character(pheno$IID)), as.character(covar$IID))
  rows <- match(id, as.character(fam$V2))
  y <- pheno$pheno[match(id, as.character(pheno$IID))]
  s <- pheno$s_pheno[match(id, as.character(pheno$IID))]
  x <- scale(as.matrix(covar[match(id, as.character(covar$IID)),
                             setdiff(names(covar), c("FID", "IID")), drop = FALSE]))
  r <- ifelse(is.na(y), 0, 1)
  w <- fun$weight(r, y, s, x)
  p_labeled <- labeled_test(y, x, r, model)

  m <- nrow(bim)
  res <- matrix(NA_real_, m, 4, dimnames = list(NULL, c("BETA", "SE", "P", "P_labeled")))
  for (idx in split(seq_len(m), ceiling(seq_len(m) / max(1, floor(5e7 / length(id)))))) {
    G <- bed[rows, idx, drop = FALSE]
    na <- which(is.na(G), arr.ind = TRUE)
    G[na] <- colMeans(G, na.rm = TRUE)[na[, 2]]

    res[idx, 4] <- p_labeled(G[r == 1, , drop = FALSE])
    run <- which(res[idx, 4] <= screen_p)
    if (length(run) == 0) next
    fit <- parallel::mclapply(run, function(j) {
      g <- scale(G[, j])
      tryCatch(c(fun$test(y, s, x, r, g, w), attr(g, "scaled:scale")),
               error = function(e) rep(NA_real_, 3))
    }, mc.cores = ncores, mc.preschedule = FALSE)
    fit <- matrix(unlist(fit), ncol = 3, byrow = TRUE)
    res[idx[run], 1:3] <- cbind(fit[, 1] / fit[, 3], fit[, 2] / fit[, 3],
                                2 * stats::pnorm(-abs(fit[, 1] / fit[, 2])))
  }

  res <- data.frame(CHR = bim$V1, SNP = bim$V2, BP = bim$V4, A1 = bim$V5, res)
  if (!is.null(out)) data.table::fwrite(res, out, sep = "\t", na = "NA")
  res
}

# P-value of each variant in the labeled samples: score test adjusted for the covariates.
labeled_test <- function(y, x, r, model) {
  lab <- r == 1
  X <- cbind(1, x[lab, , drop = FALSE])
  y <- y[lab]

  if (model == "linear") {
    fit <- stats::lm.fit(X, y)
    mu <- y - fit$residuals
    v <- rep(sum(fit$residuals^2) / (length(y) - ncol(X)), length(y))
  } else {
    family <- if (model == "logistic") stats::binomial() else stats::poisson()
    mu <- stats::glm.fit(X, y, family = family)$fitted.values
    v <- family$variance(mu)
    if (model == "poisson") v <- v * sum((y - mu)^2 / mu) / (length(y) - ncol(X))
  }
  XtVX_inv <- solve(crossprod(X, v * X))

  function(G) {
    U <- c(crossprod(G, y - mu))
    A <- crossprod(X, v * G)
    V <- colSums(G^2 * v) - colSums(A * (XtVX_inv %*% A))
    stats::pchisq(U^2 / V, 1, lower.tail = FALSE)
  }
}
