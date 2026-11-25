# Continuous phenotypes

weight_linear <- function(r, y, s, x) {
  lab <- r == 1
  y[!lab] <- 0

  Z0 <- cbind(1, x)
  Z_aug0 <- cbind(1, x, s, splines2::bSpline(s))

  beta_lab <- stats::lm(y[lab] ~ -1 + Z0[lab, ])$coefficients
  beta_aug <- stats::lm(y[lab] ~ -1 + Z_aug0[lab, ])$coefficients
  beta_aug[is.na(beta_aug)] <- 0

  p_x <- c(Z0 %*% beta_lab)
  p_s <- c(Z_aug0 %*% beta_aug)

  sum((p_s[lab] - p_x[lab]) * (y[lab] - p_x[lab])) / sum((p_s[lab] - p_x[lab])^2)
}

test_linear <- function(y, s, x, r, g, w) {
  y[is.na(y)] <- s[is.na(y)]
  n <- length(y)
  pi_h <- sum(r) / n
  lab <- r == 1

  Z <- cbind(1, g, x)
  Z_aug <- cbind(Z, s, splines2::bSpline(s))

  beta_lab <- stats::lm(y[lab] ~ -1 + Z[lab, ])$coefficients
  beta_aug <- stats::lm(y[lab] ~ -1 + Z_aug[lab, ])$coefficients
  beta_aug[is.na(beta_aug)] <- 0

  p_x <- c(Z %*% beta_lab)
  p_s <- c(Z_aug %*% beta_aug)

  y_cal <- w * p_s + (1 - w) * p_x

  moment <- (pi_h - r) / pi_h * cbind(y_cal, y_cal^2, y_cal^3)

  if (sum(abs(colMeans(moment) / apply(moment, 2, stats::sd))) < 1 / n^(1 / 4)) {
    fit <- stats::lm(y_cal ~ -1 + Z)$coefficients
    p_x <- c(Z %*% fit)
    est <- fit[2]

    IF <- c(y_cal - p_x + r / pi_h * (y - y_cal))
    se <- sqrt(stats::var(IF) / n)
  } else {
    basis <- cbind(splines2::bSpline(s), Z_aug)
    kappa <- basis %*% MASS::ginv(t(basis[lab, ]) %*% basis[lab, ]) %*% (t(basis) %*% ((1 - r / pi_h) * Z))

    fit <- stats::lm(y_cal ~ -1 + Z)$coefficients
    p_x <- c(Z %*% fit)
    M <- MASS::ginv(t(Z) %*% Z / n)
    est <- fit[2] - (M %*% colMeans(c(r * (y_cal - y)) * kappa))[2]

    IF <- c(y_cal - p_x + r / pi_h * (y - y_cal)) * Z - c(r * (y_cal - y)) * kappa
    se <- sqrt((M %*% stats::cov(IF) %*% M)[2, 2] / n)
  }

  c(est, se)
}
