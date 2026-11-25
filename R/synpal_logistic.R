# Binary phenotypes

weight_logistic <- function(r, y, s, x) {
  lab <- r == 1
  y[!lab] <- 0

  s <- 5e-2 + (1 - 1e-1) * s

  Z0 <- cbind(1, x)
  Z_aug0 <- cbind(1, x, log(s), log(1 - s), splines2::bSpline(s))

  beta_lab <- stats::glm(y[lab] ~ -1 + Z0[lab, ], family = stats::binomial)$coefficients
  beta_aug <- stats::glm(y[lab] ~ -1 + Z_aug0[lab, ], family = stats::binomial)$coefficients
  beta_aug[is.na(beta_aug)] <- 0

  p_x <- c(stats::plogis(Z0 %*% beta_lab))
  p_s <- c(stats::plogis(Z_aug0 %*% beta_aug))

  sum((p_s[lab] - p_x[lab]) * (y[lab] - p_x[lab])) / sum((p_s[lab] - p_x[lab])^2)
}

test_logistic <- function(y, s, x, r, g, w) {
  y[is.na(y)] <- s[is.na(y)]
  n <- length(y)
  pi_h <- sum(r) / n
  lab <- r == 1

  Z <- cbind(1, g, x)

  s <- 5e-2 + (1 - 1e-1) * s
  Z_aug <- cbind(Z, log(s), log(1 - s), splines2::bSpline(s))

  beta_lab <- fastglm::fastglm(x = Z[lab, ], y = y[lab], family = stats::binomial())$coefficients
  beta_aug <- fastglm::fastglm(x = Z_aug[lab, ], y = y[lab], family = stats::binomial())$coefficients
  beta_aug[is.na(beta_aug)] <- 0

  p_x <- c(stats::plogis(Z %*% beta_lab))
  p_s <- c(stats::plogis(Z_aug %*% beta_aug))

  y_cal <- w * p_s + (1 - w) * p_x
  y_cal <- pmin(pmax(y_cal, 0), 1)

  moment <- (pi_h - r) / pi_h * cbind(y_cal, y_cal^2, y_cal^3)

  if (sum(abs(colMeans(moment) / apply(moment, 2, stats::sd))) < 1 / n^(1 / 4)) {
    fit <- fastglm::fastglm(x = Z, y = y_cal, family = stats::quasibinomial())$coefficients
    p_x <- c(stats::plogis(Z %*% fit))
    p_y <- mean(y[lab])
    est <- fit[2]

    IF <- c(y_cal - p_x + r / pi_h * (y - y_cal))
    se <- sqrt(stats::var(IF) / (p_y - p_y^2)^2 / n)
  } else {
    basis <- cbind(splines2::bSpline(s), p_s - p_x, c(p_s * (1 - p_s)) * Z_aug, c(p_x * (1 - p_x)) * Z)
    kappa <- basis %*% MASS::ginv(t(basis[lab, ]) %*% basis[lab, ]) %*% (t(basis) %*% ((1 - r / pi_h) * Z))

    fit <- stats::glm(y_cal ~ -1 + Z, family = stats::quasibinomial)$coefficients
    p_x <- c(stats::plogis(Z %*% fit))
    M <- MASS::ginv(t(Z) %*% ((p_x - p_x^2) * Z) / n)
    est <- fit[2] - (M %*% colMeans(c(r * (y_cal - y)) * kappa))[2]

    IF <- c(y_cal - p_x + r / pi_h * (y - y_cal)) * Z - c(r * (y_cal - y)) * kappa
    se <- sqrt((M %*% stats::cov(IF) %*% M)[2, 2] / n)
  }

  c(est, se)
}
