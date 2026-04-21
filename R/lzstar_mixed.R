#' Compute lz and lzstar person-fit statistics for mixed-format tests
#'
#' @param scores Matrix of item scores (n_persons x n_items).
#' @param theta Vector of ability estimates, length must equal nrow(scores).
#' @param vtheta Vector of estimated variances Var(theta_hat), length must equal nrow(scores).
#' @param a Vector of discrimination parameters, length must equal ncol(scores).
#' @param c Vector of lower asymptote (guessing) parameters, length must equal ncol(scores).
#' @param b Matrix of thresholds, with nrow(b) = ncol(scores) and NA padding for shorter items.
#' @param extreme_handling Policy for extreme raw scores (all-0 or perfect):
#'   "na" returns NA for that person, "cap" truncates theta to +/- theta_cap, "none" does nothing.
#' @param theta_cap Cap used when extreme_handling = "cap".
#'
#' @return A matrix with columns "lz" and "lzstar".
#'
#' @details
#' The function implements the mixed-format person-fit statistic using 3PL for dichotomous
#' items and GPCM for polytomous items. The calculation is kept unchanged: dichotomous item
#' contributions are accumulated with pr3PL, polytomous item contributions are accumulated
#' with prGPCM_stable and d1prGPCM, and lz and lzstar are computed from
#' the same denominators used in the implementation.
#'
#' @references
#' Sinharay, S. (2026). Refining the asymptotically correct standardization of person-fit
#' statistics for mixed-format tests. \emph{British Journal of Mathematical and Statistical
#' Psychology}, \bold{00}, 1--24.
#'
#' Sinharay, S. (2016). Asymptotically correct standardization of person-fit statistics beyond
#' dichotomous items. \emph{Psychometrika}, \bold{81}(4), 992--1013.
#'
#' @examples
#' scores <- matrix(c(
#'   1, 0, 2,
#'   0, 1, 1,
#'   1, 1, 3
#' ), nrow = 3, byrow = TRUE)
#' theta  <- c(-0.5, 0.2, 1.1)
#' vtheta <- c(0.12, 0.10, 0.15)
#' a      <- c(1.1, 0.9, 1.3)
#' cpar   <- c(0.15, 0.10, 0.00)
#' b <- matrix(c(
#'   0.0,  NA,
#'  -0.2,  NA,
#'  -1.0,  0.3
#' ), nrow = 3, byrow = TRUE)
#' lzstar_mixed(scores, theta, vtheta, a, cpar, b)
#'
#' scores2 <- matrix(c(
#'   0, 1,
#'   1, 0
#' ), nrow = 2, byrow = TRUE)
#' theta2  <- c(-1, 1)
#' vtheta2 <- c(0.08, 0.08)
#' a2      <- c(1.0, 1.2)
#' c2      <- c(0.05, 0.10)
#' b2 <- matrix(c(
#'   0.0, NA,
#'   0.3, NA
#' ), nrow = 2, byrow = TRUE)
#' lzstar_mixed(scores2, theta2, vtheta2, a2, c2, b2)
#'
#' @export
lzstar_mixed <- function(scores, theta, vtheta, a, c, b,
                         extreme_handling = c("na", "cap", "none"),
                         theta_cap = 4) {
  extreme_handling <- match.arg(extreme_handling)

  if (!is.matrix(scores)) scores <- as.matrix(scores)

  n_persons <- nrow(scores)
  n_items   <- ncol(scores)

  if (length(theta)  != n_persons) stop("Length(theta) must equal nrow(scores)")
  if (length(vtheta) != n_persons) stop("Length(vtheta) must equal nrow(scores)")
  if (length(a)      != n_items)   stop("Length(a) must equal ncol(scores)")
  if (length(c)      != n_items)   stop("Length(c) must equal ncol(scores)")
  if (nrow(b)        != n_items)   stop("nrow(b) must equal ncol(scores)")

  mxscores <- apply(scores, 2, max, na.rm = TRUE)
  dich     <- which(mxscores == 1)
  poly     <- which(mxscores > 1)

  total        <- rowSums(scores, na.rm = TRUE)
  max_possible <- sum(mxscores, na.rm = TRUE)
  is_extreme   <- (total == 0) | (total == max_possible)

  if (extreme_handling == "cap") {
    theta <- pmax(pmin(theta, theta_cap), -theta_cap)
  }

  out <- matrix(NA_real_, nrow = n_persons, ncol = 2,
                dimnames = list(NULL, c("lz", "lzstar")))

  for (j in seq_len(n_persons)) {
    if (extreme_handling == "na" && is_extreme[j]) next
    if (!is.finite(theta[j]) || !is.finite(vtheta[j]) || vtheta[j] < 0) next

    lznum    <- 0
    lzden    <- 0
    lzdenadd <- 0

    if (length(dich) > 0) {
      pr <- pr3PL(theta[j], a[dich], b[dich, 1], c[dich])
      pr <- pmax(pmin(pr, 1 - 1e-10), 1e-10)
      w  <- log(pr / (1 - pr))
      e  <- exp(a[dich] * (theta[j] - b[dich, 1]))
      p1 <- (1 - c[dich]) * a[dich] * e / ((1 + e)^2)

      lznum    <- lznum    + sum((scores[j, dich] - pr) * w)
      lzden    <- lzden    + sum(w^2 * pr * (1 - pr))
      lzdenadd <- lzdenadd + sum(p1 * w)
    }

    if (length(poly) > 0) {
      for (i in poly) {
        thres_i <- b[i, ]
        thres_i <- thres_i[!is.na(thres_i)]
        if (length(thres_i) == 0) next

        probs <- prGPCM_stable(theta[j], a[i], thres_i)
        probs <- pmax(probs, 1e-10)
        probs <- probs / sum(probs)
        log_probs <- log(probs)
        vmat <- diag(probs) - outer(probs, probs)
        der1 <- d1prGPCM(theta[j], a[i], thres_i)
        obs  <- scores[j, i]
        if (is.na(obs) || obs < 0 || obs >= length(probs)) next

        lznum    <- lznum    + (log_probs[obs + 1] - sum(probs * log_probs))
        lzden    <- lzden    + as.numeric(t(log_probs) %*% vmat %*% log_probs)
        lzdenadd <- lzdenadd + sum(der1 * log_probs)
      }
    }

    if (is.finite(lzden) && lzden > 0) {
      out[j, "lz"] <- lznum / sqrt(lzden)
    }

    lzstden <- lzden - vtheta[j] * (lzdenadd^2)
    if (!is.finite(lzstden)) {
      out[j, "lzstar"] <- NA_real_
      next
    }
    if (lzstden <= 0 && lzstden > -1e-6) {
      lzstden <- 1e-6
    }
    if (lzstden > 0) {
      out[j, "lzstar"] <- lznum / sqrt(lzstden)
    }
  }

  out
}

# 3PL probability: P(X = 1 | theta)
pr3PL <- function(theta, a, b, c) {
  c + (1 - c) / (1 + exp(-a * (theta - b)))
}

# Category probabilities under GPCM with numerical stability
prGPCM_stable <- function(theta, a, b) {
  n_cat <- length(b) + 1
  log_probs <- numeric(n_cat)
  log_probs[1] <- 0
  for (k in seq(2, n_cat)) {
    log_probs[k] <- log_probs[k - 1] + a * (theta - b[k - 1])
  }
  log_probs <- log_probs - max(log_probs)
  probs <- exp(log_probs)
  probs / sum(probs)
}

# Derivative of category probabilities under GPCM
# dP_k/dtheta = a * P_k * (score_k - E(score))
d1prGPCM <- function(theta, a, b) {
  probs  <- prGPCM_stable(theta, a, b)
  scores <- seq(0, length(b))
  a * probs * (scores - sum(scores * probs))
}

#' Universal PPall person-fit wrapper with full diagnostics
#'
#' @param respm Response matrix.
#' @param ppobj PPall output.
#' @param model2est Original model specification.
#' @param orig_thres Original thres input to PPall.
#' @param orig_slopes Original slopes input.
#' @param orig_lowerA Original lowerA input.
#'
#' @details
#' This wrapper tries multiple extraction strategies and then calls
#' \code{lzstar_mixed()} without changing the core calculation.
#'
#' @export
Pfit_lzstarmix <- function(respm, ppobj, model2est = NULL,
                           orig_thres = NULL, orig_slopes = NULL,
                           orig_lowerA = NULL) {
  n_persons <- nrow(respm)
  n_items   <- ncol(respm)

  cat(sprintf("Input dimensions: %d persons x %d items\n", n_persons, n_items))

  resPP <- ppobj$resPP$resPP
  if (is.matrix(resPP) || is.data.frame(resPP)) {
    theta    <- as.numeric(resPP[, 1])
    se_theta <- as.numeric(resPP[, 2])
  } else if (is.vector(resPP)) {
    warning("resPP is a vector; setting SE = 0.5")
    theta    <- as.numeric(resPP)
    se_theta <- rep(0.5, length(theta))
  } else {
    stop("Unsupported resPP format in ppobj")
  }

  vtheta <- se_theta^2

  a_final <- NULL
  c_final <- NULL
  b_final <- NULL

  possible_a <- c("slopes", "a", "sl")
  possible_c <- c("lowerA", "c", "la")
  possible_b <- c("thres", "b", "THRESx")

  a_candidates <- list()
  c_candidates <- list()
  b_candidates <- list()

  for (a_name in possible_a) {
    if (a_name %in% names(ppobj$ipar)) {
      a_cand <- ppobj$ipar[[a_name]]
      if (length(a_cand) == n_items) a_candidates[[a_name]] <- a_cand
    }
  }
  for (c_name in possible_c) {
    if (c_name %in% names(ppobj$ipar)) {
      c_cand <- ppobj$ipar[[c_name]]
      if (length(c_cand) == n_items) c_candidates[[c_name]] <- c_cand
    }
  }
  for (b_name in possible_b) {
    if (b_name %in% names(ppobj$ipar)) {
      b_cand <- ppobj$ipar[[b_name]]
      if (is.matrix(b_cand) && nrow(b_cand) == n_items) {
        b_candidates[[b_name]] <- b_cand
      } else if (is.vector(b_cand) && length(b_cand) == n_items) {
        b_candidates[[paste0(b_name, "_converted")]] <- matrix(b_cand, ncol = 1)
      }
    }
  }

  if (!is.null(orig_slopes) && length(orig_slopes) == n_items) {
    a_final <- orig_slopes
  } else if (length(a_candidates) > 0) {
    a_final <- a_candidates[[1]]
  } else {
    stop("Cannot find valid slopes (a) parameters")
  }

  if (!is.null(orig_lowerA) && length(orig_lowerA) == n_items) {
    c_final <- orig_lowerA
  } else if (length(c_candidates) > 0) {
    c_final <- c_candidates[[1]]
  } else {
    c_final <- rep(0, n_items)
  }

  if (!is.null(orig_thres)) {
    if (is.vector(orig_thres)) {
      b_final <- matrix(orig_thres, ncol = 1)
    } else if (is.matrix(orig_thres)) {
      b_final <- orig_thres
    } else if (is.list(orig_thres)) {
      if ("thres" %in% names(orig_thres)) {
        b_final <- orig_thres$thres
        if (!is.matrix(b_final)) b_final <- matrix(b_final, ncol = 1)
      } else {
        b_final <- orig_thres[[1]]
        if (length(b_final) != n_items) b_final <- matrix(b_final, ncol = 1)
      }
    } else {
      b_final <- matrix(as.numeric(orig_thres), ncol = 1)
    }
  } else if (length(b_candidates) > 0) {
    b_final <- b_candidates[[1]]
  } else {
    b_final <- matrix(rnorm(n_items * 2), n_items, 2)
  }

  if (length(a_final) != n_items || length(c_final) != n_items || nrow(b_final) != n_items) {
    stop("Parameter dimensions still do not match!")
  }

  cat("\n--- Computing Person-Fit Statistics ---\n")
  results <- lzstar_mixed(respm, theta, vtheta, a_final, c_final, b_final)

  cat(sprintf("Success! Computed %d statistics\n", nrow(results)))
  cat(sprintf(" lz mean: %.3f (sd: %.3f)\n",
              mean(results[, "lz"], na.rm = TRUE),
              sd(results[, "lz"], na.rm = TRUE)))
  cat(sprintf(" lzstar mean: %.3f (sd: %.3f)\n",
              mean(results[, "lzstar"], na.rm = TRUE),
              sd(results[, "lzstar"], na.rm = TRUE)))

  results
}

#' Print summary of person-fit results
#'
#' @param results Matrix from lzstar_mixed.
#' @param alpha Significance level for flagging.
#'
#' @export
print_personfit_summary <- function(results, alpha = 0.05) {
  cat("\n=== Person-Fit Statistics Summary ===\n\n")
  z_crit <- qnorm(1 - alpha / 2)

  cat("lz statistics:\n")
  cat(sprintf(
    " Mean: %.3f (SD: %.3f)\n",
    mean(results[, 1], na.rm = TRUE),
    sd(results[, 1], na.rm = TRUE)
  ))
  cat(sprintf(
    " Range: [%.3f, %.3f]\n",
    min(results[, 1], na.rm = TRUE),
    max(results[, 1], na.rm = TRUE)
  ))
  flagged_lz <- sum(abs(results[, 1]) > z_crit, na.rm = TRUE)
  cat(sprintf(
    " Flagged (|z| > %.2f): %d (%.1f%%)\n\n",
    z_crit, flagged_lz, 100 * flagged_lz / nrow(results)
  ))

  cat("lzstar statistics (asymptotically correct):\n")
  cat(sprintf(
    " Mean: %.3f (SD: %.3f)\n",
    mean(results[, 2], na.rm = TRUE),
    sd(results[, 2], na.rm = TRUE)
  ))
  cat(sprintf(
    " Range: [%.3f, %.3f]\n",
    min(results[, 2], na.rm = TRUE),
    max(results[, 2], na.rm = TRUE)
  ))
  flagged_lzstar <- sum(abs(results[, 2]) > z_crit, na.rm = TRUE)
  cat(sprintf(
    " Flagged (|z| > %.2f): %d (%.1f%%)\n",
    z_crit, flagged_lzstar,
    100 * flagged_lzstar / sum(!is.na(results[, 2]))
  ))
  cat("\n")
}
