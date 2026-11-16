# Approximate moderated Wald intervals for the unshrunk log2 contrast.
add_ql_wald_ci <- function(fit, test, coef = ncol(fit$design), contrast = NULL) {
  design <- fit$design
  if (is.null(contrast)) {
    j <- if (is.character(coef)) match(coef, colnames(design)) else coef
    stopifnot(length(j) == 1L, !is.na(j), j %in% seq_len(ncol(design)))
    contrast <- diag(ncol(design))[, j]
  }
  stopifnot(is.numeric(contrast), is.null(dim(contrast)),
            length(contrast) == ncol(design), all(is.finite(contrast)),
            any(contrast != 0), all(test$df.test == 1))
  if (!is.null(names(contrast))) {
    stopifnot(!anyDuplicated(names(contrast)), setequal(names(contrast), colnames(design)))
    contrast <- contrast[colnames(design)]
  }
  result <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
  mu <- fit$fitted.values
  stopifnot(identical(rownames(result), rownames(mu)))
  n <- nrow(mu)
  expand <- function(x, byrow = FALSE) {
    if (!is.null(dim(x))) {
      stopifnot(identical(dim(x), dim(mu)))
      return(as.matrix(x))
    }
    stopifnot(length(x) %in% c(1L, if (byrow) ncol(mu) else n))
    matrix(x, n, ncol(mu), byrow = byrow)
  }
  phi <- expand(fit$dispersion)
  weights <- if (is.null(fit$weights)) 1 else expand(fit$weights, byrow = TRUE)
  posterior <- if (is.null(fit$s2.post)) fit$var.post else fit$s2.post
  df <- test$df.total
  stopifnot(length(posterior) %in% c(1L, n), length(df) %in% c(1L, n))
  posterior <- rep_len(posterior, n)
  df <- rep_len(df, n)
  beta <- if (is.null(fit$unshrunk.coefficients)) fit$coefficients else fit$unshrunk.coefficients
  center <- drop(beta %*% contrast) / log(2)
  working <- weights * mu / (1 + phi * mu)
  se <- vapply(seq_len(n), function(i) {
    if (!is.finite(posterior[i]) || posterior[i] <= 0 ||
        !isTRUE(df[i] > 0) || any(!is.finite(working[i, ]) | working[i, ] < 0)) {
      return(NA_real_)
    }
    information <- crossprod(design * sqrt(working[i, ]))
    r <- tryCatch(chol(information), error = function(e) NULL)
    if (is.null(r)) return(NA_real_)
    projection <- forwardsolve(t(r), contrast)
    sqrt(posterior[i] * sum(projection^2)) / log(2)
  }, numeric(1))
  se[!is.finite(center) | !is.finite(se) | se <= 0] <- NA_real_
  radius <- rep(NA_real_, n)
  valid <- is.finite(se)
  radius[valid] <- stats::qt(0.975, df[valid]) * se[valid]
  result$logFC_unshrunk <- center
  result$SE <- se
  result$CI_95_low <- center - radius
  result$CI_95_high <- center + radius
  result$QL_df_total <- df
  result$CI_status <- ifelse(valid, "approximate_moderated_wald", "unavailable")
  result
}
