#' @title Fitting IPW or AIPW Estimators under Nonmonotone Missing at Random Data
#' @description nmm_fit is the main function used to fit IPW or AIPW estimators under nonmonotone missing at random data
#' @param data a data.frame to fit
#' @param O missing indicator
#' @param AIPW indicator if fitting augmented IPW
#' @param formula optional formula specified to fit
#' @param func optional fitting function, currently support 'lm' and 'glm'
#' @param weights optional weights used in the estimation
#' @param ... further arguments passed to func, e.g. family = 'quasibinomial' for glm
#' @return NMMIPW returns an object of class "NMMIPW".
#' An object of class "NMMIPW" is a list containing the following components:
#' \item{coefficients}{the fitted values, only reported when formula and func are given}
#' \item{coef_sd}{the standard deviations of coefficients, only reported when formula and func are given}
#' \item{coef_IF}{the influnece function of coefficients, only reported when formula and func are given}
#' \item{gamma_para}{the first step fitted valus}
#' \item{AIPW}{an indicator of whether AIPW is fitted}
#' \item{second_step}{an indicator of whether the second step is fitted}
#' \item{second_fit}{if second step fitted, we report the fit object}
#' \item{by_prod}{a list of by products that might be useful for users, including first step IF, jacobian matrices}
#' @importFrom nloptr nloptr
#' @importFrom lava iid
#' @importFrom numDeriv jacobian
#' @export nmm_fit
#' @examples
#' n = 100
#' X = rnorm(n, 0, 1)
#' Y = rnorm(n, 1 * X, 1)
#' O1 = rbinom(n, 1, 1/(1 + exp(- 1 - 0.5 * X)))
#' O2 = rbinom(n, 1, 1/(1 + exp(+ 0.5 + 1 * Y)))
#' O = cbind(O1, O2)
#' df <- data.frame(Y = Y, X = X)
#' fit <- nmm_fit(data = df, O = O, formula = Y ~ X, func = lm)
nmm_fit <- function(data, O, AIPW = FALSE, formula = NULL, func = NULL, weights = NULL, ...) {
  N <- nrow(data)
  if (is.null(weights))
    weights <- rep(1, N)
  prop <- rep(1, N)
  if (sum(1 - O) != 0) {
    nmm_pre_result <- nmm_preprocessing(as.matrix(data), O)
    R <- nmm_pre_result$R
    pattern <- nmm_pre_result$pattern
    ## first step propensity estimation
    first_fit <- PS(data, R, pattern)
    gamma_para <- first_fit$gamma_para
    gamma_IF <- first_fit$gamma_IF
    sigma_star <- first_fit$sigma_star
    prop[R == 1] <- first_fit$prop
    nloptrresult <- first_fit$nloptrresult
    prop <- pmax(pmin(prop, 1 - 1e-8), 1e-8)
    by_prod <- list(gamma_IF = gamma_IF,
                    nloptrresult = nloptrresult)
  }
  if (!is.null(formula) && !is.null(func)) {
    ## second step IPW or AIPW
    data$weights = weights * (R == 1) / prop
    data[is.na(data)] = 0
    fit <- func(formula = formula, data = data, weights = weights, ...)
    gamma_jacobian <- jacobian(func = function(...) {
                       return(rowSums(regress_fit(...)))},
                   x = gamma_para, regress_func = func, formula = formula, data = data, weights = weights,
                   R = R, sigma_star = sigma_star, pattern = pattern, ...)
    coefficients <- fit$coefficients
    summary_fit <- summary(fit)
    bread_mat <- attr(iid(fit), "bread")
    dispersion <- summary_fit$dispersion
    if(is.null(dispersion))
      dispersion <- 1
    coef_IF <- (t(iid(fit)) + bread_mat %*% gamma_jacobian %*% gamma_IF) / dispersion
    coef_sd <- sqrt(rowSums(coef_IF^2))
    by_prod$bread_mat = bread_mat
    by_prod$gamma_jacobian = gamma_jacobian
    res <- list(coefficients = coefficients,
                coef_sd = coef_sd,
                coef_IF = coef_IF,
                gamma_para = gamma_para,
                AIPW = AIPW,
                second_step = TRUE,
                second_fit = fit,
                by_prod = by_prod)
    if (AIPW) {
      warning("AIPW is ongoing.")
    }
  } else {
    cat("Only run first step estimation.")
    res <- list(gamma_para = gamma_para,
                second_step = FALSE,
                AIPW = AIPW,
                by_prod = by_prod)
  }

  class(res) <- "NMMIPW"
  return(res)
}





#' @title Summarizing IPW or AIPW Estimators under Nonmonotone Missing at Random Data
#' @param object an object of class "NMMIPW", usually, a result of a call to NMMIPW.
#' @param x an object of class "summary.NMMIPW", usually, a result of a call to summary.NMMIPW.
#' @param ... further arguments passed to or from other methods.
#' @description summary method for class "NMMIPW".
#' @importFrom stats pnorm
#' @export summary.NMMIPW
#' @method summary NMMIPW
#' @S3method summary NMMIPW
#' @examples
#' n = 100
#' X = rnorm(n, 0, 1)
#' Y = rnorm(n, 1 * X, 1)
#' O1 = rbinom(n, 1, 1/(1 + exp(-1 - 0.5 * X)))
#' O2 = rbinom(n, 1, 1/(1 + exp(+0.5 + 1 * Y)))
#' O = cbind(O1, O2)
#' df <- data.frame(Y = Y, X = X)
#' fit <- nmm_fit(data = df, O = O, formula = Y ~ X, funct = lm)
#' summary(fit)
#' @details print.summary.NMMIPW tries to be smart about formatting coefficients, an estimated variance covariance matrix of
#' the coefficients, Z-values and the corresponding P-values.
#' @return The function summary.NMMIPW computes and returns a list of summary statistics of the fitted model given in object.
summary.NMMIPW <- function(object, ...){
  est <- object$coefficients
  se <- object$coef_sd
  zval <- est / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  res <- list(est = est,
              se = se,
              zval = zval,
              pval = pval)

  res$AIPW <- object$AIPW
  res$call <- object$second_fit$call
  res$second_step <- object$second_step

  class(res) = "summary.NMMIPW"
  res
}


#' @rdname summary.NMMIPW
#' @S3method print summary.NMMIPW
print.summary.NMMIPW <- function(x, ...){
  obj <- x
  if (obj$second_step) {
    # We print information about object:
    if (obj$AIPW) {
      cat("Fitted an AIPW:\n")
    } else {
      cat("Fitted an IPW:\n")
    }
    res <- cbind(obj$est, obj$se)
    zval <- obj$zval
    pval <- obj$pval
    res <- cbind(res, zval, pval)
    colnames(res) <- c("coef", "se(coef)", "z-value", "p-value")
    prmatrix(signif(res, 3))
    cat("   \n")
  }
  cat("Only the first step is fitted and thus there is nothing to summarize. \n")
  invisible()
}
