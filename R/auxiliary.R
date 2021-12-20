#' This function prepares necessary list of information for fitting IPW or AIPW
#' @keywords internal
#' @return returns a list containing the following components:
#' \item{gamma_para}{the fitted values for the first stage regression}
#' \item{gamma_IF}{the influence functions of the fitted values for the first stage regression}
#' \item{sigma_star}{the positivity parameter used in the link function}
#' \item{prop}{the propensities of fully observed for each R = 1 subject}
#' @export
PS <- function(data, R, pattern){
  data <- as.matrix(data)
  N <- nrow(data)
  sigma_star_set <- 10^(-c(5:9))
  M = nrow(pattern)
  work_df <- c()
  init_para <- c()
  selection <- c(0)
  for (m in 2:M) {
    work_df <- cbind(work_df, cbind(1, data[, pattern[m, ] == 1]))
    init_para <- c(init_para, log(1 / M), rep(0, sum(pattern[m, ])))
    selection <- c(selection, sum(pattern[m, ]) + 1)
  }
  work_df[is.na(work_df)] = 0
  selection <- cumsum(selection)
  para_result <- c()
  obj_value <- c()
  for (sigma_star in sigma_star_set) {
    fr <- function(x) {
      g <- 0
      tmp <- 0
      for(m in 2:M) {
        tmp_x <- rep(0, length(x))
        tmp_x[(selection[m - 1] + 1):selection[m]] <- x[(selection[m - 1] + 1):selection[m]]
        score <- as.matrix(work_df) %*% tmp_x
        g <- g + mean((R == m) * log(link_func_vec(score, sigma_star)))
        tmp <- tmp + link_func_vec(score, sigma_star)
      }
      g <- g + mean((R == 1) * log(pmax(1 - tmp, sigma_star)))
      return(-g)
    }

    gr <- function(x) {
      tmp <- 0
      for(m in 2:M) {
        design_mat <- work_df[R == 1, (selection[m - 1] + 1):selection[m]]
        score <- as.matrix(design_mat) %*% x[(selection[m - 1] + 1):selection[m]]
        tmp <- tmp + link_func_vec(score, sigma_star)
      }
      inequality <- c(tmp - 1 + sigma_star)
      return(inequality)
    }

    opts <- list( "algorithm" = "NLOPT_LN_COBYLA",
                  "xtol_rel"  = 1.0e-4,
                  "maxeval"   = 100000)
    fit <- nloptr(x0 = init_para,
                  eval_f = fr,
                  lb          = rep(-4, length(init_para)),
                  ub          = rep(4, length(init_para)),
                  eval_g_ineq = gr,
                  opts = opts)
    para_result <- rbind(para_result, fit$solution)
    obj_value <- c(obj_value, fit$objective)
  }
  minimizer <- which.min(obj_value)
  para <- para_result[minimizer, ]
  sigma_star <- sigma_star_set[minimizer]
  prop <- prop_compute(para = para, work_df = work_df, R = R, selection = selection, sigma_star = sigma_star, pattern = pattern)
  grad_fr <- function(x) {
    g <- 0
    tmp <- 0
    for(m in 2:M) {
      tmp_x <- rep(0, length(x))
      tmp_x[(selection[m - 1] + 1):selection[m]] <- x[(selection[m - 1] + 1):selection[m]]
      temp_work_df <- matrix(0, ncol = ncol(work_df), nrow = nrow(work_df))
      temp_work_df[, (selection[m - 1] + 1):selection[m]] <- work_df[, (selection[m - 1] + 1):selection[m]]
      score <- work_df %*% tmp_x
      g <- g + -(R == m)/link_func_vec(score, sigma_star) * link_func_grad_vec(score, sigma_star) * temp_work_df / N
      tmp <- tmp + link_func_vec(score, sigma_star)
    }
    for(m in 2:M) {
      tmp_x <- rep(0, length(x))
      tmp_x[(selection[m - 1] + 1):selection[m]] <- x[(selection[m - 1] + 1):selection[m]]
      temp_work_df <- matrix(0, ncol = ncol(work_df), nrow = nrow(work_df))
      temp_work_df[, (selection[m - 1] + 1):selection[m]] <- work_df[, (selection[m - 1] + 1):selection[m]]
      score <- work_df %*% tmp_x
      g <- g + (R == 1) / pmax(1 - tmp, 1 - sigma_star) * link_func_grad_vec(score, sigma_star) * temp_work_df / N
    }
    return(t(g))
  }
  res <- grad_fr(para)
  gamma_jacob <- -jacobian(func = function(x) {return(rowSums(grad_fr(x)))}, x = para)
  IF <- gamma_jacob %*% res
  fit <- list(gamma_para = para,
              gamma_IF = IF,
              sigma_star = sigma_star,
              prop = prop)
  return(fit)
}

#' This function computes propensities of entries being fully observed
#' @keywords internal
#' @return returns propensities of being fully observed
#' @export
prop_compute <- function(para, work_df, R, selection, sigma_star, pattern) {
  tmp <- 0
  M <- nrow(pattern)
  for(m in 2:M) {
    design_mat <- work_df[R == 1, (selection[m - 1] + 1):selection[m]]
    score <- as.matrix(design_mat) %*% para[(selection[m - 1] + 1):selection[m]]
    tmp <- tmp + link_func_vec(score, sigma_star)
  }
  prop <- 1 - tmp
  return(prop)
}


#' This function prepares the user-provided data into the correct format
#' @keywords internal
#' @return returns a list containing the following components:
#' \item{R}{the indicator of the missing pattern, where R = 1 is reserved for fully observed}
#' \item{pattern}{a matrix recording the pattern}
#' @export
nmm_preprocessing <- function(data, O) {
  n <- nrow(data)
  K <- ncol(data)
  O <- 1 - O
  binary <- function(x) {
    num <- 0
    for (j in 0:(K - 1)) {
      num <- num + x[j + 1] * 2^j
    }
   return(num)
  }
  R <- apply(O, 1, binary) + 1
  pattern <- c()
  for (r in sort(unique(R))) {
    tmp <- rep(0, K)
    r <- r - 1
    for (j in 1:K) {
      tmp[j] = r %% 2
      r = r %/% 2
    }
    pattern <- rbind(pattern, tmp)
  }
  pattern <- 1 - pattern
  R_sortunique <- sort(unique(R))
  for (index in 1:length(R_sortunique)) {
    R[which(R == R_sortunique[index])] <- index
  }
  nmm_pre_result <- list(R = R,
                         pattern = pattern)
  return(nmm_pre_result)
}

#' This function computes our link function that respects the positivity
#' @keywords internal
#' @return returns results our link function
#' @export
link_func <- function(x, sigma_star) {
  if(is.nan(x))
    return(1)
  if (x <= 0) {
    val <- 1/(1 + exp(-x)) * (1 - sigma_star) / 0.5
  } else {
    val <- 1 - 1/(1 + exp(x)) * sigma_star / 0.5
  }
  return(val)
}

#' This function vectorizes link_func
#' @keywords internal
#' @return returns a vector of link_func results
#' @export
link_func_vec <- Vectorize(link_func)

#' This function computes the gradient of our link function that respects the positivity
#' @keywords internal
#' @return returns the gradient of our link function
#' @export
link_func_grad <- function(x, sigma_star) {
  if(is.nan(x))
    return(0)
  if (x <= 0) {
    val <- 1/(1 + exp(-x)) * 1/(1 + exp(x)) * (1 - sigma_star) / 0.5
  } else {
    val <- 1/(1 + exp(-x)) * 1/(1 + exp(x)) * sigma_star / 0.5
  }
  return(val)
}

#' This function vectorizes link_func_grad
#' @keywords internal
#' @return returns a vector of link_func_grad results
#' @export
link_func_grad_vec <- Vectorize(link_func_grad)


#' This function returns the residuls after fiting a func regression, where func is the user-specified function
#' @keywords internal
#' @return returns the residuals of the regression by the user-specified function
#' @export
regress_fit <- function(coefficients, regress_func, formula, data, weights, R, sigma_star, pattern, ...) {
  data <- as.matrix(data)
  M = nrow(pattern)
  N <- nrow(data)
  work_df <- c()
  selection <- c(0)
  for (m in 2:M) {
    work_df <- cbind(work_df, cbind(1, data[, pattern[m, ] == 1]))
    selection <- c(selection, sum(pattern[m, ]) + 1)
  }
  work_df[is.na(work_df)] = 0
  tmp_prop <- prop_compute(para = coefficients, work_df = work_df, R = R, selection = selection, sigma_star = sigma_star, pattern = pattern)
  prop <- rep(1, N)
  prop[R == 1] <- tmp_prop
  prop <- pmax(pmin(prop, 1 - 1e-8), 1e-8)
  data <- as.data.frame(data)
  data$weights <- weights * (R == 1) / prop
  fit <- regress_func(formula = formula, data = data, weights = weights, ...)
  bread_fit <- attr(iid(fit), "bread")
  IF <- iid(fit)
  res <- t(IF %*% solve(bread_fit))
  return(res)
}

