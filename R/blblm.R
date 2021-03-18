#' @import purrr
#' @import furrr
#' @import stats
#' @import parallel
#' @importFrom magrittr %>%
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"

## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))


#' @title blblm is used to implement the bag of little bootstraps for the linear regression model.
#' @name blblm
#' @param formula an object of class formula giving relationship of variables in model
#'
#' @param data An optional data frame, list or environment containing variables for the linear model.
#' @param m An optional integer to specify the number of splits (m) in the data set.
#' @param B An optional integer for each subsample, sample observations with replacement, repeat for (B) times.
#' @param parallel logical. If TRUE the linear model with be fitted using parallel computing. Function will follow single process if FALSE.
#' @param nthreads n optional integer to specify the number of cores used only if parallel is TRUE.
#'
#' @examples
#' fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
#' fit2 <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = TRUE, nthreads = 4)
#'
#'
#' @exportbl
blblm <- function(formula, data, m = 10, B = 5000, parallel = FALSE, nthreads = detectCores()) {

  data_list <- split_data(data, m)

  set.seed(123)

  if(parallel){
    #suppressWarnings(plan(multiprocess, workers = nthreads))
    options(future.rng.onMisuse = "ignore")
    if(nthreads<2){
      nthreads = 2
    }
    else if(nthreads>detectCores()){
      nthreads = detectCores()
    }
    cl <- makeCluster(nthreads) #sets seed to be the same as the single process

    clusterSetRNGStream(cl,123)

    estimates <- future_map(data_list,
            ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))

    stopCluster(cl)
  }

  else{
    # set seed for testing purposes
    estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
    }
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}


#' split data into m parts of approximated equal sizes
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' compute the estimates
lm_each_subsample <- function(formula, data, n, B) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  m <- model.frame(formula, data)
  X <- model.matrix(formula, m)
  y <- model.response(m)
  replicate(B, lm1(X, y, n), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset
lm1 <- function(X, y, n) {
  freqs <- as.vector(rmultinom(1, n, rep(1, nrow(X))))
  fit <- lm.wfit(X, y, freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' compute the coefficients from fit
blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit
blbsigma <- function(fit) {
  p <- fit$rank
  e <- fit$residuals
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' @title print.blblm is used to print the blblm model
#' @name print.blblm
#' @param x represents the blblm model or object to be printed
#'
#' @param ... additional arguments to be passed
#'
#' @examples
#' fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
#' print.blblm(fit)
#'
#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' @title sigma.blblm is used to obtain the residual standard deviation
#' @name sigma.blblm
#'
#' @param object an R object typically from a linear model
#'
#' @param confidence logical. If TRUE will output the confidence interval for the residual standard deviation
#' @param level an alpha level for the confidence interval
#' @param ... additional arguments to be passed
#'
#'@examples
#' fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
#' sigma(fit)
#' sigma(fit, confidence = TRUE)
#'
#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' @title coef.blblm is a function to extract the model coefficients from the blblm model
#' @name coef.blblm
#' @param object an object which has relevant coefficients to be extracted
#'
#' @param ... additional arguments to be passed
#'
#' @examples
#' fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
#' coef(fit)
#'
#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}



#' @title confint.blblm is a function to compute confidence interval for one or multiple parameters in the blblm  model
#' @name confint.blblm
#'
#' @param object a blblm model
#'
#' @param parm the parameters given to calculate the confidence intervals
#' @param level an alpha level for the confidence interval
#' @param ... additional arguments to be passed
#'
#' @examples
#' fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
#' confint(fit, c("wt", "hp"))
#'
#' @export
#' @method confint blblm
#'
#'
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}


#' @title predict.blblm is a function for predictions from blblm model for new data
#' @name predict.blblm
#'
#' @param object a blblm model
#'
#'
#' @param new_data an data frame, list or env containing new data to be used for the prediction
#' @param confidence logical. If TRUE gives confidence interval of new prediction on blblm linear model/
#' @param level an alpha level for the confidence interval
#' @param ... additional arguments to be passed
#'
#' @examples
#' fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
#' predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)))
#' predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)), confidence = TRUE)
#'
#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
