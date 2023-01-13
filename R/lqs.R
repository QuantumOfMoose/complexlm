# file complexlm/R/lqs.R
# copyright (C) 2020-2023 W. L. Ryan
# copyright (C) 1998-2020 B. D. Ripley
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

#' Least Quantile Squared Regression, Sort of Complex Compatible
#' 
#' @inherit MASS::lqs description params return references seealso
#' 
#' @note This method of robust fitting relies on quantiles, which are not defined for complex numbers. While it will accept them and may return an acceptable fit, the accuracy, usefulness,
#' and rigor of it are highly questionable. Please use [rlm] instead.
#' There seems no reason other than historical to use the lms and lqs options. LMS estimation is of low efficiency (converging at rate `n^{-1/3})` whereas LTS has the same asymptotic efficiency as an M estimator with trimming at the quartiles (Marazzi, 1993, p.201). LQS and LTS have the same maximal breakdown value of `(floor((n-p)/2) + 1)/n attained if floor((n+p)/2) <= quantile <= floor((n+p+1)/2)`. The only drawback mentioned of LTS is greater computation, as a sort was thought to be required (Marazzi, 1993, p.201) but this is not true as a partial sort can be used (and is used in this implementation).
#' Adjusting the intercept for each trial fit does need the residuals to be sorted, and may be significant extra computation if n is large and p small.
#' Opinions differ over the choice of psamp. Rousseeuw and Hubert (1997) only consider p; Marazzi (1993) recommends `p+1` and suggests that more samples are better than adjustment for a given computational limit.
#' The computations are exact for a model with just an intercept and adjustment, and for LQS for a model with an intercept plus one regressor and exhaustive search with adjustment. For all other cases the minimization is only known to be approximate.
#' 
#' @export
#' @examples
#' \dontrun{
#' set.seed(4242)
#' n = 8
#' slope = complex(real = 4.23, imaginary = 2.323)
#' intercept = complex(real = 1.4, imaginary = 1.804)
#' x <- complex(real = rnorm(n), imaginary = rnorm(n))
#' y <- slope * x + intercept +complex(real= rnorm(n), imaginary= rnorm(n))/9
#' lqs(x = x, y = y) }
lqs <- function(x, ...) UseMethod("lqs")

####
### I am uncertain how rigorous / effective this is, as quantiles are not defined for complex random variables.
### It appears to work though.
####
#' @describeIn lqs S3 methad for class 'formula'
#'
#' @inherit MASS::lqs.formula params
#'
#' @export
lqs.formula <-
  function(formula, data, ...,
           method = c("lts" ,"lqs", "lms", "S", "model.frame"),
           subset, na.action,
           model = TRUE, x.ret = FALSE, y.ret = FALSE, contrasts = NULL)
  {
    trms <- terms(formula)
    respname <- as.character(attr(trms, "variables")[[attr(trms, "response") + 1]])
    cl <- match.call()
    if (is.complex(data[,respname]) == FALSE)
    {
      cl[[1]] <- MASS::lqs
      eval(cl, parent.frame())
    }
    else
    {
      method <- match.arg(method)
      mf <- match.call(expand.dots = FALSE)
      mf$method <- mf$contrasts <- mf$model <- mf$x.ret <- mf$y.ret <- mf$... <- NULL
      mf[[1L]] <- quote(stats::model.frame)
      mf <- eval.parent(mf)
      if (method == "model.frame") return(mf)
      mt <- attr(mf, "terms")
      y <- model.extract(mf, "response")
      offset <- model.offset(mf)
      if(!is.null(offset)) y <- y - offset
      x <- zmodel.matrix(mt, mf, contrasts)
      contr <- attr(x, "contrasts")
      xint <- match("(Intercept)", colnames(x), nomatch = 0L)
      if(xint) x <- x[, -xint, drop = FALSE]
      fit <- lqs.default(x, y, intercept = (xint > 0), method = method, ...)
      fit$terms <- mt
      fit$call <- match.call()
      fit$contrasts <- contr
      fit$xlevels <- .getXlevels(mt, mf)
      fit$na.action <- attr(mf, "na.action")
      if(model) fit$model <- mf
      if(x.ret) fit$x <- x
      if(y.ret) fit$y <- y
      fit
    }
  }

#' @describeIn lqs Default S3 Method
#'
#' @inherit MASS::lqs.default params
#'
#' @export
lqs.default <-
  function(x, y, intercept = TRUE, method = c("lts", "lqs", "lms", "S"),
           quantile, control = lqs.control(...), k0 = 1.548, seed, ...)
  {
    thiscall <- match.call()
    if(is.numeric(x)) # If the given data is numeric, call the lqs.default function from MASS.
    {
      thiscall[[1]] <- MASS::lqs.default
      eval(thiscall, parent.frame())
    }
    else
    {
    lqs.control <- function(psamp = NA, nsamp = "best", adjust = TRUE)
      list(psamp = psamp, nsamp = nsamp, adjust = adjust)
    
    n <- length(y)
    nmx <- deparse(substitute(x))
    if(is.null(dim(x))) {
      x <- as.matrix(x)
      colnames(x) <- nmx
    } else x <- as.matrix(x)
    p <- ncol(x)
    if(any(is.na(x)) || any(is.na(y)))
      stop("missing values are not allowed")
    nm <- colnames(x)
    if(is.null(nm))
      nm <- if(p > 1) paste("X", 1L:p, sep="") else if(p == 1) "X" else NULL
    if(intercept) {
      att <- attr(x, "contrasts")
      x <- cbind(1, x)
      nm <- c("(Intercept)", nm)
      attr(x, "contrasts") <- att
    }
    p <- ncol(x)
    if(nrow(x) != n) stop("'x' and 'y' must have the same number of rows")
    method <- match.arg(method)
    lts <- 0; beta <- 0
    if(method == "lqs" && missing(quantile)) quantile <- floor((n+p+1)/2)
    if(method == "lms") quantile <- floor((n+1)/2)
    if(method == "lts") {
      lts <- 1
      if(missing(quantile)) quantile <- floor(n/2) + floor((p+1)/2)
    }
    if(method == "S") {
      ### chi is also in lqs.c and zlqs.c, why define it twice in two different languages?
      lts <- 2
      beta <- 0.5
      quantile <- ceiling(n/2)
      chi <- function(u, k0)
      { u <- as.numeric(Conj(u/k0)*(u/k0)); ifelse(u < 1, 3*u - 3*u^2 + u^3, 1) }
    }
    if(quantile > n-1)
      stop(gettextf("'quantile' must be at most %d", n-1),
           domain = NA)
    ps <- control$psamp
    if(is.na(ps)) ps <- p
    if(ps < p) {
      ps <- p
      warning("'ps' must be at least 'p'")
    }
    adj <- control$adjust & intercept
    nsamp <- eval(control$nsamp)
    nexact <- choose(n, ps)
    if(is.character(nsamp) && nsamp == "best") {
      nsamp <- if(nexact < 5000) "exact" else "sample"
    } else if(is.numeric(nsamp) && nsamp > nexact) {
      warning(sprintf(ngettext(nexact,
                               "only %d set, so all sets will be tried",
                               "only %d sets, so all sets will be tried"),
                      nexact), domain = NA)
      nsamp <- "exact"
    }
    samp <- nsamp != "exact"
    if(samp) {
      if(nsamp == "sample") nsamp <- min(500*ps, 3000)
    } else
      nsamp <- nexact
    
    if(samp && !missing(seed)) {
      if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))  {
        seed.keep <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
        on.exit(assign(".Random.seed", seed.keep, envir=.GlobalEnv))
      }
      assign(".Random.seed", seed, envir=.GlobalEnv)
    }
    z <- .C(zlqs_fitlots, 
                 as.complex(x), as.complex(y), as.integer(n), as.integer(p),
                 as.integer(quantile), as.integer(lts), as.integer(adj),
                 as.integer(samp), as.integer(ps), as.integer(nsamp),
                 crit=complex(1), sing=integer(1L), bestone=integer(ps),
                 coefficients=complex(p), as.complex(k0), as.complex(beta)
             )[c("crit", "sing", "coefficients", "bestone")]
    if(z$sing == nsamp)
      stop("'lqs' failed: all the samples were singular", call.=FALSE)
    z$sing <- paste(z$sing, "singular samples of size", ps, "out of", nsamp)
    z$bestone <- sort(z$bestone) #If complex this will sort by real part first, then imaginary.
    names(z$coefficients) <- nm
    fitted <- drop(x %*% z$coefficients)
    z$fitted.values <- fitted
    z$residuals <- y - fitted
    c1 <- 1/qnorm((n + quantile)/(2*n))
    s <-
      if(lts == 1)
        sqrt(abs(z$crit)/quantile)/sqrt(1 - 2*n*dnorm(1/c1)/(quantile*c1)) # s will be complex if using complex data. If data isn't complex, it will be real positive. S needs to be real positive in all cases.
    else if(lts == 0) sqrt(abs(z$crit))*c1 else abs(z$crit) # Just put z$crit in abs() to make it positive real, it's worth a shot, right?
    res <- z$residuals
    ind <- abs(res) <= 2.5*abs(s) # Since s will be positive for real data, abs(s) = s in that case.
    s2 <- sum(Conj(res[ind])*res[ind])/(sum(ind) - p)
    z$scale <- c(s, sqrt(s2))
    if(method == "S") { # IWLS refinement
      psi <- function(u, k0) (1  - pmin(1, abs(u/k0))^2)^2
      resid <- z$residuals
      scale <- s
      for(i in 1L:30L) {
        w <- psi(resid/scale, k0)
        if (is.complex(x)) temp <- zlm.wfit(x, y, w, method="qr") else temp <- lm.wfit(x, y, w, method="qr")
        resid <- temp$residuals
        s2 <- scale*sqrt(abs(sum(chi(resid/scale, k0))/((n-p)*beta))) # Again, a questionable abs() to keep s2 real positive.
        if (is.complex(s2)) if(abs(abs(s2/scale) - 1) < 1e-5) break
        else if(abs(s2/scale - 1) < 1e-5) break 
        scale <- s2
      }
      z$coefficents <- temp$coefficients
      z$fitted.values <- temp$fitted.values
      z$residuals <- resid
      z$scale <- scale
    }
    class(z) <- "lqs"
    z
  }
}
# Probably don't need this.
# print.lqs <- function (x, digits = max(3, getOption("digits") - 3), ...)
# {
#   if(!is.null(cl <- x$call)) {
#     cat("Call:\n")
#     dput(cl, control=NULL)
#     cat("\n")
#   }
#   cat("Coefficients:\n")
#   print.default(format(coef(x), digits = digits), print.gap = 2,
#                 quote = FALSE)
#   cat("\nScale estimates", format(x$scale, digits = digits) ,"\n\n")
#   invisible(x)
# }

### I don't think a seperate predict.lqs is necessary for complex data....
# predict.lqs <- function (object, newdata, na.action = na.pass, ...)
# {
#   if (missing(newdata)) return(fitted(object))
#   ## work hard to predict NA for rows with missing data
#   Terms <- delete.response(terms(object))
#   m <- model.frame(Terms, newdata, na.action = na.action,
#                    xlev = object$xlevels)
#   if(!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
#   X <- model.matrix(Terms, m, contrasts = object$contrasts)
#   drop(X %*% object$coefficients)
# }

#### cov.rob as previously written is filled with calls to functions that demand numeric (real) valued inputs.
###   Here we make a new cov.rob that checks if x is real. If so, it calls the original cov.rob from MASS.
###   If x is complex, call cov.zrob, which has been modified for complex numbers.
#' Resistant Estimation of Multivariate Location and Scatter, for Numeric or Complex Variables
#' 
#' @inherit MASS::cov.rob description params details return
#'
#' @export
#' 
#' @note This function relies on interquartile range, which orders the inputs as part of the calculation. The concept of ordering does not translate to the complex numbers,
#' so this function applies interquartile range to the real and imaginary components separately. The results of this will not be rotationally invariant, so this function 
#' should be used with caution for complex variables.
#'
#' @examples
#' n <- 16
#' x <- matrix(rnorm(n), ncol = 2)
#' cov.rob(x)
cov.rob <- function(x, cor = FALSE, quantile.used = floor((n+p+1)/2),
                    method = c("mve", "mcd", "classical"), nsamp = "best", seed)
{
  thiscall <- match.call()
  method <- match.arg(method)
  x <- as.matrix(x)
  if(any(is.na(x)) || any(is.infinite(x)))
    stop("missing or infinite values are not allowed")
  n <- nrow(x); p <- ncol(x)
  if(n < p+1)
    stop(gettextf("at least %d cases are needed", p+1), domain = NA)
  if (is.numeric(x)) 
  {
    thiscall[[1]] <- MASS::cov.rob
    eval(thiscall, parent.frame())
  }
  else if (is.complex(x)) 
  {
    thiscall[[1]] <- cov.zrob
    eval(thiscall, parent.frame())
  }
  else print("Input x must be numeric or complex.")
}

### A function to calculate the unbiased sample variance of a vector of complex numbers.
### This duplicates the functionality of the masked var function. Consider removing.
### Depreciated.
#zvar <- function(x)
#{
#  sampmean <- mean(x, trim = 0)
#  return((1 / (length(x) - 1)) * sum(as.numeric((x - sampmean) * Conj(x - sampmean))))
#}

### A function for calculating the unbiased sample pseudo-variance of a vector of complex numbers.
### Can return a complex number.
### Not used in anything else at the moment.
#' Pseudo-Variance of Complex Variables
#' 
#' Calculates the pseudo-variance, also called the relational variance of a vector of complex numbers.
#' This describes the degree of covariance between the real and imaginary components.
#'
#' @param x a vector of complex numbers.
#'
#' @return complex. The pseudo-variance of `x`. If `x` is numeric, this is just the variance.
#' @export
#'
#' @examples
#' n = 6
#' z <- complex(real = rnorm(n), imaginary = rnorm(n))
#' pseuzvar(z)
pseuzvar <- function(x)
{
  sampmean <- mean(x, trim = 0)
  return((1 / (length(x) - 1)) * sum((x - sampmean) * (x - sampmean)))
}

cov.zrob <- function(x, cor = FALSE, quantile.used = floor((n+p+1)/2),
                    method = c("mve", "mcd", "classical"), nsamp = "best", seed)
{
  method <- match.arg(method)
  x <- as.matrix(x)
  if(any(is.na(x)) || any(is.infinite(x)))
    stop("missing or infinite values are not allowed")
  n <- nrow(x); p <- ncol(x)
  if(n < p+1)
    stop(gettextf("at least %d cases are needed", p+1), domain = NA)
  if(method == "classical") {
    ans <- list(center = colMeans(x), var = var(x)) # This is the var() from this package that is compatible with complex numbers.
  } else {
    if(quantile.used < p+1)
      stop(gettextf("'quantile' must be at least %d", p+1), domain = NA)
    if(quantile.used > n-1)
      stop(gettextf("'quantile' must be at most %d", n-1), domain = NA)
    ## re-scale to roughly common scale
    divisor <- complex(real = apply(Re(x), 2, IQR), imaginary = apply(Im(x), 2, IQR)) # I'm uncertain how to adapt this to complex numbers... Interquartile range doesn't make much sense for
            # complex or >= 2 dimensional numbers, it'd have to some kind of area...
    if(any(divisor == 0)) stop("at least one column has IQR 0")
    x <- x /rep(divisor, rep(n,p))
    qn <- quantile.used
    ps <- p + 1
    nexact <- choose(n, ps)
    if(is.character(nsamp) && nsamp == "best")
      nsamp <- if(nexact < 5000) "exact" else "sample"
    if(is.numeric(nsamp) && nsamp > nexact) {
      warning(sprintf(ngettext(nexact,
                               "only %d set, so all sets will be tried",
                               "only %d sets, so all sets will be tried"),
                      nexact), domain = NA)
      nsamp <- "exact"
    }
    samp <- nsamp != "exact"
    if(samp) {
      if(nsamp == "sample") nsamp <- min(500*ps, 3000)
    } else nsamp <- nexact
    if (nsamp > 2147483647) {
      if(samp)
        stop(sprintf("Too many samples (%.3g)", nsamp))
      else
        stop(sprintf('Too many combinations (%.3g) for nsamp = "exact"', nsamp))
    }
    if(samp && !missing(seed)) {
      if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))  {
        seed.keep <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
        on.exit(assign(".Random.seed", seed.keep, envir=.GlobalEnv))
      }
      assign(".Random.seed", seed, envir=.GlobalEnv)
    }
    z <-  .C(zmve_fitlots,
             as.complex(x), as.integer(n), as.integer(p),
             as.integer(qn), as.integer(method=="mcd"),
             as.integer(samp), as.integer(ps), as.integer(nsamp),
             crit=complex(1), sing=integer(1L), bestone=integer(n))
    z$sing <- paste(z$sing, "singular samples of size", ps,
                    "out of", nsamp)
    crit <- z$crit + 2*sum(log(divisor)) +
      if(method=="mcd") - p * log(qn - 1) else 0
    best <- seq(n)[z$bestone != 0]
    if(!length(best)) warning("'x' is probably collinear") # Suddenly the example always fails on this for some reason. What happens if I let it continue?
    means <- colMeans(x[best, , drop = FALSE])
    rcov <- var(x[best, , drop = FALSE]) * (1 + 15/(n - p))^2
    dist <- mahalanobis(x, means, rcov) # This maybe should return a real number. Sneaky, this would not work if inverted = T in mahalanobis().
    cut <- qchisq(0.975, p) * quantile(dist, qn/n)/qchisq(qn/n, p)
    cov <- divisor * complexlm::var(x[dist < cut, , drop = FALSE]) *
      rep(divisor, rep(p, p))
    attr(cov, "names") <- NULL
    ans <- list(center =
                  colMeans(x[dist < cut, , drop = FALSE]) * divisor,
                cov = cov, msg = z$sing, crit = crit, best = best)
  }
  if(cor) {
    sd <- sqrt(diag(ans$cov))
    ans <- c(ans, list(cor = (ans$cov/sd)/rep(sd, rep(p, p))))
  }
  ans$n.obs <- n
  ans
}

## compatibility functions for R users.
## Thought: Are these needed if I am writing the functions in R? Do I need
## compatibility functions for S users instead?

#' Comparability Wrapper for lqs
#' 
#' This is just a wrapper for [lqs]. See [lqs] for documentation.
#' 
#' @inherit MASS::ltsreg params
#' 
#' @export
lmsreg <- function(...)
{
  oc <- sys.call()
  oc$method <- "lms"
  oc[[1L]] <- quote(complexlm::lqs)
  eval.parent(oc)
}

#' Comparability Wrapper for lqs
#' 
#' @describeIn lmsreg comparability wrapper for `lqs`
#' 
#' @inherit MASS::ltsreg params
#' 
#' @export
ltsreg <- function(...)
{
  oc <- sys.call()
  oc$method <- "lts"
  oc[[1L]] <- quote(complexlm::lqs)
  eval.parent(oc)
}

#' Comparability Wrapper
#' 
#' @describeIn cov.rob comparability wrapper for `cov.rob`
#' @export
cov.mve <- function(...)
{
  oc <- sys.call()
  oc$method <- "mve"
  oc[[1L]] <- quote(complexlm::cov.rob)
  eval.parent(oc)
}

#' Comparability Wrapper
#' @describeIn cov.rob comparability wrapper for `cov.rob`
#' @export
cov.mcd <- function(...)
{
  oc <- sys.call()
  oc$method <- "mcd"
  oc[[1L]] <- quote(complexlm::cov.rob)
  eval.parent(oc)
}