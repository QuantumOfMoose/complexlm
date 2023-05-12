# file complexlm/R/rlm.R
# copyright (C) 2020-2023 W. L. Ryan
# copyright (C) 1994-2020 W. N. Venables and B. D. Ripley
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

#@useDynLib complexlm, .registration = TRUE # Not needed since lqs.R functions moved to experimental branch.

#' Robust Fitting of Linear Models, Compatible with Complex Variables
#' 
#' Uses robust M-estimation to fit a linear model to numeric or complex data. Based on [MASS::rlm()].
#' 
#' @inherit MASS::rlm details return references
#' 
#' @details M-estimation works by finding the model coefficients that minimize the sum of a function of the residuals. 
#' This function, called the objective function rho(), is a kind of statistical distance (AKA divergence), and a semimetric.
#' As a semimetric it is a function of the measured value `y` and the modeled value `Y` (residual r = y - Y) which maps from 
#' the space of the data to the positive real numbers. Semimetrics can be defined for domains of any dimensionality, including the 
#' two dimensional complex numbers, and thus so can M-estimation. 
#' What's more, all the standard algebraic operations used in the itteratively (re)weighted least-squares M-estimator robust regression
#' algorithm are defined over the set of complex numbers. While ordering is not defined for them, it is the output of rho(), a real number, that must be 
#' in M-estimation.
#' 
#' @import MASS
#' @import stats
#' @export
rlm <- function(x, ...) UseMethod("rlm")

#' @describeIn rlm S3 method for class 'formula'
#'
#' @param formula a [formula] object of the form y ~ x1 + x2. Note that algebraic expressions in formula cannot currently be used with complex data.
#' @param data a data frame containing the variables upon which a robust fit is to be applied.
#' @param weights numeric. A vector of weights to apply to the residuals.
#' @param ... additional arguments to be applied to called functions (rlm.default and the psi function).
#' @param subset an index vector specifying the cases (rows of data or x and y) to be used for fitting.
#' @param na.action a function that specifies what to do if NAs are found in the fitting data. The default is to omit them via [na.omit]. Can also be changed by [options] (na.action =).
#' @param method string. What method of robust estimation should be used. Options are "M", "MM", or "model.frame". The default is M-estimation. MM-estimation has a high breakdown point but is not compatible with complex variables or case weights. model.frame just constructs the model frame, and only works with the formula method.
#' @param wt.method string, either "inv.var" or "case". Specifies whether the weights are case weights that give the relative importance of each observation (higher weight means more important) / case, or the inverse variances of the cases (higher weight means that observation is less variable / uncertain).
#' @param model logical. Should the model frame be included in the returned object?
#' @param x.ret logical. Should the model (design) matrix be included in the returned object?
#' @param y.ret logical. Should the response be included in the returned object?
#' @param contrasts optional contrast specifications: see [stats::lm]. Not compatible with complex data.
#'
#' @export
#' 
#' @examples
#' set.seed(4242)
#' n <- 8
#' slope <- complex(real = 4.23, imaginary = 2.323)
#' interc <- complex(real = 1.4, imaginary = 1.804)
#' e <- complex(real=rnorm(n)/6, imaginary=rnorm(n)/6)
#' xx <- complex(real= rnorm(n), imaginary= rnorm(n))
#' tframe <- data.frame(x = xx, y= slope*xx + interc + e)
#' rlm(y ~ x, data = tframe, weights = rep(1,n))
rlm.formula <-
    function(formula, data, weights, ..., subset, na.action,
             method = c("M", "MM", "model.frame"),
             wt.method = c("inv.var", "case"),
             model = TRUE, x.ret = TRUE, y.ret = FALSE, contrasts = NULL)
{
    trms <- terms(formula)
    respname <- as.character(attr(trms, "variables")[[attr(trms, "response") + 1]])
    cl <- match.call()
    if (is.complex(data[,respname]) == FALSE)
    {
      cl[[1]] <- MASS::rlm
      eval(cl, parent.frame())
    }
    else
    {
      mf <- match.call(expand.dots = FALSE)
      mf$method <- mf $wt.method <- mf$model <- mf$x.ret <- mf$y.ret <- mf$contrasts <- mf$... <- NULL
      mf[[1L]] <- quote(stats::model.frame) # Works with complex numbers.
      mf <- eval.parent(mf)
      method <- match.arg(method)
      wt.method <- match.arg(wt.method)
      if(method == "model.frame") return(mf)
      mt <- attr(mf, "terms")
      y <- model.response(mf) # Also works with complex numbers.
      offset <- model.offset(mf) ## Note: offset is not the same thing as an intercept. offset is a known coefficient of the fit, for example the relationship between the response and one of the predictor variables.
      if(!is.null(offset)) y <- y - offset
      if(!is.null(contrasts)) warning("Contrasts are not supported for complex fits at the moment.")
      x <- zmodel.matrix(mt, mf)
      xvars <- as.character(attr(mt, "variables"))[-1L]
      if ((yvar <- attr(mt, "response")) > 0L)
        xvars <- xvars[-yvar]
      xlev <- if (length(xvars) > 0L) {
        xlev <- lapply(mf[xvars], levels)
        xlev[!sapply(xlev, is.null)]
      }
      weights <- model.weights(mf)
      if(!length(weights)) weights <- rep(1, nrow(x))
      fit <- rlm.default(x, y, weights, method = method,
                         wt.method = wt.method, ...)
      fit$terms <- mt
      ## fix up call to refer to the generic, but leave arg name as `formula'
      cl <- match.call()
      cl[[1L]] <- as.name("rlm")
      fit$call <- cl
      fit$contrasts <- attr(x, "contrasts")
      fit$xlevels <- .getXlevels(mt, mf)
      fit$na.action <- attr(mf, "na.action")
      if(model) fit$model <- mf
      if(!x.ret) fit$x <- NULL
      if(y.ret) fit$y <- y
      ## change in 7.3-52 suggested by Andr\'e Gillibert
      fit$offset <- offset
      if (!is.null(offset)) fit$fitted.values <- fit$fitted.values + offset
      fit
    }
}

#' @describeIn rlm Default S3 method
#'
#' @param x numeric or complex. A matrix, dataframe, or vector containing the explanatory / independent / predictor variables.
#' @param y numeric or complex. A vector containing the dependent / response variables, the same length as x.
#' @param ... additional arguments to be passed to rlm.default or to the psi function.
#' @param w (optional) initial down-weighting for each case
#' @param init (optional) initial values for the coefficients OR a method to find initial values OR the result of a fit with a coef component. Known methods are "ls" (the default) for an initial least-squares fit using weights w*weights, and "lts" for an unweighted least-trimmed squares fit with 200 samples.
#' @param psi the psi function is specified by this argument. It must give (possibly by name) a function g(x, ..., deriv) that for deriv=0 returns psi(x)/x and for deriv=1 returns psi'(x). Tuning constants will be passed in via ...
#' @param scale.est method of scale estimation: re-scaled MAD of the residuals (default) or Huber's proposal 2 (which can be selected by either "Huber" or "proposal 2"). Only MAD is implemented for complex variables.
#' @param k2 tuning constant used for Huber proposal 2 scale estimation.
#' @param maxit maximum number of IWLS iterations.
#' @param acc the accuracy for the stopping criterion.
#' @param test.vec the stopping criterion is based on changes in this vector.
#' @param lqs.control An optional list of control values for [lqs].
#' @param interc TRUE or FALSE, default is FALSE. Used with rlm.default when fitting complex valued data. If true, a y-intercept is calculated during fitting. Otherwise, the intercept is set to zero.
#'
#' @export
#' @examples
#' set.seed(4242)
#' n <- 8
#' slope <- complex(real = 4.23, imaginary = 2.323)
#' intercept <- complex(real = 1.4, imaginary = 1.804)
#' e <- complex(real=rnorm(n)/6, imaginary=rnorm(n)/6)
#' x <- complex(real = rnorm(n), imaginary = rnorm(n))
#' y <- slope * x + intercept + e
#' rlm(x = x, y = y, weights = rep(1,n), interc = TRUE)
rlm.default <-
  function(x, y, weights, ..., w = rep(1, nrow(x)),
           init = "ls", psi = psi.huber,
           scale.est = c("MAD", "Huber", "proposal 2"), k2 = 1.345,
           method = c("M", "MM"), wt.method = c("inv.var", "case"),
           maxit = 20, acc = 1e-4, test.vec = "resid", lqs.control=NULL, interc=FALSE)
    {
    thiscall <- match.call()
    if (is.numeric(x)) 
    {
      thiscall[[1]] <- MASS::rlm.default
      eval(thiscall, parent.frame())
    }
    else if (is.complex(x))
    {
    if (!is.matrix(x) && interc) {
      xx <- as.matrix(data.frame(rep(1, length(x)), x))
      attr(xx, "dimnames") <- list(as.character(1:length(x)), c("(intercept)", deparse(substitute(x))))
      attr(xx, "assign") <- c(0,1)
      x <- xx
    }
  irls.delta <- function(old, new)
    as.numeric(sqrt(sum(Conj(old - new)*(old - new))/max(1e-20, as.numeric(sum(Conj(old)*old)))))
    irls.rrxwr <- function(x, w, r)
    {
        w <- sqrt(w)
        max(abs((matrix(r * w, 1L, length(r)) %*% x)/
                sqrt(matrix(w, 1L, length(r)) %*% (x^2))))/abs(sqrt(sum(w * r^2))) 
        # What is the point of the max() here? As far as I can tell, the matrix multiplication would return
        # a single value...? The max function turns the 1x1 matrix into just a number. Unexpected use, but ok.
        # Oh, right; x can, and often will be, a matrix.
    }
    # wmad function used to be here. Moved up in rank to a user accessible function.
    method <- match.arg(method)
    wt.method <- match.arg(wt.method)
    nmx <- deparse(substitute(x))
    if(is.null(dim(x))) {
        x <- as.matrix(x)
        colnames(x) <- nmx
    } else x <- as.matrix(x)
    if(is.null(colnames(x)))
        colnames(x) <- paste("X", seq(ncol(x)), sep="")
    if(qr(x)$rank < ncol(x))
        stop("'x' is singular: singular fits are not implemented in 'rlm'")

    if(!(any(test.vec == c("resid", "coef", "w", "NULL"))
         || is.null(test.vec))) stop("invalid 'test.vec'")
    ## deal with weights
    xx <- x
    yy <- y
    if(!missing(weights)) {
        if(length(weights) != nrow(x))
            stop("length of 'weights' must equal number of observations")
        if(any(weights < 0)) stop("negative 'weights' value")
        if(wt.method == "inv.var") {
            fac <- sqrt(weights)
            y <- y*fac; x <- x* fac
            wt <- NULL
        } else {
            w <- w * weights
            wt <- weights
        }
    } else wt <- NULL

    if(method == "M") {
        scale.est <- match.arg(scale.est)
        if(!is.function(psi)) psi <- get(psi, mode="function")
        ## match any ... args to those of psi.
        arguments <- list(...)
        if(length(arguments)) {
            pm <- pmatch(names(arguments), names(formals(psi)), nomatch = 0L)
            if(any(pm == 0L)) warning("some of ... do not match")
            pm <- names(arguments)[pm> 0L]
            formals(psi)[pm] <- unlist(arguments[pm])
        }
        if(is.character(init)) {
            temp <- if(init == "ls") zlm.wfit(x, y, w, method="qr")
            else if(init == "lts") {
                if(is.null(lqs.control)) lqs.control <- list(nsamp=200L)
                do.call("lqs", c(list(x, y, intercept = FALSE), lqs.control))
            } else stop("'init' method is unknown")
            coef <- temp$coefficients
            resid <- temp$residuals
        } else {
            if(is.list(init)) coef <- init$coef
            else coef <- init
            resid <- drop(y - x %*% coef)
        }
    } else if(method == "MM") {
        scale.est <- "MM"
        temp <- do.call("lqs",
                        c(list(x, y, intercept = interc, method = "S",
                               k0 = 1.548), lqs.control))
        coef <- temp$coefficients
        resid <- temp$residuals
        psi <- psi.bisquare
        if(length(arguments <- list(...)))
            if(match("c", names(arguments), nomatch = 0L)) {
                c0 <- arguments$c
                if (c0 > 1.548) formals(psi)$c <- c0
                else
                    warning("'c' must be at least 1.548 and has been ignored")
            }
        scale <- temp$scale
    } else stop("'method' is unknown")

    done <- FALSE
    conv <- NULL
    n1 <- (if(is.null(wt)) nrow(x) else sum(wt)) - ncol(x)
    #theta <- 2*pnorm(k2) - 1 # pnorm() gives CDF of a normal distribution at k2 quantile. Do I need a complex normal pnorm()? Does such a thing exist? Maybe I should just ditch Huber's proposal 2...
  #gamma <- theta + k2^2 * (1 - theta) - 2 * k2 * dnorm(k2) # dnorm gives density of normal distribution at quantile k2. Shoot.
    ## At this point the residuals are weighted for inv.var and
    ## unweighted for case weights.  Only Huber handles case weights
    ## correctly.
    residdf <- data.frame(Re(resid), Im(resid)) 
    if(scale.est != "MM")
        scale <- if(is.null(wt)) {
        median(abs(resid - 0))/0.6745 # The median absolute deviation, centered on zero. Modeled on rlm from MASS, but I don't think that the median of the residuals can be counted on to always be the case...
          } else wmedian(abs(resid - 0), wt)/0.6745
    for(iiter in 1L:maxit) {
        if(!is.null(test.vec)) testpv <- get(test.vec)
        if(scale.est != "MM") {
            scale <- if(scale.est == "MAD")
                if(is.null(wt)) median(abs(resid))/0.6745 else wmedian(abs(resid), wt)/0.6745 ## wmad does not actually find the weighted MAD!!!! It just finds the weighted median!
            #else if(is.null(wt)) ## The two lines below are the Huber proposal 2 scale estimate. Why they didn't use the Huber function included in the package elsewhere is beyond me...
            #    sqrt(sum(pmin(Conj(resid)*resid, Conj(k2 * scale)*(k2 * scale)))/(n1*gamma))
            #else sqrt(sum(wt*pmin(Conj(resid)*resid, Conj(k2 * scale)*(k2 * scale)))/(n1*gamma))
            else 
              {
                warning("Only MAD scale is supported for complex variables. Continuing with scale.est = \"MAD.\"")
                scale <- if(scale.est == "MAD")
                  if(is.null(wt)) median(abs(resid))/0.6745 else wmedian(abs(resid), wt)/0.6745
              }
            if(scale == 0) {
                done <- TRUE
                break
            }
        }
        w <- psi(resid/scale)
        if(!is.null(wt)) w <- w * weights
        temp <- zlm.wfit(x, y, w, method="qr")
        coef <- temp$coefficients
        resid <- temp$residuals
        if(!is.null(test.vec)) convi <- irls.delta(testpv, get(test.vec))
        else convi <- irls.rrxwr(x, w, resid)
        conv <- c(conv, convi)
        done <- (convi <= acc)
        if(done) break
    }
    if(!done)
        warning(gettextf("'rlm' failed to converge in %d steps", maxit),
                domain = NA)
    fitted <- drop(xx %*% coef)
    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl <- match.call()
    cl[[1L]] <- as.name("rlm")
    fit <- list(coefficients = coef, residuals = yy - fitted, wresid = resid, # If inv.var weights are used, wresid are pre-weighted, otherwise they are the same as residuals.
                effects = temp$effects,
                rank = temp$rank, fitted.values = fitted,
                assign = temp$assign,  qr = temp$qr, df.residual = NA, w = w,
                s = scale, psi = psi, k2 = k2,
                weights = if(!missing(weights)) weights,
                conv = conv, converged = done, x = xx, call = cl)
    class(fit) <- c("rlm", "lm")
    fit
    }
    else print("ERROR: data must be numeric or complex.")
}

### I don't think I need a separate print.rlm function for complex variables...
# print.rlm <- function(x, ...)
# {
#     if(!is.null(cl <- x$call)) {
#         cat("Call:\n")
#         dput(cl, control=NULL)
#     }
#     if(x$converged)
#         cat("Converged in", length(x$conv), "iterations\n")
#     else cat("Ran", length(x$conv), "iterations without convergence\n")
#     coef <- x$coefficients
#     cat("\nCoefficients:\n")
#     print(coef, ...)
#     nobs <- length(x$residuals)
#     rdf <- nobs - length(coef)
#     cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
#     if(nzchar(mess <- naprint(x$na.action))) cat("  (", mess, ")\n", sep="")
#     cat("Scale estimate:", format(signif(x$s,3)), "\n")
#     invisible(x)
# }

#' Summary Method for Robust Linear Models
#' 
#' Summary method for objects of class "rlm", capable of processing complex variable fits. 
#' If the residuals in the passed object are numeric, this function just calls [MASS::summary.rlm()].
#'
#' @param object An object inheriting from the class rlm. That is, a fitted model that was generated by rlm. The fitted data can have been numeric or complex.
#' @param method Character string indicating if the IWLS weights should be used when calculating matrix cross-products. "XtX" does not include weights, "XtWX" does.
#' @param correlation Logical. Should the correlations be computed and printed?
#' @param ... Other arguments, passed to or from other methods.
#'
#' @details 
#' This is a method for the generic `summary()` function. It is based on the function of the same name from MASS, but has been modified to accept complex numbers.
#' In addition to the standard error of the coefficients, it calculates a "pseudo standard error" or "relational standard error" as the square root of the "pseudo-variance".
#' This is a complex number that quantifies the covariance between the real and imaginary parts. It can also be thought of as the amount and direction of anisotropy in the 
#' (presumed complex normal) probability distribution in the complex plane. The argument of this number gives the direction of the semi-major axis.
#'
#' @return 
#' Returns a list containing the following elements. If the list is printed by the function (`print.summary.rlm` or invoked as a method by `summary`), a null value is returned instead.
#' \item{`correlation`}{A numeric matrix. The computed correlation coefficient matrix for the coefficients in the model.}
#' \item{`pseudocorrelation`}{A complex matrix. The computed pseudo-correlation coefficient matrix for the coefficients in the model.}
#' \item{`cov.unscaled`}{The unscaled covariance matrix; i.e, a numeric matrix such that multiplying it by an estimate of the error variance produces an estimated covariance matrix for the coefficients.}
#' \item{`pcov.unscaled`}{The unscaled pseudo-covariance matrix; i.e, a complex matrix such that multiplying it by an estimate of the error pseudo-variance produces an estimated pseudo-covariance matrix for the coefficients.}
#' \item{`sigma`}{Numeric. The scale estimate returned by `rlm`, which was used to scale the residuals before passing them to the `psi` weight function.}
#' \item{`stddev`}{Numeric. A scale estimate used for the standard errors. Calculated from the working residuals, the `psi` weight function, and the derivative of the influence function (`psi.method(deriv = 1)`).}
#' \item{`pstddev`}{Complex. A scale estimate used for the 'pseudo' standard errors. Calculated from the working residuals, the `psi` weight function, and the derivative of the influence function (`psi.method(deriv = 1)`). See details above.}
#' \item{`df`}{The number of degrees of freedom for the model and for residuals.}
#' \item{`coefficients`}{A 4 column matrix that contains the model coefficients, their standard errors, their pseudo standard errors (see details above), and their t statistics.}
#' \item{`terms`}{The terms object used in fitting this model.}
#' @export
#'
#'@references
#' W. N. Venables and B. D. Ripley, Modern Applied Statistics with S, 4th ed (Springer, New York, 2002).
#' P. J. Huber and E. Ronchetti, Robust Statistics, 2nd ed (Wiley, Hoboken, N.J, 2009).
#'
#' @examples
#' set.seed(4242)
#' n <- 8
#' slope <- complex(real = 4.23, imaginary = 2.323)
#' intercept <- complex(real = 1.4, imaginary = 1.804)
#' e <- complex(real=rnorm(n)/6, imaginary=rnorm(n)/6)
#' x <- complex(real = rnorm(n), imaginary = rnorm(n))
#' y <- slope * x + intercept + e
#' robfit <- rlm(x = x, y = y, weights = rep(1,n), interc = TRUE)
#' summary(robfit)
summary.rlm <- function(object, method = c("XtX", "XtWX"),
                        correlation = FALSE, ...)
{
    cll <- match.call()
    if (is.numeric(object$residuals))
    {
      cll[[1]] <- MASS::summary.rlm
      eval(cll, parent.frame())
    }
    else
    {
      method <- match.arg(method)
      s <- object$s
      coef <- object$coefficients
      ptotal <- length(coef)
      wresid <- object$wresid
      res <- object$residuals
      n <- length(wresid)
      if(any(na <- is.na(coef))) coef <- coef[!na]
      cnames <- names(coef)
      p <- length(coef)
      rinv <- diag(p)
      dimnames(rinv) <- list(cnames, cnames)
      wts <- if(length(object$weights)) object$weights else rep(1, n)
      if(length(object$call$wt.method) && object$call$wt.method == "case") {
        rdf <- sum(wts) - p ## 'residual degrees of freedom'
        w <- object$psi(wresid/s)
        S <- sum(wts * (wresid*w)*Conj(wresid*w))/rdf # The estimated variance of the residuals is analogous to the area of a circle around the estimated parameters, so it is a real number.
        pS <- sum(wts * (wresid*w)*(wresid*w))/rdf # The estimated pseudo-variance is a complex number that describes the anisotropy of the distribution or set.
        # Distributions that are less circularly symmetric and more bilaterally symmetric have higher pseudovariance.
        # The direction of the pseudovariance indicates the orientation of the anisotropy; bilateral symmetry axis angle from 
        # the +real axis = pseudovariance angle divided by two. In other words, pseudovariance is the covariance between real and imaginary.
        psiprime <- object$psi(wresid/s, deriv=1) # This region of code occupies the position that finding the sum of squared deviations from the mean of the fitted values does in summary.lm().
        m1 <- sum(wts*psiprime) # What are these for? psiprime will depend on length (size) of residual, and direction. It depends on what psi function was chosen. Huber will have psiprime = 1 Exp[i (phi - pi)] for small residuals and 0 for large ones. A kind of weighted sample 1st moment of psiprime.
        m2 <- sum(wts*as.numeric(Conj(psiprime)*psiprime))
        pm2 <- sum(wts*psiprime*psiprime) # 'pseudo m2' a complex quantity.
        nn <- sum(wts) # "sample size" is sum of the weights, makes sense for case weights I suppose.
        mn <- m1/nn # 'mn' is "mean" -> This is the mean psiprime, the average derivative of the influence function for the fit in question.
        kappa <- 1 + p*(m2 - as.numeric(Conj(m1)*m1)/nn)/(nn-1)/(nn*as.numeric(Conj(mn)*mn)) #Check this. It resembles the uncertainty of a quantum operator, sigma_x = sqrt(<x^2> - <x>^2). # It almost looks like it's missing a parentheses...
        stddev <- sqrt(S)*(kappa/abs(mn))
        pkappa <- 1 + p*(pm2 - m1^2/nn)/(nn-1)/(nn*mn^2) # 'pseudo kappa', a complex thing
        pstddev <- sqrt(pS)*(pkappa/mn) # The 'pseudo standard deviation', analogous with the pseudo-variance. Might be useful, might be meaningless.
      } else {
        res <- res * sqrt(wts)
        rdf <- n - p ## 'residual degrees of freedom'
        w <- object$psi(wresid/s)
        S <- sum(as.numeric(wresid*w*Conj(wresid*w)))/rdf # The estimated variance is analogous to the area of a circle around the estimated parameters, so it is a real number.
        pS <- sum((wresid*w)*(wresid*w))/rdf # See above definition of pS.
        psiprime <- object$psi(wresid/s, deriv=1) # The derivative of the influence function.
        mn <- mean(psiprime)
        pvarpsi <- sum((psiprime - mn)^2)/(n-1)
        kappa <- 1 + p*var(psiprime)/(n*as.numeric(Conj(mn)*mn)) # This has something to do with propagation of uncertainty, I think.
        #print(var(psiprime))
        pkappa <- 1 + p*pvarpsi/(n*mn^2)
        stddev <- sqrt(S)*(kappa/abs(mn)) ## Would it be useful to do something similar with pseudo-variance? Probably.
        pstddev <- sqrt(pS)*(pkappa/mn) # The 'pseudo standard deviation', analogous with the pseudo-variance. Might be useful, might be meaningless.
      }
      X <- if(length(object$weights)) object$x * sqrt(object$weights) else object$x # The model (design) matrix, or the model matrix times the sqrt of the weights.
      if(method == "XtWX")  { # (X transposed) [weight matrix] (X) #wts are not the IWLS weights, they are the ones given by the user. w are the IWLRS weights.
        mn <- sum(wts*w)/sum(wts)
        X <- X * sqrt(w/mn)
      }
      R <- qr(X)$qr
      R <- R[1L:p, 1L:p, drop = FALSE] # Trim the bottom of R, making it a square p by p matrix.
      R[lower.tri(R)] <- 0 # Remove the lower triangular matrix, the Q from the QR decomposition.
      rinv <- solve(R, rinv) # This is efficient, we only need the diagonal matrix diag(p) to get rinv, so just set rinv <- diag(p) ahead of time. Now rinv is a different p by p matrix.
      dimnames(rinv) <- list(cnames, cnames)
      rowlen <- (abs(Conj(rinv)*rinv) %*% rep(1, p))^0.5 # Produces a real vector of length p. Elements are the sum of the rows of (X^T X)^-1, each square rooted. 
      prowlen <- (rinv^2 %*% rep(1, p))^0.5 # Produces a complex vector of length p
      names(rowlen) <- cnames # cnames are the names of the coefficients.
      if(correlation) {
        correl <- rinv * array(1/rowlen, c(p, p))
        correl <- correl %*% Conj(t(correl))
        pcorrel <- correl %*% t(correl)
      } else {
        correl <- NULL
        pcorrel <- NULL
      }
      coef <- array(coef, c(p, 4L)) # Make an array with 4 columns and p rows. put the coefficients into the first column.
      dimnames(coef) <- list(cnames, c("Value", "Std. Error", "Pseudo Std. Error", "t value"))
      #print(rinv)
      coef[, 2] <- rowlen %o% stddev # Fill the 2nd column of the coef array. These should be real numbers. Isn't stddev a single number? What is the point of the outer product? It transposes the vector into a column while multiplying all elements by stddev.
      coef[, 3] <- rowlen %o% pstddev # The 'pseudo - standard error', these things need better names...
      coef[, 4] <- coef[, 1]/coef[, 2]
      object <- object[c("call", "na.action")]
      object$residuals <- res
      object$coefficients <- coef
      object$sigma <- s
      object$stddev <- stddev
      object$pstdev <- pstddev
      object$df <- c(p, rdf, ptotal)
      object$r.squared <- NA
      object$cov.unscaled <- rinv %*% Conj(t(rinv))
      object$pcov.unscaled <- rinv %*% t(rinv) # pseudo covariance
      object$correlation <- correl
      object$pseudocorrelation <- pcorrel
      object$terms <- NA
      class(object) <- "summary.rlm"
      object
    }
}

#' @describeIn summary.rlm Print the summary of a (possibly complex) robust linear fit.
#' 
#' @param x a rlm object or an rlm summary object. Used for `print.summary.rlm` 
#' @param digits the number of digits to include in the printed report, default is three. Used for `print.summary.rlm`
#' 
#' @note For complex fits the quantiles reported by this function are based on sorting the real parts of the residuals. They should not be considered reliable..
#' 
#' @export
print.summary.rlm <- function(x, digits = max(3, .Options$digits - 3), ...)
{
  cll <- match.call()
  if (is.numeric(x$residuals))
  {
    cll[[1]] <- MASS::print.summary.rlm
    eval(cll, parent.frame())
  }
  else
  {
    cll <- match.call()
    cat("\nCall: ")
    dput(x$call, control=NULL)
    resid <- x$residuals
    df <- x$df
    rdf <- df[2L]
    print(rdf)
    cat(if(!is.null(x$weights) && diff(range(x$weights))) "Weighted ",
        "Residuals:\n", sep="")
    if(rdf > 5L) {
        if(length(dim(resid)) == 2L) {
            rq <- apply(Conj(t(resid)), 1L, quantile)
            dimnames(rq) <- list(c("Min", "1Q", "Median", "3Q", "Max"),
                                 colnames(resid))
        } else {
            rq <- quantile(resid)
            names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
        }
        print(rq, digits = digits, ...)
    } else if(rdf > 0L) {
        print(resid, digits = digits, ...)
    }
    print("NOTE: Quantiles generated by sorting real components only. Do not rely upon them.")
    if(nsingular <- df[3L] - df[1L])
        cat("\nCoefficients: (", nsingular,
            " not defined because of singularities)\n", sep = "")
    else cat("\nCoefficients:\n")
    print(format(round(x$coefficients, digits = digits)), quote = FALSE, ...)
    cat("\nResidual standard error:", format(signif(x$sigma, digits)),
        "on", rdf, "degrees of freedom\n")
    if(nzchar(mess <- naprint(x$na.action))) cat("  (", mess, ")\n", sep="")
    if(!is.null(correl <- x$correlation)) {
      p <- dim(correl)[2L]
        if(p > 1L) {
            cat("\nCorrelation of Coefficients:\n")
            ll <- lower.tri(correl)
            correl[ll] <- format(round(correl[ll], digits))
            correl[!ll] <- ""
            print(correl[-1L, -p, drop = FALSE], quote = FALSE, digits = digits, ...)
        }
    }
    invisible(x)
  }
}

### NOTE: These psi functions are actually weight functions, weight(u) = abs( influence(u) / u).
### Traditionally the influence functions are given the symbol psi.
### The derivative returned by these functions, when given deriv=1, is the (absolute value of the) derivative of the influence function,
### not the weight function.
### NOTE: The derivative of the complex influence function will be a complex number, and the direction does have meaning.
### Derivative output modified to incorporate the complex nature of the influence function.
#' Weighting functions for robust fitting
#' 
#' Weighting functions used in `rlm` for iteratively reweighted least squares. 
#' Based on the same functions from MASS, modified to accept complex variables.
#' While named 'psi', these are actually weight functions, weight(u) = abs( influence(u) / u).
#'
#' @param u Numeric or complex. When used in M-estimation, a residual.
#' @param k Numeric. A scaling constant for `psi.huber`. Default value is 1.345
#' @param a Numeric. A scaling constant for `psi.hampel`. Default value is 2
#' @param b Numeric. A scaling constant for `psi.hampel`. Default value is 4
#' @param c Numeric. A scaling constant for 'psi.hampel` or `psi.bisquare`. Default is 8 for the former and 4.685 for the later.
#' @param deriv Numeric. If `0`, return the weight at `u`. If `1` return the first derivative of the influence function at `u`.
#' 
#' @details 
#' Three weight functions for iteratively (re)weighted least squares. When used in M-estimation, `psi.huber` has
#' a unique solution. `psi.hampel` and `psi.bisquare`, which are redescending M-estimators, do not have a unique solution.
#' These are capable of completely rejecting outliers, but need a good starting point to avoid falling into unhelpful local minima.
#' If `deriv` is set to `1`, the functions return the value of the first derivative of the influence function at `u`.
#' Note that they do not return the derivative of the weight function, as might be expected.
#'
#' @return A numeric or complex that is either the value of the weight function at `u` (numeric) or the first derivative of the influence function at `u` (can be complex).
#' @export
#'
#' @examples
#' set.seed(4242)
#' z <- complex(real = rnorm(3), imaginary = rnorm(3))
#' psi.huber(z)
#' psi.hampel(z)
#' psi.bisquare(z)
#' psi.huber(z, deriv=1)
#' psi.hampel(z, deriv=1)
psi.huber <- function(u, k = 1.345, deriv=0)
{
    if(!deriv) return(pmin(1, k / abs(u)))
    ifelse(abs(u) <= k, complex(modulus = 1, argument = Arg(u) - pi), 0)
}

#' @describeIn psi.huber The weight function of the hampel objective function.
#' 
#' @export
psi.hampel <- function(u, a = 2, b = 4, c = 8, deriv=0)
{
    U <- pmin(abs(u) + 1e-50, c)
    if(!deriv) return(as.vector(ifelse(U <= a, U, ifelse(U <= b, a, a * (c - U) / (c - b) )) / U))
    ifelse(abs(u) <= c, ifelse(U <= a, complex(modulus = 1, argument = Arg(u) - pi), ifelse(U <= b, 0, complex(modulus = a/(c-b), argument = Arg(u)))), 0)
}

#' @describeIn psi.huber The weight function of Tukey's bisquare objective function.
#' 
#' @export
psi.bisquare <- function(u, c = 4.685, deriv=0)
{
  cll <- match.call()
  if (!is.complex(u)) {
    #warning("dog")
    cll[[1]] <- MASS::psi.bisquare
    eval(cll, parent.frame())
  }
  else
  {
    if(!deriv) return((1  - pmin(1, abs(u/c))^2)^2)
    t <- (u/c)
    #warning(t)
    #warning(abs(t))
    ifelse(Mod(t) < 1, complex(imaginary = -1) * (-1 + Conj(t)*t)*(-1 + 5*Conj(t)*t) * complex(modulus = 1, argument = Arg(t))^2, 0)
    #return('cat')
  }
}

# se.contrast.rlm <- # I don't understand this code. Hopefully it works fine with complex numbers as is.
#     function(object, contrast.obj,
#              coef = contr.helmert(ncol(contrast))[, 1L],
#              data = NULL, ...)
# {
#     contrast.weight.aov <- function(object, contrast)
#     {
#         asgn <- object$assign[object$qr$pivot[1L:object$rank]]
#         uasgn <- unique(asgn)
#         nterms <- length(uasgn)
#         nmeffect <- c("(Intercept)",
#                       attr(object$terms, "term.labels"))[1L + uasgn]
#         effects <- as.matrix(qr.qty(object$qr, contrast))
#         res <- matrix(0, nrow = nterms, ncol = ncol(effects),
#                       dimnames = list(nmeffect, colnames(contrast)))
#         for(i in seq(nterms)) {
#             select <- (asgn == uasgn[i])
#             res[i,] <- colSums(effects[seq_along(asgn)[select], , drop = FALSE]^2)
#         }
#         res
#     }
#     if(is.null(data)) contrast.obj <- eval(contrast.obj)
#     else contrast.obj <- eval(substitute(contrast.obj), data, parent.frame())
#     if(!is.matrix(contrast.obj)) { # so a list
#         if(!missing(coef)) {
#             if(sum(coef) != 0)
#                 stop("'coef' must define a contrast, i.e., sum to 0")
#             if(length(coef) != length(contrast.obj))
#                 stop("'coef' must have same length as 'contrast.obj'")
#         }
#         contrast <-
#             sapply(contrast.obj, function(x)
#                {
#                    if(!is.logical(x))
#                        stop(gettextf("each element of '%s' must be logical",
#                                      substitute(contrasts.list)),
#                             domain = NA)
#                    x/sum(x)
#                })
#         if(!length(contrast) || all(is.na(contrast)))
#             stop("the contrast defined is empty (has no TRUE elements)")
#         contrast <- contrast %*% coef
#     } else {
#         contrast <- contrast.obj
#         if(any(abs(colSums(contrast)) > 1e-8))
#             stop("columns of 'contrast.obj' must define a contrast (sum to zero)")
#         if(!length(colnames(contrast)))
#             colnames(contrast) <- paste("Contrast", seq(ncol(contrast)))
#     }
#     weights <- contrast.weight.aov(object, contrast)
#     summary(object)$stddev *
#         if(!is.matrix(contrast.obj)) sqrt(sum(weights)) else sqrt(colSums(weights))
# }

#' @describeIn rlm Predict new data based on the model in object. Invokes `predict.lm`.
#' @param object a `rlm` object; a fit from which you wish to predict new data.
#' @param newdata new predictor data from which to predict response data. Default is NULL.
#' @param scale this seems to be ignored. Default is NULL.
#' @param ... further arguments to be passed to NextMethod().
#' 
#' @export
predict.rlm <- function (object, newdata = NULL, scale = NULL, ...)
{
    ## problems with using predict.lm are the scale and
    ## the QR decomp which has been done on down-weighted values.
    ## Only works if explicit weights are given during the call that produced object..?
    object$qr <- qr(sqrt(object$weights) * object$x)
    NextMethod(object, scale = object$s, ...) # So this just calls predict.lm on the object.
}

#' Calculate Variance-Covariance Matrix and Pseudo Variance-Covariance Matrix for a Complex Fitted Model Object
#'
#' A version of [stats::vcov] that is compatible with complex linear models. In addition to the variance-covariance matrix,
#' the pseudo variance-covariance matrix, which is a measure of the covariance between real and imaginary components, is returned as well.
#' Can also return the "big covariance" matrix, which combines the information of the covariance matrix and the pseudo-covariance matrix, as described in (van den Bos 1995).
#' While not as compact as two seperate smaller matrices, the big covariance matrix simplifies calculations such as the [mahalanobis] distance.
#'
#' @param object a fitted model object, typically. Sometimes also a summary() object of such a fitted model.
#' @param ... Additional parameters, not currently used for anything.
#' @param merge logical. Should the covariance matrix and pseudo-covariance / relational matrix be merged into one matrix of twice the dimensions? Default is TRUE.
#'
#' @return
#' If `merge` is false, a list containing both the numeric variance-covariance matrix, and the complex pseudo variance-covariance matrix.
#' If `merge` is true, a large matrix (both dimensions twice the number of coefficients) containing both the variance-covariance matrix and the pseudo variance-covariance matrix, merged together.
#' @export
#' 
#' @references A. van den Bos, The Multivariate Complex Normal Distribution-a Generalization, IEEE Trans. Inform. Theory 41, 537 (1995).
#'
#' @examples
#' set.seed(4242)
#' n <- 8
#' slope <- complex(real = 4.23, imaginary = 2.323)
#' intercept <- complex(real = 1.4, imaginary = 1.804)
#' e <- complex(real=rnorm(n)/6, imaginary=rnorm(n)/6)
#' x <- complex(real = rnorm(n), imaginary = rnorm(n))
#' y <- slope * x + intercept + e
#' robfit <- rlm(x = x, y = y, weights = rep(1,n), interc = TRUE)
#' vcov(robfit)
vcov.rlm <- function (object, merge = TRUE, ...)
{
  cll <- match.call()
  if (is.numeric(object$residuals))
  {
    cll[[1]] <- MASS::vcov.rlm
    eval(cll, parent.frame())
  }
  else
  {
    so <- summary(object, corr = FALSE)
    varcovar <- so$stddev^2 * so$cov.unscaled
    pseudovarcovar <- so$pstddev^2 * so$pcov.unscaled
    if (merge == TRUE)
    {
      #bigcovar <- diag(rep(diag(varcovar), each = 2)) # Start by making a square diagonal matrix with two adjacent diagonal elements for each variance.
      n <- attributes(varcovar)$dim[1] # The size of the square covariance matrix.
      bigcovar <- matrix(0, ncol = 2*n, nrow = 2*n) #Start by making an empty matrix with dimensions twice those of the small covariance matrix.
      bigcovar[,seq(1, 2*n, 2)] <- matrix(as.vector(rbind(as.vector(varcovar), as.vector(Conj(pseudovarcovar)))), nrow = 2*n) # Fill the odd indexed columns of bigcovar with the entries from varcovar interleaved with the entries from pseudovarcovar conjugated.
      bigcovar[,seq(2, 2*n, 2)] <- matrix(as.vector(rbind(as.vector(pseudovarcovar), as.vector(varcovar))), nrow = 2*n) # Fill the even indexed columns of bigcovar with the entries from varcovar interleaved with the entries from pseudovarcovar, this time the later first.
      return(bigcovar)
    }
    else return(list(varcovar = as.numeric(varcovar), pseudovarcovar = pseudovarcovar))
  }  
}
