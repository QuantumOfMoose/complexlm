
####
### A function that somewhat replicates model.matrix(), but accepts complex valued data. It will probably be slower and less effecient, but mostly functional.
### It cannot handle algebraic expressions in formula.
### terms is the output of terms(formula)
####
#' Generate a Model Matrix (Design Matrix) Using Complex Variables
#' 
#' A function that somewhat replicates model.matrix(), but accepts complex valued data. It will probably be slower and less effecient, but mostly functional.
#' It cannot handle algebraic expressions in formula.
#'
#' @param trms A "terms" object as generated by [stats::terms]. 
#' @param data A data frame containing the data referenced by the symbolic formula / model in `trms`.
#' 
#'
#' @return A model matrix, AKA a design matrix for a regression, containing the relavent information from `data`.
#' @export
#'
#' @examples
#' set.seed(4242)
#' slope = complex(real = 4.23, imaginary = 2.323)
#' intercept = complex(real = 1.4, imaginary = 1.804)
#' testframe <- expand.grid(-3:3,-3:3)
#' Xtest <- complex(real = testframe[[1]], imaginary = testframe[[2]])
#' testframe <- data.frame(Xtest = Xtest, Ytest = Xtest * slope + intercept + complex(real = rnorm(length(Xtest)), imaginary = rnorm(length(Xtest))))
#' testterms <- terms(Ytest ~ Xtest)
#' zmodel.matrix(testterms, testframe)
zmodel.matrix <- function(trms, data)
{
  respname <- as.character(attr(trms, "variables")[[attr(trms, "response") + 1]])
  prednames <- attr(trms, "term.labels")
  modelframe <- data[prednames]
  if (attr(trms, "intercept") == 1) modelmatrix <- as.matrix(data.frame("(intercept)" = rep(1,length(modelframe[,1])), modelframe))
  else  modelmatrix <- as.matrix(modelframe)
  if (attr(trms, "intercept") == 1) attr(modelmatrix, "assign") <- 0:length(prednames)
  else attr(modelmatrix, "assign") <- 1:length(prednames)
  attr(modelmatrix, "dimnames") <- list(as.character(1:length(modelframe[,1])), c("(intercept)", prednames)) # Fix the dimnames to match those on standard model matrices.
  return(modelmatrix)
}

#### This function will be used instead of .Call(C_Cdqlrs, x * wts, y * wts, tol, FALSE) if x and/or y are complex.
#' An alternative to .Call(C_Cdqlrs, x * wts, y * wts, tol, FALSE)) that is compatible with complex variables
#' 
#' This serves as a wrapper for qr, replicating the behavior and output of the C++ function `C_Cdqlrs`. It is used in `zlm.wfit`,
#' and is unlikely to be needed by end users.
#'
#' @param x a complex matrix (will also accept numeric, but in that case you might as well use `C_Cdqlrs`) whose QR decomposition is to be computed.
#' @param y a complex vector or matrix of right-hand side quantities.
#' @param tol the tolerance for detecting linear dependencies in the columns of x. Not used for complex `x`.
#' @param chk not used. Included to better immitate `C_Cdqlrs`.
#'
#' @return A list that includes the qr decomposition, its coeffcionts, residuals, effects, rank, pivot information, qraux vector,
#' tolerance, and whether or not it was pivoted. This is the same output as `C_Cdqlrs`.
#'
#' @examples
Complexdqlrs <- function (x, y, tol, chk) {
  thisqr <- qr(x, tol = tol)
  coefficients <- qr.coef(thisqr, y)
  resids <- y - x %*% coefficients
  effects <- qr.qty(thisqr, y)
  xrank <- thisqr$rank
  pivot <- thisqr$pivot
  qraux <- thisqr$qraux
  pivoted <- any(pivot > (1:length(pivot)) + 1) * 1
  ans = list(qr = thisqr, coefficients = coefficients, residuals = resids, effects = effects, rank = xrank, 
             pivot = pivot, qraux = qraux, tol = tol, pivoted = pivoted)
  return(ans)
}


####
### An adaptation of lm that is compatible with complex variables. If the response is not complex, it calls the standard stats::lm
### Note: It is not capable of dealing with contrasts in the complex case. May not understand offsets either. It also can't handle algebraic expressions in formula.
### model.frame needs to be changed to allow complex variables in order to enable these features.
####
#' Linear Model Fitting for Complex or Numeric Variables
#' 
#' An adaptation of lm that is compatible with complex variables. If the response is not complex, it calls the standard `stats::lm`
#' Note: It is not capable of dealing with contrasts in the complex case. May not understand offsets either. 
#' It also can't handle algebraic expressions in formula.
#' model.frame needs to be changed to allow complex variables in order to enable these features.
#'
#' @inherit stats::lm description details params return
#'
#' @export
#'
#' @examples
#' set.seed(4242)
#' n = 8
#' slope = complex(real = 4.23, imaginary = 2.323)
#' intercept = complex(real = 1.4, imaginary = 1.804)
#' testframe <- data.frame(x = x <- complex(real = rnorm(n), imaginary = rnorm(n)), y = slope * x + intercept)
#' lm(y ~ x, data = testframe, weights = rep(1,n))
lm <- function (formula, data, subset, weights, na.action,
                method = "qr", model = TRUE, x = FALSE, y = FALSE,
                qr = TRUE, singular.ok = TRUE, contrasts = NULL,
                offset, ...)
{
  trms <- terms(formula)
  respname <- as.character(attr(trms, "variables")[[attr(trms, "response") + 1]])
  cl <- match.call()
  if (is.complex(data[,respname]) == FALSE)
  {
    cl[[1]] <- stats::lm
    eval(cl, parent.frame())
  }
  else
  {
    ret.x <- x
    ret.y <- y
    
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    ## need stats:: for non-standard evaluation
    mf[[1L]] <- quote(stats::model.frame) # It works with complex numbers :)
    mf <- eval(mf, parent.frame())
    if (method == "model.frame")
      return(mf)
    else if (method != "qr")
      warning(gettextf("method = '%s' is not supported. Using 'qr'", method),
              domain = NA)
    mt <- attr(mf, "terms") # allow model.frame to update it
    y <- model.response(mf) # Huh, this does work with complex numbers.
    ## avoid any problems with 1D or nx1 arrays by as.vector.
    w <- as.vector(model.weights(mf))
    if(!is.null(w) && !is.numeric(w))
      stop("'weights' must be a numeric vector")
    offset <- model.offset(mf) # Uncertain if this works with complex variables.
    mlm <- is.matrix(y)
    ny <- if(mlm) nrow(y) else length(y)
    if(!is.null(offset)) {
      if(!mlm) offset <- as.vector(offset)
      if(NROW(offset) != ny)
        stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                      NROW(offset), ny), domain = NA)
    }
    
    if (is.empty.model(mt)) {
      x <- NULL
      z <- list(coefficients = if(mlm) matrix(NA_real_, 0, ncol(y))
                else numeric(),
                residuals = y,
                fitted.values = 0 * y, weights = w, rank = 0L,
                df.residual = if(!is.null(w)) sum(w != 0) else ny)
      if(!is.null(offset)) {
        z$fitted.values <- offset
        z$residuals <- y - offset
      }
    }
    else {
      if (is.null(contrasts) == FALSE) warning("Contrasts are not supported for complex fits.")
      x <- zmodel.matrix(mt, mf)
      z <- if(is.null(w)) lm.fit(x, y, offset = offset,
                                 singular.ok=singular.ok, ...)
      else lm.wfit(x, y, w, offset = offset, singular.ok=singular.ok, ...)
    }
    class(z) <- c(if(mlm) "mlm", "lm")
    z$na.action <- attr(mf, "na.action")
    z$offset <- offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    if (model)
      z$model <- mf
    if (ret.x)
      z$x <- x
    if (ret.y)
      z$y <- y
    if (!qr) z$qr <- NULL
    z
  }
}

####
### Wrapper for lm.fit() If data is numeric, use lm.fit() from stats. If it is complex, use zlm.wfit().
####
#' Complex Variable Compatible Wrappers for [stats::lm.fit()] and [stats::lm.wfit()]
#' 
#' This function is just an if statement.
#' If the design matrix `x` is complex, [zlm.wfit()] is called.
#' Otherwise [stats::lm.fit()] or [stats::lm.wfit()] is called.
#' These functions are unlikely to be needed by end users, as they are called by [lm()].
#' 
#' @inherit stats::lm.wfit params return
#'
#' @export
#' 
#' @examples
#' set.seed(4242)
#' n = 6
#' p = 2
#' slope = complex(real = 4.23, imaginary = 2.323)
#' slope2 = complex(real = 2.1, imaginary = -3.9)
#' intercept = complex(real = 1.4, imaginary = 1.804)
#' designmat <- matrix(c(complex(real = rnorm(n * p), imaginary = rnorm(n * p)), rep(1, n)), n, p + 1)
#' y = designmat %*% c(slope, slope2, intercept) + complex(real = rnorm(n), imaginary = rnorm(n))
#' lm.fit(designmat, y)
lm.fit <- function(x, y, offset = NULL, method = "qr", tol = 1e-7,
       singular.ok = TRUE, ...)
{
  cll <- match.call()
  if (is.complex(x)) cll[[1]] <- zlm.wfit
  else cll[[1]] <- stats::lm.fit
  eval(cll, parent.frame())
}

####
### Wrapper for lm.wfit(). If data is numeric, use lm.wfit() from stats. If it is complex, use zlm.wfit().
####
#' describeIn lm.fit wrapper for weighted linear fitting function.
lm.wfit <- function(x, y, w, offset = NULL, method = "qr", tol = 1e-7,
        singular.ok = TRUE, ...)
{
  cll <- match.call()
  if (is.complex(x)) cll[[1]] <- zlm.wfit
  else cll[[1]] <- stats::lm.wfit
  eval(cll, parent.frame())
}

#' Least-Squares Linear Fitting for Complex Variables
#' 
#' The function eventually called by [lm()], [lm.fit()], and/or [lm.wfit()] if fed complex data.
#' Performs ordinary (least-squares) linear fitting on complex variable data.
#' Like [stats::lm.wfit()], which it is based off of, it uses qr decomposition
#' for the matrix algebra. Unlike `stats::lm.wfit()' it also handles un-weighted
#' regression by setting the weights to 1 by default.
#'
#' @param x a complex design matrix, `n` rows by `p` columns.
#' @param y a vector of observations/responses of length `n`, or a matrix with `n` rows.
#' @param w a vector of weights to be used in the fitting process. The sum of `w * r^2` is minimized, with `r` being the residuals. By default, `w` is a vector of length `n` with every element equal to 1, making this an unweighted fit.
#' @param offset optional. A complex vector of length n that will be subtracted from y prior to fitting.
#' @param method optional. a string that can be used to choose any method you would like. As long as it is "qr".
#' @param tol tolerance for the [qr] decomposition. Default is 1e-7.
#' @param singular.ok logical. If false, a singular model is an error.
#' @param ... currently disregarded.
#'
#' @inherit stats::lm.wfit return
#' @export
#'
#' @inherit lm.wfit examples
zlm.wfit <- function (x, y, w = rep(1L, ifelse(is.vector(x), length(x), nrow(x))), offset = NULL, method = "qr", tol = 1e-07, 
          singular.ok = TRUE, ...) 
{
  if (is.null(n <- nrow(x))) 
    stop("'x' must be a matrix")
  if (n == 0) 
    stop("0 (non-NA) cases")
  ny <- NCOL(y)
  if (is.matrix(y) && ny == 1L) 
    y <- drop(y)
  if (!is.null(offset)) 
    y <- y - offset
  if (NROW(y) != n | length(w) != n) 
    stop("incompatible dimensions")
  if (any(w < 0 | is.na(w))) 
    stop("missing or negative weights not allowed")
  if (method != "qr") 
    warning(gettextf("method = '%s' is not supported. Using 'qr'", 
                     method), domain = NA)
  chkDots(...)
  x.asgn <- attr(x, "assign")
  zero.weights <- any(w == 0)
  if (zero.weights) {
    save.r <- y
    save.f <- y
    save.w <- w
    ok <- w != 0
    nok <- !ok
    w <- w[ok]
    x0 <- x[!ok, , drop = FALSE]
    x <- x[ok, , drop = FALSE]
    n <- nrow(x)
    y0 <- if (ny > 1L) 
      y[!ok, , drop = FALSE]
    else y[!ok]
    y <- if (ny > 1L) 
      y[ok, , drop = FALSE]
    else y[ok]
  }
  p <- ncol(x)
  if (p == 0) {
    return(list(coefficients = numeric(), residuals = y, 
                fitted.values = 0 * y, weights = w, rank = 0L, df.residual = length(y)))
  }
  if (n == 0) {
    return(list(coefficients = rep(NA_real_, p), residuals = y, 
                fitted.values = 0 * y, weights = w, rank = 0L, df.residual = 0L))
  }
  wts <- sqrt(w)
  print(typeof(wts))
  print(typeof(x))
  z <- Complexdqlrs(x * wts, y * wts, tol, FALSE)
  if (!singular.ok && z$rank < p) 
    stop("singular fit encountered")
  coef <- z$coefficients
  pivot <- z$pivot
  r1 <- seq_len(z$rank)
  dn <- colnames(x)
  if (is.null(dn)) 
    dn <- paste0("x", 1L:p)
  nmeffects <- c(dn[pivot[r1]], rep.int("", n - z$rank))
  r2 <- if (z$rank < p) 
    (z$rank + 1L):p
  else integer()
  if (is.matrix(y)) {
    coef[r2, ] <- NA
    if (z$pivoted) 
      coef[pivot, ] <- coef
    dimnames(coef) <- list(dn, colnames(y))
    dimnames(z$effects) <- list(nmeffects, colnames(y))
  }
  else {
    coef[r2] <- NA
    if (z$pivoted) 
      coef[pivot] <- coef
    names(coef) <- dn
    names(z$effects) <- nmeffects
  }
  z$coefficients <- coef
  z$residuals <- z$residuals/wts
  z$fitted.values <- y - z$residuals
  z$weights <- w
  if (zero.weights) {
    coef[is.na(coef)] <- 0
    f0 <- x0 %*% coef
    if (ny > 1) {
      save.r[ok, ] <- z$residuals
      save.r[nok, ] <- y0 - f0
      save.f[ok, ] <- z$fitted.values
      save.f[nok, ] <- f0
    }
    else {
      save.r[ok] <- z$residuals
      save.r[nok] <- y0 - f0
      save.f[ok] <- z$fitted.values
      save.f[nok] <- f0
    }
    z$residuals <- save.r
    z$fitted.values <- save.f
    z$weights <- save.w
  }
  if (!is.null(offset)) 
    z$fitted.values <- z$fitted.values + offset
  if (z$pivoted) 
    colnames(z$qr) <- colnames(x)[z$pivot]
  qr <- z[c("qr", "qraux", "pivot", "tol", "rank")]
  c(z[c("coefficients", "residuals", "fitted.values", "effects", 
        "weights", "rank")], list(assign = x.asgn, qr = structure(qr, 
        class = "qr"), df.residual = n - z$rank))
}

####
### A wrapper for summary.lm(). If the residuals are numeric, call summary.lm() from stats.
#' Summarize Complex Linear Model Fits.
#' 
#' This extends summary.lm() to handle linear fits of complex variables.
#' If the residuals of the fit object are numeric, `stats::summary.lm()` is called.
#'
#' @param object An object of class 'lm'. Presumably returned by [lm]. May contain complex variables.
#' @param correlation Logical. If TRUE, the correlation matrix of the estimated parameters is returned and printed.
#' @param symbolic.cor Logical. If TRUE, print the correlations in a symbolic form (see [stats::symnum]) rather than as numbers. (This may not work.)
#' @param ... Further argements passed to or from other methods.
#' 
#' @details 
#' See [stats::summary.lm] for more information.
#' In addition to the output information returned by `stats::summary.lm`, this complex variable compatable version also returns
#' "pseudo standard error" or "relational standard error" which is the square root of the "pseudo-variance".
#' This is a complex number that quantifies the covariance between the real and imaginary parts. Can also be thought of as the amount and direction of anisotropy in the 
#' presumed (complex normal) probability distribution in the complex plane. The argument of this number gives the direction of the semi-major axis.
#'
#' @return
#' Returns a list containing the following elements.
#' \item{`residuals`}{Complex or numeric. The weighted residuals, that is the measured value minus the fitted value, scaled by the square root of the weights given in the call to lm.}
#' \item{`correlation`}{A numeric matrix. The computed correlation coefficient matrix for the coefficients in the model.}
#' \item{`pseudocorrelation`}{A complex matrix. The computed pseudo-correlation coefficient matrix for the coefficients in the model.}
#' \item{`cov.unscaled`}{The unscaled covariance matrix; i.e, a numeric matrix such that multiplying it by an estimate of the error variance produces an estimated covariance matrix for the coefficients.}
#' \item{`pcov.unscaled`}{The unscaled pseudo-covariance matrix; i.e, a complex matrix such that multiplying it by an estimate of the error pseudo-variance produces an estimated pseudo-covariance matrix for the coefficients.}
#' \item{`sigma`}{Numeric. The square root of the estimated variance of the random error.}
#' \item{`psigma`}{Complex. The square root of the estimated pseudo-variance of the random error. See details above.}
#' \item{`df`}{The number of degrees of freedom for the model and for residuals. A 3 element vector (p, n-p, p*), the first being the number of non-aliased coefficients, the last being the total number of coefficients.}
#' \item{`coefficients`}{A 5 column matrix that contains the model coefficients, their standard errors, their pseudo standard errors (see details above), their t statistics, and corresponding (two-sided) p-value. Aliased coefficients are omitted.}
#' \item{`aliased`}{Named logical vector showing if the original coefficients are aliased.}
#' \item{`terms`}{The terms object used in fitting this model.}
#' \item{`fstatistic`}{(for models including non-intercept terms) a 3 element numeric vector with the value of the F-statistic with its numerator and denominator degrees of freedom.}
#' \item{`r.squared`}{Numeric. The ???fraction of variance explained by the model???.}
#' \item{`adj.r.squared`}{the above R^2 statistic ???adjusted???, penalizing for higher p.}
#' \item{`symbolic.cor`}{(only if correlation is true.) The value of the argument symbolic.cor.}
#' \item{`na.action`}{from `object`, if present there.}
#' 
#' @export
#'
#' @examples
#' set.seed(4242)
#' n = 8
#' slope = complex(real = 4.23, imaginary = 2.323)
#' intercept = complex(real = 1.4, imaginary = 1.804)
#' testframe <- data.frame(x = x <- complex(real = rnorm(n), imaginary = rnorm(n)), y = slope * x + intercept)
#' fit <- lm(y ~ x, data = testframe, weights = rep(1,n))
#' summary.lm(fit)
summary.lm <- function (object, correlation = FALSE, symbolic.cor = FALSE, ...)
{
  cll <- match.call()
  if (is.numeric(object$residuals)) 
  {
    cll[[1]] <- stats::summary.lm
    eval(cll, parent.frame())
  }
  else
  {
    z <- object
    p <- z$rank
    rdf <- z$df.residual
    if (p == 0) {
      r <- z$residuals
      n <- length(r)
      w <- z$weights
      if (is.null(w)) {
        rss <- sum(as.numeric(Conj(r)*r))
        prss <- sum(r^2) # 'pseudo' sum of squared residuals.
      } else {
        rss <- sum(w * as.numeric(Conj(r)*r))
        prss <- sum(r^2) # 'pseudo' sum of squared residuals.
        r <- sqrt(w) * r
      }
      resvar <- rss/rdf # Variance of residuals.
      respvar <- prss / rdf # Pseudo-variance of residuals.
      ans <- z[c("call", "terms", if(!is.null(z$weights)) "weights")]
      class(ans) <- "summary.lm"
      ans$aliased <- is.na(coef(object))  # used in print method
      ans$residuals <- r
      ans$df <- c(0L, n, length(ans$aliased))
      ans$coefficients <- matrix(NA_real_, 0L, 5L, dimnames =
                                   list(NULL, c("Estimate", "Std. Error", "Pseudo Std. Error", "t value", "Pr(>|t|)")))
      ans$sigma <- sqrt(resvar)
      ans$psigma <- sqrt(respvar) # Pseudo standard deviation of residuals.
      ans$r.squared <- ans$adj.r.squared <- 0
      ans$cov.unscaled <- matrix(NA_real_, 0L, 0L)
      if (correlation) ans$correlation <- ans$cov.unscaled
      return(ans)
    }
    if (is.null(z$terms))
      stop("invalid 'lm' object:  no 'terms' component")
    if(!inherits(object, "lm"))
      warning("calling summary.lm(<fake-lm-object>) ...")
    Qr <- stats:::qr.lm(object) # Internal function that just returns the thing in z$qr. Unless it's not there, in which case it gives an error.
    n <- NROW(Qr$qr)
    if(is.na(z$df.residual) || n - p != z$df.residual)
      warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
    ## do not want missing values substituted here
    r <- z$residuals
    f <- z$fitted.values
    w <- z$weights
    if (is.null(w)) {
      mss <- if (attr(z$terms, "intercept"))
        sum(as.numeric(Conj(f - mean(f))*(f - mean(f)))) else sum(as.numeric(Conj(f)*f)) # Sum of conjugate squared deviations from the mean of the fitted values.
      pmss <- if (attr(z$terms, "intercept"))
        sum((f - mean(f))^2) else sum(f^2) # Sum of squared deviations from the mean of the fitted values. A complex number.
      rss <- sum(as.numeric(Conj(r)*r))
      prss <- sum(r^2)
    } else {
      mss <- if (attr(z$terms, "intercept")) {
        m <- sum(w * f /sum(w))
        sum(w * as.numeric(Conj(f - m)*(f - m)))
      } else sum(w * as.numeric(Conj(f)*f))
      pmss <- if (attr(z$terms, "intercept")) {
        m <- sum(w * f /sum(w))
        sum(w * (f - m)^2)
      } else sum(w * f^2)
      rss <- sum(w * as.numeric(Conj(r)*r))
      prss <- sum(w * r^2)
      r <- sqrt(w) * r
    }
    resvar <- rss/rdf # Variance of the residuals.
    respvar <- prss/rdf # Pseudo variance of the residuals.
    ## see thread at https://stat.ethz.ch/pipermail/r-help/2014-March/367585.html
    if (is.finite(resvar) &&
        resvar < (as.numeric(Conj(mean(f))*mean(f)) + var(c(f))) * 1e-30)  # a few times .Machine$double.eps^2
      warning("essentially perfect fit: summary may be unreliable")
    p1 <- 1L:p
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    se <- sqrt(diag(R) * resvar)
    pse <- sqrt(diag(R) * respvar)
    est <- z$coefficients[Qr$pivot[p1]]
    tval <- est/se
    ans <- z[c("call", "terms", if(!is.null(z$weights)) "weights")]
    ans$residuals <- r
    ans$coefficients <-
      cbind(Estimate = est, "Std. Error" = se, "Pseudo Std. Error" = pse, "t value" = tval,
            "Pr(>|t|)" = 2*pt(abs(tval), rdf, lower.tail = FALSE))
    ans$aliased <- is.na(z$coefficients)  # used in print method
    ans$sigma <- sqrt(resvar)
    ans$psigma <- sqrt(respvar)
    ans$df <- c(p, rdf, NCOL(Qr$qr))
    if (p != attr(z$terms, "intercept")) {
      df.int <- if (attr(z$terms, "intercept")) 1L else 0L
      ans$r.squared <- mss/(mss + rss)
      ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int)/rdf)
      ans$fstatistic <- c(value = (mss/(p - df.int))/resvar,
                          numdf = p - df.int, dendf = rdf)
    } else ans$r.squared <- ans$adj.r.squared <- 0
    ans$cov.unscaled <- R
    dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1,1)]
    if (correlation) {
      ans$correlation <- (R * resvar)/outer(se, se)
      dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
      ans$symbolic.cor <- symbolic.cor
    }
    if(!is.null(z$na.action)) ans$na.action <- z$na.action
    class(ans) <- "summary.lm"
    ans
  }
}

#' Calculate Variance-Covariance Matrix and Pseudo Variance-Covariance Matrix for a Complex Fitted Model Object
#'
#' A version of [stats::vcov] that is compatible with complex linear models. In addition to variance-covariance matrix,
#' the pseudo variance-covariance matrix, which is a measure of the covariance between real and imaginary components, is returned as well.
#' 
#' @param object a fitted model object, typically. Sometimes also a summary() object of such a fitted model.
#' @param ... Additional parameters, not currently used for anything.
#'
#' @return A list containing both the numeric variance-covariance matrix, and the complex pseudo variance-covariance matrix.
#' @export
#'
#' @examples
#' set.seed(4242)
#' n = 8
#' slope = complex(real = 4.23, imaginary = 2.323)
#' intercept = complex(real = 1.4, imaginary = 1.804)
#' testframe <- data.frame(x = x <- complex(real = rnorm(n), imaginary = rnorm(n)), y = slope * x + intercept)
#' fit <- lm(y ~ x, data = testframe, weights = rep(1,n))
#' vcov.lm(fit)
vcov.lm <- function (object, ...)
{
  cll <- match.call()
  if (is.numeric(object$residuals))
  {
    cll[[1]] <- stats:::vcov.lm
    eval(cll, parent.frame())
  }
  else
  {
    so <- summary(object, corr = FALSE)
    varcovar <- so$stddev^2 * so$cov.unscaled
    pseudovarcovar <- so$pstddev^2 * so$pcov.unscaled
    return(list(varcovar = varcovar, pseudovarcovar = pseudovarcovar))
  }  
}

