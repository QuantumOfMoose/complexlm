
####
### A function that somewhat replicates model.matrix(), but accepts complex valued data. It will probably be slower and less effecient, but mostly functional.
### It cannot handle algebraic expressions in formula.
### terms is the output of terms(formula)
####
#' zmodel.matrix
#'
#' @param trms 
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
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
lm.wfit <- function(x, y, w, offset = NULL, method = "qr", tol = 1e-7,
        singular.ok = TRUE, ...)
{
  cll <- match.call()
  if (is.complex(x)) cll[[1]] <- zlm.wfit
  else cll[[1]] <- stats::lm.wfit
  eval(cll, parent.frame())
}

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
    Qr <- qr.lm(object) # Internal function that just returns the thing in z$qr. Unless it's not there, in which case it gives an error.
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
