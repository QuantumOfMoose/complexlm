% File MASS/man/zlm.wfit.Rd
% Original Copyright 1995-2014 R Core Team
% Modifications Copyright 2021 William Ryan
% Distributed under GPL 2 or later

\name{zlm.wfit}
\title{Perform a weighted linear fit on complex valued data}
\usage{
zlm.wfit(x, y, w, offset = NULL, method = "qr", tol = 1e-7,
        singular.ok = TRUE, \dots)
}
\description{
  Works like \code{\link{lm.wfit}} but for complex valued data.
}
\arguments{
  \item{x}{design matrix of dimension \code{n * p}.}
  \item{y}{vector of observations of length \code{n}, or a matrix with
    \code{n} rows.}
  \item{w}{vector of weights (length \code{n}) to be used in the fitting
    process for the \code{wfit} functions.  Weighted least squares is
    used with weights \code{w}, i.e., \code{sum(|w * e^2|)} is minimized.}
  \item{offset}{(numeric of length \code{n}).  This can be used to
    specify an \emph{a priori} known component to be included in the
    linear predictor during fitting.}

  \item{method}{currently, only \code{method = "qr"} is supported.}

  \item{tol}{tolerance for the \code{\link{qr}} decomposition.  Default
    is 1e-7. Since \code{\link{qr}} uses LAPACK ZGEQP3 for complex matrices, tol is disregarded.}

  \item{singular.ok}{logical. If \code{FALSE}, a singular model is an
    error.}

  \item{\dots}{currently disregarded.}
}
\value{
  a \code{\link{list}} with components (for \code{lm.fit} and \code{lm.wfit})
  \item{coefficients}{\code{p} vector}
  \item{residuals}{\code{n} vector or matrix}
  \item{fitted.values}{\code{n} vector or matrix}
  \item{effects}{\code{n} vector of orthogonal single-df
    effects.  The first \code{rank} of them correspond to non-aliased
    coefficients, and are named accordingly.}
  \item{weights}{\code{n} vector --- \emph{only} for the \code{*wfit*}
    functions.}
  \item{rank}{integer, giving the rank}
  \item{df.residual}{degrees of freedom of residuals}
  \item{qr}{the QR decomposition, see \code{\link{qr}}.}

  Fits without any columns or non-zero weights do not have the
  \code{effects} and \code{qr} components.
}
\Note{
  \code{\link{lm}} does not work with complex valued data. To find a unweighted fit on complex data, use zlm.wfit with uniform weights.
}
\examples{
%% TODO: Add some examples.

%require(utils)
%% FIXME: Do something more sensible (non-random data) !!
%set.seed(129)

%n <- 7 ; p <- 2
%X <- matrix(rnorm(n * p), n, p) # no intercept!
%y <- rnorm(n)
%w <- rnorm(n)^2

%str(lmw <- lm.wfit(x = X, y = y, w = w))

%str(lm. <- lm.fit (x = X, y = y))
%\dontshow{
%  ## These are the same calculations at C level, but a parallel BLAS
%  ## might not do them the same way twice, and if seems serial MKL does not.
%  lm.. <- .lm.fit(X,y)
%  lm.w <- .lm.fit(X*sqrt(w), y*sqrt(w))
%  id <- function(x, y) all.equal(x, y, tolerance = 1e-15, scale = 1)
%  stopifnot(id(unname(lm.$coef), lm..$coef),
%	    id(unname(lmw$coef), lm.w$coef))
%}
%\donttest{
%if(require("microbenchmark")) {
%  mb <- microbenchmark(lm(y~X), lm.fit(X,y), .lm.fit(X,y))
%  print(mb)
%  boxplot(mb, notch=TRUE)
%}
%}
%% do an example which sets 'tol' and gives a difference!
}
\keyword{regression}
\keyword{array}
