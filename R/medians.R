# file complexlm/R/medians.R
# copyright (C) 2022 W. L. Ryan
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

####
### This function combines the median from stats and the geo_median from pracma. 
### The former is used when given numeric data, the later for complex input.
### Returns a numeric or complex.
### Maybe turn this into a method for median, to do so, use: setMethod("median", "complex", this function without the outer)
####
#' Median, for Numeric or Complex
#' 
#' Extends [stats::median] to understand complex variable input. If x is complex, the geometric median
#' is calculated by [pracma::geo_median], and returned as a complex number. Otherwise, [stats::median] is called.
#'
#' @inherit median references
#' @inherit pracma::geo_median references
#'
#' @param x a numeric or complex vector, of which the median is to be calculated.
#' @param na.rm logical. Should NA values be removed before finding the median?
#' @param tol numeric. Relative tolerance to be passed to [pracma::geo_median]. Default is 1e-07.
#' @param maxiter maximum number of iterations for calculating geometric median. Not used if x is numeric.
#'
#' @return The median of x. If x is complex, the geometric median as calculated by Weiszfeld's algorithm.
#' @export
#' 
#' @seealso [stats::median] and [pracma::geo_median]
#'
#' @examples
#' set.seed(4242)
#' n <- 7
#' foo <- complex(real = rnorm(n), imaginary = rnorm(n))
#' median(foo)
median <- function(x, na.rm = FALSE, tol = 1e-07, maxiter = 200)
{
  matchcall <- match.call()
  matchcall[[1]] <- stats::median
  if (is.numeric(x)) eval(matchcall, parent.frame())
  else
  {
    Zxmatrix <- as.matrix(data.frame(re = Re(x), im = Im(x)))
    gmed <- pracma::geo_median(Zxmatrix, tol, maxiter)
    return(complex(real = gmed['p']$p[1], imaginary =  gmed['p']$p[2]))
  }
}

####
### Median absolute deviation, adapted to operate on complex data as well as numeric.
### In the later case it simply calls the mad from stats.
### For complex x it uses the geometric median, geo_median(), from pracma as the center,
### then returns the median absolute difference between center and each element of x.
####
#' Median Absolute Deviation, compatible with complex variables
#' 
#' Median absolute deviation, adapted to operate on complex data as well as numeric.
#' In the later case it simply calls [stats::mad()].
#' For complex x it uses the geometric median, [pracma::geo_median()], as the center,
#' then returns the median absolute difference between center and each element of x, multiplied by `constant`.
#'
#' @param x a numeric or complex vector.
#' @param center optional, numeric or complex. The center about which to calculate MAD. Defaults to median for numeric, and geo_median for complex.
#' @param constant a constant by which to multiply the median absolute deviation from center. Default is 1.4826, which is one over the quantile of 3/4 for the normal distribution.
#' @param na.rm logical. Should NAs be removed from x before calculating.
#' @param low logical. If TRUE, compute the "lo-median", i.e., for even sample size, do not average the two middle values, but take the smaller one. Not used if x is complex.
#' @param high logical. If TRUE, compute the "hi-median", i.e., take the larger of the two middle values for even sample size. Not used if x is complex.
#'
#' @note The concept of quantile requires ordering to be defined, which the complex numbers lack. 
#' The usefulness of multiplying by `constant` is thus called into question. However, for no more rigorous
#' reason than consistency, the default behavior of this function is to do so.
#'
#' @return numeric. The median absolute deviation (MAD) from center.
#' @export
#' 
#' @seealso [median], [stats::mad]
#'
#' @examples
#' set.seed(4242)
#' n <- 8
#' foo <- complex(real = rnorm(n), imaginary = rnorm(n))
#' mad(foo)
mad <- function(x, center = median(x), constant = 1.4826, na.rm = FALSE, low = FALSE, high = FALSE)
{
  cll <- match.call()
  cll[[1]] <- stats::mad
  if (is.numeric(x)) eval(cll, parent.frame())
  else 
  {
    if (na.rm == TRUE) x <- x[!is.na(x)]
    distances <- Mod(x - center)
    zmad <- median(distances)
    # Note: the standard mad() function from the stats package scales the mad by the constant, which is set to 
    # one over the quantile of 3/4 for the normal distribution. The Quantile function requires ordering
    # to be defined though, which the complex numbers lack. For lack of a more rigorous idea, I'll just multiply
    # the complex mad by constant as well.
    return(constant * zmad)
  }
}

####
### This actually just finds the weighted absolute median. Used in rlm.default, 
### it is passed the residuals as x, since these residuals are (kind of) equal 
### to measurement - median it behaves like the WMAD in that context. Name changed to reflect this.
### TO DO: Come up with some kind of complex weighted median.
####
#' Weighted Median
#' 
#' This calculates the weighted median of a vector `x` using the weights in `w`. Weights are re-scaled based on their sum.
#'
#' @param x numeric, a vector containing the data from which to calculate the weighted median.
#' @param w numeric, a vector of weights to give the data in x.
#'
#' @details Sorts `x` and `w` by size of the elements of `x`. Then re-scales the elements of `w` to be between 0 and 1.
#' Then sets n equal to the sum of all scaled weights with values less than 0.5. If the (n+1)-th element of the 
#' rescaled weights is greater than 0.5, the weighted median is the (n+1)-th element of the sorted `x`. Otherwise
#' it is the average of the (n+1)-th and (n+2)-th elements of the sorted `x`.
#'
#' @note This is not compatible with complex data.
#'
#' @return numeric. The weighted median of `x` using `w` as the weights.
#' @export
#' 
#' @references F. Y. Edgeworth, XXII. On a New Method of Reducing Observations Relating to Several Quantities, (1888).
#' Also see the Wikipedia article on weighted median for a very good explanation and a model algorithm.
#' 
#' @seealso [median]
#'
#' @examples
#' xx <- rnorm(10, 4L, 1.5)
#' ww <- runif(10)
#' wmedian(xx, ww)
#' 
wmedian <- function(x, w = rep(1, length(x)))
{
  if (is.numeric(x)) { ## Works fine for complex regression because we take the absolute value of the residual in finding median absolute deviation.
    o <- sort.list(x); x <- x[o]; w <- w[o] # Removed abs() from around x, so that the function is more general.
    p <- cumsum(w)/sum(w)
    n <- sum(p < 0.5) ## Count how many of elements of p are greater than .5
    if (p[n + 1L] > 0.5) return(x[n + 1L]) else return((x[n + 1L] + x[n + 2L])/2) ## For a normal distribution, the standard deviation about equals MAD/0.6745 #Removed the /0.6745, thus making this just a weighted median function.
  }
  else { ## Could be nice to have a weighted median function for complex variables, but not strictly necessary.
    ## Zxmatrix <- as.matrix(data.frame(re = Re(x), im = Im(x)))
    ## geomed <- pracma::geo_median(Zxmatrix) # From the pracma-package
    ## distances <- apply(X = Zxmatrix, MARGIN = 1, FUN = function(z) sum((z - geomed$p)^2)^0.5)
    print("Sorry, this function only works with numeric at the moment. If you have an idea for a complex weighted median algorithm, please contact the maintainer.")
  }
}

####
### A wrapper for var from the stats package that will accept (and use) complex numbers.
### var is used in summary.rlm to find variance of a set of complex numbers (psiprime).
### Perhaps add cor and cov functionality in the future.
####
#' Variance, Covariance, and Correlation for Complex Data
#'
#'  Wrappers of [stats::var], [stats::cov], and [stats::cor] that are capable of handling complex input.
#'
#' @param x a numeric or complex vector, matrix, or dataframe.
#' @param y NULL (default) or a numeric vector, matrix, or dataframe with dimensions compatible with x.
#' @param na.rm logical. Should missing values be removed? Only considered when `x` is a vector.
#' @param use character string giving the desired method of computing covariances in the presence of missing values. Options are "everything" (default),
#' "all.obs", "complete.obs", or "na.or.complete". See [stats::cov] for explanation of what each one does. Note that "pairwise.complete.obs" is not available for this complex method.
#' @param method The method for calculating correlation coefficient. Only `"pearson"` is supported for complex variables, so this parameter is ignored.
#' @param ... Other parameters, ignored.
#' 
#' @details For vector input, the sample variance is calculated as,\cr
#'  \eqn{sum(Conj( mean(x) - x ) * ( mean(x) - x )) / (length(x) - 1)}\cr
#'  And the sample covariance is calculated as, \cr
#'  \eqn{sum(Conj( mean(x) - x ) * ( mean(y) - y )) / (length(x) - 1)}\cr
#'  The Pearson correlation coefficient, which is the only kind available for complex data, is the covariance divided by the product of the standard deviations of all variables.
#'
#' @return numeric or complex the sample variance, covariance, or correlation of the input data.
#' @export
#'
#' @examples
#' set.seed(4242)
#' n <- 9
#' foo <- complex(real = rnorm(n), imaginary = rnorm(n))
#' var(foo)
#' bar <- complex(real = rnorm(n), imaginary = rnorm(n))
#' var(x = foo, y = bar)
#' foobar <- data.frame(foo, bar)
#' cov(foobar)
#' cor(foobar)
cov <- function(x, y = NULL, na.rm = FALSE, method = "pearson", use = "everything", ...)
  {
  matdf <- is.matrix(x) || is.data.frame(x) # Is x a matrix or dataframe?
  cll <- match.call()
  if (matdf) {
    cll[[1]] <- stats::cov
    if (is.numeric(x[[1]])) eval(cll, parent.frame())
  }
  cll[[1]] <- stats::var
  if (is.numeric(x)) eval(cll, parent.frame())
  else
  {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.data.frame(y)) y <- as.matrix(y)
    if (!is.matrix(x)) # Deal with the vector case first, it shares little with the matrix / dataframe code.
    {
      if (na.rm == TRUE) x <- x[!is.na(x)]
      if (length(x) == 1) return(NA)
      mn <- mean(x)
      if (is.null(y)) return(sum(as.numeric(Conj(x - mn)*(x - mn))) / (length(x) - 1)) # as.numeric() needed to convert type to numeric.
      if (na.rm == TRUE) {
        y <- y[!is.na(x)]
        x <- x[!is.na(y)]
        y <- y[!is.na(y)]
      }
      mny <- mean(y)
      return(sum((x - mn) * Conj(y - mny)) / (length(x) - 1))
    }
    # Code for dealing with dataframe or matrix input. Uses lots of if statements to comprehend the 'use' parameter.
    if (grepl("complete", use)) {
      drops <- is.na(rowSums(x)) # Clever way to find rows with NAs without explicitly using apply.
      if (!is.null(y)) drops <- drops | is.na(rowSums(y))
      x <- x[!drops,]
      if (!is.null(y)) y <- y[!drops,]
      if (length(x) == 0) {
        if (grepl("na.or", use)) return(NA)
        stop("Error, no complete cases (rows) of observations (input data).")
      }
    }
    if (grepl("all", use)) {
      if (!is.null(y)) nays <- any(is.na(y)) else nays <- FALSE
      if (any(is.na(x)) || nays) stop("Error. Missing values not accepted with use = \"all.obs\"")
    }
    if (grepl("pairwise.complete.obs", use)) warning("use option pairwise.complete.obs is not available for complex input. Defaulting to use = everything.")
    # Now we actually calculate stuff. Since use = "everything" is default, we don't have an if statement for it.
    mns <- colMeans(x)
    x <- x - mns[col(x)]
    if (is.null(y)) return((Conj(t(x)) %*% x)  / (length(x[,1] - 1)))
    y <- y - colMeans(y)[col(y)]
    return((Conj(t(x)) * y)  / (length(x[,1] - 1)))
  }
  # matdf <- is.matrix(x) || is.data.frame(x) # Is x a matrix or dataframe?
  # cll <- match.call()
  # cll[[1]] <- stats::var
  # if (matdf) {
  #   if (is.numeric(x[[1]])) eval(cll, parent.frame())
  # }
  # if (is.numeric(x)) eval(cll, parent.frame())
  # else 
  # {
  #   if (na.rm == TRUE) x <- x[!is.na(x)]
  #   mn <- mean(x)
  #   if (length(x) == 1) return(NA)
  #   else return(vvar <- sum(as.numeric(Conj(x - mn)*(x - mn))) / (length(x) - 1)) # as.numeric() needed to convert type to numeric.j
  # }
}

#' @describeIn cov Correlation coefficient of complex variables.
#' @export
cor <- function(x, y = NULL, na.rm = FALSE, use = "everything", method = "pearson", ...)
{
  matdf <- is.matrix(x) || is.data.frame(x) # Is x a matrix or dataframe?
  cll <- match.call()
  cll[[1]] <- stats::cor
  if (matdf) {
    if (is.numeric(x[[1]])) eval(cll, parent.frame())
  }
  if (is.numeric(x)) eval(cll, parent.frame())
  else
  {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.data.frame(y)) y <- as.matrix(y)
    if (!is.matrix(x)) # Deal with the vector case first, it shares little with the matrix / dataframe code.
    {
      if (na.rm == TRUE) x <- x[!is.na(x)]
      if (length(x) == 1) return(NA)
      mn <- mean(x)
      sdx <- sqrt(sum((x - mn) * Conj(x - mn)) / (length(x) - 1))
      if (is.null(y)) return(sum(as.numeric(Conj(x - mn)*(x - mn))) / (length(x) - 1)) # as.numeric() needed to convert type to numeric.
      if (na.rm == TRUE) {
        y <- y[!is.na(x)]
        x <- x[!is.na(y)]
        y <- y[!is.na(y)]
      }
      mny <- mean(y)
      sdy <- sqrt(sum((y - mny) * Conj(y - mny)) / (length(y) - 1))
      covv <- sum((x - mn) * Conj(y - mny)) / (length(x) - 1)
      return(covv / (sdx * sdy))
    }
    # Code for dealing with dataframe or matrix input. Uses lots of if statements to comprehend the 'use' parameter.
    if (grepl("complete", use)) {
      drops <- is.na(rowSums(x)) # Clever way to find rows with NAs without explicitly using apply.
      if (!is.null(y)) drops <- drops | is.na(rowSums(y))
      x <- x[!drops,]
      if (!is.null(y)) y <- y[!drops,]
      if (length(x) == 0) {
        if (grepl("na.or", use)) return(NA)
        stop("Error, no complete cases (rows) of observations (input data).")
      }
    }
    if (grepl("all", use)) {
      if (!is.null(y)) nays <- any(is.na(y)) else nays <- FALSE
      if (any(is.na(x)) || nays) stop("Error. Missing values not accepted with use = \"all.obs\"")
    }
    if (grepl("pairwise.complete.obs", use)) warning("use = \"pairwise.complete.obs\" is not available for complex input. Defaulting to use = everything.")
    # Now we actually calculate stuff. Since use = "everything" is default, we don't have an if statement for it.
    mns <- colMeans(x)
    x <- x - mns[col(x)]
    sdx <- colSums(x * Conj(x)) / (length(x[,1] - 1))
    if (is.null(y)) {
      covv <- (Conj(t(x)) %*% x) / (length(x[,1] - 1))
      sdmat <- outer(sdx, sdx)
    return(covv / sdmat)
    }
    y <- y - colMeans(y)[col(y)]
    sdy <- colSums(y * Conj(y)) / (length(x[,1] - 1))
    sdmat <- outer(sdx, sdy)
    covv <- (Conj(t(x)) * y)  / (length(x[,1] - 1))
    return(covv / sdmat)
  }
}

#' @describeIn cov S3 Variance of complex variables, a synonym for [complexlm::cov].
#' @export
var <- function(x, y = NULL, na.rm = FALSE, use = "everything", ...)
{
  matdf <- is.matrix(x) || is.data.frame(x) # Is x a matrix or dataframe?
  cll <- match.call()
  if (matdf) {
    cll[[1]] <- stats::cov
    if (is.numeric(x[[1]])) eval(cll, parent.frame())
  }
  cll[[1]] <- stats::var
  if (is.numeric(x)) eval(cll, parent.frame())
  else
  {
    cll <- match.call()
    cll[[1]] <- cov
    eval(cll, parent.frame())
  }
}


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

#' Combine covariance matrix and pseudo covariance matrix into a "double covariance matrix"
#' 
#' Interleaves the elements of a \eqn{(p x p)} matrix with those of a different \eqn{(p x p)} matrix to form a \eqn{(2p x 2p)} matrix. 
#' This function was originally made to combine the covariance and pseudo covariance matrices of a 
#' complex random vector into a single "double covariance matrix", as described in (van den Bos 1995). However, it will accept
#' and operate on matrices of any storage mode.
#'
#' @param cov A square matrix, such as one describing the covariance between two complex random vectors.
#' @param pcov A square matrix with the same size as cov. Perhaps a pseudo covariance matrix.
#' 
#' @return A square matrix with dimension twice that of the input matrices. Each element is an element from one of the inputs, and its adjacent elements are from the other input.
#' @export
#' 
#' @references A. van den Bos, The Multivariate Complex Normal Distribution-a Generalization, IEEE Trans. Inform. Theory 41, 537 (1995).
#' 
#' @seealso [mahalanobis], [vcov.zlm], [vcov.rzlm]
#'
#' @examples
#' set.seed(4242)
#' mata <- matrix(rnorm(9), nrow = 3)
#' matb <- matrix(rnorm(9), nrow = 3)
#' matrixweave(mata, matb)
matrixweave <- function(cov, pcov)
{
  n <- attributes(cov)$dim[1] # The size of the square covariance matrix.
  if (attributes(pcov)$dim[1] != n) 
  {
    #warning("cov and pcov must be the same size.")
    stop("cov and pcov must be the same size.")
  }
  # print(attributes(cov))
  bigcovar <- matrix(0, ncol = 2*n, nrow = 2*n) #Start by making an empty matrix with dimensions twice those of the small covariance matrix.
  bigcovar[,seq(1, 2*n, 2)] <- matrix(as.vector(rbind(as.vector(cov), as.vector(Conj(pcov)))), nrow = 2*n) # Fill the odd indexed columns of bigcovar with the entries from cov interleaved with the entries from pcov conjugated.
  bigcovar[,seq(2, 2*n, 2)] <- matrix(as.vector(rbind(as.vector(pcov), as.vector(cov))), nrow = 2*n) # Fill the even indexed columns of bigcovar with the entries from cov interleaved with the entries from pcov, this time the later first.
  return(bigcovar)
}

#' Mahalanobis Distance, with better complex behavior
#' 
#' The Mahalanobis distance function included in the `stats` package returns a complex number when given complex values of `x`.
#' But a distance (and thus its square) is always positive real. This function calculates the Mahalanobis distance using
#' the conjugate transpose if given complex data, otherwise it calls [stats::mahalanobis].
#'
#' @inherit stats::mahalanobis return examples
#' 
#' @param x A length \eqn{p} vector or matrix with row length \eqn{p}. Or, a length \eqn{2p} vector or matrix with row length \eqn{2p}.
#' @param center A vector of length equal to that of `x`.
#' @param cov The covariance matrix \eqn{(p x p)} of the distribution. Or, the "double covariance matrix" of the distribution, which contains the information from `cov` and `pcov` in a single \eqn{(2p x 2p)} matrix. 
#' Can be generated by [matrixweave], [vcov.zlm], or [vcov.rzlm].
#' vcov.rzlm].
#' @param pcov The pseudo covariance matrix \eqn{(p x p)} of the distribution. Optional.
#' @param inverted Boolean, if TRUE, `cov` and `pcov` are not taken to be the \emph{inverse} covariance and pseudo covariance matrices.
#' @param ... Optional arguments to be passed to [solve], which is used for computing the inverse of `cov`. If `inverted = TRUE`, unused.
#' 
#' @details Depending on the relative sizes of `x`, `cov`, and `pcov`, the function will perform slightly different calculations. If `pcov` is not included, 
#' the Mahalanobis distance is calculated using only `cov`. In this case if the dimension of `cov` is twice that of `x`, `x` is interleaved with its complex conjugate 
#' so that it becomes the same length as `cov`. Note that in this case the resulting Mahalanobis distance will only incorporate information about the interactions between
#'  the real and imaginary components if the "double covariance matrix is given as `cov` . If `pcov` is included in the input, `pcov` and `cov` are interleaved to form the "double covariance", and this is used to 
#' calculate the Mahalanobis distance, interleaving `x` if necessary. This gives the user a great deal of flexibility when it comes to input. 
#' 
#' 
#' @references D. Dai and Y. Liang, High-Dimensional Mahalanobis Distances of Complex Random Vectors, Mathematics 9, 1877 (2021).
#' 
#' @export
#' 
#' @seealso [matrixweave]
#'
#' @examples
#' set.seed(4242)
#' n <- 8
#' x <- matrix(complex(real = rnorm(n), imaginary = rnorm(n)), ncol = 2)
#' mu <- complex(real = 1.4, imaginary = 0.4)
#' sigma <- 3.4
#' mahalanobis(x, mu, sigma * diag(2))
mahalanobis <- function(x, center, cov, pcov = NULL, inverted=FALSE, ...)
{
  cll <- match.call()
  cll[[1]] <- stats::mahalanobis
  if (!is.complex(x)) eval(cll, parent.frame())
  else 
  {
    x <- if(is.vector(x)) matrix(x, ncol = length(x)) else as.matrix(x)
    p <- ncol(x)
    ## save speed in customary case
    if(!isFALSE(center))
      x <- sweep(x, 2L, center)# = "x - center"
    if(!is.null(pcov)) cov <- matrixweave(cov, pcov) # if pcov is specified, assume that the user would like to calculate the mahalanobis distance with the double covariance matrix.
    
    if(ncol(cov) > p) # If the dimension of cov is greater than the length of x (or the length of its rows), it probably means that x and center have not been interleaved with their complex conjugates yet.
    {
      star <- vector(mode = mode(x), 2 * p)
      star[,seq(1, 2 * p, 2)] <- x
      star[,seq(2, 2 * p, 2)] <- Conj(x)
      x <- star
    }
    
    if(!inverted)
      cov <- solve(cov, ...) # Returns inverse of cov.
    return(setNames(rowSums(Conj(x) %*% cov * x), rownames(x)))
  }
}

#' summary method for complex objects
#' 
#' The base summary method for complex objects only reports their length and that they are complex..
#' This improved method returns the mean, median, variance, and pseudo variance of the given complex object.
#'
#' @param object a complex vector or scalar.
#' @param ... additional arguments, not used.
#' @param digits integer specifying the number of digits to include in the summary values.
#' 
#' @return A complex vector containing in order the length, median, mean, variance, and pseudo variance of the object. 
#' The length element will be a positive integer, despite being in the complex mode.
#' @export
#'
#' @examples
#' set.seed(4242)
#' n <- 8
#' foo <- complex(real = rnorm(n), imaginary = rnorm(n))
#' summary(foo)
summary.complex <- function(object, ..., digits)
{
  value <- c(length(object), median(object), mean(object), var(object), pseuzvar(object))
  if(!missing(digits)) value <- signif(value, digits)
  names(value) <- c("length", "median", "mean", "var.", "pvar.")
  class(value) <- c("summaryDefault", "table")
  return(value)
}

#' Range For Complex Objects
#'
#' This function extends [base::range] to the field of complex numbers.
#' It returns a vector containing two complex numbers that are the diagonal points of a rectangle,
#' with sides parallel to the real and imaginary axes, that just contains all the complex numbers 
#' given as arguments. If given non complex input it calls [base::range], please see the documentation
#' for that function for an explanation of its behavior with other input.
#'
#' @param ... Any complex, numeric, or character object
#' @param na.rm logical, indicates if `NA`'s should be removed.
#' @param finite logical, indicates if non-finite elements should be omitted.
#'
#' @return A complex vector describing a rectangle that all input values fall within.
#' @export
#' 
#' @seealso [base::range]
#'
#' @examples
#' set.seed(4242)
#' n <- 8
#' foo <- complex(real = rnorm(n), imaginary = rnorm(n))
#' range(foo)
range <- function(..., na.rm = FALSE, finite = FALSE)
{
  cll <- match.call()
  cll[[1]] <- base::range
  x <- c(..., recursive = TRUE)
  if (!is.complex(x)) eval(cll, parent.frame())
  else 
  {
    real.range <- range(Re(x), na.rm, finite)
    imag.range <- range(Im(x), na.rm, finite)
    return(c(complex(real = min(real.range), imaginary = min(imag.range)), complex(real = max(real.range), imaginary = max(imag.range))))
  }
}
