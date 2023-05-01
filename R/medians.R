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
#' @note The concept of Quantile requires ordering to be defined, which the complex numbers lack. 
#' The usefulness of multiplying by `constant` is thus called into question. However, for no more rigorous
#' reason than consistency, the default behavior of this function is to do so.
#'
#' @return numeric or complex. The median absolute deviation (MAD) from center.
#' @export
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
    if (p[n + 1L] > 0.5) x[n + 1L] else (x[n + 1L] + x[n + 2L])/2 ## For a normal distribution, the standard deviation about equals MAD/0.6745 #Removed the /0.6745, thus making this just a weighted median function.
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
#' Variance for Complex Data
#'
#' A wrapper for [stats::var] that can handle complex variables.
#'
#' @param x a numeric or complex vector, matrix, or dataframe.
#' @param y NULL (default) or a numeric vector, matrix, or dataframe with dimensions compatible with x. For complex x,
#' this parameter is not used as `cov()` and `cor()` are not implemented for complex variables.
#' @param na.rm logical. Should missing values be removed?
#' @param use character string giving the desired method of computing covariances in the presense of missing values. Not used.
#' 
#' @details The sample variance is calculated as `sum(Conj( mean(x) - x ) * ( mean(x) - x )) / (length(x) - 1)`.
#'
#' @return numeric, the sample variance of the input data.
#' @export
#'
#' @examples
#' set.seed(4242)
#' n <- 9
#' foo <- complex(real = rnorm(n), imaginary = rnorm(n))
#' var(foo)
var <- function(x, y = NULL, na.rm = FALSE, use)
  {
  cll <- match.call()
  cll[[1]] <- stats::var
  if (is.numeric(x)) eval(cll, parent.frame())
  else 
  {
    if (na.rm == TRUE) x <- x[!is.na(x)]
    if (length(x) == 1) return(NA)
    else return(vvar <- sum(as.numeric(Conj(mean(x) - x )*(mean(x) - x ))) / (length(x) - 1)) # as.numeric() needed to convert type to numeric.j
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

#' Combine covariance matrix and pseudo covariance matrix into a "big covariance matrix"
#' 
#' Interleaves the elements of a pxp matrix with those of a different pxp matrix to form a 2px2p matrix. 
#' This function was originally made to combine the covariance and pseudo covariance matrices of a 
#' complex random vector into a single "big covariance matrix", as described in [1]. However, it will accept
#' and operate on matrices of any mode.
#'
#' @param cov A square matrix, such as one describing the covariance between two complex random vectors.
#' @param pcov A square matrix with the same size as cov. Perhaps a pseudo covariance matrix.
#' 
#' @return
#' @export
#' 
#' @references [1] A. van den Bos, The Multivariate Complex Normal Distribution-a Generalization, IEEE Trans. Inform. Theory 41, 537 (1995).
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
#' the conjugate transpose if given complex data, otherwise it calls [stats::mahalanobis()].
#'
#' @inherit stats::mahalanobis return examples
#' 
#' @param x A length \eqn{p} vector or matrix with row length \eqn{p}. Or, a length \eqn{2p} vector or matrix with row length \eqn{2p}.
#' @param center A vector of length equal to that of `x`.
#' @param cov The covariance matrix \eqn{(p x p)} of the distribution. Or, the "big covariance matrix" of the distribution, which contains the information from `cov` and `pcov` in a single \eqn{(2p x 2p)} matrix. Can be generated by [matrixweave], [lm.vcov], or [rlm.vcov].
#' @param pcov The pseudo covariance matrix \eqn{(p x p)} of the distribution. Optional.
#' @param inverted Boolean, if TRUE, `cov` and `pcov` are not taken to be the \emph{inverse} covariance and pseudo covariance matrices.
#' 
#' @details Depending on the relative sizes of `x`, `cov`, and `pcov`, the function will perform slightly different calculations. If `pcov` is not included, 
#' the Mahalanobis distance is calculated using only `cov`. In this case if the dimension of `cov` is twice that of `x`, `x` is interleaved with its complex conjugate 
#' so that it becomes the same length as `cov`. Note that in this case the resulting Mahalanobis distance will only incorporate information about the interactions between
#'  the real and imaginary components if the "big covariance matrix is given as `cov` . If `pcov` is included in the input, `pcov` and `cov` are interleaved to form the "big covariance", and this is used to 
#' calculate the Mahalanobis distance, interleaving `x` if necessary. This gives the user a great deal of flexibility when it comes to input. 
#' 
#' 
#' @references D. Dai and Y. Liang, High-Dimensional Mahalanobis Distances of Complex Random Vectors, Mathematics 9, 1877 (2021).
#' 
#' @export
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
    if(!is.null(pcov)) cov <- matrixweave(cov, pcov) # if pcov is specified, assume that the user would like to calculate the mahalanobis distance with the big covariance matrix.
    
    if(ncol(cov) > p) # If the dimension of cov is greater than the length of x (or the length of its rows), it probably means that x and center have not been interleaved with their complex conjugates yet.
    {
      star <- vector(mode = mode(x), 2 * p)
      star[,seq(1, 2 * p, 2)] <- x
      star[,seq(2, 2 * p, 2)] <- Conj(x)
      x <- star
    }
    
    if(!inverted)
      cov <- solve(cov, ...) # Returns inverse of cov.
    setNames(rowSums(Conj(x) %*% cov * x), rownames(x))
  }
}