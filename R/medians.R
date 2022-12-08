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
####
median <- function(x, na.rm = FALSE, tol = 1e-07, maxiter = 200)
{
  matchcall <- match.call()
  matchcall[[1]] <- stats::median
  if (is.numeric(x)) eval(matchcall, parent.frame())
  else
  {
    Zxmatrix <- as.matrix(data.frame(re = Re(x), im = Im(x)))
    gmed <- geo_median(Zxmatrix, tol, maxiter)
    return(complex(real = gmed['p']$p[1], imaginary =  gmed['p']$p[2]))
  }
}

####
### Median absolute deviation, adapted to operate on complex data as well as numeric.
### In the later case it simply calls the mad from stats.
### For complex x it uses the geometric median, geo_median(), from pracma as the center,
### then returns the median absolute difference between center and each element of x.
####
mad <- function(x, center = median(x), constant = 1.4826, na.rm = FALSE, low = FALSE, high = FALSE)
{
  cll <- match.call()
  cll[[1]] <- stats::mad
  if (is.numeric(x)) eval(cll, parent.frame())
  else 
  {
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
wmed <- function(x, w) 
{
  if (is.numeric(x)) { # Works fine for complex regression because we take the absolute value of the residual in finding median absolute deviation.
    o <- sort.list(x); x <- x[o]; w <- w[o] # Removed abs() from around x, so that the function is more general.
    p <- cumsum(w)/sum(w)
    n <- sum(p < 0.5) # Count how many of elements of p are greater than .5
    if (p[n + 1L] > 0.5) x[n + 1L] else (x[n + 1L] + x[n + 2L])/2 # For a normal distribution, the standard deviation ~equals MAD/0.6745 #Removed the /0.6745, thus making this just a weighted median function.
  }
  else { # Could be nice to have a weighted median function for complex variables, but not strictly necessary.
    Zxmatrix <- as.matriix(data.frame(re = Re(x), im = Im(x)))
    geomed <- geo_median(Zxmatrix) # From the pracma-package
    #distances <- apply(X = Zxmatrix, MARGIN = 1, FUN = function(z) sum((z - geomed$p)^2)^0.5)
    
  }
}

#####
#### A wrapper for var from the stats package that will accept (and use) complex numbers.
#### var is used in summary.rlm to find variance of a set of complex numbers (psiprime).
#### Perhaps add cor and var functionality in the future.
#####
var <- function(x, y = NULL, na.rm = FALSE, use)
{
  cll <- match.call()
  cll[[1]] <- stats::var
  if (is.numeric(x)) eval(cll, parent.frame())
  else 
  {
    if (length(x) == 1) return(NA)
    else return(vvar <- Sum(Conj(mean(x) - x )*(mean(x) - x )) / (length(x) - 1))
  }
}
