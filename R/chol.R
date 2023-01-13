#  File complexlm/R/chol.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 2023 William L. Ryan
#  Copyright (C) 1995-2020 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

#' @inherit base::chol
#' @export
#'
#' @note Now compatible with complex variables. Uses the `zpotri` LAPACK routine
chol <- function(x, ...) UseMethod("chol")

#' @describeIn chol Default S3 method:
#' @export
#'
chol.default <- function(x, pivot = FALSE, LINPACK = FALSE, tol = -1, ...)
{
    if (!missing(LINPACK))
        stop("the LINPACK argument has been defunct since R 3.1.0")
    if (is.complex(x))
    {
      #stop("complex matrices not permitted at present")
      .Internal(La_chol_cmlpx(as.matrix(x), pivot, tol))
    }
    else
      .Internal(La_chol(as.matrix(x), pivot, tol))
}

#' @inherit base::chol2inv
#' @export
#'
#' @note Now compatible with complex variables. Uses the `zpotri` LAPACK routine
chol2inv <- function(x, size = NCOL(x), LINPACK = FALSE)
{
    if (!missing(LINPACK))
        stop("the LINPACK argument has been defunct since R 3.1.0")
    if (is.complex(x)) .Internal(La_chol2inv_cmplx(x, size))
    else .Internal(La_chol2inv(x, size))
}
