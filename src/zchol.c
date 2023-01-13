/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2023 	      William L. Ryan
 *  Copyright (C) 2001--2019  The R Core Team.
 *  Copyright (C) 2003--2010  The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/*
 *  These functions are complex versions of functions found in r-source/src/modules/lapack/Lapack.c
 *
 * Exports
 *	La_chol_cmplx(...)
 *	La_chol2inv_cmplx(...)
 *
 * to be called as  .C(.) or maybe .Internal(.) in ../R/chol.R
 */

#include <math.h>  // currently in R.h
//#include "Defn.h"
#include <Rmath.h>
#include <stddef.h>  // currently in R.h
#include <complex.h> // This contains functions to act on complex numbers.
#include <Rinternals.h>

#include <R_ext/Complex.h> // This contains the Rcomplex type. (It's a struct.)
//#include "Rinlinedfuns.h"
#include <R.h>
#include "Rcomplex.h"
//#include "arithmetic.h"
#include <Rconfig.h>
#include <Rdefines.h>

#include <ctype.h>  /* for toupper */
#include <limits.h> /* for PATH_MAX */
#include <stdlib.h> /* for realpath */

#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>	/* for the QR	  routines */
#include <R_ext/Utils.h>	/* for the *sort() routines */
#include <R_ext/Arith.h>
#include <R_ext/Altrep.h>	/* for R_xlen_t */
#define BIG DBL_MAX


/* Make a zero Rcomplex number. NVM, doesn't work.
const Rcomplex rcplxzero;
rcplxzero->r = 0;
rcplxzero->i = 0;
*/

/* ---------------------------------------------------------- */

/* Complex case for Choleski Decomposition - Don't need La_chol_cmplx and can't get it to compile since I can't seem to make a Rcomplex zero. Skip it.*/
/* Bizzare formatting in these Choleski functions...*/
/*
static SEXP La_chol_cmplx(SEXP A, SEXP pivot, SEXP stol)
{
    if (!(isMatrix(A) && isComplex(A))) error(_("'a' must be a complex matrix"));
  
    SEXP ans = PROTECT(isComplex(A) ? duplicate(A): coerceVector(A, CPLXSXP));
    SEXP adims = getAttrib(A, R_DimSymbol);
    if (TYPEOF(adims) != INTSXP) error("non-integer dims");
    int m = INTEGER(adims)[0], n = INTEGER(adims)[1];
  
    if (m != n) error(_("'a' must be a square matrix"));
    if (m <= 0) error(_("'a' must have dims > 0"));
    size_t N = n;
    for (int j = 0; j < n; j++) */	/* zero the lower triangle *//*
  for (int i = j+1; i < n; i++) COMPLEX(ans)[i + N * j] = rcplxzero;
  
    int piv = asInteger(pivot);
    if (piv != 0 && piv != 1) error("invalid '%s' value", "pivot");
    if(!piv) {
  int info;
  F77_CALL(zpotrf)("U", &m, COMPLEX(ans), &m, &info FCONE);
  if (info != 0) {
       if (info > 0)
    error(_("the leading minor of order %d is not positive definite"),
          info);
      error(_("argument %d of Lapack routine %s had invalid value"),
      -info, "dpotrf");
  }
    } else {
  double tol = asReal(stol);
  SEXP piv = PROTECT(allocVector(INTSXP, m));
  int *ip = INTEGER(piv);
  double *work = (double *) R_alloc(2 * (size_t)m, sizeof(double));
  int rank, info;
  F77_CALL(zpstrf)("U", &m, COMPLEX(ans), &m, ip, &rank, &tol, work, &info
       FCONE);
  if (info != 0) {
      if (info > 0)
    warning(_("the matrix is either rank-deficient or indefinite"));
      else
    error(_("argument %d of Lapack routine %s had invalid value"),
          -info, "dpstrf");
  }
  setAttrib(ans, install("pivot"), piv);
  SEXP s_rank = install("rank");
  setAttrib(ans, s_rank, ScalarInteger(rank));
  SEXP cn, dn = getAttrib(ans, R_DimNamesSymbol);
  if (!isNull(dn) && !isNull(cn = VECTOR_ELT(dn, 1))) {
      // need to pivot the colnames
      SEXP dn2 = PROTECT(duplicate(dn));
      SEXP cn2 = VECTOR_ELT(dn2, 1);
      for(int i = 0; i < m; i++)
    SET_STRING_ELT(cn2, i, STRING_ELT(cn, ip[i] - 1)); // base 1
      setAttrib(ans, R_DimNamesSymbol, dn2);
      UNPROTECT(1);
  }
  UNPROTECT(1); // piv
    }
    UNPROTECT(1); // ans
    return ans;
}
*/

/* Complex case of changing Choleski decomposition to invese - The errror cause problems*/
static SEXP La_chol2inv_cmplx(SEXP A, SEXP size)
{
    int sz = asInteger(size);
    if (sz == NA_INTEGER || sz < 1) {
  error(_("'size' argument must be a positive integer"));
  return R_NilValue; /* -Wall */
    } else {
  SEXP ans, Amat = A; /* -Wall: we initialize here as for the 1x1 case */
  int m = 1, n = 1, nprot = 0;

  if (sz == 1 && !isMatrix(A) && isComplex(A)) {
    /* nothing to do; m = n = 1; ... */
  } else if (isMatrix(A)) {
    SEXP adims = getAttrib(A, R_DimSymbol);
    if (TYPEOF(adims) != INTSXP) error("non-integer dims");
    Amat = PROTECT(coerceVector(A, CPLXSXP)); nprot++;
    m = INTEGER(adims)[0]; n = INTEGER(adims)[1];
  } else error(_("'a' must be a complex matrix"));

  if (sz > n) { UNPROTECT(nprot); error(_("'size' cannot exceed ncol(x) = %d"), n); }
  if (sz > m) { UNPROTECT(nprot); error(_("'size' cannot exceed nrow(x) = %d"), m); }
  ans = PROTECT(allocMatrix(CPLXSXP, sz, sz)); nprot++;
  size_t M = m, SZ = sz;
  for (int j = 0; j < sz; j++) {
      for (int i = 0; i <= j; i++)
    COMPLEX(ans)[i + j * SZ] = COMPLEX(Amat)[i + j * M];
  }
  int info;
  F77_CALL(zpotri)("U", &sz, COMPLEX(ans), &sz, &info FCONE);
  if (info != 0) {
    UNPROTECT(nprot);
    if (info > 0)
      error(_("element (%d, %d) is zero, so the inverse cannot be computed"),
            info, info);
    error(_("argument %d of Lapack routine %s had invalid value"),
          -info, "dpotri");
  }
  for (int j = 0; j < sz; j++)
      for (int i = j+1; i < sz; i++)
    COMPLEX(ans)[i + j * SZ] = COMPLEX(ans)[j + i * SZ];
  UNPROTECT(nprot);
  return ans;
    }
}
