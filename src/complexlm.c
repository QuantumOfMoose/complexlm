/*
 *  MASS/src/MASS.c by W. N. Venables and B. D. Ripley  Copyright (C) 1994-2016
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
 *  A copy of the GNU General Public License is available at
 *  http://www.r-project.org/Licenses/
 */

#include <stdlib.h>
#include <math.h>
#include <stddef.h>

#include <complex.h> // This contains functions to act on complex numbers.
#include <R.h>
#include <R_ext/Applic.h>
#include <R_ext/Complex.h> // This contains the Rcomplex type. (It's a struct.)
#include <R_ext/Lapack.h>

#ifndef max
#  define max(a,b) ((a) > (b) ? (a) : (b))
#  define min(a,b) ((a) < (b) ? (a) : (b))
#endif

#define abs9(a) (a > 0 ? a:-a)

#include "R_ext/Rdynload.h"

void
zlqs_fitlots(Rcomplex *Rx, Rcomplex *Ry, int *n, int *p, int *qn,
	    int *lts, int *adj, int *sample, int *nwhich,
	    int *ntrials, Rcomplex *Rcrit, int *sing, int *bestone,
	    Rcomplex *Rbestcoef, Rcomplex *Rpk0, Rcomplex *Rbeta);   

void
zmve_fitlots(Rcomplex *x, int *n, int *p, int *qn, int *mcd,
	    int *sample, int *nwhich, int *ntrials,
	    Rcomplex *crit, int *sing, int *bestone);

static const R_CMethodDef CEntries[] = {
    {"zlqs_fitlots", (DL_FUNC) &zlqs_fitlots, 16},
    {"zmve_fitlots", (DL_FUNC) &zmve_fitlots, 11},
    {NULL, NULL, 0}
};


void R_init_complexlm(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
