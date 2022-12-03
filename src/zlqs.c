/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998-2016	B. D. Ripley
 *  Copyright (C) 1999          R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 or 3 of the License
 *  (at your option).
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available at
 *  http://www.r-project.org/Licenses/
 *
 *
 * Exports
 *	zlqs_fitlots(...)
 *	zmve_fitlots(...)
 *
 * to be called as  .C(.)  in ../R/zlqs.R
 *
 * Modified by William Ryan to handle complex number data.
 * Note: In many cases, I didn't know if it would be better to define variables or pointers, so I tried to imitate the original code.
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

#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>	/* for the QR	  routines */
#include <R_ext/Utils.h>	/* for the *sort() routines */
#include <R_ext/Arith.h>
#include <R_ext/Altrep.h>	/* for R_xlen_t */
#define BIG DBL_MAX

/* GLOBAL Variables, explicitly allocated and freed: They would obviously need to be changed to complex.*/
static Rcomplex *coef, *qraux, *work, *res, *yr, *xr, *means, *d2, *d2copy, *rwork, *lwork, *rb ;
static double complex *ccoef, *cqraux, *cres, *y, *x, *cyr, *cxr, *cmeans, *cd2, *cd2copy;
static int *pivot, *which, *which2;
static int *ind;
static double *s;

/* This is probably not the correct way to do it, but try defining duplicate global variables for the Rcomplex ones above, which are the same, but represented as C99 double complex. */
static double complex *ccoef, *cqraux, *cwork, *cres, *cyr, *cxr, *cd2, *cd2copy, *cb;

/* A function that allocates memoy and */

/* A function to create a matching array of Rcomplex from a double complex array of length n. It does not allocate memory, or create a variable, it just copies over values. */
static void to_Rcmplx_array(double complex *carray, int n, Rcomplex *Rarray)
{
   int i;
   for (i = 0; i <= n; i++) {
       Rarray[i].r = creal(carray[i]);
       Rarray[i].i = cimag(carray[i]);
   }
}

/* A function to creat ea matcing array of double complex from a Rcomplex array of length n. It does not create the array or allocate memory, it just copies the values. */
static void to_dblcmplx_array (Rcomplex *Rarray, int n, double complex *carray)
{
    int i;
    for (i = 0; i <= n; i++) {
        carray[i] = toC99(&Rarray[i]);
    }
}

static void lqs_setup(int *n, int *p, int *nwhich)
{
    coef = (Rcomplex *) R_alloc(*p, sizeof(Rcomplex));
    qraux = (Rcomplex *) R_alloc(*p, sizeof(Rcomplex));
    //work = (Rcomplex *) R_alloc(2*(*p), sizeof(Rcomplex));
    res = (Rcomplex *) R_alloc(*n, sizeof(Rcomplex));
    yr = (Rcomplex *) R_alloc(*n, sizeof(Rcomplex));
    xr = (Rcomplex *) R_alloc((*n)*(*p), sizeof(Rcomplex));
    rb = (Rcomplex *) R_alloc((*n)*(*p), sizeof(Rcomplex));
    pivot = (int *) R_alloc(*p, sizeof(int));
    ind = (int *) R_alloc(*n, sizeof(int));
    which = (int *) R_alloc(*nwhich, sizeof(int));
    s = (double *) R_alloc((*nwhich > *p) ? *p : *nwhich, sizeof(double));
    /*bestone = R_alloc(*nwhich, int);*/
    ccoef = (double complex *) R_alloc(*p, sizeof(double complex));
    cqraux = (double complex *) R_alloc(*p, sizeof(double complex));
    //cwork = (double complex *) R_alloc(*p, sizeof(double complex));
    cres = (double complex *) R_alloc(*p, sizeof(double complex));
    y = (double complex *) R_alloc(*p, sizeof(double complex));
    x = (double complex *) R_alloc((*n)*(*p), sizeof(double complex));
    cyr = (double complex *) R_alloc(*p, sizeof(double complex));
    cxr = (double complex *) R_alloc((*n)*(*p), sizeof(double complex));
    cd2 = (double complex *) R_alloc(*p, sizeof(double complex));
    cd2copy = (double complex *) R_alloc(*p, sizeof(double complex));
}


/*
   Sampling k from 0:n-1 without replacement. -> No change needed.
 */
static void sample_noreplace(int *x, int n, int k)
{
    int i, j, nn=n;

    for (i = 0; i < n; i++) ind[i] = i;
    for (i = 0; i < k; i++) {
	j = (int)(nn * unif_rand());
	x[i] = ind[j];
	ind[j] = ind[--nn];
    }
}

/*
   Find all subsets of size k in order: this gets a new one each call. -> No change needed.
 */
static void next_set(int *x, int n, int k)
{
    int i, j, tmp;

    j = k - 1;
    tmp = x[j]++;
    while(j > 0 && x[j] >= n - (k - 1 -j)) tmp = ++x[--j];
    for(i = j+1; i < k; i++)  x[i] =  ++tmp;
}


/*
   Adjust the constant for an LMS fit. This is the midpoint of the
   qn contiguous observations of shortest length. -> Do this whole thing with double complex.
 */
static double complex lmsadj(double complex *x, int n, int qn, double complex *ssbest)
{
    int i, k = qn - 1;
    //Rcomplex adj, RRssbest; Leftover from trying to juggle double complex and Rcomplex.
    double complex len, best, cadj;
    //double complex x = toC99(Rx);
    //double complex ssbest = toC99(Rssbest);

    best = x[k] - x[0];
    cadj = 0.5*(x[k] + x[0]);
    for(i = 1; i < n-k; i++){
	len = x[i+k] - x[i];
	if(cabs(len) < cabs(best)) {
	    best = len;
	    cadj = 0.5*(x[i+k] + x[i]);
	}
    }
    *ssbest = 0.25*best*best;
    return(cadj);
}

/*
   Adjust the constant for an LTS fit. This is the mean of the
   qn contiguous observations of smallest variance. -> Convert Rcomplex to C99 double complex and back. -> Didn't work. Do the whole thing with double complex.
 */
static double complex ltsadj(double complex *x, int n, int qn, double complex *ssbest)
{
    int i, k = qn - 1;
    //Rcomplex adj;
    double complex ss, best, m1, m2, cadj;
    //double complex *x = toC99(Rx);
    //double complex *ssbest = toC99(Rssbest);

    /*printf("qn = %d\n", qn);*/
    m1 = m2 = 0.0;
    for(i=0; i < qn; i++) {
	m1 += x[i];
	m2 += x[i]*x[i];
    }
    cadj = m1/qn;
    best = m2 - m1*m1/qn;

    for(i = 1; i < n-k; i++){
	m1 += x[i+k] - x[i-1];
	m2 += x[i+k]*x[i+k] - x[i-1]*x[i-1];
	ss = m2 - m1*m1/qn;
	/* if(ss < best) is undefined for complex numbers. x is a sorted array of (squared? No, they don't appear to be.) residuals, so ss will be positive as well. But sorting complex numbers unambiguously is impossible...*/
	if(cabs(ss) < cabs(best)) {
	    best = ss;
	    cadj = m1/qn;
	}
    }
    *ssbest = best;
    //*Rssbest->r = creal(ssbest);
    //*Rssbest->i = cimag(ssbest);
    //adj->r = creal(cadj);
    //adj->i = cimag(cadj);
    return(cadj);
}

/* the chi function for the S estimator: the integral of biweight. -> Convert Rcomplex to C99 double complex and back. -> Decided to do it all in double complex */
static double complex chi(double complex x, double complex a)
{
    //double complex x = toC99(Rx);
    //double complex a = toC99(Ra);
    double complex chi;
    //Rcomplex Rchi;
    
    x /= a; x *= x;
    if(cabs(x) > 1) return(1.0); /* Use the real part of x. After all, if the real part is > 1, the abs will be too. Will return(1.0) work? 1.0 is a float... Changed from return(1.0), I have no idea if this will work... Maybe check if cabs(x) > 1, and if it is set chi = 1/sqrt(2) + 1/sqrt(2)*I... Since chi() will usually be dealing with doubles, I'll leave this as return(1.0).*/
    else return(x*(3 + x*(-3 + x)));
    
    //Rchi->r = creal(chi);
    //Rchi->i = cimag(chi);
    //return(chi)
}

/*
   For lots of subsets of size *nwhich, compute the exact fit to those
   data points and the residuals from all the data points. -> Convert lots of Rcomplex to C99 double complex, then back for returning.
   -> Finish changing to zgelsd. -> Based on La_qr_cmplx() from Lapack.c in r-source.
 */
void
zlqs_fitlots(Rcomplex *Rx, Rcomplex *Ry, int *n, int *p, int *qn,
	    int *lts, int *adj, int *sample, int *nwhich,
	    int *ntrials, Rcomplex *Rcrit, int *sing, int *bestone,
	    Rcomplex *Rbestcoef, Rcomplex *Rpk0, Rcomplex *Rbeta)
{
    int nnew = *nwhich, pp = *p;
    int i, iter, j, k,  nn = *n, this, trial, lwork, *iwork, itmp;
    int rank, info, n100 = 100;
    int firsttrial = 1;
    int nrhs = 1;
    double tol = 1.0e-7, *rwork, rtmp;
    Rcomplex *work, tmp;
    double complex a = 0.0, sum, target,
	old, new, dummy, k0 = toC99(Rpk0), beta = toC99(Rbeta), best = BIG + BIG*I, *bestcoef,
	crit = toC99(Rcrit), thiscrit; /* I hope that works for setting the Rcomplex values. */

    lqs_setup(n, p, nwhich);
    /* Copy values from Rcomplx arrays to double complex arrays. */
    to_dblcmplx_array(Rx, nn*pp, x);
    to_dblcmplx_array(Ry, nn, y);    

    *sing = 0;
    target = (nn - pp)* (beta); // This is fine, beta is not an array.

    if(!(*sample)) {
	for(i = 0; i < nnew; i++) which[i] = i;
    } else GetRNGstate();

    for(trial = 0; trial < *ntrials; trial++) {

	R_CheckUserInterrupt();

	if(!(*sample)) {if(trial > 0) next_set(which, nn, nnew);}
	else sample_noreplace(which, nn, nnew);

	/* The following loop formats b (aka y) and A (aka x) before passing them to zgelsd. It's a relic from using dqrdc2 and dqrsl. Would it be better to pass them directly? Probably not. I think this takes a random sample of the data in b (aka y) and A (aka x). Copy yr into rb, since zgels will replace it with the solution x (aka b) ((That's some confusing notation...)).*/
	for(j = 0; j < nnew; j++) {
	    this = which[j];
	    yr[j] = Ry[this];
	    rb[j] = yr[j];
	    for(k = 0; k < pp; k++) xr[j + nnew*k] = Rx[this + nn*k];
	}
	
	/* Workspace query. Get optimal value of lwork, and minimum values of rwork and iwork. Then allocate space for work, rwork, and iwork. */
	lwork = -1;
	F77_CALL(zgelsd)(&nnew, &pp, &nrhs, xr, &nnew, rb, &nnew, s, &tol, &rank, &tmp, &lwork, &rtmp, &itmp, &info);
	//if (info != 0)
	//	error(_("error code %d from Lapack routine '%s'"), &info, "zgelsd");
	lwork = (int) tmp.r;
	work = (Rcomplex *) R_alloc(lwork, sizeof(Rcomplex));
	rwork = (double *) R_alloc(rtmp, sizeof(double));
	iwork = (int *) R_alloc(itmp, sizeof(int));
	/* compute fit, find residuals */
	F77_CALL(zgelsd)(&nnew, &pp, &nrhs, xr, &nnew, rb, &nnew, s, &tol, &rank, work, &lwork, rwork, iwork, &info);
	//if(rank < pp) { (*sing)++; continue; }
	/* F77_CALL(zqrsl)(xr, &nnew, &nnew, &rank, qraux, yr, &dummy, yr,
			coef, &dummy, &dummy, &n100, &info); zgelsd does the work of both dqrdc2 and dqrsl, though in a different way. */
/* the only thing that gets used later from thes FORTRAN calls is coef, as far as I can tell. */
	coef = rb;
	to_dblcmplx_array(coef, pp, ccoef);
	/* This loop calculates the residuals.*/
	for(i = 0; i < nn; i++) {
	    sum = y[i];
	    for(j = 0; j < rank; j++) sum -= ccoef[j] * x[i + nn*j];
	    cres[i] = sum;
	    res[i].r = creal(sum);
	    res[i].i = cimag(sum);
	}

	if(*lts < 2) {/* lqs or lts estimation */
	    /* find the constant subtracted from the residuals that minimizes
	       the criterion. As this is a univariate problem, has an exact
	       solution.  */
	    if(*adj) {
		R_csort(res, nn); /* R sorts complex numbers by the real component first, then the imaginary. I'm not sure how well that will work... There is a R_csort and cPsort in R_ext/Utils.h! Let's use those! */
		to_dblcmplx_array(res, nn, cres);
		if(*lts) a = ltsadj(cres, nn, *qn, &thiscrit);
		else a = lmsadj(cres, nn, *qn, &thiscrit);
	    } else {
		for(i = 0; i < nn; i++) {
		    sum = cres[i] - a;
		    cres[i] = sum*sum; /* Should this be changed to conjf(sum)*sum? I don't think so, since res are complex.*/
		    res[i].r = creal(cres[i]);
		    res[i].i = cimag(cres[i]);
		}
		cPsort(res, nn, *qn-1); /* partial sort - Unknown if it can handle complex numbers. -> cPsort can. Since it's an R function, can I just pass Rcomplex?*/
		to_dblcmplx_array(res, nn, cres);
		if(!(*lts)) thiscrit = cres[*qn-1];
		else {
		    sum = 0.0;
		    for(i = 0; i < *qn; i++) sum += cres[i];
		    thiscrit = sum;
		}
	    }
	} else { /* S estimation */
	    if(firsttrial) {
		for(i = 0; i < nn; i ++) {
		    cres[i] = cabs(cres[i]);
		    res[i].r = creal(cres[i]);
		    res[i].i = cimag(cres[i]);
		}
		cPsort(res, nn, nn/2); /* Another rPsort, changed to cPsort. */
		to_dblcmplx_array(res, nn, cres);
		old = cres[nn/2]/0.6745;	 /* MAD provides the initial scale */
		firsttrial = 0;
	    } else {
		/* only find optimal scale if it will be better than
		   existing best solution */
		sum = 0.0;
		for(i = 0; i < nn; i ++) sum += chi(cres[i], k0 * best);
		if(cabs(sum) > cabs(target)) continue; /* A comparison. since complex numbers are not ordered, this is not defined. We'll compare the absolute values of sum and target. I'm not sure if cabs() is the right function here, hard to tell from the R source.*/
		old = best;
	    } /* now solve for scale S by re-substitution */
	    for(iter = 0; iter < 30; iter++) {
		/*printf("iter %d, s = %f sum = %f %f\n", iter, old, sum, target);*/
		sum = 0.0;
		for(i = 0; i < nn; i ++) sum += chi(cres[i], k0 * old);
		new = csqrt(sum/target) * old; /* Not sure if it's necessary to change sqrt to csqrt. */
		if(fabs(cabs(sum/target) - 1.) < 1e-4) break; //The specific modifications to this line are unjustified. Someting needed changing, but this may not be the proper way to do it.
		old = new;
	    }
	    thiscrit = new;
	}

	if(cabs(thiscrit) < cabs(best)) {  /* first trial might be singular, so use fence. Again, we'll compare absolute values, since complex numbers are not an ordered field. */
	    sum = 0.0;
	    for(i = 0; i < nn; i ++) sum += chi(cres[i], k0 * best);
	    best = thiscrit;
	    /* printf("trial %d, best = %f sum = %f %f\n", trial, best, sum, target);*/
	    for(i = 0; i < nnew; i++) bestone[i] = which[i] + 1;
	    bestcoef = (double complex *) R_alloc(pp, sizeof(double complex)); //Allocate memory for double complex bestcoef.
	    for(i = 0; i < pp; i++) bestcoef[i] = toC99(&coef[i]);
	    bestcoef[0] += a;
	    for(i = 0; i < pp; i++) {
	    	Rbestcoef[i].r = creal(bestcoef[i]);
	    	Rbestcoef[i].i = cimag(bestcoef[i]);
	    }
	}
    }/* for(trial in 0:ntrials) */

    crit = (creal(best) < 0.0 && cimag(best) < 0.0) ? 0.0 : best; /* I'm not sure about this line, but it seemed to be the best way to adapt the real version to complex. */
    Rcrit->r = creal(crit);
    Rcrit->i = cimag(crit);
    if(*sample) PutRNGstate();
    /* lqs_free(); */
}


static void mve_setup(int *n, int *p, int *ps)
{
    xr = (Rcomplex *) R_alloc((*ps)*(*p), sizeof(Rcomplex));
    qraux = (Rcomplex *) R_alloc(*p, sizeof(Rcomplex));
    pivot = (int *) R_alloc(*p, sizeof(int));
    //work = (Rcomplex *) R_alloc(2*(*p), sizeof(Rcomplex)); Create when needed.
    d2 = (Rcomplex *) R_alloc(*n, sizeof(Rcomplex));
    d2copy = (Rcomplex *) R_alloc(*n, sizeof(Rcomplex));
    means = (Rcomplex *) R_alloc((*p), sizeof(Rcomplex));
    ind = (int *) R_alloc(*n, sizeof(int));
    which = (int *) R_alloc(*ps, sizeof(int));
    which2 = (int *) R_alloc(*ps, sizeof(int));
    cxr = (double complex *) R_alloc((*ps)*(*p), sizeof(double complex));
    x = (double complex *) R_alloc((*ps)*(*p), sizeof(double complex));
    cqraux = (double complex *) R_alloc(*p, sizeof(double complex));
    //cwork = (double complex *) R_alloc(2*(*p), sizeof(double complex));
    cd2 = (double complex *) R_alloc(*n, sizeof(double complex));
    cd2copy = (double complex *) R_alloc(*n, sizeof(double complex));
    cmeans = (double complex *) R_alloc((*p), sizeof(double complex));
}


/* find the squared Mahalanobis distance to x via QR decomposition in xr. */
static Rcomplex mah(Rcomplex *xr, int nnew, int p, Rcomplex *x)
{
    int i, j;
    //Rcomplex s, ss = {.r = 0.0, .i = 0.0};
    Rcomplex mahh;
    double complex s, ss = 0.0;

    for(j = 0; j < p; j++) {
	s = toC99(&x[j]);
	if(j > 0) for(i = 0; i < j; i++) s -= toC99(&work[i]) * toC99(&xr[i + nnew*j]);
	cwork[j] = s / toC99(&xr[j + nnew*j]);
	work[j].r = creal(cwork[j]);
	work[i].i = cimag(cwork[j]);
	ss += cwork[j] * cwork[j];
    }
    mahh.r = creal(ss*(nnew-1));
    mahh.i = cimag(ss*(nnew-1));
    return(mahh);
}

/*
   Compute the squared Mahalanobis distances, in d2, to all points in x
   from the mean of the subset in which using the covariance of that
   subset. -> I think I could just do this in R...
*/
static int do_one(Rcomplex *Rx, int *which, int n, int nnew, int p,
       double complex *det, Rcomplex *d2)
{
    int i, j, k, lwork;
    int rank;
    int info;
    double tol = 1.0e-7, *rwork, rtmp;
    Rcomplex *work, tmp;
    double complex sum;
    //cxr = toC99(xr);
    to_dblcmplx_array(Rx, n*p, x);

    for(j = 0; j < nnew; j++)
	for(k = 0; k < p; k++) cxr[j + nnew*k] = x[which[j] + n*k];
    for(k = 0; k < p; k++) {
	sum = 0.0;
	for(j = 0; j < nnew; j++) sum += cxr[j + nnew*k];
	sum /= nnew;
	cmeans[k] = sum;
	for(j = 0; j < nnew; j++) cxr[j + nnew*k] -= sum;
    }
    to_Rcmplx_array(cxr, n*p, xr);

    /* First carry out a workspace query, then use the info from that to call zgeqp3 optimally. */
    lwork = -1;
    F77_CALL(zgeqp3)(&nnew, &p, xr, &nnew, pivot, qraux, &tmp, &lwork, &rtmp, &info);
    //if (info != 0)
	//	error(_("error code %d from Lapack routine '%s'"), &info, "zgeqp3");
    lwork = (int) tmp.r;
    rwork = (double *) R_alloc(2*p, sizeof(double));
    work = (Rcomplex *) R_alloc(lwork, sizeof(Rcomplex));   
    F77_CALL(zgeqp3)(&nnew, &p, xr, &nnew, pivot, qraux, work, &lwork, rwork, &info);
    to_dblcmplx_array(qraux, p, cqraux);
    to_dblcmplx_array(xr, n*p, cxr);
    if(rank < p) return(1);

    sum = 0.0;
    for(k = 0; k < p; k++)
	sum += log(fabs(cxr[k + nnew*k]));
    *det = sum;

    /* now solve R^T b = (x[i, ] - means) and find squared length of b */
    for(i = 0; i < n; i++) {
	for(j = 0; j < p; j++) cqraux[j] = x[i + n*j] - cmeans[j];
	to_Rcmplx_array(cqraux, p, qraux);
	d2[i] = mah(xr, nnew, p, qraux);
    }
    return(0);
}


void
zmve_fitlots(Rcomplex *x, int *n, int *p, int *qn, int *mcd,
	    int *sample, int *nwhich, int *ntrials,
	    Rcomplex *crit, int *sing, int *bestone)
{
    int i, iter, j, nn = *n, quan = *qn, trial, this_sing;
    int nnew = *nwhich;
    Rcomplex thiscrit, lim, Rdet;
    double complex det, best = BIG + BIG*I, *cx;

    if(*mcd != 1)
	mve_setup(n, p, nwhich);
    else
	mve_setup(n, p, n); /* could get ties */

    *sing = 0;
    if(!*sample) {
	for(i = 0; i < nnew; i++) which[i] = i;
    } else GetRNGstate();

    thiscrit.r = 0.0;
    thiscrit.i = 0.0;		/* -Wall */

    for(trial = 0; trial < *ntrials; trial++) {

	R_CheckUserInterrupt();

	if(!(*sample)) {if(trial > 0) next_set(which, nn, nnew);}
	else sample_noreplace(which, nn, nnew);

	/* for(i = 0; i < nnew; i++) printf("%d ", which[i]); printf("\n");
	   fflush(stdout);*/


	/* Find the mean and covariance matrix of the sample. Check if singular.
	   Compute Mahalanobis distances of all points to the means using
	   this covariance matrix V, and find quantile largest. Volume is
	   then monotone in determinant(V * dist2). */

	this_sing = do_one(x, which, nn, nnew, *p, &det, d2);
	if(this_sing)  {(*sing)++; continue;}

	/*for(i = 0; i < nnew; i++) printf(" %d", which[i]); printf("\n");*/

	for(i = 0; i < nn; i++) d2copy[i] = d2[i];
	cPsort(d2copy, nn, quan-1);
	lim = d2copy[*qn-1];
	if(!*mcd) {
		thiscrit.r = creal((*p) * log(toC99(&lim)) + 2*det);
		thiscrit.i = cimag((*p) * log(toC99(&lim)) + 2*det);
	}
	else {
	    for(iter = 0; iter < 4; iter++) {
		/*for(i = 0; i < nn; i++) printf(" %f", d2[i]); printf("\n");*/
		if(iter > 0) {
		    for(i = 0; i < nn; i++) d2copy[i] = d2[i];
		    cPsort(d2copy, nn, quan-1);
		    lim = d2copy[*qn-1];
		}
		j = 0;
		for(i = 0; i < nn; i++)
		    if(cabs(toC99(&d2[i])) <= cabs(toC99(&lim))) which2[j++] = i;
		/* note: we take all points that meet this limit:
		   there could be more than quan. */
		(void) do_one(x, which2, nn, quan, *p, &det, d2);
		if(iter > 0 && 2 * cabs(det) >= 0.999*cabs(toC99(&thiscrit))) break;
		thiscrit.r = creal(2 * det);
		thiscrit.i = cimag(2 * det);
		/* printf("iter %d %f", iter, thiscrit);
		   for(i = 0; i < quan; i++) printf(" %d", which2[i]);
		   printf("\n"); fflush(stdout);*/
	    }

	}
	/*   printf("this %f\n", thiscrit);*/


	if(cabs(toC99(&thiscrit)) < cabs(best)) { /* warning: first might be singular */
	    best = toC99(&thiscrit);
	    for(i = 0; i < nn; i++) bestone[i] = (cabs(toC99(&d2[i])) <= cabs(toC99(&lim)));
	}
    }
    crit->r = creal(best);
    crit->i = cimag(best);
    if(*sample) PutRNGstate();
    /* mve_free(); */
}

