/*
 *
 *
 *        _/       _/   _/       _/    _/_/_/_/_/   _/_/_/          _/
 *       _/       _/   _/_/     _/    _/           _/    _/       _/  _/
 *      _/       _/   _/ _/    _/    _/           _/      _/    _/    _/
 *     _/       _/   _/  _/   _/    _/_/_/_/     _/      _/   _/      _/ 
 *    _/       _/   _/   _/  _/    _/           _/      _/   _/_/_/_/_/  
 *   _/       _/   _/    _/ _/    _/           _/      _/   _/      _/          
 *   _/     _/    _/     _/_/    _/           _/     _/    _/      _/          
 *    _/_/_/     _/       _/    _/_/_/_/_/   _/_/_/_/     _/      _/   
 *
 *
 *   UNEDA - The Universal Engine for Decision Analysis
 *
 *   Website: https://people.dsv.su.se/~mad/UNEDA
 *   GitHub:  https://github.com/uneda-cda/UNEDA
 *
 *   Licensed under CC BY 4.0: https://creativecommons.org/licenses/by/4.0/.
 *   Provided "as is", without warranty of any kind, express or implied.
 *   Reuse and modifications are encouraged, with proper attribution.
 *
 *
 *
 *                   UNEDA Decision Tree Layer (DTL)
 *                   -------------------------------
 *
 *    +----- o o o ------------------------------------------------+
 *    |    o       o              Prof. Mats Danielson             |
 *    |   o  STHLM  o             DECIDE Research Group            |
 *    |   o         o    Dept. of Computer and Systems Sciences    |
 *    |   o   UNI   o             Stockholm University             |
 *    |    o       o      PO Box 1203, SE-164 25 Kista, SWEDEN     |
 *    +----- o o o ------------------------------------------------+
 *
 *                Copyright (c) 2012-2025 Mats Danielson
 *                     Email: mats.danielson@su.se
 *
 */

/*
 *   File: DTLbnormal.c
 *
 *   Purpose: Support routines for B-normal calculations
 *
 *
 *   Functions exported outside DTL
 *   ------------------------------
 *   NONE
 *
 *   Functions outside of module, inside DTL
 *   ---------------------------------------
 *   sgn
 *   log1px
 *   n_cdf
 *   owens_t
 *   bn_cdf
 *   bn_inv_cdf
 *   b_delta
 *
 *   Functions internal to module
 *   ----------------------------
 *   f_abs
 *   ra
 *   inv_n_cdf
 *   inv_bn_cdf
 *
 */

#include "DTL.h"
#include "DTLinternal.h"


 /*********************************************************
  *
  *  Business normal (B-normal) function handling
  *
  *  DTL layer 2: all functions are below DTL API level
  *
  *********************************************************/

/* f_abs exists due to problems with the library function fabs
 * when implementing the B-normal distribution. In later times,
 * it might be possible to unify f_abs with fabs. But thanks to
 * MS VC++ we are not there at the moment. */

static double f_abs(double x) {
	double value;

	/* Why not simply define it as "return x>=0.0?x:-x;"?
	 * Does not compile correctly under MS VC++ 6.0 :-( */
	if (0.0 <= x) {
		value = x;
		}
	else {
		value = -x;
		}
	return value;
	}


double sgn(double x) {

//return x>0.0?1.0:(x<0.0?-1.0:0.0);
	return x>=0.0?1.0:-1.0;
	}


/* log1p() is missing in MS VC++ 6.0 which was released before the C99 standard */

#ifndef LOG1P
double log1px(double x) {

	// x is large enough to do standard calculation (incl. error handling for x <= -1.0)
	if (f_abs(x) > 1.0E-8)
		return log(1.0+x);
	// otherwise use a Taylor approximation: log(1+x) = x-x^2/2 with error roughly x^3/3
	// -> for |x| < 1.0E-8 we have |x|^3 < 1.0E-24 -> relative error less than 1.0E-16
	return (-0.5*x+1.0)*x;
	}
#else
/* Use library log1p */
#define log1px log1p
#endif


/* CDF of the standard normal N(0,1) distribution.
 * erf() is not available in C89 so HMF is used. */

double n_cdf(double x) {
	double a1 =  0.254829592;
	double a2 = -0.284496736;
	double a3 =  1.421413741;
	double a4 = -1.453152027;
	double a5 =  1.061405429;
	double p  =  0.3275911;
	double t,y;
	int sign = 1;

	// Handbook of Mathematical Functions by
	// Abramowitz and Stegun, formula 7.1.26
	if (x < 0)
		sign = -1;
	x = f_abs(x)/sqrt(2.0);
	t = 1.0/(1.0 + p*x);
	y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
	return 0.5*(1.0 + sign*y);
	}


/* Inverse CDF of the standard normal N(0,1) distribution */

static double ra(double t) {
	double c[] = {2.515517, 0.802853, 0.010328};
	double d[] = {1.432788, 0.189269, 0.001308};

	// Handbook of Mathematical Functions by
	// Abramowitz and Stegun, formula 26.2.23
	// |error| is < 5.0E-4
	return t - ((c[2]*t + c[1])*t + c[0]) /
						 (((d[2]*t + d[1])*t + d[0])*t + 1.0);
	}


static double inv_n_cdf(double x) {

	/* Handle round-off issues */
	if (x <= 0.0)
		return 0.0; // out of range
	if (x >= 1.0)
		return 1.0; // out of range
	/* Legit input -> process */
	if (x < 0.5)
		// F^-1(x) is -G^-1(x)
		return -ra(sqrt(-2.0*log(x)));
	else
		// F^-1(x) is G^-1(1-x)
		return ra(sqrt(-2.0*log(1-x)));
	}


/* CDF of Owen's standard T-distribution */

double owens_t(double x, double alpha) {

	double alphasq;
	int i;
	double r[5] = {
		0.1477621,
		0.1346334,
		0.1095432,
		0.0747257,
		0.0333357 };
	double r1;
	double r2;
	double rt;
	double tp = 0.159155;
	double tv1 = 1.0E-35;
	double tv2 = 15.0;
	double tv3 = 15.0;
	double tv4 = 1.0E-05;
	double u[5] = {
		0.0744372,
		0.2166977,
		0.3397048,
		0.4325317,
		0.4869533 };
	double value;
	double x1;
	double x2;
	double xs;

	// Handbook of Mathematical Functions by
	// Abramowitz and Stegun, reference 26.22
	if (f_abs(x) < tv1) {
		value = tp*atan(alpha);
		return value;
		}
	if (tv2 < f_abs(x)) {
		value = 0.0;
		return value;
		}
	if (f_abs(alpha) < tv1) {
		value = 0.0;
		return value;
		}
	xs = -0.5*x*x;
	x2 = alpha;
	alphasq = alpha*alpha;
	/* Newton iteration */
#ifdef LOG1P
	if (tv3 <= log1p(alphasq) - xs*alphasq) {
#else
	if (tv3 <= log1px(alphasq) - xs*alphasq) {
#endif
		x1 = 0.5*alpha;
		alphasq = 0.25*alphasq;
		for (;;) {
			rt = alphasq+1.0;
			x2 = x1 + (xs*alphasq+tv3-log(rt)) / (2.0*x1*(1.0/rt-xs));
			alphasq = x2*x2;
			if (f_abs(x2-x1) < tv4)
				break;
			x1 = x2;
			}
		}
	/* Gaussian quadrature */
	rt = 0.0;
	for (i=0; i<5; i++) {
		r1 = 1.0 + alphasq*pow(0.5+u[i],2.0);
		r2 = 1.0 + alphasq*pow(0.5-u[i],2.0);
		rt = rt + r[i]*(exp(xs*r1)/r1 + exp(xs*r2)/r2);
		}
	value = rt*x2*tp;
	return value;
	}


/* Inverse CDF of Owen's standard T-distribution */

#define INV_LOOPS 100

static double inv_bn_cdf(double cdf, double alpha) {
	int loop,loop2;
	double val,new_val,diff,new_cdf;

	/* Find inverse to skewed B-normal value from cdf by Newton-Raphson */
	new_val = 0.0; // start at midpoint for unskewed N(0,1)
	loop = loop2 = 0;
	do {
		val = new_val;
		do { // seek next cdf test value
			new_cdf = cdf+2.0*owens_t(val,alpha);
			if (new_cdf < 0.0) // value is too big
				val -= 0.1;
			else if (new_cdf > 1.0) // value is too small
				val += 0.1;
			loop2++; // oscillating risk
			} while (((new_cdf < 0.0) || (new_cdf > 1.0)) && (loop2 < INV_LOOPS));
		if (loop2 >= INV_LOOPS)
			goto inv_bn_cdf_error;
		new_val = inv_n_cdf(new_cdf);
		diff = new_val-val;
		new_val = (new_val+val)/2.0;
		loop++;
		} while ((fabs(diff) > DTL_EPS) && (loop < INV_LOOPS));
	if (loop < INV_LOOPS)
		return new_val;
	/* Hopeless: if no value found within INV_LOOPS, there is none to be found */
inv_bn_cdf_error:
	return inv_n_cdf(cdf); // return unskewed value as an approximation
	}


/* CDF and inverse for the B-normal (business normal) distribution. Derived from the
   skew-normal distribution using moderated skew and linear interpolation for truncated
   tails. Adapted to the corporate management case of finite business-wide risk handling. */

double bn_cdf(double val, double mean, double var, double alpha) {
	double owen,cdf;

	/* Catch pointwise mass */
	if (var < DTL_EPS)
		if (val < mean-DTL_EPS)
			return 0.0;
		else if (val > mean+DTL_EPS)
			return 1.0;
		else // at the mean
			return 0.5;
	/* Normal interval case */
	owen = owens_t((val-mean)/sqrt(var),alpha);
	cdf = n_cdf((val-mean)/sqrt(var))-2.0*owen;
	/* Catch round-off errors */
	if (cdf < 1.0E-6)
		cdf = 0.0;
	else if (cdf > 1.0-1.0E-6)
		cdf = 1.0;
	return cdf;
	}


/* Find the inverse of an unskewed or skewed B-normal CDF */

double bn_inv_cdf(double cdf, double mean, double var, double alpha) {
	double ival;

	if (alpha)
		/* Get inverse to skewed CDF */
		ival = inv_bn_cdf(cdf,alpha);
	else
		/* Get inverse to unskewed CDF */
		ival = inv_n_cdf(cdf);
	return max(min(ival*sqrt(var)+mean,1.0),-1.0);
	}


/* b_delta is the unsigned and moderated delta for the B-normal distribution.
   Code from 2012, vindicated by formulas 2.24 and 2.28 pp.30-32 in Azzalini &
   Capitanio, 2014. See the math documentation for details and constants. */

double b_delta(double skew) {
	double tau,b_skew;

	// case A: unmoderated skew
	// below 0.9
	b_skew = f_abs(skew);
	if (b_skew > 2.0)
		// case C: max skew
		// (larger than this blows alpha)
		b_skew = 0.955;
	else if (b_skew > 0.9)
		// case B: moderated skew 0.9+(x-0.9)/20.0
		// maps [0.9,2.0] -> [0.9,0.955]
		b_skew = (17.1+b_skew)/20.0;
	tau = pow(b_skew,2.0/3.0);
	return sqrt(PI*tau/(2.0*tau+DELTAPI));
	}
