/*
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
 *   Copyright (c) 2012-2025  Prof. Mats Danielson, Stockholm University
 *
 *   Website: https://people.dsv.su.se/~mad/UNEDA
 *   GitHub:  https://github.com/uneda-cda/UNEDA
 *
 *   Licensed under CC BY 4.0: https://creativecommons.org/licenses/by/4.0/.
 *   Provided "as is", without warranty of any kind, express or implied.
 *   Reuse and modifications are encouraged, with proper attribution.
 *
 */

/*
 *   File: random.c
 *
 *   Purpose: UCT random procedures
 *
 */

#include "uct.h"

/* Random generators - better than built-in */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836

static long iseed;


void init_random() {

	srand((unsigned int)time(NULL));
	iseed = (long)time(NULL);
	}


double xrandom() {
	long k;
	double rnd;

	/* Returns uniform [0,1] but never exactly 0 or 1 */
	k = iseed/IQ;
	iseed = IA*(iseed-k*IQ) - IR*k;
	if (iseed < 0)
		iseed += IM;
	rnd = AM *iseed;
	return rnd;
	}


double drandom(double lobo, double upbo) {

	/* Returns equally distributed double btw lobo and upbo exclusive */
	return xrandom()*(upbo-lobo)+lobo;
	}


int irandom(int lobo, int upbo) {

	/* Returns equally distributed int btw lobo and upbo inclusive */
	return (int)(xrandom()*(upbo-lobo+1)+lobo);
	}


double rand_td(double lobo, double top, double upbo) {
	double rnd,x;

	/* Triangle distributed random numbers */
	if (upbo-lobo < EPS)
		rnd = upbo;
	else {
		/* rand() is a lousy linear congruential generator */
		x = lobo + rand()/(RAND_MAX+1.0)*(upbo-lobo);
		if (x < top)
			rnd = lobo + sqrt((double)((x-lobo)*(top-lobo)));
		else
			rnd = upbo - sqrt((double)((upbo-x)*(upbo-top)));
		}
	return rnd;
	}


double random_td(double lobo, double top, double upbo) {
	double rnd,x;

	/* Triangle distributed random numbers */
	if (upbo-lobo < EPS)
		rnd = upbo;
	else {
		x = lobo + xrandom()*(upbo-lobo);
		if (x < top)
			rnd = lobo + sqrt((double)((x-lobo)*(top-lobo)));
		else
			rnd = upbo - sqrt((double)((upbo-x)*(upbo-top)));
		}
	return rnd;
	}
