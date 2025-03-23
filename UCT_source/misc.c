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
 *   File: misc.c
 *
 *   Purpose: support functions for UCT - UNEDA core tester
 *
 */

#include "uct.h"
#include <signal.h>


int call(int rc) {

	/* The call is actually made in evaluating the argument */
	/* This procedure is the common post processing */
	if (rc) {
		printf("Command failed: %s\n",TCL_get_errtxt(rc));
		user_error++;
		}
	return rc;
	}


static int first=TRUE;

void prep_input() {
#ifndef _MSC_VER
	int c;
#endif

#ifdef _MSC_VER
	fflush(stdin);
#else
	if (first) // input not started
		first = FALSE;
	else
		while ((c = getchar()) != '\n' && c != EOF);
#endif
#ifdef MPW
	printf("\n");
#endif
	}


void read_cac(struct d_frame *df, int *alt, int *node, int mode) {

/* call     PS mode
 * ----     ----- -- 
 * P     ?alt  -> alt
 *       ?node -> node
 * V     ?alt  -> alt
 *       ?node -> node
 * rma   ?alt  -> alt
 * rmk         --    
 * adc   ?alt  -> alt
 * rmc   ?alt  -> alt
 *       ?node -> node
 */

	do {
		printf("Alternative number (1..%d): ",df->n_alts);
		prep_input();
		scanf(" %d",alt);
		} while ((*alt < 1) || (*alt > df->n_alts));
	if (mode > ALT_MODE) {
		do {
			printf("Node number (1..%d): ",df->tot_cons[*alt]);
			prep_input();
			scanf(" %d",node);
			} while ((*node < 1) || (*node > df->tot_cons[*alt]));
		}
	else
		*node = 1;
	}


void read_P_interval(double *lobo, double *upbo) {

	printf("Lower bound: ");
	prep_input();
	scanf(" %lf",lobo);
	printf("Upper bound: ");
	prep_input();
	scanf(" %lf",upbo);
	}


/* x is external value range */

double x1v(double x) {
	/* Convert x -> v for range */
	return (x - v_min) / v_mm;
	}

double x2v(double x) {
	/* Convert x -> v for diff */
	return x / v_mm;
	}

/* v is internal value range */

double v0x(double v, int n_terms) {
	/* Convert v -> x for an arbitrary number of terms (scales) 
	 * where a positive term counts as 1 and a negative as -1. */
	return v * v_mm + n_terms * v_min;
	}

double v1x(double v) {
	/* Convert v -> x for range */
	return v0x(v,1);
	}

double v2x(double v) {
	/* Convert v -> x for diff */
	return v0x(v,0);
	}


void read_V_interval(double xv(), double *lobo, double *upbo);

void read_V1_interval(double *lobo, double *upbo) {

	read_V_interval(x1v, lobo, upbo);
	}

void read_V2_interval(double *lobo, double *upbo) {

	read_V_interval(x2v, lobo, upbo);
	}

void read_V_interval(double xv(), double *lobo, double *upbo) {
	double lower, upper;

	printf("Lower bound: ");
	prep_input();
	scanf(" %lf",&lower);
	printf("Upper bound: ");
	prep_input();
	scanf(" %lf",&upper);
	*lobo = xv(lower);
	*upbo = xv(upper);
	}
