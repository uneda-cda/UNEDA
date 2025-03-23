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
 *                     UNEDA Tree Core Layer (TCL)
 *                     --------------------------- 
 *
 *    +----- o o o -----------------------------------------------+
 *    |    o       o             Prof. Mats Danielson             |
 *    |   o  STHLM  o            DECIDE Research Group            |
 *    |   o         o   Dept. of Computer and Systems Sciences    |
 *    |   o   UNI   o            Stockholm University             |
 *    |    o       o     PO Box 1203, SE-164 25 Kista, SWEDEN     |
 *    +----- o o o -----------------------------------------------+
 *
 *                Copyright (c) 2012-2025 Mats Danielson
 *                     Email: mats.danielson@su.se
 *
 */

/*
 *   File: TCLevaluate.c
 *
 *   Purpose: compare expected utilities for the alternatives
 *
 *
 *   Functions exported outside TCL
 *   ------------------------------
 *   TCL_evaluate
 *   TCL_evaluate_omega
 *
 *   Functions outside module, inside TCL
 *   ------------------------------------
 *   sort_dom2
 *
 *   Functions internal to module
 *   ----------------------------
 *   mm_cmp
 *   mm_cmp_rev
 *   omega
 *   calc_omega
 *   calc_psi
 *   calc_delta
 *   calc_gamma
 *   calc_digamma
 *
 */

#include "TCLinternal.h"


 /*********************************************************
  *
  *  Configuration
  *
  *********************************************************/

/* Allowing PSI to be calculated using the digamma function */

#define DIGAMMA_PSI


 /*********************************************************
  *
  *  Sorting function for optimisation substitute function
  *
  *********************************************************/

static double *mm;

static int mm_cmp(const void *k1, const void *k2) {

	if (*(mm+(*(int*)k1)) > *(mm+(*(int*)k2)))
		return 1;
	else if (*(mm+(*(int*)k1)) < *(mm+(*(int*)k2)))
		return -1;
	else
		return 0;
	}


static int mm_cmp_rev(const void *k1, const void *k2) {

	if (*(mm+(*(int*)k1)) > *(mm+(*(int*)k2)))
		return -1;
	else if (*(mm+(*(int*)k1)) < *(mm+(*(int*)k2)))
		return 1;
	else
		return 0;
	}


void sort_dom2(i_row lin_order, d_row maxmin, int start, int stop, bool rev) {

	mm = maxmin;
	qsort(lin_order+start,stop-start+1,sizeof(int),(rev?mm_cmp_rev:mm_cmp));
	}


 /*********************************************************
  *
  *  Evaluation procedures
  *
  *********************************************************/

static d_row V_lobo,V_upbo;
static d_row P_point,im_P_point;
static d_row P_mid,V_mid;

static double omega(int Ai) {
	int i,Ai_begin,Ai_end;
	double omega;

	/* Evaluate Ai */
	mpoint_P(P_mid);
	mpoint_V(V_mid);
	Ai_begin = get_V_start(Ai);
	Ai_end = get_V_end(Ai);
	omega = 0.0;
	for (i=Ai_begin; i<=Ai_end; i++) {
		omega += P_mid[i]*V_mid[i];
		}
	return omega;
	}


static void calc_omega(int Ai, a_result result) {

	result[Ai][E_MIN] = result[Ai][E_MID] = result[Ai][E_MAX] = omega(Ai);
	}


static void calc_psi(int Ai, a_result result) {

	/* Evaluate Ai */
	fhull_V(V_lobo,V_upbo);
	result[Ai][E_MIN] = eval_P_min(Ai,0,1,V_lobo,P_point,im_P_point,TRUE);
	result[Ai][E_MID] = omega(Ai);
	result[Ai][E_MAX] = eval_P_max(Ai,0,1,V_upbo,P_point,im_P_point,TRUE);
	}


static void calc_delta(int Ai, int Aj, a_result result) {

	/* Compare Ai to Aj */
	fhull_V(V_lobo,V_upbo);
	result[Ai][E_MIN]  = eval_P_min(Ai,0,1,V_lobo,P_point,im_P_point,TRUE);
	result[Ai][E_MIN] -= eval_P_max(Aj,0,1,V_upbo,P_point,im_P_point,TRUE);
	result[Ai][E_MID]  = omega(Ai)-omega(Aj);
	result[Ai][E_MAX]  = eval_P_max(Ai,0,1,V_upbo,P_point,im_P_point,TRUE);
	result[Ai][E_MAX] -= eval_P_min(Aj,0,1,V_lobo,P_point,im_P_point,TRUE);
	}


static void calc_gamma(struct d_frame *df, int Ai, a_result result) {
	int Aj;
	double scale;

	/* Compare Ai to all other alternatives */
	scale = df->n_alts - 1.0;
	fhull_V(V_lobo,V_upbo);
	result[Ai][E_MIN] = eval_P_min(Ai,0,1,V_lobo,P_point,im_P_point,TRUE);
	result[Ai][E_MID] = omega(Ai);
	result[Ai][E_MAX] = eval_P_max(Ai,0,1,V_upbo,P_point,im_P_point,TRUE);
	for (Aj=1; Aj<=df->n_alts; Aj++)
		if (Aj != Ai) {
			result[Ai][E_MIN] -= eval_P_max(Aj,0,1,V_upbo,P_point,im_P_point,TRUE)/scale;
			result[Ai][E_MID] -= omega(Aj)/scale;
			result[Ai][E_MAX] -= eval_P_min(Aj,0,1,V_lobo,P_point,im_P_point,TRUE)/scale;
			}
	}


static void calc_digamma(struct d_frame *df, int Ai, int alts, a_result result) {
	int Aj,n_active;
	double scale;

	/* Find nbr of active alts */
	for (n_active=0, Aj=1; Aj<=df->n_alts; Aj++)
		if ((Aj!=Ai) && (alts&(0x01<<(Aj-1))))
			n_active++;
	scale = n_active;
	fhull_V(V_lobo,V_upbo);
	/* Compare Ai to subset of other alternatives */
	result[Ai][E_MIN] = eval_P_min(Ai,0,1,V_lobo,P_point,im_P_point,TRUE);
	result[Ai][E_MID] = omega(Ai);
	result[Ai][E_MAX] = eval_P_max(Ai,0,1,V_upbo,P_point,im_P_point,TRUE);
#ifdef DIGAMMA_PSI
	if (!n_active) // if called with empty alts, will return psi
		return;
#endif
	for (Aj=1; Aj<=df->n_alts; Aj++)
		if ((Aj!=Ai) && (alts&(0x01<<(Aj-1)))) {
			result[Ai][E_MIN] -= eval_P_max(Aj,0,1,V_upbo,P_point,im_P_point,TRUE)/scale;
			result[Ai][E_MID] -= omega(Aj)/scale;
			result[Ai][E_MAX] -= eval_P_min(Aj,0,1,V_lobo,P_point,im_P_point,TRUE)/scale;
			}
	}


 /*************************************************************
  *
  *  Evaluation entry points (omega, delta, gamma, digamma, and
  *  psi in one function, omega also in a separate function).
  *
  *  In reality, all functions could be calls to calc_digamma.
  *
  *************************************************************/

rcode TCL_evaluate(struct d_frame *df, int Ai, int Aj, int eval_method, a_result result) {

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;
	if ((Ai < 1) || (Ai > df->n_alts))
		return TCL_INPUT_ERROR;
	if (eval_method == DELTA) {
		if ((Aj < 1) || (Aj > df->n_alts) || (Ai == Aj))
			return TCL_INPUT_ERROR;
			}
	else if (eval_method < DIGAMMA)
		if (Aj)
			return TCL_INPUT_ERROR;

	/* OMEGA evaluation */
	if (eval_method == OMEGA)
		calc_omega(Ai,result);
	/* PSI evaluation */
	else if (eval_method == PSI)
		calc_psi(Ai,result);
	/* DELTA evaluation (Aj is an int) */
	else if (eval_method == DELTA)
		calc_delta(Ai,Aj,result);
	/* GAMMA evaluation */
	else if (eval_method == GAMMA)
		calc_gamma(df,Ai,result);
	/* DIGAMMA evaluation (Aj is a bitmask) */
	else if (eval_method == DIGAMMA) {
		/* Not ok to compare with self */
		if (0x01<<(Ai-1) & Aj)
			return TCL_INPUT_ERROR;
#ifndef DIGAMMA_PSI
		/* Not ok to compare with nothing */
		if (!((0x01<<df->n_alts)-1 & Aj))
			return TCL_INPUT_ERROR;
#endif
		calc_digamma(df,Ai,Aj,result);
		}
	else // unknown method
		return TCL_INPUT_ERROR;
	return TCL_OK;
	}


rcode TCL_evaluate_omega(struct d_frame *df, int Ai, double *result) {

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;
	if ((Ai < 1) || (Ai > df->n_alts))
		return TCL_INPUT_ERROR;
	/* Get omega EV */
	*result = omega(Ai);
	return TCL_OK;
	}
