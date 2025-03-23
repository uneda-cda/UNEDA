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
 *   File: eval.c
 *
 *   Purpose: evaluation for UCT - UNEDA core tester
 *
 */

#include "uct.h"
#include <signal.h>
#include <time.h>

void rescue() {
	char yn;

	if (uf->dirty && ask_to_save) {
		if (folder[0])
			printf("Save frame '%s' in folder '%s' (y/n): ",uf->frame_name,folder);
		else
			printf("Save frame '%s' in home folder (y/n): ",uf->frame_name);
		prep_input();
		scanf(" %c",&yn);
		if (tolower(yn) == 'y') {
			if (call(write_ufile(uf->frame_name,u_folder,uf)))
				return;
			uf->dirty = FALSE;
			}
		ask_to_save = FALSE;
		}
	}


static a_result cube0[MAX_ALTS+1],cube1;

void compare_delta() {
	int Ai,Aj;
	struct d_frame *df;

	df = uf->df;
	for (Ai=1; Ai<df->n_alts; Ai++)
		for (Aj=Ai+1; Aj<=df->n_alts; Aj++)
			if (call(TCL_evaluate(df,Ai,Aj,DELTA,cube0[Aj])))
				return;
	printf("         min      mid      max\n");
	for (Ai=1; Ai<df->n_alts; Ai++)
		for (Aj=Ai+1; Aj<=df->n_alts; Aj++)
			printf("E%d-E%d%8.3lf %8.3lf %8.3lf\n",Ai,Aj,
					v2x(cube0[Aj][Ai][E_MIN]),v2x(cube0[Aj][Ai][E_MID]),v2x(cube0[Aj][Ai][E_MAX]));
	}


void compare_gamma() {
	int Ai;
	struct d_frame *df;

	df = uf->df;
	for (Ai=1; Ai<=df->n_alts; Ai++)
		if (call(TCL_evaluate(df,Ai,0,GAMMA,cube1)))
			return;
	printf("         min      mid      max\n");
	for (Ai=1; Ai<=df->n_alts; Ai++)
		printf("E%d   %8.3lf %8.3lf %8.3lf\n",Ai,
				v2x(cube1[Ai][E_MIN]),v2x(cube1[Ai][E_MID]),v2x(cube1[Ai][E_MAX]));
	}


void compare_psi() {
	int Ai;
	struct d_frame *df;

	df = uf->df;
	for (Ai=1; Ai<=df->n_alts; Ai++)
		if (call(TCL_evaluate(df,Ai,0,PSI,cube1)))
			return;
	printf("         min      mid      max\n");
	for (Ai=1; Ai<=df->n_alts; Ai++)
		printf("E%d   %8.3lf %8.3lf %8.3lf\n",Ai,
				v2x(cube1[Ai][E_MIN]),v2x(cube1[Ai][E_MID]),v2x(cube1[Ai][E_MAX]));
	}


void compare_omega() {
	int Ai;
	struct d_frame *df;

	df = uf->df;
	for (Ai=1; Ai<=df->n_alts; Ai++)
		if (call(TCL_evaluate_omega(df,Ai,cube1[Ai]+E_MID)))
			return;
	for (Ai=1; Ai<=df->n_alts; Ai++)
		printf("E%d: %8.3lf\n",Ai,v1x(cube1[Ai][E_MID]));
	}


/* See (Danielson, 1997) for explanations of strong, marked, weak */

void security_level() {
	int Ai;
	struct d_frame *df;
	double min_value;
	a_vector strong,marked,weak;
	double xmin;

	df = uf->df;
	printf("Minimum acceptable value: ");
	prep_input();
	scanf("%lf",&xmin);
	min_value = x1v(xmin);
	if (call(TCL_security_level(uf->df,min_value,strong,marked,weak)))
		return;
	printf("         min      mid      max\n");
	for (Ai=1; Ai<=df->n_alts; Ai++) {
		printf("A%d: ",Ai);
		if (strong[Ai] > EPS)
			printf(" %8.3lf",strong[Ai]);
		else
			printf("   -ok-");
		if (marked[Ai] > EPS)
			printf(" %8.3lf",marked[Ai]);
		else
			printf("   -ok-");
		if (weak[Ai] > EPS)
			printf(" %8.3lf\n",weak[Ai]);
		else
			printf("   -ok-\n");
		}
	}

