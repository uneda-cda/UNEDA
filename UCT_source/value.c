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
 *   File: value.c
 *
 *   Purpose: value handling for UCT - UNEDA core tester
 *
 */

#include "uct.h"
#include <signal.h>


void print_V_stmt(struct base *V, int snbr) {

	printf("%2d: V%d.%d%s= [%.3lf %.3lf]\n",
					snbr,V->stmt[snbr].alt[1],V->stmt[snbr].cons[1],
					V->stmt[snbr].cons[1]<10?" ":"",
					v1x(V->stmt[snbr].lobo),v1x(V->stmt[snbr].upbo));
	}


void make_V_midpoint() {
	struct stmt_rec stmt;

	stmt.n_terms = 1;
	read_cac(uf->df,&(stmt.alt[1]),&(stmt.cons[1]),V_MODE);
	stmt.sign[1] = 1;
	read_V1_interval(&(stmt.lobo),&(stmt.upbo));
	if (call(TCL_add_V_mstatement(uf->df,&stmt)))
		return;
	uf->dirty = TRUE;
	ask_to_save = TRUE;
	}


void remove_V_midpoint() {
	struct stmt_rec stmt;

	stmt.n_terms = 1;
	read_cac(uf->df,&(stmt.alt[1]),&(stmt.cons[1]),V_MODE);
	stmt.sign[1] = 1;
	if (call(TCL_delete_V_mstatement(uf->df,&stmt)))
		return;
	uf->dirty = TRUE;
	ask_to_save = TRUE;
	}


void add_V_constraint() {
	struct stmt_rec stmt;

	stmt.n_terms = 1;
	read_cac(uf->df,&(stmt.alt[1]),&(stmt.cons[1]),V_MODE);
	stmt.sign[1] = 1;
	read_V1_interval(&(stmt.lobo),&(stmt.upbo));
	if (call(TCL_add_V_constraint(uf->df,&stmt)))
		return;
	uf->dirty = TRUE;
	ask_to_save = TRUE;
	}


void remove_V_constraint() {
	int vstmt;

	do {
		printf("Constraint number (1..%d): ",uf->df->V_base->n_stmts);
		prep_input();
		scanf(" %d",&vstmt);
		} while ((vstmt<0) || (vstmt>uf->df->V_base->n_stmts));
	if (!vstmt)
		return;
	print_V_stmt(uf->df->V_base,vstmt);
	if (call(TCL_delete_V_constraint(uf->df,vstmt)))
		return;
	printf("    *DELETED*\n");
	uf->dirty = TRUE;
	ask_to_save = TRUE;
	}


void replace_V_constraint() {
	struct stmt_rec stmt;
	int stmt_number;
	struct base *V;

	V = uf->df->V_base;
	do {
		printf("Constraint number (1..%d): ",V->n_stmts);
		prep_input();
		scanf(" %d",&stmt_number);
		} while ((stmt_number<0) || (stmt_number>V->n_stmts));
	if (!stmt_number)
		return;
	if (V->stmt[stmt_number].n_terms != 1) {
		printf("Error in constraint\n");
		return;
		}
	print_V_stmt(V,stmt_number);
	stmt.n_terms = 1;
	read_cac(uf->df,&(stmt.alt[1]),&(stmt.cons[1]),V_MODE);
	stmt.sign[1] = 1;
	read_V1_interval(&(stmt.lobo),&(stmt.upbo));
	if (call(TCL_replace_V_constraint(uf->df,stmt_number,&stmt)))
		return;
	uf->dirty = TRUE;
	ask_to_save = TRUE;
	}


void change_V_constraint() {
	int stmt_number;
	struct base *V;
	double lobo,upbo;

	V = uf->df->V_base;
	do {
		printf("Constraint number (1..%d): ",V->n_stmts);
		prep_input();
		scanf(" %d",&stmt_number);
		} while ((stmt_number<0) || (stmt_number>V->n_stmts));
	if (!stmt_number)
		return;
	print_V_stmt(V,stmt_number);
	if (V->stmt[stmt_number].n_terms == 1)
		read_V1_interval(&lobo,&upbo);
	else
		read_V2_interval(&lobo,&upbo);
	if (call(TCL_change_V_constraint(uf->df,stmt_number,lobo,upbo)))
		return;
	uf->dirty = TRUE;
	ask_to_save = TRUE;
	}


void show_V_base() {
	int i;
	struct base *V;

	V = uf->df->V_base;
	printf("The value base contains %d constraint%s\n",V->n_stmts,(V->n_stmts==1?"":"s"));
	for (i=1; i<=V->n_stmts; i++)
		print_V_stmt(V,i);
	}


static d_row lo_midpoint,up_midpoint;

void show_V_midpoints() {
	int i,j,k;
	bool show = FALSE;
	struct base *V;

	V = uf->df->V_base;
	if (call(TCL_get_V_mbox(uf->df,lo_midpoint,up_midpoint)))
		return;
	k = 1;
	for (i=1; i<=uf->df->n_alts; i++)
		for (j=1; j<=uf->df->tot_cons[i]; j++, k++)
			if (lo_midpoint[k] > -1.0) {
				if (!show)
					printf("Value midpoints\n");
				printf("V%d.%d = [%.3lf %.3lf]\n",i,j,
								v1x(lo_midpoint[k]),v1x(up_midpoint[k]));
				show = TRUE;
				}
	if (!show)
		printf("No value midpoints\n");
	}


void show_V_hull() {
	int i,j,k,kr;
	d_row lobo,upbo;
	d_row lobo2,upbo2;
	struct d_frame *df;

	df = uf->df;
	if (call(TCL_get_V_hull(df,lobo,upbo)))
		return;
	if (call(TCL_get_V_hull(df,lobo2,upbo2)))
		return;
	k = 1;
	kr = 1;
	for (i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++) {
			if (lobo[k] > -1.0) {
				printf("V%d.%d%s= [%.3lf %.3lf]\n",i,j,j<10?" ":"",v1x(lobo[k]),v1x(upbo[k]));
				kr++;
				}
			else
				printf("V%d.%d%s= -IM-\n",i,j,j<10?" ":"");
			k++;
			}
	}


void show_V_core() {
	int i,j,k;
	d_row V_core;
	struct d_frame *df;

	df = uf->df;
	if (call(TCL_get_V_masspoint(df,V_core)))
		return;
	k = 1;
	for (i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++) {
			if (V_core[k] > -1.0)
				printf("V%d.%d%s= %.3lf\n",i,j,j<10?" ":"",v1x(V_core[k]));
			else
				printf("V%d.%d%s= -IM-\n",i,j,j<10?" ":"");
			k++;
			}
	}
