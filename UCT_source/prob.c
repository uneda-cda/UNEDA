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
 *   File: ddt.c
 *
 *   Purpose: main loop for UCT - UNEDA core tester
 *
 */

#include "uct.h"
#include <signal.h>


void print_P_stmt(struct base *P, int snbr) {

	printf("%2d: P%d.%d%s= [%.3lf %.3lf]\n",
					snbr,P->stmt[snbr].alt[1],P->stmt[snbr].cons[1],
					P->stmt[snbr].cons[1]<10?" ":"",
					P->stmt[snbr].lobo,P->stmt[snbr].upbo);
	}


void make_P_midpoint() {
	struct stmt_rec stmt;

	stmt.n_terms = 1;
	read_cac(uf->df,&(stmt.alt[1]),&(stmt.cons[1]),PW_MODE);
	stmt.sign[1] = 1;
	read_P_interval(&(stmt.lobo),&(stmt.upbo));
	if (call(TCL_add_P_mstatement(uf->df,&stmt)))
		return;
	uf->dirty = TRUE;
	ask_to_save = TRUE;
	}


void remove_P_midpoint() {
	struct stmt_rec stmt;

	stmt.n_terms = 1;
	read_cac(uf->df,&(stmt.alt[1]),&(stmt.cons[1]),PW_MODE);
	stmt.sign[1] = 1;
	if (call(TCL_delete_P_mstatement(uf->df,&stmt)))
		return;
	uf->dirty = TRUE;
	ask_to_save = TRUE;
	}


void add_P_constraint() {
	struct stmt_rec stmt;

	stmt.n_terms = 1;
	read_cac(uf->df,&(stmt.alt[1]),&(stmt.cons[1]),PW_MODE);
	stmt.sign[1] = 1;
	read_P_interval(&(stmt.lobo),&(stmt.upbo));
	if (call(TCL_add_P_constraint(uf->df,&stmt)))
		return;
	uf->dirty = TRUE;
	ask_to_save = TRUE;
	}


void remove_P_constraint() {
	int pstmt;

	do {
		printf("Constraint number (1..%d): ",uf->df->P_base->n_stmts);
		prep_input();
		scanf("%d",&pstmt);
		} while ((pstmt<0) || (pstmt>uf->df->P_base->n_stmts));
	if (!pstmt)
		return;
	print_P_stmt(uf->df->P_base,pstmt);
	if (call(TCL_delete_P_constraint(uf->df,pstmt)))
		return;
	printf("    *DELETED*\n");
	uf->dirty = TRUE;
	ask_to_save = TRUE;
	}


void replace_P_constraint() {
	struct stmt_rec stmt;
	int stmt_number;
	struct base *P;

	P = uf->df->P_base;
	do {
		printf("Constraint number (1..%d): ",P->n_stmts);
		prep_input();
		scanf(" %d",&stmt_number);
		} while ((stmt_number<0) || (stmt_number>P->n_stmts));
	if (!stmt_number)
		return;
	if (P->stmt[stmt_number].n_terms != 1) {
		printf("Error in constraint\n");
		return;
		}
	print_P_stmt(P,stmt_number);
	stmt.n_terms = 1;
	read_cac(uf->df,&(stmt.alt[1]),&(stmt.cons[1]),PW_MODE);
	stmt.sign[1] = 1;
	read_P_interval(&(stmt.lobo),&(stmt.upbo));
	if (call(TCL_replace_P_constraint(uf->df,stmt_number,&stmt)))
		return;
	uf->dirty = TRUE;
	ask_to_save = TRUE;
	}


void change_P_constraint() {
	int stmt_number;
	struct base *P;
	double lobo,upbo;

	P = uf->df->P_base;
	do {
		printf("Constraint number (1..%d): ",P->n_stmts);
		prep_input();
		scanf(" %d",&stmt_number);
		} while ((stmt_number<0) || (stmt_number>P->n_stmts));
	if (!stmt_number)
		return;
	print_P_stmt(P,stmt_number);
	read_P_interval(&lobo,&upbo);
	if (call(TCL_change_P_constraint(uf->df,stmt_number,lobo,upbo)))
		return;
	uf->dirty = TRUE;
	ask_to_save = TRUE;
	}


void show_P_base() {
	int i;
	struct base *P;

	P = uf->df->P_base;
	printf("The probability base contains %d constraint%s\n",P->n_stmts,(P->n_stmts==1?"":"s"));
	for (i=1; i<=P->n_stmts; i++)
		print_P_stmt(P,i);
	}


static d_row lo_midpoint,up_midpoint;

void show_P_midpoints() {
	int i,j,k;
	bool show = FALSE;
	struct base *P;

	P = uf->df->P_base;
	if (call(TCL_get_P_mbox(uf->df,lo_midpoint,up_midpoint)))
		return;
	k = 1;
	for (i=1; i<=uf->df->n_alts; i++)
		for (j=1; j<=uf->df->tot_cons[i]; j++, k++)
			if (lo_midpoint[k] >= 0.0) {
				if (!show)
					printf("Probability midpoints\n");
				printf("P%d.%d = [%.3lf %.3lf]\n",i,j,lo_midpoint[k],up_midpoint[k]);
				show = TRUE;
				}
	if (!show)
		printf("No probability midpoints\n");
	}


void show_P_hull() {
	int i,j,k;
	d_row hlobo,hupbo,llobo,lupbo;
	struct d_frame *df;

	df = uf->df;
	if (call(TCL_get_P_hull(df,hlobo,hupbo,llobo,lupbo)))
		return;
	printf("\t  GLOBAL\t  LOCAL\n");
	for (i=1,k=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++)
			printf("P%d.%d%s= [%.3lf %.3lf]  [%.3lf,%.3lf]\n",i,j,j<10?" ":"",
						hlobo[k],hupbo[k],llobo[k],lupbo[k]);
	}


void show_P_core() {
	int i,j,k;
	d_row P_core,P_lcore;
	struct d_frame *df;

	df = uf->df;
	if (call(TCL_get_P_masspoint(df,P_core,P_lcore)))
		return;
	printf("     GLOBAL   LOCAL\n");
	k = 1;
	for (i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++) {
			printf("P%d.%d%s= %.3lf   %.3lf\n",i,j,j<10?" ":"",P_core[k],P_lcore[k]);
			k++;
			}
	}
