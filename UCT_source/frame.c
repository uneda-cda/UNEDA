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
 *   File: frame.c
 *
 *   Purpose: frame handling for UCT - UNEDA core tester
 *
 */

#include "uct.h"
#include <signal.h>

static t_row tnext,tdown;

void show_index() {
	int i;
	struct d_frame *df;

	df = uf->df;
	printf("PS-%s '%s' ",uf->multilevel?"tree":"frame",uf->frame_name);
	if (folder[0])
		printf("in folder '%s' ",folder);
	else
		printf("in home folder ");
	printf("has %d alternatives\n",df->n_alts);
	for (i=1; i<=df->n_alts; i++)	{
		printf("A%d ('%s') with %d node%s ",i,uf->alt_name[i],
						df->tot_cons[i],(df->n_cons[i]==1?"":"s"));
		if (df->im_cons[i])
			printf("(%d double + %d intermediate)",
							df->n_cons[i],df->im_cons[i]);
		printf("\n");
		}
	}


void show_all() {

	/* void valued, can't stop'em */
	show_index();
	show_P_base();
	printf("Probability hull\n");
	show_P_hull();
	show_P_midpoints();
	show_V_base();
	printf("Value hull\n");
	show_V_hull();
	show_V_midpoints();
	}


void show_frame_info() {

	printf("Frame version %d.%d.%d\n",uf_ver_main,uf_ver_func,uf_ver_tech);
	show_index();
	printf("Probability base has %d constraint%s\n",
					uf->df->P_base->n_stmts,uf->df->P_base->n_stmts==1?"":"s");
	printf("Value base has %d constraint%s\n",
					uf->df->V_base->n_stmts,uf->df->V_base->n_stmts==1?"":"s");
	printf("Value range is [%.3lf %.3lf]\n",uf->v_lo,uf->v_up);
	}


void draw_tree(int alt, int snode, int level) {
	struct d_frame *df = uf->df;
	int i,tnode,is_at=level;

	for (tnode=df->down[alt][snode]; tnode; tnode=df->next[alt][tnode]) {
		for (i=is_at; i<level; i++)
			printf("     ");
		if (!is_at)
			printf(" ");
		if (df->down[alt][tnode]) {
			printf("%2d---",tnode);
			if (level<12)
				draw_tree(alt,tnode,level+1);
			else
				printf("*\n");
			}
		else
			printf("%2d\n",tnode);
		is_at = 0;
		}
	}


void show_tree_structure() {
	int i;

	if (uf->multilevel)
		for (i=1; i<=uf->df->n_alts; i++) {
			printf("Alternative %d\n",i);
			draw_tree(i,0,0);
			}
	else
		printf("Flat structure\n");
	}


void draw_drow(int alt, int snode, d_row drow, int level) {
	struct d_frame *df = uf->df;
	int i,tnode,is_at=level;

	for (tnode=df->down[alt][snode]; tnode; tnode=df->next[alt][tnode]) {
		for (i=is_at; i<level; i++)
			printf("       ");
		if (df->down[alt][tnode]) {
			printf("%.3lf--",drow[TCL_get_real_index(alt,tnode)]);
			draw_drow(alt,tnode,drow,level+1);
			}
		else
			printf("%.3lf\n",drow[TCL_get_real_index(alt,tnode)]);
		is_at = 0;
		}
	}


void show_tree_fp() {
	int i;
	d_row P_ccp,LP_ccp,V_ccp;
	struct d_frame *df;

	df = uf->df;
	if (call(TCL_get_P_masspoint(df,P_ccp,LP_ccp)))
		return;
	if (call(TCL_get_V_masspoint(df,V_ccp)))
		return;
	if (uf->multilevel) {
		printf("GLOBAL\n");
		for (i=1; i<=df->n_alts; i++) {
			printf("Alt %d\n",i);
			draw_drow(i,0,P_ccp,0);
			}
		printf("\nLOCAL\n");
		for (i=1; i<=df->n_alts; i++) {
			printf("Alt %d\n",i);
			draw_drow(i,0,LP_ccp,0);
			}
		}
	else /* flat */
		printf("Not applicable\n");
	}


void show_fp() {
	int i,j,k;
	d_row P_ccp,P_lccp,V_ccp;
	struct d_frame *df;

	df = uf->df;
	if (call(TCL_get_P_masspoint(df,P_ccp,P_lccp)))
		return;
	if (call(TCL_get_V_masspoint(df,V_ccp)))
		return;
		printf("Cons.    GP      LP       V\n");
		k = 1;
		for (i=1; i<=df->n_alts; i++) {
			for (j=1; j<=df->tot_cons[i]; j++) {
				if (V_ccp[k] > -1.0)
					printf("C%d.%d: %s%.3lf   %.3lf   %.3lf\n",i,j,j<10?" ":"",P_ccp[k],P_lccp[k],
									v1x(V_ccp[k]));
				else
					printf("C%d.%d: %s%.3lf   %.3lf   -IM-\n",i,j,j<10?" ":"",P_ccp[k],P_lccp[k]);
				k++;
				}
			}
	}


static d_row P_sd, V_sd;

void show_moments() {
	int i,j,k;
	a_row rm1,cm2,cm3;
	struct d_frame *df;

	df = uf->df;
	if (call(TCL_get_moments(df,rm1,cm2,cm3)))
		return;
	printf("Alt    RM1    CM2\n");
	for (i=1; i<=uf->df->n_alts; i++) {
		printf("A%-2d   %.3lf  %.3lf\n",i,rm1[i],cm2[i]);
		}
	printf("\nCons  Psdev  Vsdev\n");
	for (i=1,k=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			if (call(TCL_get_P_sd(df,k,P_sd+k)))
				return;
			if (call(TCL_get_V_sd(df,k,FALSE,V_sd+k)))
				return;
			if (V_sd[k] > -1.0)
				printf("C%d.%d%s %.3lf  %.3lf\n",i,j,j<10?" ":"",P_sd[k],V_sd[k]);
			else
				printf("C%d.%d%s %.3lf  -IM-\n",i,j,j<10?" ":"",P_sd[k]);
			}
	}
