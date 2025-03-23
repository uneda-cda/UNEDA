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
 *   File: TCLinternal.h
 *
 *   Purpose: the internal header file for TCL. Not for export.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <setjmp.h>
#include "DTLparameters.h"
#include "TCL.h"

/* Compatible with DMC/TDL or not */
#define TDL_COMPAT

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#ifndef max
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#endif

#define EPS 1.0E-8
#define EPS100 1.0E-6

/* Watermarks */
#define D_MARK 0xC572
#define P_MARK 0x6A1E
#define V_MARK 0x94BD

/* Max folder name */
#define FOLDER_SIZE 224

/* TCLpbase.c */
rcode create_P(struct d_frame *df);
rcode dispose_P(struct d_frame *df);
rcode load_P(struct d_frame *df);
int get_P_index(int alt, int cons);
int get_P_im_index(int alt, int cons);
void hull_P(d_row hlobo, d_row hupbo);
void l_hull_P(d_row hlobo, d_row hupbo);
void cpoint_P(d_row mid);
void l_cpoint_P(d_row mid);
void mpoint_P(d_row masspt);
int get_P_max(d_row objective, d_row maxpoint);
int get_TP_max(d_row objective, d_row im_objective, d_row maxpoint, d_row im_maxpoint);
double eval_P_max(int alt, int snode, int k_start, 
				d_row V_pt, d_row P_pt, d_row im_P_pt, bool positive);
double eval_P_min(int alt, int snode, int k_start, 
				d_row V_pt, d_row P_pt, d_row im_P_pt, bool positive);

/* TCLvbase.c */
rcode create_V(struct d_frame *df);
rcode dispose_V(struct d_frame *df);
rcode load_V(struct d_frame *df);
int get_V_index(int alt, int cons);
int get_V_start(int alt);
int get_V_end(int alt);
int get_im_V_start(int alt);
int get_im_V_end(int alt);
void hull_V(d_row lobo, d_row upbo);
void fhull_V(d_row lobo, d_row upbo);
void cpoint_V(d_row mid);
void mpoint_V(d_row masspt);

/* TCLevaluate.c */
void sort_dom2(i_row lin_order, d_row maxmin, int start, int stop, bool rev);

/* TCLevalp.c */
double ixset_P_max(int alt, int snode, i_row ixset);
double ixset_P_min(int alt, int snode, i_row ixset);

/* Globals */
extern t_matrix t2f,t2r,t2i,r2t,i2t;
extern tn_row f2r,f2i,r2f,i2f,i2end;
extern int n_alts;
extern int alt_inx[];
extern int n_vars;
extern int im_alt_inx[];
extern int im_vars;
extern int tot_alt_inx[];
extern int tot_vars;

/* Index conversions between modes A1(t), A2(r&i), B1(f), B2(r&i) */

/* Within mode A (=within alt) */
/* A1 -> A2: t2r[a][t], t2i[a][t] */
/* A2 -> A1: r2t[a][r], i2t[a][i] */

/* Within mode B (=over all alts) */
/* B1 -> B2: f2r[f], f2i[f] */
/* B2 -> B1: r2f[r], i2f[i] */

/* From mode A to B */
/* A1 -> B1: at2f(a,t) */
#define at2f(alt,tc) tot_alt_inx[alt-1]+tc
/* A1 -> B2: at2r(a,t), at2i(a,t) */
#define at2r(alt,tc) alt_inx[alt-1]+t2r[alt][tc]
#define at2i(alt,tc) im_alt_inx[alt-1]+t2i[alt][tc]

/* From mode B to A */
/* B2 -> A1: r2at(a,r), i2at(a,i) */
#define r2at(alt,rc) r2t[alt][rc-alt_inx[alt-1]]
#define i2at(alt,ic) i2t[alt][ic-im_alt_inx[alt-1]]
