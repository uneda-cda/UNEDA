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
 *   File: TCLvbase.c
 *
 *   Purpose: maintain the value base
 *
 *
 *   Functions exported outside TCL
 *   ------------------------------
 *   NONE
 *
 *   Functions outside module, inside TCL
 *   ------------------------------------
 *   create_V
 *   dispose_V
 *   load_V
 *   get_V_start
 *   get_V_end
 *   get_V_index
 *   hull_V
 *   fhull_V
 *   cpoint_V
 *   mpoint_V
 *
 *   Functions internal to module
 *   ----------------------------
 *   calc_V_hull
 *
 */

#include "TCLinternal.h"


 /*********************************************************
  *
  *  Data structures
  *
  *********************************************************/

/* Local structures */
static d_row box_lobo;
static d_row box_upbo;
static d_row hull_lobo;
static d_row hull_upbo;
static d_row mbox_lobo;
static d_row mbox_upbo;
static d_row mass_point;

static rcode calc_V_hull(struct base *V);


 /*********************************************************
  *
  *  Create/dispose value base
  *
  *********************************************************/

rcode create_V(struct d_frame *df) {
	int i;

	/* Get new memory chunk */
	df->V_base = (struct base *)mem_alloc(sizeof(struct base),"struct base","create_V");
	if (!df->V_base)
		return TCL_OUT_OF_MEMORY;
	/* Pre-fill entries */
	df->V_base->watermark = V_MARK;
	df->V_base->n_stmts = 0;
	df->V_base->box = FALSE;
	for (i=0; i<=MAX_CONS; i++) {
		df->V_base->lo_midbox[i] = -1.0;
		df->V_base->up_midbox[i] = -1.0;
		}
	return TCL_OK;
	}


rcode dispose_V(struct d_frame *df) {

	/* Check input parameters */
	if (df->V_base->watermark != V_MARK)
		return TCL_CORRUPTED;
	/* Prevent accidental reuse */
	df->V_base->watermark = 0;
	return mem_free((void *)df->V_base);
	}


 /*********************************************************
  *
  *  Load and verify value base
  *
  *********************************************************/

rcode load_V(struct d_frame *df) {
	rcode rc;
	int i,alt,cons,tcons,var_nbr;
	struct base *V;

	/* Check input parameters */
	if (df->V_base->watermark != V_MARK)
		return TCL_CORRUPTED;
	V = df->V_base;
	if (V->box) {
		/* User supplied ranges */
		memcpy(box_lobo,V->box_lobo,(n_vars+1)*sizeof(double));
		memcpy(box_upbo,V->box_upbo,(n_vars+1)*sizeof(double));
		}
	else {
		/* Init all ranges to [0.0,1.0] */
		for (i=1; i<=n_vars; i++) {
			box_lobo[i] = 0.0;
			box_upbo[i] = 1.0;
			}
		}
	/* Check input parameters for each statement */
	for (i=1; i<=V->n_stmts; i++) {
		if (V->stmt[i].n_terms != 1)
			return TCL_INPUT_ERROR;
		if (V->stmt[i].lobo < 0.0)
			return TCL_INPUT_ERROR;
		if (V->stmt[i].upbo < V->stmt[i].lobo)
			return TCL_INPUT_ERROR;
		if (V->stmt[i].upbo > 1.0)
			return TCL_INPUT_ERROR;
		alt = V->stmt[i].alt[1];
		tcons = V->stmt[i].cons[1];
		if ((alt < 1) || (alt > n_alts))
			return TCL_INPUT_ERROR;
		if ((tcons < 1) || (tcons > df->tot_cons[alt]))
			return TCL_INPUT_ERROR;
		if (V->stmt[i].sign[1] != 1)
			return TCL_INPUT_ERROR;
		cons = t2r[alt][V->stmt[i].cons[1]];
		/* Im-node not allowed */
		if (!cons)
			return TCL_ILLEGAL_NODE;
		/* Enter into box */
		var_nbr = alt_inx[alt-1] + cons;
		box_lobo[var_nbr] = max(box_lobo[var_nbr],V->stmt[i].lobo);
		box_upbo[var_nbr] = min(box_upbo[var_nbr],V->stmt[i].upbo);
		if (box_lobo[var_nbr] > box_upbo[var_nbr])
			return TCL_INCONSISTENT;
#ifdef NO_ZERO_INTERVALS
		if (box_upbo[var_nbr]-box_lobo[var_nbr] < MIN_WIDTH)
			return TCL_TOO_NARROW_STMT;
#endif
		}
	rc = calc_V_hull(V);
	return rc;
	}


static rcode calc_V_hull(struct base *V) {
	int i,j;

	/* 1. Initialize hull_upbo & hull_lobo */
	for (j=1; j<=n_vars; j++) {
		hull_upbo[j] = box_upbo[j];
		hull_lobo[j] = box_lobo[j];
		if (hull_lobo[j] > hull_upbo[j])
			return TCL_INCONSISTENT;
		}

	/* 2. Check midbox consistency */
	for (i=1; i<=n_alts; i++)
		for (j=alt_inx[i-1]+1; j<=alt_inx[i]; j++)
			if (V->lo_midbox[j] >= 0.0) {
				if ((V->lo_midbox[j] < hull_lobo[j]-EPS) ||
						(V->up_midbox[j] > hull_upbo[j]+EPS) ||
						(V->lo_midbox[j] > V->up_midbox[j]))
					return TCL_INCONSISTENT;
				mbox_lobo[j] = V->lo_midbox[j];
				mbox_upbo[j] = V->up_midbox[j];
				}
			else /* No midbox for this variable, use hull */ {
				mbox_lobo[j] = hull_lobo[j];
				mbox_upbo[j] = hull_upbo[j];
				}
	/* 3. Allocate mass point */
	for (j=1; j<=n_vars; j++) {
		/* Symmetric trapezoid/triangle */
		mass_point[j] = (mbox_lobo[j] + mbox_upbo[j]) / 2.0;
		}
	return TCL_OK;
	}


 /*********************************************************
  *
  *  Access operations for evaluation
  *
  *********************************************************/

int get_V_start(int alt) {

	/* Check input parameters */
	if ((alt < 1) || (alt > n_alts))
		return 0;
	/* Return start index */
	return alt_inx[alt-1]+1;
	}


int get_V_end(int alt) {

	/* Check input parameters */
	if ((alt < 1) || (alt > n_alts))
		return 0;
	/* Return end index */
	return alt_inx[alt];
	}


 /*********************************************************
  *
  *  Access operations for internal TCL use
  *
  *********************************************************/

/* Get B2(r&i) index from A1(a,t) */

int get_V_index(int alt, int cons) {

	/* Check input parameters */
	if ((alt<1) || (alt>n_alts))
		return 0;
	if ((cons<1) || (cons>tot_alt_inx[alt]-tot_alt_inx[alt-1]))
		return 0;
	/* Return B2 index if real value node */
	if (t2r[alt][cons])
		return at2r(alt,cons);
	else
		return 0;
	}


/* Deliver B1 indexed lobo&upbo from B2 indexed data */

void hull_V(d_row lobo, d_row upbo) {
	int i;

	for (i=1; i<=tot_vars; i++) {
		if (f2r[i])	{
			lobo[i] = hull_lobo[f2r[i]];
			upbo[i] = hull_upbo[f2r[i]];
			}
		else {
			lobo[i] = -1.0;
			upbo[i] = -1.0;
			}
		}
	}


/* Internal B2 indexing only */

void fhull_V(d_row lobo, d_row upbo) {
	int i;

	for (i=1; i<=n_vars; i++) {
		lobo[i] = hull_lobo[i];
		upbo[i] = hull_upbo[i];
		}
	}


/* B2 -> B1 external indexing */

void cpoint_V(d_row mid) {
	int i;

	for (i=1; i<=tot_vars; i++)
		if (f2r[i])
			mid[i] = mass_point[f2r[i]];
		else
			mid[i] = -1.0;
	}


/* Internal B2 indexing */

void mpoint_V(d_row masspt) {
	int i;

	for (i=1; i<=n_vars; i++)
		masspt[i] = mass_point[i];
	}
