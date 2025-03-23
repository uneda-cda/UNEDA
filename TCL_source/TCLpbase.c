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
 *   File: TCLpbase.c
 *
 *   Purpose: maintain the probability base
 *
 *
 *   Functions exported outside TCL
 *   ------------------------------
 *   NONE
 *
 *   Functions outside module, inside TCL
 *   ------------------------------------
 *   create_P
 *   dispose_P
 *   load_P
 *   get_P_index
 *   get_P_im_index
 *   l_hull_P
 *   hull_P
 *   l_cpoint_P
 *   cpoint_P
 *   mpoint_P (but no l_mpoint_P)
 *
 *   Functions internal to module
 *   ----------------------------
 *   check_norm
 *   loc_2_glob
 *   calc_tree_hull
 *   calc_tree_mhull
 *   renorm_mp
 *   save_mp
 *   check_mp
 *   N_DoF_mp
 *
 */

#include "TCLinternal.h"


 /*********************************************************
  *
  *  Data structures
  *
  *********************************************************/

/* Local data structures */
static d_row box_lobo;
static d_row box_upbo;
static d_row im_box_lobo;
static d_row im_box_upbo;
static d_row hull_lobo;
static d_row hull_upbo;
static d_row im_hull_lobo;
static d_row im_hull_upbo;
static d_row L_hull_lobo;
static d_row L_hull_upbo;
static d_row im_L_hull_lobo;
static d_row im_L_hull_upbo;
static d_row mass_point;
static d_row im_mass_point;
static d_row L_mass_point;
static d_row im_L_mass_point;
static d_row mbox_lobo;
static d_row mbox_upbo;
static d_row im_mbox_lobo;
static d_row im_mbox_upbo;
static d_row mhull_lobo;
static d_row mhull_upbo;
static d_row im_mhull_lobo;
static d_row im_mhull_upbo;
static d_row L_mhull_lobo;
static d_row L_mhull_upbo;
static d_row im_L_mhull_lobo;
static d_row im_L_mhull_upbo;

/* Forward declarations */
static t_matrix tnext,tprev,tdown,tup;


 /*********************************************************
  *
  *  Create/dispose P-base
  *
  *********************************************************/

rcode create_P(struct d_frame *df) {
	int i;

	/* Get new memory chunk */
	df->P_base = (struct base *)mem_alloc(sizeof(struct base),"struct base","create_P");
	if (!df->P_base)
		return TCL_OUT_OF_MEMORY;
	/* Pre-fill entries */
	df->P_base->watermark = P_MARK;
	df->P_base->n_stmts = 0;
	df->P_base->box = FALSE;
	for (i=0; i<=MAX_CONS; i++) {
		df->P_base->lo_midbox[i] = -1.0;
		df->P_base->up_midbox[i] = -1.0;
		df->P_base->lo_im_midbox[i] = -1.0;
		df->P_base->up_im_midbox[i] = -1.0;
		}
	return TCL_OK;
	}


rcode dispose_P(struct d_frame *df) {

	/* Check input parameters */
	if (df->P_base->watermark != P_MARK)
		return TCL_CORRUPTED;
	/* Prevent accidental reuse */
	df->P_base->watermark = 0;
	return mem_free((void *)df->P_base);
	}


 /*********************************************************
  *
  *  Support functions for P validation
  *
  *********************************************************/

static int check_norm(int alt) {
	double norm_sum;
	int j;

	/* Check global normalisation */
	norm_sum = 0.0;
	for (j=alt_inx[alt-1]+1; j<=alt_inx[alt]; j++)
		norm_sum += mass_point[j];
	if ((norm_sum < 1.0-EPS100) || (norm_sum > 1.0+EPS100)) {
		return TCL_INCONSISTENT;
		}
	return TCL_OK;
	}


static void loc_2_glob(int alt, int snode, double norm, d_row P_loc, d_row im_P_loc, d_row P_glob, d_row im_P_glob) {
	int inx,tnode;

	for (tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode])
		if (tdown[alt][tnode]) {
			inx = at2i(alt,tnode);
			/* Im-node */
			im_P_glob[inx] = norm * im_P_loc[inx];
			/* Globalise children */
			loc_2_glob(alt,tnode,im_P_glob[inx],P_loc,im_P_loc,P_glob,im_P_glob);
			}
		else {
			inx = at2r(alt,tnode);
			/* Re-node */
			P_glob[inx] = norm * P_loc[inx];
			}
	}


static rcode calc_tree_hull(int alt, int snode, double p_lobo, double p_upbo) {
	int tnode,inx;
	double pmin=0.0,pmax=0.0;

	/* Step 1: Calculate min and max prob for this level */
	for (tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode]) {
		if (tdown[alt][tnode]) {
			/* Im-node */
			inx = at2i(alt,tnode);
			pmin += im_box_lobo[inx];
			pmax += im_box_upbo[inx];
			}
		else {
			inx = at2r(alt,tnode);
			pmin += box_lobo[inx];
			pmax += box_upbo[inx];
			}
		}
	/* Step 2: Check consistency */
	if ((pmin > 1.0+EPS) || (pmax < 1.0-EPS))
		return TCL_INCONSISTENT;
	pmin = min(pmin,1.0);
	pmax = max(pmax,1.0);
	/* Step 3: Calculate hull */
	for (tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode]) {
		if (tdown[alt][tnode]) {
			/* Im-node */
			inx = at2i(alt,tnode);
			im_L_hull_lobo[inx] = max(im_box_lobo[inx],im_box_upbo[inx]+1.0-pmax);
			im_L_hull_upbo[inx] = min(im_box_upbo[inx],im_box_lobo[inx]+1.0-pmin);
			im_hull_lobo[inx] = im_L_hull_lobo[inx]*p_lobo;
			im_hull_upbo[inx] = im_L_hull_upbo[inx]*p_upbo;
			if (calc_tree_hull(alt,tnode,im_hull_lobo[inx],im_hull_upbo[inx]))
				return TCL_INCONSISTENT;
			}
		else	{
			/* Re-node */
			inx = at2r(alt,tnode);
			L_hull_lobo[inx] = max(box_lobo[inx],box_upbo[inx]+1.0-pmax);
			L_hull_upbo[inx] = min(box_upbo[inx],box_lobo[inx]+1.0-pmin);
			hull_lobo[inx] = L_hull_lobo[inx]*p_lobo;
			hull_upbo[inx] = L_hull_upbo[inx]*p_upbo;
			}
		}
	return TCL_OK;
	}


static rcode calc_tree_mhull(int alt, int snode, double p_lobo, double p_upbo) {
	int tnode,inx;
	double pmin=0.0,pmax=0.0;

	/* Step 1: Calculate min and max prob for this level */
	for (tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode]) {
		if (tdown[alt][tnode]) {
			/* Im-node */
			inx = at2i(alt,tnode);
			pmin += im_mbox_lobo[inx];
			pmax += im_mbox_upbo[inx];
			}
		else {
			inx = at2r(alt,tnode);
			pmin += mbox_lobo[inx];
			pmax += mbox_upbo[inx];
			}
		}
	/* Step 2: Check consistency */
	if ((pmin > 1.0+EPS) || (pmax < 1.0-EPS))
		return TCL_INCONSISTENT;
	pmin = min(pmin,1.0);
	pmax = max(pmax,1.0);
	/* Step 3: Calculate hull */
	for (tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode]) {
		if (tdown[alt][tnode]) {
			/* Im-node */
			inx = at2i(alt,tnode);
			im_L_mhull_lobo[inx] = max(im_mbox_lobo[inx],im_mbox_upbo[inx]+1.0-pmax);
			im_L_mhull_upbo[inx] = min(im_mbox_upbo[inx],im_mbox_lobo[inx]+1.0-pmin);
			im_mhull_lobo[inx] = im_L_mhull_lobo[inx]*p_lobo;
			im_mhull_upbo[inx] = im_L_mhull_upbo[inx]*p_upbo;
			if (calc_tree_mhull(alt,tnode,im_mhull_lobo[inx],im_mhull_upbo[inx]))
				return TCL_INCONSISTENT;
			}
		else	{
			/* Re-node */
			inx = at2r(alt,tnode);
			L_mhull_lobo[inx] = max(mbox_lobo[inx],mbox_upbo[inx]+1.0-pmax);
			L_mhull_upbo[inx] = min(mbox_upbo[inx],mbox_lobo[inx]+1.0-pmin);
			mhull_lobo[inx] = L_mhull_lobo[inx]*p_lobo;
			mhull_upbo[inx] = L_mhull_upbo[inx]*p_upbo;
			}
		}
	return TCL_OK;
	}


 /*********************************************************
  *
  *  Mass point handling
  *
  *********************************************************/

static void renorm_mp(int alt, int snode) {
	int tnode,inx;
	double cc_sum,cc_dist,frac;

	cc_sum = 0.0;
	for (tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode])
		cc_sum += (tdown[alt][tnode]?im_L_mass_point[at2i(alt,tnode)]:L_mass_point[at2r(alt,tnode)]);
	cc_dist = 0.0;
	if (cc_sum < 1.0-EPS) {
		for (tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode]) {
			if (tdown[alt][tnode]) {
				inx = at2i(alt,tnode);
				cc_dist += im_L_mhull_upbo[inx]-im_L_mass_point[inx];
				}
			else {
				inx = at2r(alt,tnode);
				cc_dist += L_mhull_upbo[inx]-L_mass_point[inx];
				}
			}
		frac = (1.0-cc_sum) / cc_dist;
		for (tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode]) {
			if (tdown[alt][tnode]) {
				inx = at2i(alt,tnode);
				im_L_mass_point[inx] += frac * (im_L_mhull_upbo[inx]-im_L_mass_point[inx]);
				}
			else {
				inx = at2r(alt,tnode);
				L_mass_point[inx] += frac * (L_mhull_upbo[inx]-L_mass_point[inx]);
				}
			}
		}
	else if (cc_sum > 1.0+EPS) {
		for (tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode]) {
			if (tdown[alt][tnode]) {
				inx = at2i(alt,tnode);
				cc_dist += im_L_mass_point[inx]-im_L_mhull_lobo[inx];
				}
			else {
				inx = at2r(alt,tnode);
				cc_dist += L_mass_point[inx]-L_mhull_lobo[inx];
				}
			}
		frac = (1.0-cc_sum) / cc_dist;
		for (tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode]) {
			if (tdown[alt][tnode]) {
				inx = at2i(alt,tnode);
				im_L_mass_point[inx] += frac * (im_L_mass_point[inx]-im_L_mhull_lobo[inx]);
				}
			else {
				inx = at2r(alt,tnode);
				L_mass_point[inx] += frac * (L_mass_point[inx]-L_mhull_lobo[inx]);
				}
			}
		}
	for (tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode])
		if (tdown[alt][tnode]) {
			/* Distribute mass over children */
			renorm_mp(alt,tnode);
			}
	}


#include "TCLwarp.c"

/* The basic algorithm has N degrees of freedom (DoF). Latest research shows that
   it is impossible to know whether the decision-maker has an N or N-1 DoF model
   in his/her head -> function adjust_vx combines it with an N-1 DoF model. */

static void N_DoF_mp(int alt, int snode, double norm) {
	int tnode,inx;
	double pmin,pmax,lofrac,upfrac;

	/* Collect at current level */
	pmin = pmax = 0.0;
	for (tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode]) {
		pmin += (tdown[alt][tnode]?im_L_mhull_lobo[at2i(alt,tnode)]:L_mhull_lobo[at2r(alt,tnode)]);
		pmax += (tdown[alt][tnode]?im_L_mhull_upbo[at2i(alt,tnode)]:L_mhull_upbo[at2r(alt,tnode)]);
		}
	if (pmin >= 1.0) {
		lofrac = 1.0;
		upfrac = 0.0;
		}
	else if (pmax <= 1.0) {
		lofrac = 0.0;
		upfrac = 1.0;
		}
	else if (pmax > pmin+EPS) {
		lofrac = (pmax-1.0) / (pmax-pmin);
		upfrac = 1.0 - lofrac;
		}
	else
		lofrac = upfrac = 0.5;

	/* Distribute over current level */
	for (tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode])
		if (tdown[alt][tnode]) {
			inx = at2i(alt,tnode);
			/* Im-node: Convex combination of upper and lower midpoint bounds */
			im_L_mass_point[inx] = lofrac * im_L_mhull_lobo[inx] + upfrac * im_L_mhull_upbo[inx];
			}
		else {
			inx = at2r(alt,tnode);
			/* Re-node: Convex combination of upper and lower midpoint bounds */
			L_mass_point[inx] = lofrac * L_mhull_lobo[inx] + upfrac * L_mhull_upbo[inx];
			}
	/* Combine N-1 DoF model with N DoF model */
	adjust_vx(alt,snode,lofrac);

	/* Distribute over children */
	for (tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode])
		if (tdown[alt][tnode]) {
			inx = at2i(alt,tnode);
			im_mass_point[inx] = norm * im_L_mass_point[inx];
			N_DoF_mp(alt,tnode,im_mass_point[inx]);
			}
		else {
			inx = at2r(alt,tnode);
			mass_point[inx] = norm * L_mass_point[inx];
			}
	}


 /*********************************************************
  *
  *  Load and verify P-base
  *
  *********************************************************/

rcode load_P(struct d_frame *df) {
	int i,j,alt,tcons,var_nbr;
	struct base *P;

	/* Stage 1: Box formation. Local data only. */

	/* Check input parameters */
	if (df->P_base->watermark != P_MARK)
		return TCL_CORRUPTED;
	P = df->P_base;
	if (P->box) {
		/* User supplied ranges */
		memcpy(box_lobo,P->box_lobo,(n_vars+1)*sizeof(double));
		memcpy(box_upbo,P->box_upbo,(n_vars+1)*sizeof(double));
		memcpy(im_box_lobo,P->im_box_lobo,(im_vars+1)*sizeof(double));
		memcpy(im_box_upbo,P->im_box_upbo,(im_vars+1)*sizeof(double));
		}
	else {
		/* Init all ranges to [0.0,1.0] */
		for (i=1; i<=n_vars; i++) {
			box_lobo[i] = 0.0;
			box_upbo[i] = 1.0;
			}
		for (i=1; i<=im_vars; i++) {
			im_box_lobo[i] = 0.0;
			im_box_upbo[i] = 1.0;
			}
		}
	/* Copy df data into local store */
	for (i=1; i<=n_alts; i++)
		for (j=0; j<=df->tot_cons[i]; j++) {
			tnext[i][j] = df->next[i][j];
			tprev[i][j] = df->prev[i][j];
			tdown[i][j] = df->down[i][j];
			tup[i][j] = df->up[i][j];
			}
	/* Check input parameters for each statement */
	for (i=1; i<=P->n_stmts; i++) {
		if (P->stmt[i].n_terms != 1)
			return TCL_INPUT_ERROR;
		if (P->stmt[i].lobo < 0.0)
			return TCL_INPUT_ERROR;
		if (P->stmt[i].upbo < P->stmt[i].lobo)
			return TCL_INPUT_ERROR;
		if (P->stmt[i].upbo > 1.0)
			return TCL_INPUT_ERROR;
		alt = P->stmt[i].alt[1];
		tcons = P->stmt[i].cons[1];
		if ((alt < 1) || (alt > n_alts) || (tcons < 1) || (tcons > df->tot_cons[alt]))
			return TCL_INPUT_ERROR;
		if (P->stmt[i].sign[1] != 1)
			return TCL_INPUT_ERROR;
		if (t2r[alt][tcons]) {
			/* Real consequence */
			var_nbr = at2r(alt,tcons);
			box_lobo[var_nbr] = max(box_lobo[var_nbr],P->stmt[i].lobo);
			box_upbo[var_nbr] = min(box_upbo[var_nbr],P->stmt[i].upbo);
			if (box_upbo[var_nbr]-box_lobo[var_nbr] < 0.0) // -EPS ??
				return TCL_INCONSISTENT;
#ifdef NO_ZERO_INTERVALS
			if (box_upbo[var_nbr]-box_lobo[var_nbr] < MIN_WIDTH)
				return TCL_TOO_NARROW_STMT;
#endif
			}
		else {
			/* Intermediate consequence */
			var_nbr = at2i(alt,tcons);
			im_box_lobo[var_nbr] = max(im_box_lobo[var_nbr],P->stmt[i].lobo);
			im_box_upbo[var_nbr] = min(im_box_upbo[var_nbr],P->stmt[i].upbo);
			if (im_box_upbo[var_nbr]-im_box_lobo[var_nbr] < 0.0)
				return TCL_INCONSISTENT;
#ifdef NO_ZERO_INTERVALS
			if (im_box_upbo[var_nbr]-im_box_lobo[var_nbr] < MIN_WIDTH)
				return TCL_TOO_NARROW_STMT;
#endif
			}
		}

	/* Stage 2: Consistency checks and hull formation.
		 Local input transformed into local (and global) hull. */

	/* Calculate tree hull for each alternative */
	for (i=1; i<=n_alts; i++)
		if (calc_tree_hull(i,0,1.0,1.0))
			return TCL_INCONSISTENT;

	/* Load midbox */
	for (i=1; i<=n_alts; i++) {
		/* Load real (end node) midbox */
		for (j=alt_inx[i-1]+1; j<=alt_inx[i]; j++) {
			if (P->lo_midbox[j] >= 0.0) {
				/* Check midbox consistency */
				if ((P->lo_midbox[j] < L_hull_lobo[j]-EPS) ||
						(P->up_midbox[j] > L_hull_upbo[j]+EPS) ||
						(P->lo_midbox[j] > P->up_midbox[j])) {
					return TCL_INCONSISTENT;
					}
				mbox_lobo[j] = P->lo_midbox[j];
				mbox_upbo[j] = P->up_midbox[j];
				}
			else /* No midbox for this variable, use hull */ {
				mbox_lobo[j] = L_hull_lobo[j];
				mbox_upbo[j] = L_hull_upbo[j];
				}
			}
		/* Load intermediate midbox */
		for (j=im_alt_inx[i-1]+1; j<=im_alt_inx[i]; j++) {
			if (P->lo_im_midbox[j] >= 0.0) {
				/* Check midpoint consistency */
				if ((P->lo_im_midbox[j] < im_L_hull_lobo[j]-EPS) ||
						(P->up_im_midbox[j] > im_L_hull_upbo[j]+EPS) ||
						(P->lo_im_midbox[j] > P->up_im_midbox[j])) {
					return TCL_INCONSISTENT;
					}
				im_mbox_lobo[j] = P->lo_im_midbox[j];
				im_mbox_upbo[j] = P->up_im_midbox[j];
				}
			else /* No midbox for this variable, use hull */ {
				im_mbox_lobo[j] = im_L_hull_lobo[j];
				im_mbox_upbo[j] = im_L_hull_upbo[j];
				}
			}
		}

	/* Calculate tree mhull for each alternative */
	for (i=1; i<=n_alts; i++)
		if (calc_tree_mhull(i,0,1.0,1.0))
			return TCL_INCONSISTENT;

	/* Stage 3: Mass point distribution (all input is local) */

	/* Distibute mass point */
	for (i=1; i<=n_alts; i++) {
		N_DoF_mp(i,0,1.0);
		/* Check global normalisation */
		if (check_norm(i))
			return TCL_INCONSISTENT;
		}

	return TCL_OK;
	}


 /*********************************************************
  *
  *  Access operations for user interface
  *
  *********************************************************/

/* Get B2(r&i) index from A1(a,t) */

int get_P_index(int alt, int cons) {

	/* Check input parameters */
	if ((alt<1) || (alt>n_alts))
		return 0;
	if ((cons<1) || (cons>tot_alt_inx[alt]-tot_alt_inx[alt-1]))
		return 0;
	/* Return B2 real index (or 0 if intermediate node) */
	if (t2r[alt][cons])
		return at2r(alt,cons);
	else
		return 0;
	}


int get_P_im_index(int alt, int cons) {

	/* Check input parameters */
	if ((alt<1) || (alt>n_alts))
		return 0;
	if ((cons<1) || (cons>tot_alt_inx[alt]-tot_alt_inx[alt-1]))
		return 0;
	/* Return B2 intermediate index (or 0 if real node) */
	if (t2i[alt][cons])
		return at2i(alt,cons);
	else
		return 0;
	}


/* Deliver B1 indexed lobo & upbo from B2 indexed data */

void l_hull_P(d_row hlobo, d_row hupbo) {
	int i;

	for (i=1; i<=tot_vars; i++)
		if (f2r[i]) {
			hlobo[i] = L_hull_lobo[f2r[i]];
			hupbo[i] = L_hull_upbo[f2r[i]];
			}
		else {
			hlobo[i] = im_L_hull_lobo[f2i[i]];
			hupbo[i] = im_L_hull_upbo[f2i[i]];
			}
	}


void hull_P(d_row hlobo, d_row hupbo) {
	int i;

	for (i=1; i<=tot_vars; i++)
		if (f2r[i]) {
			hlobo[i] = hull_lobo[f2r[i]];
			hupbo[i] = hull_upbo[f2r[i]];
			}
		else {
			hlobo[i] = im_hull_lobo[f2i[i]];
			hupbo[i] = im_hull_upbo[f2i[i]];
			}
	}


/* Local prob (B2) to external (B1) indexing */

void l_cpoint_P(d_row mid) {
	int i;

	for (i=1; i<=tot_vars; i++)
		if (f2r[i])
			mid[i] = L_mass_point[f2r[i]];
		else
			mid[i] = im_L_mass_point[f2i[i]];
	}


/* Global prob (B2) to external (B1) indexing */

void cpoint_P(d_row mid) {
	int i;

	for (i=1; i<=tot_vars; i++)
		if (f2r[i])
			mid[i] = mass_point[f2r[i]];
		else
			mid[i] = im_mass_point[f2i[i]];
	}


/* Global prob (B2), keep internal indexing (B2) */

void mpoint_P(d_row masspt) {
	int i;

	for (i=1; i<=n_vars; i++)
		masspt[i] = mass_point[i];
	}

#include "TCLevalp.c"
