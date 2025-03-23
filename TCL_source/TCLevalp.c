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
 *   File: TCLevalp.c
 *
 *   Purpose: evaluation in P-base
 *
 *   Algorithms described in US patent 7257566 (expired June 6, 2025)
 *
 *
 *   Functions exported outside TCL
 *   ------------------------------
 *   TCL_get_P_max
 *   TCL_get_P_min
 *
 *   Functions outside module, inside TCL
 *   ------------------------------------
 *   eval_P_max
 *   eval_P_min
 *
 *   Functions internal to module
 *   ----------------------------
 *   calc_EV
 *   eval_P1_max
 *   eval_P1_min
 *
 */


 /*********************************************************
  *
  *  Access operations for evaluation
  *
  *********************************************************/

static d_row local_p_lobo,local_p_upbo,local_v,p_max;
static d_row im_local_v;
static d_row im_V_pt;
static i_row order;

static double calc_EV(int alt, int snode, d_row V_pt, d_row im_V_pt, d_row P_pt, d_row im_P_pt) {
	int inx,tnode;
	double EV;

	/* Calculate EV from local tree node data */
	EV = 0.0;
	/* Sum all nodes on this level */
	for (tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode]) {
		/* Create local v-order */
		if (tdown[alt][tnode]) {
			/* Im-node */
			inx = at2i(alt,tnode);
			im_V_pt[inx] = calc_EV(alt,tnode,V_pt,im_V_pt,P_pt,im_P_pt);
			EV += im_P_pt[inx] * im_V_pt[inx];
			}
		else {
			/* Re-node */
			inx = at2r(alt,tnode);
			EV += P_pt[inx] * V_pt[inx];
			}
		}
	return EV;
	}


static double eval_P1_max(int alt, int snode, int k_start, d_row V_pt, d_row P_pt, d_row im_P_pt, bool positive) {
	int inx,i,j,k,tnode;
	double pmin,pmass,add_on,EV;

	/* Collect all nodes on this level */
	for (k=k_start, tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode], k++) {
		/* Create local v-order */
		if (tdown[alt][tnode]) {
			/* Im-node */
			inx = at2i(alt,tnode);
			local_p_lobo[k] = im_L_hull_lobo[inx];
			local_p_upbo[k] = im_L_hull_upbo[inx];
			local_v[k] = eval_P1_max(alt,tnode,k+1,V_pt,P_pt,im_P_pt,positive);
			if (!positive)
				local_v[k] = -local_v[k];
			}
		else {
			/* Re-node */
			inx = at2r(alt,tnode);
			local_p_lobo[k] = L_hull_lobo[inx];
			local_p_upbo[k] = L_hull_upbo[inx];
			local_v[k] = V_pt[inx];
			}
		}
	for (i=k_start; i<k; i++)
		order[i] = i;
	sort_dom2(order,local_v,k_start,k-1,TRUE);
	EV = 0.0;
	pmin = 0.0;
	for (j=k_start; j<k; j++)
		pmin += local_p_lobo[j];
	pmass = 1.0 - pmin;
	/* Pick max in max-order. */
	for (j=k_start; j<k; j++) {
		add_on = min(local_p_upbo[order[j]]-local_p_lobo[order[j]],pmass);
		p_max[order[j]] = local_p_lobo[order[j]] + add_on;
		if (!positive)
			p_max[order[j]] = -p_max[order[j]];
		pmass -= add_on;
		EV += p_max[order[j]] * local_v[order[j]];
		}
	/* Update all nodes on this level */
	for (k=k_start, tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode], k++) {
		/* Update global nodes */
		if (tdown[alt][tnode]) {
			/* Im-node */
			inx = at2i(alt,tnode);
			im_P_pt[inx] = p_max[k];
			}
		else {
			/* Re-node */
			inx = at2r(alt,tnode);
			P_pt[inx] = p_max[k];
			}
		}
	return EV;
	}


static double eval_P1_min(int alt, int snode, int k_start, d_row V_pt, d_row P_pt, d_row im_P_pt, bool positive) {
	int inx,i,j,k,tnode;
	double pmin,pmass,add_on,EV;

	/* Collect all nodes on this level */
	for (k=k_start, tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode], k++) {
		/* Create local v-order */
		if (tdown[alt][tnode]) {
			/* Im-node */
			inx = at2i(alt,tnode);
			local_p_lobo[k] = im_L_hull_lobo[inx];
			local_p_upbo[k] = im_L_hull_upbo[inx];
			local_v[k] = eval_P1_min(alt,tnode,k+1,V_pt,P_pt,im_P_pt,positive);
			if (!positive)
				local_v[k] = -local_v[k];
			}
		else {
			/* Re-node */
			inx = at2r(alt,tnode);
			local_p_lobo[k] = L_hull_lobo[inx];
			local_p_upbo[k] = L_hull_upbo[inx];
			local_v[k] = V_pt[inx];
			}
		}
	for (i=k_start; i<k; i++)
		order[i] = i;
	sort_dom2(order,local_v,k_start,k-1,FALSE);
	EV = 0.0;
	pmin = 0.0;
	for (j=k_start; j<k; j++)
		pmin += local_p_lobo[j];
	pmass = 1.0 - pmin;
	/* Pick max in min-order. */
	for (j=k_start; j<k; j++) {
		add_on = min(local_p_upbo[order[j]]-local_p_lobo[order[j]],pmass);
		p_max[order[j]] = local_p_lobo[order[j]] + add_on;
		if (!positive)
			p_max[order[j]] = -p_max[order[j]];
		pmass -= add_on;
		EV += p_max[order[j]] * local_v[order[j]];
		}
	/* Update all nodes on this level */
	for (k=k_start, tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode], k++) {
		/* Update global nodes */
		if (tdown[alt][tnode]) {
			/* Im-node */
			inx = at2i(alt,tnode);
			im_P_pt[inx] = p_max[k];
			}
		else {
			/* Re-node */
			inx = at2r(alt,tnode);
			P_pt[inx] = p_max[k];
			}
		}
	return EV;
	}


double eval_P_max(int alt, int snode, int k_start, d_row V_pt, d_row P_pt, d_row im_P_pt, bool positive) {
	double EV;

	EV = eval_P1_max(alt,snode,k_start,V_pt,P_pt,im_P_pt,TRUE);
	if (!positive)
		EV = -EV;
	return EV;
	}


double eval_P_min(int alt, int snode, int k_start, d_row V_pt, d_row P_pt, d_row im_P_pt, bool positive) {
	double EV;

	EV = eval_P1_min(alt,snode,k_start,V_pt,P_pt,im_P_pt,TRUE);
	if (!positive)
		EV = -EV;
	return EV;
	}


/* Find max manually for P-base */

rcode TCL_get_P_max(struct d_frame *df, int alt, int snode,
				d_row V_pt, d_row P_pt, d_row im_P_pt, bool positive, double *maxval) {

	/* Input checks */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;

	*maxval = eval_P_max(alt,snode,1,V_pt,P_pt,im_P_pt,positive);
	return TCL_OK;
	}


/* Find min manually for P-base */

rcode TCL_get_P_min(struct d_frame *df, int alt, int snode,
				d_row V_pt, d_row P_pt, d_row im_P_pt, bool positive, double *minval) {

	/* Input checks */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;

	*minval = eval_P_min(alt,snode,1,V_pt,P_pt,im_P_pt,positive);
	return TCL_OK;
	}


 /*********************************************************
  *
  *  Security levels
  *
  *********************************************************/

double ixset_P_max(int alt, int snode, i_row ixset) {
	int inx,tnode;
	double LP_sum,LP_csum,sprob;

	/* Collect the maximum probability for the index set
	   and the minimum probability for the complementary set */
	LP_sum = 0.0;
	LP_csum = 1.0;
	/* Visit all nodes at this level */
	for (tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode]) {
		/* Collect dangerous set and complement set */
		if (tdown[alt][tnode]) {
			/* Im-node */
			inx = at2i(alt,tnode);
			sprob = ixset_P_max(alt,tnode,ixset);
			/* Distribute fraction over both sets */
			LP_sum += sprob * im_L_hull_upbo[inx];
			LP_csum -= (1.0-sprob) * im_L_hull_lobo[inx];
			}
		else {
			/* Re-node */
			inx = at2r(alt,tnode);
			if (ixset[inx])
				LP_sum += L_hull_upbo[inx];
			else
				LP_csum -= L_hull_lobo[inx];
			}
		}
	return min(LP_sum,LP_csum);
	}


double ixset_P_min(int alt, int snode, i_row ixset) {
	int inx,tnode;
	double LP_sum,LP_csum,sprob;

	/* Collect the minimum probability for the index set
	   and the maximum probability for the complementary set */
	LP_sum = 0.0;
	LP_csum = 1.0;
	/* Visit all nodes at this level */
	for (tnode=tdown[alt][snode]; tnode; tnode=tnext[alt][tnode]) {
		/* Collect dangerous set and complement set */
		if (tdown[alt][tnode]) {
			/* Im-node */
			inx = at2i(alt,tnode);
			sprob = ixset_P_min(alt,tnode,ixset);
			/* Distribute fraction over both sets */
			LP_sum += sprob * im_L_hull_lobo[inx];
			LP_csum -= (1.0-sprob) * im_L_hull_upbo[inx];
			}
		else {
			/* Re-node */
			inx = at2r(alt,tnode);
			if (ixset[inx])
				LP_sum += L_hull_lobo[inx];
			else
				LP_csum -= L_hull_upbo[inx];
			}
		}
	return max(LP_sum,LP_csum);
	}
