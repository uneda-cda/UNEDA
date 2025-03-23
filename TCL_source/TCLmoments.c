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
 *   File: TCLmoments.c
 *
 *   Purpose: compute moments for the alternatives
 *            using the NEMO (net moment) calculus
 *
 *   Algorithms described in M. Danielson: The BEDA Method
 *
 *
 *   Functions exported outside TCL
 *   ------------------------------
 *   TCL_get_moments
 *   TCL_get_mc_moments
 *   TCL_get_P_sd
 *   TCL_get_V_sd
 *
 *   Functions outside module, inside TCL
 *   ------------------------------------
 *   NONE
 *
 *   Functions internal to module
 *   ----------------------------
 *   calc_nemo_Pnode
 *   calc_nemo_Vnode
 *   mult_moments
 *   calc_nemo_tree
 *   calc_nemo_mc_tree
 *
 */

#include "TCLinternal.h"


 /*********************************************************
  *
  *  Calculate central moments for single variables
  *
  *********************************************************/

/* Calculate the NEMO parameters for a P-node */

static void calc_nemo_Pnode(double lobo, double mid, double upbo, 
				double lambda, double *var, double *cov) {
	double t;

	/* This is a generalised (bounded) joint m-variate Dirichlet (beta marginals)
	   distribution with non-standard parameters due to lower and upper bounds.
	   Mid interpreted as mean. */
	t = upbo-lobo;
	if (t > EPS) {
		*var = t*t*mid*(1.0-mid)/(lambda+1.0);
		*cov = t*t*mid*mid/(lambda+1.0);
		}
	else { // prevent breakdown in Dirac delta-pulse singular points
		*var = 0.0;
		*cov = 0.0;
		}
	}


/* Calculate NEMO parameters for a V-node. The configuration conditional flag
 * V_MID_SNAP controls whether a non-physical mean value is corrected or not. */

#define V_MID_SNAP      // do not use mean literally, snap to physical triangle
#ifdef V_MID_SNAP
#define V_MID_SNAP_HALF // use mean semi-literally, i.e. snap halfway
#endif

static void calc_nemo_Vnode(double lobo, double mid, double upbo, double *var, double *tcm) {
	double t,q;
	double mean,mode;

	/* Employing a triangular distribution */
	t = upbo-lobo;
	if (t > EPS) {
		/* Mid interpreted as mean -> convert to mode */
#ifdef V_MID_SNAP
		/* Snap to geometrically compliant triangle ('M' in BRS).
		 * Note: this entails that unphysical midpoints have no influence. */
		mean = max(mid,(2.0*lobo+upbo)/3.0);
		mean = min(mean,(lobo+2.0*upbo)/3.0);
#ifdef V_MID_SNAP_HALF
		mean += (mid-mean)/2.0;
#endif
#else
		/* Non-snap: triangle top point can extend to be at 2 times the base. This
		 * entails that the variance can be up to 3 times a std physical triangle. */
		mean = mid;
#endif
		mode = 3.0*mean-lobo-upbo;
		q = (mode-lobo)/t;
		*var = t*t*(1.0-q+q*q)/18.0;
		*tcm = t*t*t*(2.0-3.0*q-3.0*q*q+2.0*q*q*q)/270.0;
		if (*tcm < 1E-18) // EPS^3
			*tcm = 0.0;
		}
	else { // prevent breakdown in Dirac-delta-pulse singular points
		*var = 0.0;
		*tcm = 0.0;
		}
	}


 /*********************************************************
  *
  *  Multiply P and V moments -> NEMO product moments
  *
  *********************************************************/

static void mult_moments(double Pmean, double Pvar, double Pcov, 
				double Vmean, double Vvar, double Vtcm, 
				double *PVmean, double *PVvar, double *PVcov, double *PVtcm) {

	/* NEMO mean, variance, covariance and tcm for PV products */
	*PVmean= Pmean*Vmean;
	*PVvar = Pvar*Vvar + Pvar*Vmean*Vmean + Pmean*Pmean*Vvar;
	*PVcov = sqrt(Pcov)*Vmean;
	*PVtcm = Pmean*Vtcm;
	}


 /*********************************************************
  *
  *  Multiply and add PV products in a tree, node by node
  *
  *********************************************************/

static d_row P_lobo,P_mid,P_upbo,V_lobo,V_mid,V_upbo,P_sd,V_sd;
// note: covar matrix entries used selectively in the recursion
static double covar[MAX_NODES+1][MAX_NODES+1];

static void calc_nemo_tree(struct d_frame *df, int alt, int snode, double *tot_mean,
		double *tot_var, double *tot_tcm) {
	int tnode,t2node,inx,inx1,inx2,n_nodes;
	double P_var,P_cov,V_mean,V_var,V_tcm;
	double PV_mean,PV_var,PV_cov,PV_tcm;
	double node_cov,lambda,lobo_s;

	/* Calculate lambda scaling */
	lambda = lobo_s = 0.0;
	for (tnode=df->down[alt][snode]; tnode; tnode=df->next[alt][tnode]) {
		inx = at2f(alt,tnode);
		lambda += P_upbo[inx]-P_lobo[inx];
		lobo_s += P_lobo[inx];
		}
	if (lobo_s < 1.0-EPS)
		lambda /= (1.0-lobo_s);
	else // 0-div-by-0 -> l'Hospital -> converges to 1
		lambda = 1.0;
	/* Calculate moments from local tree node data */
	*tot_mean = 0.0;
	*tot_var  = 0.0;
	*tot_tcm  = 0.0;
	/* Sum all nodes at this level */
	for (tnode=df->down[alt][snode]; tnode; tnode=df->next[alt][tnode]) {
		inx = at2f(alt,tnode);
		if (df->down[alt][tnode])	{
			/* Im-node */
			calc_nemo_Pnode(P_lobo[inx],P_mid[inx],P_upbo[inx],lambda,&P_var,&P_cov);
			calc_nemo_tree(df,alt,tnode,&V_mean,&V_var,&V_tcm);
			mult_moments(P_mid[inx],P_var,P_cov,V_mean,V_var,V_tcm,
					&PV_mean,&PV_var,&PV_cov,&PV_tcm);
			}
		else {
			/* Re-node */
			calc_nemo_Pnode(P_lobo[inx],P_mid[inx],P_upbo[inx],lambda,&P_var,&P_cov);
			calc_nemo_Vnode(V_lobo[inx],V_mid[inx],V_upbo[inx],&V_var,&V_tcm);
			mult_moments(P_mid[inx],P_var,P_cov,V_mid[inx],V_var,V_tcm,
					&PV_mean,&PV_var,&PV_cov,&PV_tcm);
			}
		*tot_mean += PV_mean;
		*tot_var  += PV_var;
		*tot_tcm  += PV_tcm;
		covar[inx][0] = PV_cov; // store separable covariance in column 0
		P_sd[inx] = sqrt(P_var);
		V_sd[inx] = sqrt(V_var);
		}
	/* Generate separable covariance matrix (upper right triangle, above diagonal)
	 * The lower left triangle is a mirror image, so the sum is doubled later. */
	for (tnode=df->down[alt][snode]; tnode; tnode=df->next[alt][tnode]) {
		inx1 = at2f(alt,tnode);
		for (t2node=df->next[alt][tnode]; t2node; t2node=df->next[alt][t2node]) {
			inx2 = at2f(alt,t2node);
			covar[inx1][inx2] = -covar[inx1][0];
			}
		for (t2node=df->down[alt][snode]; t2node!=tnode; t2node=df->next[alt][t2node]) {
			inx2 = at2f(alt,t2node);
			covar[inx2][inx1] *= covar[inx1][0];
			}
		}
	/* Calculate variance from separable covariance (sum of upper right triangle entries) */
	n_nodes = 0;
	node_cov = 0.0;
	for (tnode=df->down[alt][snode]; tnode; tnode=df->next[alt][tnode]) {
		inx1 = at2f(alt,tnode);
		for (t2node=df->next[alt][tnode]; t2node; t2node=df->next[alt][t2node]) {
			inx2 = at2f(alt,t2node);
			node_cov += covar[inx1][inx2];
			}
		n_nodes++;
		}
	*tot_var += 2.0*node_cov; // upper + lower triangle
	*tot_tcm /= (double)n_nodes;
	if (*tot_var < EPS) // catch roundoff errors
		*tot_var = 0.0;
	if (*tot_tcm < EPS) // catch roundoff errors
		*tot_tcm = 0.0;
	}


rcode TCL_get_moments(struct d_frame *df, a_row rm1, a_row cm2, a_row cm3) {
	int Ai;
	double A_mean,A_var,A_tcm;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;

	/* Pick up start vectors */
	l_hull_P(P_lobo,P_upbo);
	l_cpoint_P(P_mid);
	hull_V(V_lobo,V_upbo);
	cpoint_V(V_mid);
	/* Calculate moments for each alternative */
	for (Ai=1; Ai<=df->n_alts; Ai++) {
		calc_nemo_tree(df,Ai,0,&A_mean,&A_var,&A_tcm);
		rm1[Ai] = A_mean;
		cm2[Ai] = A_var;
		cm3[Ai] = A_tcm;
		}
	return TCL_OK;
	}


 /**************************************************************
  *
  *  Multiply and add WV products in an MC tree, node by node.
  *  Differences compared to the standard calculation above:
  *  1. Criteria V parameters are already supplied in the call.
  *  2. Delta/gamma range [-1,1] -> Pcov mult can be negative.
  *
  **************************************************************/

static void calc_nemo_mc_tree(struct d_frame *df, int snode, d_row Vx_rm1, d_row Vx_cm2,
			d_row Vx_cm3, double *tot_mean, double *tot_var, double *tot_tcm) {
	int tnode,t2node,inx,inx1,inx2,inxr,n_nodes;
	double P_var,P_cov,V_mean,V_var,V_tcm;
	double PV_mean,PV_var,PV_cov,PV_tcm;
	double node_cov,lambda,lobo_s;

	/* Calculate lambda scaling */
	lambda = lobo_s = 0.0;
	for (tnode=df->down[1][snode]; tnode; tnode=df->next[1][tnode]) {
		inx = at2f(1,tnode);
		lambda += P_upbo[inx]-P_lobo[inx];
		lobo_s += P_lobo[inx];
		}
	if (lobo_s < 1.0-EPS)
		lambda /= (1.0-lobo_s);
	else // 0-div-by-0 -> l'Hospital -> converges to 1
		lambda = 1.0;
	/* Calculate moments from local W tree node data and external V moments */
	*tot_mean = 0.0;
	*tot_var = 0.0;
	*tot_tcm = 0.0;
	/* Sum all nodes at this level */
	for (tnode=df->down[1][snode]; tnode; tnode=df->next[1][tnode]) {
		inx = at2f(1,tnode);
		if (df->down[1][tnode]) {
			/* Im-node */
			calc_nemo_Pnode(P_lobo[inx],P_mid[inx],P_upbo[inx],lambda,&P_var,&P_cov);
			calc_nemo_mc_tree(df,tnode,Vx_rm1,Vx_cm2,Vx_cm3,&V_mean,&V_var,&V_tcm);
			mult_moments(P_mid[inx],P_var,P_cov,V_mean,V_var,V_tcm,
					&PV_mean,&PV_var,&PV_cov,&PV_tcm);
			V_sd[inx] = sqrt(V_var);
			}
		else {
			/* Re-node */
			calc_nemo_Pnode(P_lobo[inx],P_mid[inx],P_upbo[inx],lambda,&P_var,&P_cov);
			inxr = at2r(1,tnode);
			mult_moments(P_mid[inx],P_var,P_cov,Vx_rm1[inxr],Vx_cm2[inxr],Vx_cm3[inxr],
					&PV_mean,&PV_var,&PV_cov,&PV_tcm);
			V_sd[inx] = sqrt(Vx_cm2[inxr]);
			}
		*tot_mean += PV_mean;
		*tot_var += PV_var;
		*tot_tcm += PV_tcm;
		covar[inx][0] = PV_cov; // store separable covariance in column 0
		P_sd[inx] = sqrt(P_var);
		}
	/* Generate separable covariance matrix (upper right triangle, above diagonal) */
	for (tnode=df->down[1][snode]; tnode; tnode=df->next[1][tnode]) {
		inx1 = at2f(1,tnode);
		for (t2node=df->next[1][tnode]; t2node; t2node=df->next[1][t2node]) {
			inx2 = at2f(1,t2node);
			covar[inx1][inx2] = -covar[inx1][0];
			}
		for (t2node=df->down[1][snode]; t2node!=tnode; t2node=df->next[1][t2node]) {
			inx2 = at2f(1,t2node);
			covar[inx2][inx1] *= covar[inx1][0];
			}
		}
	/* Calculate variance from separable covariance (sum of upper right triangle entries).
	 * Symmetric along the diagonal -> calculate sum of upper right and double it. */
	n_nodes = 0;
	node_cov = 0.0;
	for (tnode=df->down[1][snode]; tnode; tnode=df->next[1][tnode]) {
		inx1 = at2f(1,tnode);
		for (t2node=df->next[1][tnode]; t2node; t2node=df->next[1][t2node]) {
			inx2 = at2f(1,t2node);
			node_cov += covar[inx1][inx2];
			}
		n_nodes++;
		}
	*tot_var += 2.0*node_cov; // upper + lower triangle
	*tot_tcm /= (double)n_nodes;
	if (*tot_var < EPS) // catch roundoff errors
		*tot_var = 0.0;
	if (*tot_tcm < EPS) // catch roundoff errors
		*tot_tcm = 0.0;
	}


rcode TCL_get_mc_moments(struct d_frame *df, int snode, d_row Vx_rm1, d_row Vx_cm2, 
			d_row Vx_cm3, double *rm1, double *cm2, double *cm3) {

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;

	/* Pick up start vectors only for W since V is supplied in the call. */
	l_hull_P(P_lobo,P_upbo);
	l_cpoint_P(P_mid);
	/* Calculate moments for MC */
	calc_nemo_mc_tree(df,snode,Vx_rm1,Vx_cm2,Vx_cm3,rm1,cm2,cm3);
	return TCL_OK;
	}


 /****************************************************************
  *
  *  Get the standard deviation for individual P and V variables
  *
  ****************************************************************/

rcode TCL_get_P_sd(struct d_frame *df, int inx, double *sd) {

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;
	if ((inx < 1) || (inx > df->tot_cons[0]))
		return TCL_INPUT_ERROR;
	/* Deliver standard deviation for re+im-nodes */
	*sd = P_sd[inx];
	return TCL_OK;
	}


/* im=0: return only nodes that can have user input -> not im-nodes
 * im=1: return also nodes that cannot have user input -> all nodes */

rcode TCL_get_V_sd(struct d_frame *df, int inx, int im, double *sd) {

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;
	if ((inx < 1) || (inx > df->tot_cons[0]))
		return TCL_INPUT_ERROR;
	/* Deliver standard deviation for re(+im)-nodes */
	if (f2r[inx] || im)
		*sd = V_sd[inx]; // re-node or im-flag
	else
		*sd = -1.0; // im-node and no im-flag
	return TCL_OK;
	}
