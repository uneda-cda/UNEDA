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
 *                   UNEDA Decision Tree Layer (DTL)
 *                   -------------------------------
 *
 *    +----- o o o ------------------------------------------------+
 *    |    o       o              Prof. Mats Danielson             |
 *    |   o  STHLM  o             DECIDE Research Group            |
 *    |   o         o    Dept. of Computer and Systems Sciences    |
 *    |   o   UNI   o             Stockholm University             |
 *    |    o       o      PO Box 1203, SE-164 25 Kista, SWEDEN     |
 *    +----- o o o ------------------------------------------------+
 *
 *                Copyright (c) 2012-2025 Mats Danielson
 *                     Email: mats.danielson@su.se
 *
 */

/*
 *   File: DTLtornado.c
 *
 *   Purpose: evaluation of tornados
 *
 *
 *   Functions exported outside DTL
 *   ------------------------------
 *   DTL_get_P_tornado
 *   DTL_get_MCP_tornado
 *   DTL_get_V_tornado
 *   DTL_get_MCV_tornado
 *   DTL_get_W_tornado
 *   DTL_get_W_tornado_alt
 *   DTL_get_cons_influence
 *
 *   Functions outside of module, inside DTL
 *   ---------------------------------------
 *   NONE
 *
 *   Functions internal to module
 *   ----------------------------
 *   rollback_PW_base
 *   dtl_mass_PW_tornado
 *   dtl_get_PW_tornado
 *   dtl_get_MCP_tornado
 *   rollback_V_base
 *   dtl_mass_V_tornado
 *   dtl_get_V_tornado
 *   dtl_get_MCV_tornado
 *   dtl_get_W_tornado
 *   dtl_get_cons_influence
 *
 */


 /***********************************************************
  *
  *  Tornados for sensitivity analyses
  *  ---------------------------------
  *
  *  Tornados describe how much a single variable influences
  *  the evaluation result, either in the form of EV or mass.
  *
  *  NOTE: If CAR data has been generated with CAR_light=OFF,
  *  P/W-tornados with odd mode numbers will be used (no mid),
  *  since a full set of midpoints cannot move around and
  *  would only result in empty tornados. This is, however,
  *  not preserved on file since mboxes are stored as stmts.
  *  This way, entered and stored frames can yield differing
  *  results when it comes to P/W-tornados. The debugging
  *  remedy is to explicitly run with floating midpoints.
  *
  *  NOTE2: In the very unusual case of split midpoints, i.e.
  *  the mbox contains different upper and lower points, the
  *  tornados will not always be consistent. Split midpoints
  *  are not encouraged, but if they are still being used,
  *  refrain from using tornado functions with such frames.
  *  Especially W/P-tornados are susceptible because of the
  *  normalisation, i.e. summing the constituents to one.
  *
  ***********************************************************/

#include "DTL.h"
#include "DTLinternal.h"


 /********************************************************
  *
  *  Configuration parameters
  *
  ********************************************************/

#define T_EPS 4.0E-6 // must be larger than 2*CAR_EPS


 /*********************************************************
  *
  *  Local data areas
  *
  *********************************************************/

static e_matrix t_result;
static d_row h_lobo,h_upbo,m_lobo,m_upbo,ms_lobo,ms_upbo;
static d_row W_mid,LW_mid,V_mid;


 /*********************************************************
  *
  *  Probability/weight tornados
  *
  *********************************************************/

static void rollback_PW_base(int pstmt, d_row lobo, d_row upbo) {

	/* To preserve offending rc, the rollback may not use call() */
	if (pstmt)
		TCL_delete_P_constraint(uf->df,pstmt);
	TCL_set_P_mbox(uf->df,lobo,upbo);
	}


static rcode dtl_mass_PW_tornado(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo) {
	rcode rc;
	int i,j,k,alt,sym,global_cst;
	double baseline_ev,baseline_mass;
	struct d_frame *df;

	alt = max(-crit,0); // W, not P
	crit = max(crit,0); // P, not W
	if (load_df0(crit))
		return DTL_CRIT_UNKNOWN;
	df = uf->df;
	/* Clean mode field */
	mode &= 0x01;
	if (call(TCL_get_P_mbox(df,ms_lobo,ms_upbo),"TCL_get_P_mbox"))
		return dtl_kernel_error();
	if (!mode) {
		for (k=1,i=1; i<=df->n_alts; i++)
			for (j=1; j<=df->tot_cons[i]; j++,k++) {
				m_lobo[k] = -1.0;
				m_upbo[k] = -1.0;
				}
		/* Clear mhull */
		if (call(TCL_set_P_mbox(df,m_lobo,m_upbo),"TCL_set_P_mbox"))
			return dtl_kernel_error();
		}
	/* Stop cst reporting */
	global_cst = cst_on;
	cst_on = FALSE;
	for (k=1,i=1; i<=df->n_alts; i++) {
		if (rc = crit?evaluate_frame(crit,E_PSI,i,0,t_result):
						 evaluate_frameset(0,E_PSI,alt,0,t_result)) {
			if (!mode)
				rollback_PW_base(0,ms_lobo,ms_upbo); // restore
			cst_on = global_cst;
			return rc;
			}
		baseline_ev = t_result[E_MID][0];
		if (rc = dtl_ev_to_cdf(crit,baseline_ev,&baseline_mass)) {
			if (!mode)
				rollback_PW_base(0,ms_lobo,ms_upbo); // restore
			cst_on = global_cst;
			return rc;
			}
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			/* Convert EV into mass in situ */
			sym = fabs(t_upbo[i][j]+t_lobo[i][j]) < 2.0E-4;
			if (t_lobo[i][j] < -DTL_EPS) {
				if (rc = dtl_ev_to_cdf(crit,max(baseline_ev+t_lobo[i][j],0.0),t_lobo[i]+j)) {
					if (!mode)
						rollback_PW_base(0,ms_lobo,ms_upbo); // restore
					cst_on = global_cst;
					return rc;
					}
				t_lobo[i][j] = baseline_mass-t_lobo[i][j];
				}
			else
				t_lobo[i][j] = 0.0;
			if (t_upbo[i][j] > DTL_EPS) {
				if (rc = dtl_ev_to_cdf(crit,min(baseline_ev+t_upbo[i][j],1.0),t_upbo[i]+j)) {
					if (!mode)
						rollback_PW_base(0,ms_lobo,ms_upbo); // restore
					cst_on = global_cst;
					return rc;
					}
				t_upbo[i][j] = baseline_mass-t_upbo[i][j];
				}
			else
				t_upbo[i][j] = 0.0;
			/* Rebalance to counter the Ericsson effect */
			if (sym) {
				t_upbo[i][j] = (t_upbo[i][j]-t_lobo[i][j])/2.0;
				t_lobo[i][j] = -t_upbo[i][j];
				}
			}
		}
	if (!mode)
		rollback_PW_base(0,ms_lobo,ms_upbo); // restore
	cst_on = global_cst;
	eval_cache_invalidate();
	return DTL_OK;
	}


static rcode dtl_get_PW_tornado(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo) {
	rcode rc;
	int i,j,k,jj,kk,ts_nbr,global_cst;
	int h_start;
	double bound1,bound2,baseline_ev,h_sum,th_lo_bound,th_up_bound;
	struct stmt_rec stmt;
	struct d_frame *df;

	/* Check input parameters */
	if (load_df0(crit))
		return DTL_CRIT_UNKNOWN;
	df = uf->df;
	/* Clean mode field */
	mode &= 0x01;
	if (uf->WP_autogen[crit])
		mode = 0; // autogen needs floating midpoint
	/* Collect starting point */
	dtl_abort_init();
	if (call(TCL_get_P_hull(df,m_lobo,m_upbo,h_lobo,h_upbo),"TCL_get_P_hull"))
		return dtl_kernel_error();
	if (call(TCL_get_P_mbox(df,m_lobo,m_upbo),"TCL_get_P_mbox"))
		return dtl_kernel_error();
	if (call(TCL_get_P_mbox(df,ms_lobo,ms_upbo),"TCL_get_P_mbox"))
		return dtl_kernel_error();
	/* Lower cst reporting level */
	global_cst = cst_on;
	cst_on = cst_ext;
	if (!mode) {
		/* Mode 0: midpoint removed */
		for (k=1,i=1; i<=df->n_alts; i++)
			for (j=1; j<=df->tot_cons[i]; j++,k++) {
				m_lobo[k] = -1.0;
				m_upbo[k] = -1.0;
				}
		/* Clear mhull */
		if (call(TCL_set_P_mbox(df,m_lobo,m_upbo),"TCL_set_P_mbox")) {
			cst_on = global_cst;
			return dtl_kernel_error();
			}
		}
	/* Insert dummy P-statement */
	stmt.n_terms = 1;
	stmt.alt[1] = 1;
	stmt.cons[1] = 1;
	stmt.sign[1] = 1;
	stmt.lobo = 0.0;
	stmt.upbo = 1.0;
	if (call(TCL_add_P_constraint(df,&stmt),"TCL_add_P_constraint")) {
		rollback_PW_base(0,ms_lobo,ms_upbo);
		cst_on = global_cst;
		return dtl_kernel_error();
		}
	ts_nbr = df->P_base->n_stmts;
	for (k=1,i=1; i<=df->n_alts; i++) {
		dtl_abort_check();
		h_start = k;
		if (rc = evaluate_frame(crit,E_PSI,i,0,t_result)) {
			rollback_PW_base(ts_nbr,ms_lobo,ms_upbo);
			cst_on = global_cst;
			return rc;
			}
		baseline_ev = t_result[E_MID][0];
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			/* Find lower boundary */
			if (mode) {
				/* Mode 1: explicit midpoint kept */
				/* Calc mass point movability limits */
				h_sum = 0.0;
				for (jj=1,kk=h_start; jj<=df->tot_cons[i]; jj++,kk++)
					if ((jj != j) && !TCL_different_parents(df,i,j,jj)) {
						/* Limit is midpoint if it exists, else hull */
						if (ms_upbo[kk] >= 0.0)
							h_sum += ms_upbo[kk];
						else
							h_sum += min(h_upbo[kk],1.0);
						}
				if (ms_lobo[k] >= 0.0)
					th_lo_bound = max(ms_lobo[k],1.0-h_sum);
				else
					th_lo_bound = max(h_lobo[k],1.0-h_sum);
				}
			else
				/* Mode 0: midpoint removed */
				th_lo_bound = max(h_lobo[k],0.0);
			/* Find upper boundary */
			if (mode) {
				/* Mode 1: explicit midpoint kept */
				/* Calc mass point movability limits */
				h_sum = 0.0;
				for (jj=1,kk=h_start; jj<=df->tot_cons[i]; jj++,kk++)
					if ((jj != j) && !TCL_different_parents(df,i,j,jj)) {
						/* Limit is midpoint if it exists, else hull */
						if (ms_lobo[kk] >= 0.0)
							h_sum += ms_lobo[kk];
						else
							h_sum += max(h_lobo[kk],0.0);
						}
				/* Catch roundoff errors */
				h_sum = min(h_sum,1.0);
				if (ms_upbo[k] >= 0.0)
					th_up_bound = min(ms_upbo[k],1.0-h_sum);
				else
					th_up_bound = min(h_upbo[k],1.0-h_sum);
				}
			else
				/* Mode 0: midpoint removed */
				th_up_bound = min(h_upbo[k],1.0);
			/* Remove midpoint for this consequence */
			if (mode) {
				m_lobo[j] = -1.0;
				m_upbo[j] = -1.0;
				/* Set mhull */
				if (call(TCL_set_P_mbox(df,m_lobo,m_upbo),"TCL_set_P_mbox")) {
					rollback_PW_base(ts_nbr,ms_lobo,ms_upbo);
					cst_on = global_cst;
					return dtl_kernel_error();
					}
				}
			/* Find movement in mass point */
			if (th_up_bound-th_lo_bound < 5.0*T_EPS)
				/* Narrow trap */
				t_lobo[i][j] = t_upbo[i][j] = 0.0;
			else {
				/* Prepare statement */
				stmt.alt[1] = i;
				stmt.cons[1] = j;
				stmt.lobo = max(th_lo_bound,0.0);
				stmt.upbo = min(th_lo_bound+T_EPS,1.0);
				/* Explore lower boundary */
				if (call(TCL_replace_P_constraint(df,ts_nbr,&stmt),"TCL_replace_P_constraint")) {
					rollback_PW_base(ts_nbr,ms_lobo,ms_upbo);
					cst_on = global_cst;
					return dtl_kernel_error();
					}
				if (dtl_abort_request) {
					rollback_PW_base(ts_nbr,ms_lobo,ms_upbo);
					cst_on = global_cst;
					}
				dtl_abort_check();
				if (rc = evaluate_frame(crit,E_PSI,i,0,t_result)) {
					rollback_PW_base(ts_nbr,ms_lobo,ms_upbo);
					cst_on = global_cst;
					return rc;
					}
				bound1 = t_result[E_MID][0];
				/* Explore upper boundary */
				if (call(TCL_change_P_constraint(df,ts_nbr,max(th_up_bound-T_EPS,0.0),
						min(th_up_bound,1.0)),"TCL_change_P_constraint")) {
					rollback_PW_base(ts_nbr,ms_lobo,ms_upbo);
					cst_on = global_cst;
					return dtl_kernel_error();
					}
				if (dtl_abort_request) {
					rollback_PW_base(ts_nbr,ms_lobo,ms_upbo);
					cst_on = global_cst;
					}
				dtl_abort_check();
				if (rc = evaluate_frame(crit,E_PSI,i,0,t_result)) {
					rollback_PW_base(ts_nbr,ms_lobo,ms_upbo);
					cst_on = global_cst;
					return rc;
					}
				bound2 = t_result[E_MID][0];
				/* An increase in P can decrease EV and v.v. -> must sort boundaries */
				if (bound1 < bound2) {
					t_lobo[i][j] = bound1 - baseline_ev;
					t_upbo[i][j] = bound2 - baseline_ev;
					}
				else {
					t_lobo[i][j] = bound2 - baseline_ev;
					t_upbo[i][j] = bound1 - baseline_ev;
					}
				/* Catch roundoff errors */
				if (t_lobo[i][j] > -T_EPS)
					t_lobo[i][j] = 0.0;
				if (t_upbo[i][j] < +T_EPS)
					t_upbo[i][j] = 0.0;
				}
			/* Clean-up for next variable */
			if (call(TCL_change_P_constraint(df,ts_nbr,0.0,1.0),"TCL_change_P_constraint")) {
				rollback_PW_base(ts_nbr,ms_lobo,ms_upbo);
				cst_on = global_cst;
				return dtl_kernel_error();
				}
			if (mode) {
				/* Restore mhull */
				m_lobo[k] = ms_lobo[k];
				m_upbo[k] = ms_upbo[k];
				if (call(TCL_set_P_mbox(df,m_lobo,m_upbo),"TCL_set_P_mbox")) {
					rollback_PW_base(ts_nbr,ms_lobo,ms_upbo);
					cst_on = global_cst;
					return dtl_kernel_error();
					}
				}
			}
		}
	rollback_PW_base(ts_nbr,ms_lobo,ms_upbo);
	eval_cache_invalidate();
	cst_on = global_cst;
	return DTL_OK;
	}


/* Mode: 0 = Midpoint kept (default)
 *       1 = Force floating midpoint
 *      +2 = Belief mass output */

rcode DTLAPI DTL_get_P_tornado(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo) {
	rcode rc;
	int i,j;

	/* Begin single thread semaphore */
	_smx_begin("TOP");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_P_tornado(%d,%d)\n",crit,mode);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(t_lobo,1);
	_certify_ptr(t_upbo,2);
	_dtl_assert(t_lobo!=t_upbo,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (!crit)
		return dtl_error(DTL_CRIT_UNKNOWN);
	if ((mode < 0) || (mode > 3))
		return dtl_error(DTL_INPUT_ERROR);
	mode ^= 0x01; // flip lowest mode bit (compat reasons)
	/* Get tornado */
	if (rc = dtl_get_PW_tornado(crit,mode,t_lobo,t_upbo))
		return dtl_error(rc);
	if (mode&0x02) {
		if (rc = dtl_mass_PW_tornado(crit,mode,t_lobo,t_upbo))
			return dtl_error(rc);
		}
	/* Log function result */
	if (cst_ext)
		for (i=1; i<=uf->df->n_alts; i++)
			for (j=1; j<=uf->df->tot_cons[i]; j++) {
				sprintf(msg,"    P%d.%d.%-2d [%.3lf %.3lf]\n",crit,i,j,t_lobo[i][j],t_upbo[i][j]);
				cst_log(msg);
				}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


static h_matrix x_lobo,x_upbo;

static rcode dtl_get_MCP_tornado(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo) {
	rcode rc;
	int i,j,node;

	/* Check input parameters */
	if (!crit)
		return DTL_CRIT_UNKNOWN;
	/* Get local tornado in criterion */
	if (rc = dtl_get_PW_tornado(crit,mode,t_lobo,t_upbo))
		return rc;
	if (mode&0x02) {
		if (rc = dtl_mass_PW_tornado(crit,mode,t_lobo,t_upbo))
			return rc;
		}

	/* Get global weight for this criterion */
	if (load_df0(0))
		return DTL_SYS_CORRUPT;
	if (call(TCL_get_P_masspoint(uf->df,W_mid,LW_mid),"TCL_get_P_masspoint"))
		return dtl_kernel_error();
	/* Get global weight node index */
	node = r2t[1][crit];

	/* Apply to local criterion tornado */
	if (load_df1(crit))
		return DTL_CRIT_UNKNOWN;
	for (i=1; i<=uf->df->n_alts; i++)
		for (j=1; j<=uf->df->tot_cons[i]; j++) {
			t_lobo[i][j] *= W_mid[node];
			t_upbo[i][j] *= W_mid[node];
			}
	return DTL_OK;
	}


/* Mode: 0 = Midpoint kept (default)
 *       1 = Force floating midpoint
 *      +2 = Belief mass output */

rcode DTLAPI DTL_get_MCP_tornado(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo) {
	rcode rc;
	int i,j;

	/* Begin single thread semaphore */
	_smx_begin("TMCP");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_MCP_tornado(%d,%d)\n",crit,mode);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(t_lobo,1);
	_certify_ptr(t_upbo,2);
	_dtl_assert(t_lobo!=t_upbo,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (!PM)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if ((mode < 0) || (mode > 3))
		return dtl_error(DTL_INPUT_ERROR);
	mode ^= 0x01; // flip lowest mode bit (compat reasons)
	/* Get global tornado in criterion */
	if (rc = dtl_get_MCP_tornado(crit,mode,t_lobo,t_upbo))
		return dtl_error(rc);
	/* Log function result */
	if (cst_ext)
		for (i=1; i<=uf->df->n_alts; i++)
			for (j=1; j<=uf->df->tot_cons[i]; j++) {
				sprintf(msg,"    P%d.%d.%-2d [%.3lf %.3lf]\n",crit,i,j,t_lobo[i][j],t_upbo[i][j]);
				cst_log(msg);
				}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /************************************************************
  *
  *  Value tornados
  *
  *  NOTE: due to there currently being no dependencies btw
  *  value variables, the 'rough' binary tree tornados work
  *  equally well. But staying prepared for the future, it's
  *  a safer bet to keep using the more advanced V-tornados.
  *
  ************************************************************/

static void rollback_V_base(int vstmt, d_row lobo, d_row upbo) {

	/* To preserve offending rc, the rollback may not use call() */
	if (vstmt)
		TCL_delete_V_constraint(uf->df,vstmt);
	TCL_set_V_mbox(uf->df,lobo,upbo);
	}


static rcode dtl_mass_V_tornado(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo) {
	rcode rc;
	int i,j,k,sym,global_cst;
	double baseline_ev,baseline_mass;
	struct d_frame *df;

	if (load_df1(crit))
		return DTL_CRIT_UNKNOWN;
	df = uf->df;
	/* Clean mode field */
	mode &= 0x01;
	if (call(TCL_get_V_mbox(df,ms_lobo,ms_upbo),"TCL_get_V_mbox"))
		return dtl_kernel_error();
	if (!mode) {
		for (k=1,i=1; i<=df->n_alts; i++)
			for (j=1; j<=df->tot_cons[i]; j++,k++) {
				m_lobo[k] = -1.0;
				m_upbo[k] = -1.0;
				}
		/* Clear mhull */
		if (call(TCL_set_V_mbox(df,m_lobo,m_upbo),"TCL_set_V_mbox"))
			return dtl_kernel_error();
		}
	/* Stop cst reporting */
	global_cst = cst_on;
	cst_on = FALSE;
	for (k=1,i=1; i<=df->n_alts; i++) {
		if (rc = evaluate_frame(crit,E_PSI,i,0,t_result)) {
			if (!mode)
				rollback_V_base(0,ms_lobo,ms_upbo); // restore
			cst_on = global_cst;
			return rc;
			}
		baseline_ev = t_result[E_MID][0];
		if (rc = dtl_ev_to_cdf(crit,baseline_ev,&baseline_mass)) {
			if (!mode)
				rollback_V_base(0,ms_lobo,ms_upbo); // restore
			cst_on = global_cst;
			return rc;
			}
		for (j=1; j<=df->tot_cons[i]; j++,k++)
			if (TCL_get_V_index(df,i,j)) { // is a real-cons
				/* Convert EV into mass in situ */
				sym = fabs(t_upbo[i][j]+t_lobo[i][j]) < 2.0E-4;
				if (t_lobo[i][j] < -DTL_EPS) {
					if (rc = dtl_ev_to_cdf(crit,max(baseline_ev+t_lobo[i][j],0.0),t_lobo[i]+j)) {
						if (!mode)
							rollback_V_base(0,ms_lobo,ms_upbo); // restore
						cst_on = global_cst;
						return rc;
						}
					t_lobo[i][j] = baseline_mass-t_lobo[i][j];
					}
				else
					t_lobo[i][j] = 0.0;
				if (t_upbo[i][j] > DTL_EPS) {
					if (rc = dtl_ev_to_cdf(crit,min(baseline_ev+t_upbo[i][j],1.0),t_upbo[i]+j)) {
						if (!mode)
							rollback_V_base(0,ms_lobo,ms_upbo); // restore
						cst_on = global_cst;
						return rc;
						}
					t_upbo[i][j] = baseline_mass-t_upbo[i][j];
					}
				else
					t_upbo[i][j] = 0.0;
				/* Rebalance to counter the Ericsson effect */
				if (sym) {
					t_upbo[i][j] = (t_upbo[i][j]-t_lobo[i][j])/2.0;
					t_lobo[i][j] = -t_upbo[i][j];
					}
				}
			else
				t_lobo[i][j] = t_upbo[i][j] = -1.0;
		}
	if (!mode)
		rollback_V_base(0,ms_lobo,ms_upbo); // restore
	cst_on = global_cst;
	eval_cache_invalidate();
	return DTL_OK;
	}


static rcode dtl_get_V_tornado(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo) {
	rcode rc;
	int i,j,k,ts_nbr,global_cst;
	double baseline_ev;
	struct stmt_rec stmt;
	struct d_frame *df;

	if (load_df1(crit))
		return DTL_CRIT_UNKNOWN;
	df = uf->df;
	/* Set parameters */
	mode &= 0x01;
	/* Collect starting point */
	dtl_abort_init();
	if (call(TCL_get_V_hull(df,h_lobo,h_upbo),"TCL_get_V_hull"))
		return dtl_kernel_error();
	if (call(TCL_get_V_mbox(df,m_lobo,m_upbo),"TCL_get_V_mbox"))
		return dtl_kernel_error();
	if (call(TCL_get_V_mbox(df,ms_lobo,ms_upbo),"TCL_get_V_mbox"))
		return dtl_kernel_error();
	/* Lower cst reporting level */
	global_cst = cst_on;
	cst_on = cst_ext;
	if (!mode) {
		/* Mode 0: midpoint removed */
		for (k=1,i=1; i<=df->n_alts; i++)
			for (j=1; j<=df->tot_cons[i]; j++,k++) {
				m_lobo[k] = -1.0;
				m_upbo[k] = -1.0;
				}
		/* Clear mhull */
		if (call(TCL_set_V_mbox(df,m_lobo,m_upbo),"TCL_set_V_mbox")) {
			cst_on = global_cst;
			return dtl_kernel_error();
			}
		}
	/* Insert dummy V-statement */
	stmt.n_terms = 1;
	stmt.alt[1] = 1;
	stmt.cons[1] = 1;
	/* Find first re-node */
	while (!TCL_get_V_index(df,1,stmt.cons[1]))
		stmt.cons[1]++;
	stmt.sign[1] = 1;
	stmt.lobo = 0.0;
	stmt.upbo = 1.0;
	if (call(TCL_add_V_constraint(df,&stmt),"TCL_add_V_constraint")) {
		rollback_V_base(0,ms_lobo,ms_upbo);
		cst_on = global_cst;
		return dtl_kernel_error();
		}
	ts_nbr = df->V_base->n_stmts;
	for (k=1,i=1; i<=df->n_alts; i++) {
		dtl_abort_check();
		if (rc = evaluate_frame(crit,E_PSI,i,0,t_result)) {
			rollback_V_base(ts_nbr,ms_lobo,ms_upbo);
			cst_on = global_cst;
			return rc;
			}
		baseline_ev = t_result[E_MID][0];
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			if (TCL_get_V_index(df,i,j)) {
				/* Remove midpoint for real consequence */
				if (mode) {
					m_lobo[k] = -1.0;
					m_upbo[k] = -1.0;
					/* Set mhull */
					if (call(TCL_set_V_mbox(df,m_lobo,m_upbo),"TCL_set_V_mbox")) {
						rollback_V_base(ts_nbr,ms_lobo,ms_upbo);
						cst_on = global_cst;
						return dtl_kernel_error();
						}
					}
				if (h_upbo[k]-h_lobo[k] < 5.0*T_EPS)
					/* Narrow trap */
					t_lobo[i][j] = t_upbo[i][j] = 0.0;
				else {
					/* Enter statement */
					stmt.alt[1] = i;
					stmt.cons[1] = j;
					stmt.lobo = max(h_lobo[k],0.0);
					stmt.upbo = h_lobo[k]+T_EPS;
					if (call(TCL_replace_V_constraint(df,ts_nbr,&stmt),"TCL_replace_V_constraint")) {
						rollback_V_base(ts_nbr,ms_lobo,ms_upbo);
						cst_on = global_cst;
						return dtl_kernel_error();
						}
					if (dtl_abort_request) {
						rollback_V_base(ts_nbr,ms_lobo,ms_upbo);
						cst_on = global_cst;
						}
					dtl_abort_check();
					if (rc = evaluate_frame(crit,E_PSI,i,0,t_result)) {
						rollback_V_base(ts_nbr,ms_lobo,ms_upbo);
						cst_on = global_cst;
						return rc;
						}
					t_lobo[i][j] = t_result[E_MID][0] - baseline_ev;
					/* Change statement */
					if (call(TCL_change_V_constraint(df,ts_nbr,h_upbo[k]-T_EPS,min(h_upbo[k],1.0)),"TCL_change_V_constraint")) {
						rollback_V_base(ts_nbr,ms_lobo,ms_upbo);
						cst_on = global_cst;
						return dtl_kernel_error();
						}
					if (dtl_abort_request) {
						rollback_V_base(ts_nbr,ms_lobo,ms_upbo);
						cst_on = global_cst;
						}
					dtl_abort_check();
					if (rc = evaluate_frame(crit,E_PSI,i,0,t_result)) {
						rollback_V_base(ts_nbr,ms_lobo,ms_upbo);
						cst_on = global_cst;
						return rc;
						}
					t_upbo[i][j] = t_result[E_MID][0] - baseline_ev;
					/* Catch roundoff errors */
					if (t_lobo[i][j] > -T_EPS)
						t_lobo[i][j] = 0.0;
					if (t_upbo[i][j] < +T_EPS)
						t_upbo[i][j] = 0.0;
					}
				/* Prepare for next consequence */
				if (call(TCL_change_V_constraint(df,ts_nbr,0.0,1.0),"TCL_change_V_constraint")) {
					rollback_V_base(ts_nbr,ms_lobo,ms_upbo);
					cst_on = global_cst;
					return dtl_kernel_error();
					}
				if (mode) {
					/* Restore mhull */
					m_lobo[k] = ms_lobo[k];
					m_upbo[k] = ms_upbo[k];
					if (call(TCL_set_V_mbox(df,m_lobo,m_upbo),"TCL_set_V_mbox")) {
						rollback_V_base(ts_nbr,ms_lobo,ms_upbo);
						cst_on = global_cst;
						return dtl_kernel_error();
						}
					}
				}
			else {
				t_lobo[i][j] = t_upbo[i][j] = -1.0;
				}
			}
		}
	rollback_V_base(ts_nbr,ms_lobo,ms_upbo);
	eval_cache_invalidate();
	cst_on = global_cst;
	return DTL_OK;
	}


/* Mode: 0 = Midpoint kept (default)
 *       1 = Force floating midpoint
 *      +2 = Belief mass output */

rcode DTLAPI DTL_get_V_tornado(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo) {
	rcode rc;
	int i,j;

	/* Begin single thread semaphore */
	_smx_begin("TOV");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_V_tornado(%d,%d)\n",crit,mode);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(t_lobo,1);
	_certify_ptr(t_upbo,2);
	_dtl_assert(t_lobo!=t_upbo,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if ((mode < 0) || (mode > 3))
		return dtl_error(DTL_INPUT_ERROR);
	mode ^= 0x01; // flip lowest mode bit (compat reasons)
	/* Get tornado */
	if (rc = dtl_get_V_tornado(crit,mode,t_lobo,t_upbo))
		return dtl_error(rc);
	if (mode&0x02) {
		if (rc = dtl_mass_V_tornado(crit,mode,t_lobo,t_upbo))
			return dtl_error(rc);
		}
	/* Log function result */
	if (cst_ext)
		for (i=1; i<=uf->df->n_alts; i++)
			for (j=1; j<=uf->df->tot_cons[i]; j++) {
				sprintf(msg,"    V%d.%d.%-2d [%.3lf %.3lf]\n",crit,i,j,t_lobo[i][j],t_upbo[i][j]);
				cst_log(msg);
				}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


static rcode dtl_get_MCV_tornado(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo) {
	rcode rc;
	int i,j,node;

	/* Get local tornado in criterion */
	if (rc = dtl_get_V_tornado(crit,mode,t_lobo,t_upbo))
		return rc;
	if (mode&0x02) {
		if (rc = dtl_mass_V_tornado(crit,mode,t_lobo,t_upbo))
			return rc;
		}

	/* Get global weight for this criterion */
	if (load_df0(0))
		return DTL_SYS_CORRUPT;
	if (call(TCL_get_P_masspoint(uf->df,W_mid,LW_mid),"TCL_get_P_masspoint"))
		return dtl_kernel_error();
	/* Get global weight node index */
	node = r2t[1][crit];

	/* Apply to local criterion tornado */
	if (load_df1(crit))
		return DTL_CRIT_UNKNOWN;
	for (i=1; i<=uf->df->n_alts; i++)
		for (j=1; j<=uf->df->tot_cons[i]; j++)
			if (t_lobo[i][j] > -1.0) {
				/* Real node */
				t_lobo[i][j] *= W_mid[node];
				t_upbo[i][j] *= W_mid[node];
				}
	return DTL_OK;
	}


/* Mode: 0 = Midpoint kept (default)
 *       1 = Force floating midpoint
 *      +2 = Belief mass output */

rcode DTLAPI DTL_get_MCV_tornado(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo) {
	rcode rc;
	int i,j;

	/* Begin single thread semaphore */
	_smx_begin("TMCV");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_MCV_tornado(%d,%d)\n",crit,mode);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(t_lobo,1);
	_certify_ptr(t_upbo,2);
	_dtl_assert(t_lobo!=t_upbo,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (!PM)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if ((mode < 0) || (mode > 3))
		return dtl_error(DTL_INPUT_ERROR);
	mode ^= 0x01; // flip lowest mode bit (compat reasons)
	/* Get global tornado in criterion */
	if (rc = dtl_get_MCV_tornado(crit,mode,t_lobo,t_upbo))
		return dtl_error(rc);
	/* Log function result */
	if (cst_ext)
		for (i=1; i<=uf->df->n_alts; i++)
			for (j=1; j<=uf->df->tot_cons[i]; j++) {
				sprintf(msg,"    V%d.%d.%-2d [%.3lf %.3lf]\n",crit,i,j,t_lobo[i][j],t_upbo[i][j]);
				cst_log(msg);
				}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*********************************************************
  *
  *  Weight tornado
  *
  *********************************************************/

static int w_map[2*MAX_CRIT+1];
static d_row omega_ev;
/* WARNING: omega_ev needs 2*MAX_CRIT+MAX_ALTS+1 elements ->
 *          MAX_ALTS cannot be bloated to many hundreds
 * static double omega_ev[2*MAX_CRIT+MAX_ALTS+1]; */
static e_matrix ev_result;

static rcode dtl_get_W_tornado(int alt, int mode, h_matrix t_lobo, h_matrix t_upbo) {
	rcode rc;
	int i,j,k;

	if ((alt < 1) || (alt > uf->n_alts))
		return DTL_ALT_UNKNOWN;
	/* Collect all mass point expected values in sub-frames */
	for (j=1; j<=w_map[0]; j++) {
		if (k = w_map[j]) {
			/* Re-wt-node */
			if (load_df1(k)) {
				/* No criterion, do stand-in eval */
				omega_ev[j] = 0.5;
				}
			else {
				if (rc = evaluate_frame(k,E_PSI,alt,0,ev_result))
					return rc;
				omega_ev[j] = ev_result[E_MID][0];
				}
			}
		else
			/* Im-wt-node */
			omega_ev[j] = 0.0;
		}
	for (i=2; i<=uf->n_alts; i++)
		/* Dummy alternatives */
		omega_ev[j++] = 0.0;
	/* Load EV-box into V-base */
	if (load_df0(0))
		return DTL_SYS_CORRUPT;
	if (call(TCL_reset_V_base(uf->df),"TCL_reset_V_base"))
		return dtl_kernel_error();
	if (call(TCL_set_V_box(uf->df,omega_ev,omega_ev),"TCL_set_V_box"))
		return dtl_kernel_error();

	if (rc = dtl_get_PW_tornado(0,mode,t_lobo,t_upbo)) {
		TCL_unset_V_box(uf->df);
		return rc;
		}
	if (mode&0x02) {
		if (rc = dtl_mass_PW_tornado(-alt,mode,t_lobo,t_upbo)) {
			TCL_unset_V_box(uf->df);
			return rc;
			}
		}
	/* Unload EV-box from V-base */
	TCL_unset_V_box(uf->df);
	return DTL_OK;
	}


/* Mode: 0 = Midpoint kept (default)
 *       1 = Force floating midpoint
 *      +2 = Belief mass output */

rcode DTLAPI DTL_get_W_tornado(int mode, h_matrix t_lobo, h_matrix t_upbo) {
	rcode rc;
	int i,j,cst_global;

	/* Begin single thread semaphore */
	_smx_begin("TOW");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_W_tornado(%d)\n",mode);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(t_lobo,1);
	_certify_ptr(t_upbo,2);
	_dtl_assert(t_lobo!=t_upbo,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (!PM)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if (load_df0(0))
		return dtl_error(DTL_SYS_CORRUPT);
	if ((mode < 0) || (mode > 3))
		return dtl_error(DTL_INPUT_ERROR);
	mode ^= 0x01; // flip lowest mode bit (compat reasons)
	/* Initialise weight map vector */
	w_map[0] = uf->df->tot_cons[1];
	for (j=1; j<=w_map[0]; j++)
		w_map[j] = TCL_get_V_index(uf->df,1,j);
	/* Get weight tornado for all alternatives */
	for (j=1; j<=w_map[0]; j++) {
		t_lobo[0][j] = 0.0;
		t_upbo[0][j] = 0.0;
		}
	cst_global = cst_on;
	cst_on = FALSE;
	for (i=1; i<=uf->n_alts; i++) {
		if (rc = dtl_get_W_tornado(i,mode,x_lobo,x_upbo)) {
			cst_on = cst_global;
			return dtl_error(rc);
			}
		for (j=1; j<=w_map[0]; j++) {
			t_lobo[i][j]  = x_lobo[1][j];
			t_upbo[i][j]  = x_upbo[1][j];
			}
		}
	cst_on = cst_global;
	/* Log function result */
	if (cst_ext)
		for (i=1; i<=uf->n_alts; i++)
			for (j=1; j<=w_map[0]; j++) {
				sprintf(msg,"    A%d W%-2d [%.3lf %.3lf]\n",i,j,t_lobo[i][j],t_upbo[i][j]);
				cst_log(msg);
				}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


/* Mode: 0 = Midpoint kept (default)
 *       1 = Force floating midpoint
 *      +2 = Belief mass output
 *
 * Alt: >0 = Weight tornado for one alternative
 *      <0 = Weight tornado for single criteria node */

rcode DTLAPI DTL_get_W_tornado_alt(int alt, int mode, h_vector t_lobo, h_vector t_upbo) {
	rcode rc;
	int i,j,cst_global;

	/* This function is for compatibility */
	_smx_begin("TOWA");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_W_tornado_alt(%d,%d)\n",alt,mode);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(t_lobo,1);
	_certify_ptr(t_upbo,2);
	_dtl_assert(t_lobo!=t_upbo,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (!PM)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if (load_df0(0))
	  return dtl_error(DTL_SYS_CORRUPT);
	if ((mode < 0) || (mode > 3))
		return dtl_error(DTL_INPUT_ERROR);
	mode ^= 0x01; // flip lowest mode bit (compat reasons)
	/* Initialise weight map vector */
	w_map[0] = uf->df->tot_cons[1];
	for (j=1; j<=w_map[0]; j++)
		w_map[j] = TCL_get_V_index(uf->df,1,j);
	if (alt >= 0) {
		/* Get weight tornado for single alternative */
		if (rc = dtl_get_W_tornado(alt,mode,x_lobo,x_upbo))
			return dtl_error(rc);
		for (j=1; j<=w_map[0]; j++) {
			t_lobo[j] = x_lobo[1][j];
			t_upbo[j] = x_upbo[1][j];
			}
		}
	else {
		/* Get weight tornado for single tree node (-alt is node nbr) */
		if (-alt > w_map[0])
			return dtl_error(DTL_INPUT_ERROR);
		cst_global = cst_on;
		cst_on = FALSE;
		for (i=1; i<=uf->n_alts; i++) {
			if (rc = dtl_get_W_tornado(i,mode,x_lobo,x_upbo)) {
				cst_on = cst_global;
				return dtl_error(rc);
				}
			t_lobo[i] = x_lobo[1][-alt];
			t_upbo[i] = x_upbo[1][-alt];
			}
		cst_on = cst_global;
		}
	/* Log function result */
	if (cst_ext)
		for (j=1; j<=w_map[0]; j++) {
			if (w_map[j])
				sprintf(msg,"    W%-2d [%6.3lf %.3lf]  K%d\n",j,t_lobo[j],t_upbo[j],w_map[j]);
			else
				sprintf(msg,"    W%-2d [%6.3lf %.3lf]  *\n",j,t_lobo[j],t_upbo[j]);
			cst_log(msg);
			}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*********************************************************
  *
  *  Consequence influence
  *
  *********************************************************/

static rcode dtl_get_cons_influence(int crit, double mult, h_matrix result) {
	int i,j,k;
	struct d_frame *df;

	/* Check input parameters */
	if (load_df1(crit))
		return DTL_CRIT_UNKNOWN;
	df = uf->df;
	/* Get omega components */
	if (call(TCL_get_P_masspoint(df,W_mid,LW_mid),"TCL_get_P_masspoint"))
		return dtl_kernel_error();
	if (call(TCL_get_V_masspoint(df,V_mid),"TCL_get_V_masspoint"))
		return dtl_kernel_error();
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			if (V_mid[k] > -DTL_EPS)
				result[i][j] = mult*max(W_mid[k],0.0)*max(V_mid[k],0.0); // p*v
			else // im-node
				result[i][j] = -1.0;
			/* Log function result */
			if (cst_ext) {
				sprintf(msg,"    V%d.%d.%-2d %.3lf\n",crit,i,j,result[i][j]);
				cst_log(msg);
				}
			}
	eval_cache_invalidate();
	return DTL_OK;
	}


/* Crit: partial weight tree not supported (same as for tornados)
 * Output mode: 0 = local, 1 = global
 * DTL_OUTPUT_ERROR not supported due to lack of information */

rcode DTLAPI DTL_get_cons_influence(int crit, int mode, h_matrix result) {
	rcode rc;
	int node;
	double mult;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("CINF");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_cons_influence(%d,%d)\n",crit,mode);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(result,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if ((mode < 0) || (mode > 1))
		return dtl_error(DTL_INPUT_ERROR);
	if (PM && mode) {
		if (load_df0(0))
			return dtl_error(DTL_SYS_CORRUPT);
		df = uf->df;
		if (call(TCL_get_P_masspoint(df,W_mid,LW_mid),"TCL_get_P_masspoint"))
			return dtl_kernel_error();
		node = r2t[1][crit];
		mult = W_mid[node];
		}
	else // either PS or PM+local
		mult = 1.0;
	if (rc = dtl_get_cons_influence(crit,mult,result))
		return dtl_error(rc);
	/* End single thread semaphore */
	_smx_end();
	return rc;
	}
