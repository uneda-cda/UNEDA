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
 *   File: DTLeval.c
 *
 *   Purpose: evaluation of alternatives
 *
 *
 *   Functions exported outside DTL
 *   ------------------------------
 *   DTL_evaluate_frame
 *   DTL_evaluate_full
 *   DTL_evaluate_omega
 *   DTL_evaluate_omega1
 *   DTL_get_mass_above
 *   DTL_get_mass_below
 *   DTL_get_mass_range
 *   DTL_get_mass_density
 *   DTL_get_support_mass
 *   DTL_get_support_lower
 *   DTL_get_support_upper
 *   DTL_get_aversion_value
 *   DTL_compare_alternatives
 *   DTL_delta_mass
 *   DTL_rank_alternatives
 *   DTL_daisy_chain/1/2
 *   DTL_pie_chart/1/2
 *
 *   Functions outside of module, debug use
 *   --------------------------------------
 *   DTI_get_lo_inflexion
 *   DTI_get_up_inflexion
 *   DTI_bn_to_dtl_cdf
 *   DTI_get_support_mid
 *   DTI_get_mass_moments
 *   DTI_get_psi_moments
 *   DTI_get_bn_params
 *
 *   Functions outside of module, inside DTL
 *   ---------------------------------------
 *   eval_cache_invalidate
 *   evaluate_frame
 *   evaluate_frameset
 *   dtl_ev_to_cdf
 *   dtl_cdf_to_ev
 *
 *   Functions internal to module
 *   ----------------------------
 *   sort_b
 *   sq
 *   eval_cache_mass_init
 *   eval_cache_mass
 *   eval_cache_mc_mass
 *   get_cdf_ev
 *   expand_eval_result1/3
 *   dtl_evaluate_omega
 *   dtl_mass_validity
 *   dti_cdf_to_ev
 *   bp4corr_lo/up
 *   bn2dtl
 *   bn4dtl
 *   dtl2bn
 *   dtl4bn
 *   add_dom
 *
 */

#include "DTL.h"
#include "DTLinternal.h"


 /*******************************************************
  *
  *  Configuration parameters
  *  ------------------------
  *  DOMINANCE   include dominance-based evaluations
  *  DOM_2024    include 2024 dominance extensions
  *  Q_SORT      faster sort for large nbr of alts
  *
  *******************************************************/

#define noQ_SORT // for configs with up to 10000 alternatives

#define DOMINANCE // include dominance-based evaluations

#ifdef DOMINANCE
#define DOM_2024  // include 2024 dominance extension
#endif

#ifdef _MSC_VER
/* VC++ 6.0 compiler mistake: warns for non-existent problem numerous times:
 * "warning: local variable may be used without having been initialized" */
#pragma warning(disable:4701)
#endif


 /*********************************************************
  *
  *  Local data areas
  *
  *********************************************************/

static e_matrix e_cache[MAX_CRIT+1];
static int dtl_latest_mc_eval;
static a_result eval_result;


 /*********************************************************
  *
  *  Local functions
  *
  *********************************************************/

#ifdef Q_SORT

#define sort_b sort_dom2

#else // no Q_SORT

void sort_b(int order[], double maxmin[], int start, int stop, bool max) {
	int i,tmp;
	bool done,nok;

	/* Bubble sort, which is fast for few entries. But more importantly, the tolerance for
	 * equality is DTL_EPS which is the limit/horizon for round-off errors in calculations. */
	do {
		done = TRUE;
		for (i=start; i<=stop-1; i++) {
			nok = (max ? maxmin[order[i]] < maxmin[order[i+1]]-DTL_EPS :
									 maxmin[order[i]] > maxmin[order[i+1]]+DTL_EPS);
			if (nok) {
				/* Switch the two elements */
				done = FALSE;
				tmp = order[i];
				order[i] = order[i+1];
				order[i+1] = tmp;
				}
			}
		} while (!done);
	}

#endif // Q_SORT

/* Not a macro, to protect against double evaluation
 * of arguments - a stack push is most often cheaper */

static double sq(double a) {

	return a*a;
	}


 /*********************************************************
  *
  *  Evaluation cache
  *
  *********************************************************/

static struct bn_rec ec[MAX_CRIT+1];
static cr_col ecache_rm1;
static cr_col ecache_cm2;
static cr_col ecache_cm3;

void eval_cache_invalidate() {
	int j;

	if (frame_loaded)
		for (j=0; j<=uf->n_crit; j++)
			ec[j].valid = FALSE;
	}


static void eval_cache_mass_init() {
	int j;

	for (j=0; j<=uf->n_crit; j++) {
		ecache_rm1[j] = 0.0;
		ecache_cm2[j] = 0.0;
		ecache_cm3[j] = 0.0;
		ec[j].valid = FALSE;
		}
	}


static rcode eval_cache_mass(int crit, int method, int Ai, int Aj) {
	rcode rc;
	int m_field;
	int j,n_alts,n_active;
	struct d_frame *df;
	a_row rm1,cm2,cm3;
	double m1,m2,m3,skew=0.0,delta;

	df = uf->df;
	n_alts = df->n_alts;
	if (rc = TCL_get_moments(df,rm1,cm2,cm3))
		return rc;
	/* Process current evaluation rule */
	m_field = method & M_EVAL;
	switch (m_field) {
		case E_DELTA:
			m1 = rm1[Ai]-rm1[Aj];
			m2 = cm2[Ai]+cm2[Aj];
			m3 = cm3[Ai]-cm3[Aj];
			break;
		case E_GAMMA:
			m1 = rm1[Ai];
			m2 = cm2[Ai];
			m3 = cm3[Ai];
			for (j=1; j<=n_alts; j++)
				if (Ai!=j) {
					m1 -= rm1[j]/(double)(n_alts-1);
					m2 += cm2[j]/(double)(n_alts-1);
					m3 -= cm3[j]/(double)(n_alts-1);
					}
			break;
		case E_PSI:
			m1 = rm1[Ai];
			m2 = cm2[Ai];
			m3 = cm3[Ai];
			break;
		case E_DIGAMMA:
			m1 = m2 = m3 = 0.0;
			n_active = 0;
			for (j=1; j<=n_alts; j++)
				if ((Ai!=j) && (Aj&(0x01<<(j-1)))) {
					m1 -= rm1[j];
					m2 += cm2[j];
					m3 -= cm3[j];
					n_active++;
					}
			if (n_active) {
				m1 /= (double)n_active;
				m2 /= (double)n_active;
				m3 /= (double)n_active;
				}
			m1 += rm1[Ai];
			m2 += cm2[Ai];
			m3 += cm3[Ai];
			break;
		default:
			return DTL_WRONG_METHOD;
		}
	ecache_rm1[crit] = m1;
	ecache_cm2[crit] = m2;
	ecache_cm3[crit] = m3;
	if (m2 > DTL_EPS)
		skew = m3/pow(m2,1.5);
	else
		skew = 0.0;
	/* See math documentation for B-normal parameters */
	if (m2 > DTL_EPS)
		delta = sgn(skew)*b_delta(skew);
	else
		delta = 0.0;
	ec[crit].alpha = delta/sqrt(1.0-delta*delta);
	ec[crit].scale2 = m2/(1.0-2.0*delta*delta/PI);
	ec[crit].location = m1-sqrt(ec[crit].scale2)*delta*sqrt(2.0/PI);
	return DTL_OK;
	}


static rcode eval_cache_mc_mass(int snode, double V_rm1[], double V_cm2[], double V_cm3[]) {
	rcode rc;
	struct d_frame *df;
	double m1,m2,m3,skew=0.0,delta;

	/* NOTE: Implicit type cast cr_col -> d_row is ok since d_row is larger
	 *       and only MAX_CRIT entries are used by TCL_get_mc_moments */
	df = uf->df;
	if (rc = TCL_get_mc_moments(df,snode,V_rm1,V_cm2,V_cm3,&m1,&m2,&m3))
		return rc;
	ecache_rm1[0] = m1;
	ecache_cm2[0] = m2;
	ecache_cm3[0] = m3;
	/* See math documentation for B-normal parameters */
	if (m2 > DTL_EPS)
		delta = sgn(skew)*b_delta(skew);
	else
		delta = 0.0;
	ec[0].alpha =  delta/sqrt(1.0-delta*delta);
	ec[0].scale2 = m2/(1.0-2.0*delta*delta/PI);
	ec[0].location = m1-sqrt(ec[0].scale2)*delta*sqrt(2.0/PI);
	return DTL_OK;
	}


 /*********************************************************
  *
  *  Standard evaluation
  *
  *********************************************************/

 /*
  * Call semantics: All alternatives are evaluated using the delta, gamma,
  * psi, or digamma rules. For the requested alternative(s) Ai (and Aj), the
  * result is stored in eval_result. Note that for digamma, Aj is a bitmap.
  *
  * DTL layer 1: at DTL API level
  */

rcode evaluate_frame(int crit, int method, int Ai, int Aj, e_matrix e_result) {
	rcode trc;
	int m_field,eval_rule;
	int i;
	struct d_frame *df;

#ifdef INTERNAL_TEST
	// for testing error exit paths
	return dtl_error(DTL_INTERNAL_ERROR);
#endif
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (dtl_error_count)
		return dtl_error(DTL_OUTPUT_ERROR);
	df = uf->df;
	/* Process eval rule */
	m_field = method & M_EVAL;
	switch (m_field) {
		case E_DELTA:
			eval_rule = DELTA;   break;
		case E_GAMMA:
			eval_rule = GAMMA;   break;
		case E_PSI:
			eval_rule = PSI;     break;
		case E_DIGAMMA:
			eval_rule = DIGAMMA; break;
		default:
			return dtl_error(DTL_WRONG_METHOD);
		}
	/* Process alternatives */
	if ((Ai < 1) || (Ai > df->n_alts))
		return dtl_error(DTL_ALT_UNKNOWN);
	if (eval_rule == DELTA) {
		/* Check Aj also */
		if ((Aj < 1) || (Aj > df->n_alts))
			return dtl_error(DTL_ALT_UNKNOWN);
		if (Ai == Aj)
			return dtl_error(DTL_INPUT_ERROR);
		}
	else if (eval_rule < DIGAMMA)
		Aj = 0;
	/* Get TCL results */
	if (call(TCL_evaluate(df,Ai,Aj,eval_rule,eval_result),"TCL_evaluate")) {
		return dtl_kernel_error();
		}
	/* Transfer result to caller */
	e_cache[crit][E_MIN][0] = eval_result[Ai][E_MIN];
	e_cache[crit][E_MID][0] = eval_result[Ai][E_MID];
	e_cache[crit][E_MAX][0] = eval_result[Ai][E_MAX];
	trc = eval_cache_mass(crit,method,Ai,Aj);
	if (trc) {
		ec[crit].valid = FALSE;
		return dtl_error(DTL_INTERNAL_ERROR);
		}
	else
		ec[crit].valid = TRUE;
	/* Transfer results */
	e_result[E_MIN][0] = e_cache[crit][E_MIN][0];
	e_result[E_MID][0] = e_cache[crit][E_MID][0];
	e_result[E_MAX][0] = e_cache[crit][E_MAX][0];
	if (cst_ext) {
		for (i=E_MIN; i<=E_MAX; i++) {
			sprintf(msg," %6.3lf",e_result[i][0]);
			cst_log(msg);
			}
		cst_log("\n");
		}
	return DTL_OK;
	}


// VERY dangerous name resolution bug in the VC++ 6.0 compiler if variables begin with V_ and W_
// Resolved correctly by gcc - but _varies_ in different Microsoft VC++ versions! Horrible! :-(
// Names changed to Vc_ and Wc_ beginnings to circumvent the compiler bug (c for compiler bug).

static d_row Vc_lobo,Vc_upbo,Wc_point,im_Wc_point;

// DTL layer 1: at DTL API level

rcode evaluate_frameset(int crit, int method, int Ai, int Aj, e_matrix e_result) {
	rcode rc,drc,trc;
	int i,c,m_field;
	double minval,maxval;

	/* Criterion is loaded */
	eval_cache_mass_init();
	if (crit > 0)
		rc = evaluate_frame(crit,method,Ai,Aj,e_result);
	else {
		/* Full PM frame eval */
		dtl_abort_init();
		/* PM-0 MC evaluation */
		if ((Ai < 1) || (Ai > uf->df->n_alts))
			return dtl_error(DTL_ALT_UNKNOWN);
		if (method == E_DELTA) {
			/* Check Aj also */
			if ((Aj < 1) || (Aj > uf->df->n_alts))
				return dtl_error(DTL_ALT_UNKNOWN);
			if (Ai == Aj)
				return dtl_error(DTL_INPUT_ERROR);
			}
		m_field = method & M_EVAL;
		for (c=1; c<=uf->n_crit; c++) {
			rc = load_df1(c);
			if (rc == DTL_CRIT_UNKNOWN) {
				/* Stand-in evaluation for criterion with empty frame */
				if ((m_field == E_PSI) || ((m_field == E_DIGAMMA) && !Aj)) {
					Vc_upbo[c] = 1.0;
					Vc_lobo[c] = 0.0;
					ecache_rm1[c] = 0.5;
					ecache_cm2[c] = 1.0/24.0;
					}
				else {
					Vc_upbo[c] = 1.0;
					Vc_lobo[c] = -1.0;
					ecache_rm1[c] = 0.0;
					ecache_cm2[c] = 1.0/12.0;
					}
				if (cst_on)
					cst_log(" dtl_standin_eval: ok\n");
				}
			else if (rc)
				return dtl_error(rc);
			else {
				if (rc = evaluate_frame(c,method,Ai,Aj,e_result))
					return dtl_error(rc);
				Vc_upbo[c] = e_result[E_MAX][0];
				Vc_lobo[c] = e_result[E_MIN][0];
				dtl_abort_check();
				}
			}
		/* Find MC result */
		if (load_df0(0))
			return dtl_error(DTL_SYS_CORRUPT);
		/* PM tree */
		if (drc = call(TCL_get_P_min(uf->df,1,-crit,
				Vc_lobo,Wc_point,im_Wc_point,FALSE,&minval),"TCL_get_TP_min"))
			return dtl_error(DTL_KERNEL_ERROR+drc);
		if (drc = call(TCL_get_P_max(uf->df,1,-crit,
				Vc_upbo,Wc_point,im_Wc_point,TRUE,&maxval),"TCL_get_TP_max"))
			return dtl_error(DTL_KERNEL_ERROR+drc);
		e_result[E_MIN][0] = e_cache[0][E_MIN][0] = -minval;
		e_result[E_MAX][0] = e_cache[0][E_MAX][0] =  maxval;
		trc = eval_cache_mc_mass(-crit,ecache_rm1,ecache_cm2,ecache_cm3);
		if (trc) {
			ec[0].valid = FALSE;
			e_result[E_MID][0] = e_cache[0][E_MID][0] = (maxval-minval)/2.0;
			return dtl_error(DTL_INTERNAL_ERROR); // mass not ok
			}
		else {
			ec[0].valid = TRUE;
			e_result[E_MID][0] = e_cache[0][E_MID][0] = ecache_rm1[0];
			dtl_latest_mc_eval = crit;
			rc = DTL_OK; // mass ok
			}
		if (cst_on)
			cst_log(" dtl_evaluate_mc: ok\n");
		if (cst_ext) {
			for (i=E_MIN; i<=E_MAX; i++) {
				sprintf(msg," %6.3lf",e_result[i][0]);
				cst_log(msg);
				}
			cst_log("\n");
			}
		}
	return rc;
	}


 /* Each expansion has the form of a matrix {min,mid,max} x {steps},
  * with values from increasing support level steps. There are 21
  * values corresponding to belief levels of 0-100% (in 5% steps).
  *
  * expand_eval_result1 (0x0040): cone converging to 50% cdf
  * expand_eval_result2 (0x0080): cone converging to 50% cdf + swap midpoints
  * expand_eval_result3 (0x00C0): cone converging to mass point (default in eval)
  * expand_eval_result4 (0x0100): cone converging to mass point + interpolate
  *
  * NOTE: expand_eval_result1 is also used internally by the dominance
  * functions below and expand_eval_result3 is used by DTL_evaluate_rpf. */

static void expand_eval_result1(int crit, int swap, e_matrix e_result) {
	rcode rc;
	int i;
	double level,lobo,upbo;

	for (i=1; i<MAX_RESULTSTEPS; i++) {
		level = 1.0-(double)i/(double)(MAX_RESULTSTEPS-1);
		rc = dtl_cdf_to_ev(crit,max(level,1.0E-5),&lobo,&upbo);
		if (!rc) {
			e_result[E_MIN][i] = lobo;
			e_result[E_MID][i] = (lobo+upbo)/2.0;
			e_result[E_MAX][i] = upbo;
			}
		else {
			e_result[E_MIN][i] = -1.0;
			e_result[E_MID][i] = -1.0;
			e_result[E_MAX][i] = -1.0;
			}
		}
	if (swap) {
		/* Swap midpoints */
		e_result[E_MIN][MAX_RESULTSTEPS-1] = e_result[E_MID][0];
		e_result[E_MAX][MAX_RESULTSTEPS-1] = e_result[E_MID][0];
		e_result[E_MID][0] = e_result[E_MID][MAX_RESULTSTEPS-1];
		e_result[E_MID][MAX_RESULTSTEPS-1] = e_result[E_MAX][MAX_RESULTSTEPS-1];
		}
	}


static rcode get_cdf_ev(int crit, double cdf, double *ev) {
	double ev2;

	if (cdf < 0.5)
		return dtl_cdf_to_ev(crit,min(max(1.0-2.0*cdf,1.0E-5),0.999),ev,&ev2);
	else
		return dtl_cdf_to_ev(crit,min(max(2.0*cdf-1.0,1.0E-5),0.999),&ev2,ev);
	}


static rcode expand_eval_result3(int crit, int ip, e_matrix e_result) {
	int i;
	double level,lobo,upbo,mid_pdf,shift,step;

	if (dtl_ev_to_cdf(crit,e_result[E_MID][0],&mid_pdf))
		return DTL_INTERNAL_ERROR;
	e_result[E_MID][0] = (e_result[E_MIN][0]+e_result[E_MAX][0])/2.0;
	for (i=1; i<MAX_RESULTSTEPS; i++) {
		level = 1.0-(double)i/(double)(MAX_RESULTSTEPS-1);
		shift = (1.0-level)*(mid_pdf-0.5);
		if (get_cdf_ev(crit,(1.0-level)/2.0-shift,&lobo))
			return DTL_INTERNAL_ERROR;
		e_result[E_MIN][i] = lobo;
		if (get_cdf_ev(crit,(1.0+level)/2.0-shift,&upbo))
			return DTL_INTERNAL_ERROR;
		e_result[E_MAX][i] = upbo;
		e_result[E_MID][i] = (lobo+upbo)/2.0;
		}
	if (ip) { // interpolate first (base) entry to constant derivative
		step = sq(e_result[E_MIN][2]-e_result[E_MIN][1])/(e_result[E_MIN][3]-e_result[E_MIN][2]);
		e_result[E_MIN][0] = max(e_result[E_MIN][0],e_result[E_MIN][1]-step);
		step = sq(e_result[E_MAX][1]-e_result[E_MAX][2])/(e_result[E_MAX][2]-e_result[E_MAX][3]);
		e_result[E_MAX][0] = min(e_result[E_MAX][0],e_result[E_MAX][1]+step);
		e_result[E_MID][0] = (e_result[E_MIN][0]+e_result[E_MAX][0])/2.0;
		}
	return DTL_OK;
	}


rcode DTLAPI DTL_evaluate_frame(int crit, int method, int Ai, int Aj, e_matrix e_result) {
	rcode rc;

	/* Begin single thread semaphore */
	_smx_begin("EVAL");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_evaluate_frame(%d,%d,%d,%d)\n",crit,method,Ai,Aj);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(e_result,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df00(crit)) // must validate input here
		return dtl_error(DTL_CRIT_UNKNOWN);
	/* Evaluate */
	rc = evaluate_frameset(crit,method,Ai,Aj,e_result);
	/* End single thread semaphore */
	if (!rc)
		_smx_end();
	return rc;
	}


// DTL layer 0: above DTL proper

rcode DTLAPI DTL_evaluate_full(int crit, int method, int Ai, int Aj, e_matrix e_result) {
	rcode rc;
	int i,j,eval_method,exp_mode;

	eval_method = method & M_EVAL;
	/* Note that due to piggybacking, cst_log will not show the exp_mode bits in
	   the function call. They will instead be shown in the log result below. */
	exp_mode = (method-eval_method)>>6; // stored from 0x40 up
	/* Evaluate */
	if (rc = DTL_evaluate_frame(crit,eval_method,Ai,Aj,e_result))
		return rc;
	/* Expand to many support levels */
	switch (exp_mode) {
		case 1:
		case 2:
			expand_eval_result1(crit,exp_mode-1,e_result);
			break;
		case 0: // mode 3 is the default
			exp_mode = 3;
		case 3:
		case 4:
			rc = expand_eval_result3(crit,exp_mode-3,e_result);
			break;
		default:
			rc = DTL_WRONG_METHOD;
		}
	/* Log function result */
	if (cst_ext) {
		sprintf(msg," expand_eval_result%d: %s\n",exp_mode,rc<0?DTL_get_errtxt(rc):"ok");
		cst_log(msg);
		if (!rc)
			for (i=E_MIN; i<=E_MAX; i++) {
				for (j=0; j<MAX_RESULTSTEPS; j+=4) {
					sprintf(msg," %6.3lf",e_result[i][j]);
					cst_log(msg);
					}
				cst_log("\n");
				}
		}
	return rc;
	}


static d_row W_mid,LW_mid;
static i_row t_inx;

// DTL layer 1: at DTL API level

static rcode dtl_evaluate_omega(int Ai, cr_col o_result) {
	rcode rc;
	int c;
	double omega;

	if (load_df0(0))
		return dtl_error(DTL_SYS_CORRUPT);
	if (dtl_error_count)
		return dtl_error(DTL_OUTPUT_ERROR);
	/* Check input parameter */
	if ((Ai < 1) || (Ai > uf->n_alts))
		return dtl_error(DTL_ALT_UNKNOWN);
	/* Get MC weights */
	if (call(TCL_get_P_masspoint(uf->df,W_mid,LW_mid),"TCL_get_P_masspoint"))
		return dtl_kernel_error();
	/* Collect index type A2 from A1 */
	for (c=1; c<=uf->n_crit; c++)
		if (!(t_inx[c] = TCL_get_tot_index(1,c)))
			return dtl_error(DTL_INTERNAL_ERROR);
	o_result[0] = 0.0;
	/* Get omega from each criterion */
	for (c=1; c<=uf->n_crit; c++) {
		rc = load_df1(c);
		if (rc == DTL_CRIT_UNKNOWN)
			/* Stand-in evaluation for empty frame */
			omega = 0.5;
		else if (rc)
			return dtl_error(rc);
		else
			if (call(TCL_evaluate_omega(uf->df,Ai,&omega),"TCL_evaluate_omega"))
				return dtl_kernel_error();
		o_result[c] = W_mid[t_inx[c]] * omega;
		o_result[0] += o_result[c];
		/* Log function result */
		if (cst_ext) {
			sprintf(msg," W%-2d %.3lf -> %.3lf\n",c,W_mid[c],o_result[c]);
			cst_log(msg);
			}
		}
	/* Log function total result */
	if (cst_ext) {
		sprintf(msg," Tot 1.000 -> %.3lf\n",o_result[0]);
		cst_log(msg);
		}
	return DTL_OK;
	}


/* Output mode: 0=order
 *              1=olympic rank
 *              2=strict rank
 *              3=group rank
 *             +4=percent of omega EV (default: percent of entire scale) */

static cr_col o2_result,o3_result;
static ci_col o_order;

rcode DTLAPI DTL_evaluate_omega(int Ai, int mode, cr_col o_result, ci_col o_rank) {
	rcode rc;
	int i,j,level,renorm,cst_global;

	/* Begin single thread semaphore */
	_smx_begin("OMEGA");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_evaluate_omega(%d,%d)\n",Ai,mode);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(o_result,1);
	_certify_ptr(o_rank,2);
	_dtl_assert(o_result!=o2_result,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (PS)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check and split input parameter */
	if ((mode < 0) || (mode > 7))
		return dtl_error(DTL_INPUT_ERROR);
	renorm = mode&0x04;
	mode &= 0xFB;
	/* Evaluate one alt or all */
	if (Ai)
		rc = dtl_evaluate_omega(Ai,o_result);
	else {
		cst_global = cst_on;
		cst_on = FALSE;
		rc = dtl_evaluate_omega(1,o_result);
		if (!rc) {
			for (i=2; i<=uf->n_alts; i++) {
				if (rc = dtl_evaluate_omega(i,o2_result))
					break;
				for (j=0; j<=uf->n_crit; j++)
					o_result[j] += o2_result[j];
				}
			for (j=0; j<=uf->n_crit; j++)
				o_result[j] /= (double)uf->n_alts;
			}
		cst_on = cst_global;
		}
	if (rc) // catch eval error
		return dtl_error(rc);
	/* Sort omega in descending order */
	for (j=1; j<=uf->n_crit; j++)
		o_order[j] = j;
	sort_b(o_order,o_result,1,uf->n_crit,TRUE);
	if (mode) {
		/* Convert order to ranking. Mode 1 is olympic weak ranking:
		 * if two share the gold medal, no silver medal is awarded. */
		o_rank[o_order[1]] = level = 1;
		for (j=2; j<=uf->n_crit; j++)
			if ((o_result[o_order[j-1]]-o_result[o_order[j]]) > DTL_EPS)
				o_rank[o_order[j]] = mode>2?++level:j;
			else
				o_rank[o_order[j]] = mode>1?(mode>2?level:j):o_rank[o_order[j-1]];
		}
	else
		/* Transfer order verbatim (mode 0) */
		for (j=1; j<=uf->n_crit; j++)
			o_rank[j] = o_order[j];
	o_rank[0] = uf->n_crit; // nbr of entries
	/* Renormalise if percent of EV */
	if (renorm)
		for (j=1; j<=uf->n_crit; j++)
			o_result[j] /= o_result[0];
	/* Log function result (per real crit/end node)
	 *
	 * For mode=0 (order)
	 * First row: T0 Nm 0.nnn, m=nbr of crit, n=total omega
	 * Next rows: Rp Ks 0.nnn, p=order pos, n=partial omega, s=crit nbr
	 *
	 * For mode>0 (ranking)
	 * First row: T0 Nm 0.nnn, m=nbr of crit, n=total omega
	 * Next rows: Ks Rp 0.nnn, p=order pos, n=partial omega, s=crit nbr
	 */
	if (cst_ext)
		for (j=0; j<=uf->n_crit; j++) {
			sprintf(msg," %c%-2d %c%-2d  %.3lf\n",j?(mode?'K':'R'):'T',j,
					j?(mode?'R':'K'):'N',o_rank[j],o_result[j?(mode?j:o_rank[j]):0]);
			cst_log(msg);
			}
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


/* Omega1 is the total contribution from each node at
 * the first weight tree level (akin to part worth)
 *
 * Output mode: 0 = percent of entire scale
 *              4 = percent of omega EV (renormalise) */

rcode DTLAPI DTL_evaluate_omega1(int Ai, int mode, cr_col o_result, ci_col o_node) {
	rcode rc;
	int i,j,k,cst_global;
	int pos,wnode,wnode1,state;
	double ev_sum;

	/* Begin single thread semaphore */
	_smx_begin("OMEGA1");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_evaluate_omega1(%d,%d)\n",Ai,mode);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(o_result,1);
	_certify_ptr(o_node,2);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (PS)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameter */
	if (mode & 0xFB)
		return dtl_error(DTL_INPUT_ERROR);
	/* Evaluate */
	if (Ai)
		rc = dtl_evaluate_omega(Ai,o2_result);
	else {
		cst_global = cst_on;
		cst_on = FALSE;
		rc = dtl_evaluate_omega(1,o2_result);
		if (!rc) {
			for (i=2; i<=uf->n_alts; i++) {
				if (rc = dtl_evaluate_omega(i,o3_result))
					break;
				for (j=0; j<=uf->n_crit; j++)
					o2_result[j] += o3_result[j];
				}
			for (j=0; j<=uf->n_crit; j++)
				o2_result[j] /= (double)uf->n_alts;
			}
		cst_on = cst_global;
		}
	if (rc) // catch eval error
		return dtl_error(rc);
	/* Collect omega parts at first tree level, adding at deeper levels */
	pos = 1;
	ev_sum = 0.0;
	/* Traverse the tree and collect omega parts in all re-nodes encountered.
	 * State machine with two states: vacuuming at level 1 (L1) or deeper (L2+). */
	for (wnode=state=1; wnode<=uf->df->tot_cons[1]; wnode++) {
		k = dtl_node2crit(wnode);
		if (!dtl_W_node_parents(1,wnode)) {
			// L1
			if (state > 1) {
				o_result[pos] = ev_sum;
				o_node[pos++] = wnode1;
				ev_sum = 0.0;
				state = 1; // state change -> S1
				}
			if (k) {
				// L1 & re-node
				o_result[pos] = o2_result[k];
				o_node[pos++] = wnode;
				}
			else {
				// L1 & im-node
				wnode1 = wnode;
				state = 2; // state change -> S2
				}
			}
		else // L2+
			if (k)
				// L2+ & re-node
				ev_sum += o2_result[k];
		}
	if (state > 1) {
		// collect the last L2+ entry
		o_result[pos] = ev_sum;
		o_node[pos++] = wnode1; // bug in VC++ 6.0 compiler "used w/o init"
		}
	// result entry 0 contains total omega (also for Ai=0 average)
	o_result[0] = o2_result[0];
	// node entry 0 contains nbr of entries (nodes at level 1)
	o_node[0] = pos-1;
	/* Renormalise if percent of EV */
	if (mode)
		for (j=1; j<=o_node[0]; j++)
			o_result[j] /= o_result[0];
	/* Log function result (per level 1 node - whether crit or im node)
	 * First row: Tm 0.nnn, m=nbr of crit, n=total omega
	 * Next rows: Ws 0.nnn, s=wt node nbr, n=partial omega */
	if (cst_ext)
		for (j=0; j<pos; j++) {
			sprintf(msg," %c%-2d %.3lf\n",j?'W':'T',o_node[j],o_result[j]);
			cst_log(msg);
			}
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*************************************************************************
  *
  *  Mass calculation functions are of two different types:
  *
  *          DTL API                  Application output
  *          -------                  ------------------
  *  Type 1: have EV level (range) -> want to know belief mass in %
  *          DTL_get_mass_above       mass % above EV level (inv. cdf)
  *          DTL_get_mass_below       mass % below EV level (inv. cdf)
  *          DTL_get_mass_range       mass % between 2 EV levels (inv. cdf)
  *          DTL_get_mass_density     density (integrand function)
  *
  *  Type 2: have belief mass in % -> want to know EV interval (range)
  *          DTL_get_support_mass     EV interval centred around midpoint
  *          DTL_get_support_lower    EV interval from lower scale endpoint
  *          DTL_get_support_upper    EV interval from upper scale endpoint
  *
  *************************************************************************/

/* Definition of infinite mass point boundaries */

#define INF_MASS_VAR 1.0E-8
#define MIN_SUPPORT_LEVEL 1.0E-5
#define MAX_SUPPORT_LEVEL 0.9990234375 // risk aversion = 10


static double bn2dtl(double ref_lo, double ref_up, double bn) {

	if (ref_up > ref_lo+DTL_EPS)
		return (bn-ref_lo)/(ref_up-ref_lo);
	else
		return 0.5; // degenerated (Dirac) case
	}

static double dtl2bn(double ref_lo, double ref_up, double dtl) {

	if (ref_up > ref_lo+DTL_EPS)
		return ref_lo+dtl*(ref_up-ref_lo);
	else
		return 0.5; // degenerated
	}


 /**********************************************************************
  *
  *  Internal mass distribution validity and verification
  *
  **********************************************************************/

static rcode dtl_mass_validity(int crit) {
	double ref_lo,ref_up;

	if (e_cache[crit][E_MAX][0]-e_cache[crit][E_MIN][0] < DTL_EPS)
		return DTL_INFINITE_MASS; // almost infinite mass
	/* "Valid" if skew does not make interpolation do too much to the function */
	ref_lo = bn_cdf(e_cache[crit][E_MIN][0],ec[crit].location,ec[crit].scale2,ec[crit].alpha);
	ref_up = bn_cdf(e_cache[crit][E_MAX][0],ec[crit].location,ec[crit].scale2,ec[crit].alpha);
	if (ref_up-ref_lo < 0.9)
		return DTL_WEAK_MASS_DISTR; // spans less than 90%
	return DTL_OK;
	}


 /**********************************************************************
  *
  *  Type 0: have b-normal cdf (mass), want to know support (mass) in %.
  *  How many % of the belief mass on the DTL truncated cdf do we have?
  *
  **********************************************************************/

rcode DTLAPI DTI_bn_to_dtl_cdf(int crit, double cdf_bn, double *cdf_dtl) {
	double ref_lo,ref_up;

	/* Begin single thread semaphore */
	_smx_begin("BNDTL");
	/* Log function call */
	if (cst_ext) {
		sprintf(msg,"DTI_bn_to_dtl_cdf(%d,%.3lf)\n",crit,cdf_bn);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(cdf_dtl,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df00(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if ((crit < 1) && (dtl_latest_mc_eval-crit))
		return dtl_error(DTL_OUTPUT_ERROR);
	crit = max(crit,0); // partial/subcrit trees in slot 0
	if (!ec[crit].valid)
		return dtl_error(DTL_OUTPUT_ERROR);
	if ((cdf_bn < 0.0) || (cdf_bn > 1.0))
		return dtl_error(DTL_INPUT_ERROR);
	/* Calculate DTL cdf for given b-normal cdf */
	ref_lo = bn_cdf(e_cache[crit][E_MIN][0],ec[crit].location,ec[crit].scale2,ec[crit].alpha);
	ref_up = bn_cdf(e_cache[crit][E_MAX][0],ec[crit].location,ec[crit].scale2,ec[crit].alpha);
	if (cdf_bn < ref_lo)
		*cdf_dtl = 0.0;
	else if (cdf_bn > ref_up)
		*cdf_dtl = 1.0;
	else
		/* Convert from b-normal cdf scale to DTL cdf scale */
		*cdf_dtl = bn2dtl(ref_lo,ref_up,cdf_bn);
	/* Log function result */
	if (cst_ext) {
		sprintf(msg," cdf_dtl = %6.3lf\n",*cdf_dtl);
		cst_log(msg);
		}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*******************************************************************
  *
  *  Type 1: have interval (range), want to know support (mass) in %.
  *  How many % of the belief mass fall _above_ the EV supplied, i.e
  *  where on the cdf does this EV reside? The cdf originally tells
  *  how much has already been passed or is to the left of the curve,
  *  i.e. the mass _below_ the point.
  *
  *******************************************************************/

/* Dirac splitting yields the semantics "above", otherwise "at or above".
 * Range mass function requires the "above" definition to work properly. */

#define DIRAC_SPLIT
#define DENS_EPS 1.0E-6

rcode dtl_ev_to_cdf(int crit, double ev_level, double *mass) {
	double ref_lo,ref_up,cur;

	/* Check input parameters */
	if (load_df00(crit))
		return DTL_CRIT_UNKNOWN;
	if ((crit < 1) && (dtl_latest_mc_eval-crit))
		return DTL_OUTPUT_ERROR;
	crit = max(crit,0); // partial/subcrit trees in slot 0
	if (!ec[crit].valid)
		return DTL_OUTPUT_ERROR;
	if ((ev_level < -1.0) || (ev_level > 1.0))
		return DTL_INPUT_ERROR;
	/* Calculate mass above the given level */
	if (ev_level < e_cache[crit][E_MIN][0]-DTL_EPS)
		*mass = 1.0;
	else if (ev_level > e_cache[crit][E_MAX][0]+DTL_EPS)
		*mass = 0.0;
	else if ((e_cache[crit][E_MAX][0]-e_cache[crit][E_MIN][0] > DTL_EPS) && (ecache_cm2[crit] > INF_MASS_VAR)) {
		ref_lo = bn_cdf(e_cache[crit][E_MIN][0],ec[crit].location,ec[crit].scale2,ec[crit].alpha);
		ref_up = bn_cdf(e_cache[crit][E_MAX][0],ec[crit].location,ec[crit].scale2,ec[crit].alpha);
		cur = bn_cdf(ev_level,ec[crit].location,ec[crit].scale2,ec[crit].alpha);
		/* Catch roundoff errors */
		cur = min(max(cur,ref_lo),ref_up);
		/* Convert from b-normal cdf scale to DTL cdf scale */
		*mass = max(1.0-bn2dtl(ref_lo,ref_up,cur),0.0);
		}
	/* If the total area is (close to) zero, return position relative to mass point instead */
	else if (ev_level < ecache_rm1[crit]-DTL_EPS)
		*mass = 1.0;
	else if (ev_level > ecache_rm1[crit]+DTL_EPS)
		*mass = 0.0;
#ifdef DIRAC_SPLIT
	else { // (almost) at mass point, split above and below
		if (ev_level < -1.0+DTL_EPS)
			*mass = 1.0;
		else if (ev_level > 1.0-DTL_EPS)
			*mass = 0.0;
		else
			*mass = 0.5;
		}
#else
	else // at or above is always full mass, no split
		*mass = 1.0;
#endif
	return DTL_OK; // dtl_mass_validity at next level up
	}


rcode DTLAPI DTL_get_mass_above(int crit, double lo_level, double *mass) {
	rcode rc;

	/* Begin single thread semaphore */
	_smx_begin("AMASS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_mass_above(%d,%.3lf)\n",crit,lo_level);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(mass,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Fetch the "reverse cdf" of the evaluation */
	if (rc = dtl_ev_to_cdf(crit,lo_level,mass))
		return dtl_error(rc);
	/* Log function result */
	if (cst_ext) {
		sprintf(msg," mass above = %6.3lf\n",*mass);
		cst_log(msg);
		}
	/* End single thread semaphore */
	_smx_end();
	return dtl_mass_validity(crit);
	}


rcode DTLAPI DTL_get_mass_below(int crit, double up_level, double *mass) {
	rcode rc;

	/* Begin single thread semaphore */
	_smx_begin("BMASS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_mass_below(%d,%.3lf)\n",crit,up_level);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(mass,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Fetch the "reverse cdf" of the evaluation */
	if (rc = dtl_ev_to_cdf(crit,up_level,mass))
		return dtl_error(rc);
	/* Revert to "true cdf" = mass below (up to) specified level */
	*mass = 1.0-*mass;
	/* Log function result */
	if (cst_ext) {
		sprintf(msg," mass below = %6.3lf\n",*mass);
		cst_log(msg);
		}
	/* End single thread semaphore */
	_smx_end();
	return dtl_mass_validity(crit);
	}


rcode DTLAPI DTL_get_mass_range(int crit, double lo_level, double up_level, double *mass) {
	rcode rc;
	double lo_mass,up_mass;

	/* Begin single thread semaphore */
	_smx_begin("RMASS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_mass_range(%d,%.3lf,%.3lf)\n",crit,lo_level,up_level);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(mass,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if ((lo_level < -1.0) || (up_level > 1.0))
		return dtl_error(DTL_INPUT_ERROR);
	if (lo_level > up_level)
		return dtl_error(DTL_INPUT_ERROR);
	/* Fetch the lower and upper "cdf" of the evaluation */
	if (rc = dtl_ev_to_cdf(crit,max(lo_level-2.0*DTL_EPS,-1.0),&lo_mass))
		return dtl_error(rc);
	if (rc = dtl_ev_to_cdf(crit,min(up_level+2.0*DTL_EPS,1.0),&up_mass))
		return dtl_error(rc);
	*mass = lo_mass-up_mass;
	/* Log function result */
	if (cst_ext) {
		sprintf(msg," range mass = %6.3lf\n",*mass);
		cst_log(msg);
		}
	/* End single thread semaphore */
	_smx_end();
	return dtl_mass_validity(crit);
	}


rcode DTLAPI DTL_get_mass_density(int crit, double ev_level, double *density) {
	rcode rc;
	int crit0;
	double level1,level2,mass1,mass2;

	/* Begin single thread semaphore */
	_smx_begin("MDENS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_mass_density(%d,%.3lf)\n",crit,ev_level);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(density,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df00(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if ((crit < 1) && (dtl_latest_mc_eval-crit))
		return dtl_error(DTL_OUTPUT_ERROR);
	crit0 = max(crit,0); // partial/subcrit trees in slot 0
	if (!ec[crit0].valid)
		return dtl_error(DTL_OUTPUT_ERROR);
	if ((ev_level < -1.0) || (ev_level > 1.0))
		return dtl_error(DTL_INPUT_ERROR);
	/* Fetch the "pdf" of the evaluation */
	level1 = min(ev_level+DENS_EPS,1.0);
	if (rc = dtl_ev_to_cdf(crit,level1,&mass1))
		return dtl_error(rc);
	level2 = max(ev_level-DENS_EPS,-1.0);
	if (rc = dtl_ev_to_cdf(crit,level2,&mass2))
		return dtl_error(rc);
	if ((level2 > e_cache[crit0][E_MIN][0]) && 
			(level1 < e_cache[crit0][E_MAX][0]) && 
			(level1 > level2)) {
		*density = (mass2-mass1)/(level1-level2);
		if (*density < 1.0E-10) // catch roúnd-off
			*density = 0.0;
		}
	else if ((fabs(e_cache[crit0][E_MIN][0]-e_cache[crit0][E_MAX][0]) < DTL_EPS) && 
					 (fabs(e_cache[crit0][E_MIN][0]-ev_level) < DTL_EPS))
#ifdef INFINITY
		*density = INFINITY; // infinite Dirac delta density
#elif defined(DBL_MAX)
		*density = DBL_MAX;
#else
		*density = 1.79E308; // "almost infinite" (max is around 1.7977E308 in IEEE754 binary64)
#endif
	else
		*density = 0.0; // outside of range
	/* Log function result */
	if (cst_ext) {
		sprintf(msg," density = %6.3lf\n",*density);
		cst_log(msg);
		}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*********************************************************************
  *
  *  Type 2: have support (mass) in %, want to know interval (EV range).
  *  The interval covering the central 'belief_level' % on the EV axis.
  *  bn_cdf is a one-way function, loop through EV values ever closer.
  *
  *  Unfortunately, there exist no closed analytical PDF expression.
  *  In DTLbnormal.c, there is an analytical inverse to the normal
  *  distribution CDF but not to Owen's T-function distribution CDF.
  *  Thus, the Owen's T inverse call is a loop as well -> no real gain
  *  calling that compared to performing the calculations in DTLeval.c.
  *  For historical reasons, the DTL_ version uses a local loop while
  *  the DTI_ version instead makes use of the DTLbnormal.c function.
  *
  *********************************************************************/

rcode dtl_cdf_to_ev(int crit, double belief_level, double *lobo, double *upbo) {
	double cdf,ref_lo,ref_up,step,val,new_val,lo_target,up_target;

	/* Check input parameters */
	if (load_df00(crit))
		return DTL_CRIT_UNKNOWN;
	if ((crit < 1) && (dtl_latest_mc_eval-crit))
		return DTL_OUTPUT_ERROR;
	crit = max(crit,0); // subcrit trees in slot 0
	if (!ec[crit].valid)
		return DTL_OUTPUT_ERROR;
	if ((belief_level < MIN_SUPPORT_LEVEL) || (belief_level > MAX_SUPPORT_LEVEL))
		return DTL_INPUT_ERROR; // out-of-range, round-off errors will start interfering
	/* Calculate support mass */
	if ((e_cache[crit][E_MAX][0]-e_cache[crit][E_MIN][0] > DTL_EPS) && (ecache_cm2[crit] > INF_MASS_VAR)) {
		ref_lo = bn_cdf(e_cache[crit][E_MIN][0],ec[crit].location,ec[crit].scale2,ec[crit].alpha);
		ref_up = bn_cdf(e_cache[crit][E_MAX][0],ec[crit].location,ec[crit].scale2,ec[crit].alpha);
		/* Lower bound */
		step = (e_cache[crit][E_MAX][0]-e_cache[crit][E_MIN][0])/2.0;
		new_val = (e_cache[crit][E_MAX][0]+e_cache[crit][E_MIN][0])/2.0;
		/* Convert from DTL cdf scale to b-normal cdf scale */
		lo_target = dtl2bn(ref_lo,ref_up,(1.0-belief_level)/2.0);
		do {
			val = new_val;
			cdf = bn_cdf(val,ec[crit].location,ec[crit].scale2,ec[crit].alpha);
			step /= 2.0;
			if (cdf > lo_target)
				new_val = val-step;
			else
				new_val = val+step;
			} while ((fabs(cdf-lo_target) > 0.000001) && (step > 1.0E-7));
		*lobo = val;
		/* Upper bound */
		step = (e_cache[crit][E_MAX][0]-e_cache[crit][E_MIN][0])/2.0;
		new_val = (e_cache[crit][E_MAX][0]+e_cache[crit][E_MIN][0])/2.0;
		/* Convert from DTL cdf scale to b-normal cdf scale */
		up_target = dtl2bn(ref_lo,ref_up,(1.0+belief_level)/2.0);
		do {
			val = new_val;
			cdf = bn_cdf(val,ec[crit].location,ec[crit].scale2,ec[crit].alpha);
			step /= 2.0;
			if (cdf > up_target)
				new_val = val-step;
			else
				new_val = val+step;
			} while ((fabs(cdf-up_target) > 0.000001) && (step > 1.0E-7));
		*upbo = val;
		}
	else {
		*lobo = ecache_rm1[crit];
		*upbo = ecache_rm1[crit];
		}
	return DTL_OK;
	}


/* Collector functions: collect the mass from a previous evaluation.
 * Valid belief levels are [50%,99.9%]. (100% is the entire EV range) */

/* DTL_get_support_mass:
 * Get the centred interval around the 50% mass point
 * that contains 'belief_level' percent of the mass */

rcode DTLAPI DTL_get_support_mass(int crit, double belief_level, double *lobo, double *upbo) {
	rcode rc;

	/* Begin single thread semaphore */
	_smx_begin("SMASS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_support_mass(%d,%.3lf)\n",crit,belief_level);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(lobo,1);
	_certify_ptr(upbo,2);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (belief_level < 0.5)
		return dtl_error(DTL_INPUT_ERROR);
	/* Fetch the support mass of the evaluation */
	if (rc = dtl_cdf_to_ev(crit,belief_level,lobo,upbo))
		return dtl_error(rc);
	/* Log function result */
	if (cst_ext) {
		sprintf(msg," [%6.3lf %6.3lf]\n",*lobo,*upbo);
		cst_log(msg);
		}
	/* End single thread semaphore */
	_smx_end();
	return dtl_mass_validity(crit);
	}


/* DTL_get_support_lower:
 * Get the interval from 0% mass (lower hull) up to 'belief_level' percent of the mass */

rcode DTLAPI DTL_get_support_lower(int crit, double belief_level, double *lobo, double *upbo) {
	rcode rc;

	/* Begin single thread semaphore */
	_smx_begin("SMASL");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_support_lower(%d,%.3lf)\n",crit,belief_level);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(lobo,1);
	_certify_ptr(upbo,2);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (belief_level < 0.5)
		return dtl_error(DTL_INPUT_ERROR);
	/* Fetch the support mass of the evaluation */
	if (rc = dtl_cdf_to_ev(crit,2.0*belief_level-1.0,lobo,upbo))
		return dtl_error(rc);
	crit = max(crit,0); // partial/subcrit trees in slot 0
	*lobo = e_cache[crit][E_MIN][0];
	/* Log function result */
	if (cst_ext) {
		sprintf(msg," [%6.3lf %6.3lf]\n",*lobo,*upbo);
		cst_log(msg);
		}
	/* End single thread semaphore */
	_smx_end();
	return dtl_mass_validity(crit);
	}


/* DTL_get_support_upper:
 * Get the interval from 'belief_level' percent of the mass up to 100% mass (upper hull) */

rcode DTLAPI DTL_get_support_upper(int crit, double belief_level, double *lobo, double *upbo) {
	rcode rc;

	/* Begin single thread semaphore */
	_smx_begin("SMASU");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_support_upper(%d,%.3lf)\n",crit,belief_level);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(lobo,1);
	_certify_ptr(upbo,2);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (belief_level < 0.5)
		return dtl_error(DTL_INPUT_ERROR);
	/* Fetch the support mass of the evaluation */
	if (rc = dtl_cdf_to_ev(crit,2.0*belief_level-1.0,lobo,upbo))
		return dtl_error(rc);
	crit = max(crit,0); // partial/subcrit trees in slot 0
	*upbo = e_cache[crit][E_MAX][0];
	/* Log function result */
	if (cst_ext) {
		sprintf(msg," [%6.3lf %6.3lf]\n",*lobo,*upbo);
		cst_log(msg);
		}
	/* End single thread semaphore */
	_smx_end();
	return dtl_mass_validity(crit);
	}


/* DTI_get_support_mid:
 * Get the location of the cdf 50% midpoint */

rcode DTLAPI DTI_get_support_mid(int crit, double *cdf) {
	double ref_lo,ref_up,step,val,new_val,target;

	/* Begin single thread semaphore */
	_smx_begin("SMASM");
	/* Log function call */
	if (cst_ext) {
		sprintf(msg,"DTI_get_support_mid(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(cdf,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df00(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if ((crit < 1) && (dtl_latest_mc_eval-crit))
		return dtl_error(DTL_OUTPUT_ERROR);
	crit = max(crit,0); // subcrit trees in slot 0
	if (!ec[crit].valid)
		return dtl_error(DTL_OUTPUT_ERROR);
	/* Calculate support mass */
	if ((e_cache[crit][E_MAX][0]-e_cache[crit][E_MIN][0] > DTL_EPS) && (ecache_cm2[crit] > INF_MASS_VAR)) {
		ref_lo = bn_cdf(e_cache[crit][E_MIN][0],ec[crit].location,ec[crit].scale2,ec[crit].alpha);
		ref_up = bn_cdf(e_cache[crit][E_MAX][0],ec[crit].location,ec[crit].scale2,ec[crit].alpha);
		/* Lower bound */
		step = (e_cache[crit][E_MAX][0]-e_cache[crit][E_MIN][0])/2.0;
		new_val = (e_cache[crit][E_MAX][0]+e_cache[crit][E_MIN][0])/2.0;
		/* Convert from DTL cdf scale to b-normal cdf scale */
		target = dtl2bn(ref_lo,ref_up,0.5);
		do {
			val = new_val;
			*cdf = bn_cdf(val,ec[crit].location,ec[crit].scale2,ec[crit].alpha);
			step /= 2.0;
			if (*cdf > target)
				new_val = val-step;
			else
				new_val = val+step;
			} while ((fabs(*cdf-target) > 0.000001) && (step > 1.0E-7));
		*cdf = val;
		}
	else
		*cdf = ecache_rm1[crit];
	/* Log function result */
	if (cst_ext) {
		sprintf(msg," cdf = %6.3lf\n",*cdf);
		cst_log(msg);
		}
	/* End single thread semaphore */
	_smx_end();
	return dtl_mass_validity(crit);
	}


/* Obtain EV range by way of the DTLbnormal.c function for the CDF inverse */

static rcode dti_cdf_to_ev(int crit, double belief_level, double *lobo, double *upbo) {
	double ref_lo,ref_up,lo_target,up_target;

	/* Check input parameters */
	if (load_df00(crit))
		return DTL_CRIT_UNKNOWN;
	if ((crit < 1) && (dtl_latest_mc_eval-crit))
		return DTL_OUTPUT_ERROR;
	crit = max(crit,0); // subcrit trees in slot 0
	if (!ec[crit].valid)
		return DTL_OUTPUT_ERROR;
	if ((belief_level < MIN_SUPPORT_LEVEL) || (belief_level > MAX_SUPPORT_LEVEL))
		return DTL_INPUT_ERROR; // out-of-range, round-off errors will start interfering
	/* Calculate support mass */
	if ((e_cache[crit][E_MAX][0]-e_cache[crit][E_MIN][0] > DTL_EPS) && (ecache_cm2[crit] > INF_MASS_VAR)) {
		ref_lo = bn_cdf(e_cache[crit][E_MIN][0],ec[crit].location,ec[crit].scale2,ec[crit].alpha);
		ref_up = bn_cdf(e_cache[crit][E_MAX][0],ec[crit].location,ec[crit].scale2,ec[crit].alpha);
		/* Lower bound  - convert from DTL cdf scale to b-normal cdf scale */
		lo_target = dtl2bn(ref_lo,ref_up,(1.0-belief_level)/2.0);
		*lobo = bn_inv_cdf(lo_target,ec[crit].location,ec[crit].scale2,ec[crit].alpha);
		/* Upper bound  - convert from DTL cdf scale to b-normal cdf scale */
		up_target = dtl2bn(ref_lo,ref_up,(1.0+belief_level)/2.0);
		*upbo = bn_inv_cdf(up_target,ec[crit].location,ec[crit].scale2,ec[crit].alpha);
		}
	else { // Dirac point
		*lobo = ecache_rm1[crit];
		*upbo = ecache_rm1[crit];
		}
	return DTL_OK;
	}


rcode DTLAPI DTI_get_support_mass(int crit, double belief_level, double *lobo, double *upbo) {
	rcode rc;

	/* Begin single thread semaphore */
	_smx_begin("IMASS");
	/* Log function call */
	if (cst_ext) {
		sprintf(msg,"DTI_get_support_mass(%d,%.3lf)\n",crit,belief_level);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(lobo,1);
	_certify_ptr(upbo,2);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (belief_level < 0.5)
		return dtl_error(DTL_INPUT_ERROR);
	/* Fetch the support mass of the evaluation */
	if (rc = dti_cdf_to_ev(crit,belief_level,lobo,upbo))
		return dtl_error(rc);
	/* Log function result */
	if (cst_ext) {
		sprintf(msg," [%6.3lf %6.3lf]\n",*lobo,*upbo);
		cst_log(msg);
		}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


/* DTL_get_aversion_value
 * Get the risk aversion EV value from an expressed risk aversion/utility curve
 * (where =0 is risk neutral, >0 is risk averse, and <0 is risk prone).
 * Aversion parameter range: [-10,10] */

rcode DTLAPI DTL_get_aversion_value(int crit, double risk_aversion, double *ra_value) {
	rcode rc;
	double ra_level,ra_lo,ra_up,ip;

	/* Begin single thread semaphore */
	_smx_begin("AVERS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_aversion_value(%d,%.3lf)\n",crit,risk_aversion);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(ra_value,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	/* Risk_aversion in range |([0.00,0.01]+)[0.01,1.00]+[1.00,3/lg(2)]| */
	ra_level = 1.0-pow(2.0,-fabs(risk_aversion));
	if (ra_level > MIN_SUPPORT_LEVEL) {
		/* Negative risk_aversion means risk prone user */
		rc = dtl_cdf_to_ev(crit,ra_level,&ra_lo,&ra_up);
		*ra_value = risk_aversion>0.0?ra_lo:ra_up;
		}
	else {
		/* Low risk_aversion -> interpolate toward 50% cdf (not midpoint) */
		rc = dtl_cdf_to_ev(crit,MIN_SUPPORT_LEVEL,&ra_lo,&ra_up);
		ip = (MIN_SUPPORT_LEVEL-ra_level)/(2.0*MIN_SUPPORT_LEVEL);
		if (risk_aversion > 0.0)
			*ra_value = (1.0-ip)*ra_lo + ip*ra_up;
		else
			*ra_value = ip*ra_lo + (1.0-ip)*ra_up;
		}
	if (rc)
		return dtl_error(rc);
	/* Log function result */
	if (cst_ext) {
		sprintf(msg," aversion value = %6.3lf\n",*ra_value);
		cst_log(msg);
		}
	/* End single thread semaphore */
	_smx_end();
	return dtl_mass_validity(crit);
	}


 /**********************************************************
  *
  *  Mass debug information (note: partial tree rides on 0)
  *
  **********************************************************/

	/* Obtain DTL evaluation cache moments as calculated in this module.
	 * Undocumented - only for debugging.
	 * NOTE: different prerequisites if moved to DMC, should be checked. */

rcode DTLAPI DTI_get_mass_moments(int crit, double *rm1, double *cm2, double *cm3) {

	/* Begin single thread semaphore */
	_smx_begin("GMMOM");
	/* Log function call */
	if (cst_ext) {
		sprintf(msg,"DTI_get_mass_moments(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(rm1,1);
	_certify_ptr(cm2,2);
	_certify_ptr(cm3,3);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df00(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if ((crit < 1) && (dtl_latest_mc_eval-crit))
		return dtl_error(DTL_OUTPUT_ERROR);
	crit = max(crit,0); // subcrit trees in slot 0
	if (!ec[crit].valid)
		return dtl_error(DTL_OUTPUT_ERROR);
	/* Get result from cache */
	*rm1 = ecache_rm1[crit];
	*cm2 = ecache_cm2[crit];
	*cm3 = ecache_cm3[crit];
	/* Log function result */
	if (cst_ext) {
		sprintf(msg," rm1 = %.6le  cm2 = %.6le  cm3 = %.6le\n",*rm1,*cm2,*cm3);
		cst_log(msg);
		}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


static a_row rm1dum,cm2dum,cm3dum;

	/* Obtain the alternatives' TCL moments as calculated by the moment function.
	 * Undocumented - only for debugging (inefficient, only one alt at a time). */

rcode DTLAPI DTI_get_psi_moments(int crit, int alt, double *rm1, double *cm2, double *cm3) {

	/* Begin single thread semaphore */
	_smx_begin("GPMOM");
	/* Log function call */
	if (cst_ext) {
		sprintf(msg,"DTI_get_psi_moments(%d,%d)\n",crit,alt);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(rm1,1);
	_certify_ptr(cm2,2);
	_certify_ptr(cm3,3);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit)) // only PV-tree, not weights
		return dtl_error(DTL_CRIT_UNKNOWN);
	if ((alt < 1) || (alt > uf->df->n_alts))
		return dtl_error(DTL_ALT_UNKNOWN);
	/* Fetch moments for all alts from TCL */
	if (call(TCL_get_moments(uf->df,rm1dum,cm2dum,cm3dum),"TCL_get_moments"))
		return dtl_kernel_error();
	/* Transfer moments for only one alt */
	*rm1 = rm1dum[alt];
	*cm2 = cm2dum[alt];
	*cm3 = cm3dum[alt];
	/* Log function result */
	if (cst_ext) {
		sprintf(msg," rm1 = %.6le  cm2 = %.6le  cm3 = %.6le\n",*rm1,*cm2,*cm3);
		cst_log(msg);
		}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


rcode DTLAPI DTI_get_bn_params(int crit, double *loc, double *scale, double *alpha) {

	/* Undocumented - only for debugging */
	/* Begin single thread semaphore */
	_smx_begin("BNPAR");
	/* Log function call */
	if (cst_ext) {
		sprintf(msg,"DTI_get_bn_params(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(loc,1);
	_certify_ptr(scale,2);
	_certify_ptr(alpha,3);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df00(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if ((crit < 1) && (dtl_latest_mc_eval-crit))
		return dtl_error(DTL_OUTPUT_ERROR);
	crit = max(crit,0); // partial/subcrit trees in slot 0
	if (!ec[crit].valid)
		return dtl_error(DTL_OUTPUT_ERROR);
	/* Get result from cache */
	*loc = ec[crit].location;
	*scale = ec[crit].scale2;
	*alpha = ec[crit].alpha;
	/* Log function result */
	if (cst_ext) {
		sprintf(msg," loc = %.6le  scale = %.6le  alpha = %.6le\n",*loc,*scale,*alpha);
		cst_log(msg);
		}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*********************************************************
  *
  *  Compound evaluations
  *
  *********************************************************/

static e_matrix rank_result;
static ai_col gamma_order,omega_order;
static ar_col dummy_value;

rcode DTLAPI DTL_compare_alternatives(int crit, int method, double belief_level, ar_col lo_value, ar_col up_value) {
	rcode rc;
	int Ai,cst_global;
	struct d_frame *df;
	double lo_val,up_val,ip;

	/* Begin single thread semaphore */
	_smx_begin("COMP");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_compare_alternatives(%d,%d,%.3lf)\n",crit,method,belief_level);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(lo_value,1);
	_certify_ptr(up_value,2);
	_dtl_assert(lo_value!=up_value,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if ((belief_level < 0.0) || (belief_level > 1.0))
		return dtl_error(DTL_INPUT_ERROR);
	/* Handle partial MC eval */
	if (load_df00(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	df = uf->df;
	/* Evaluate frameset and collect support range */
	cst_global = cst_on;
	cst_on = FALSE;
	dtl_abort_init();
	for (Ai=1; Ai<=df->n_alts; Ai++) {
		dtl_abort_check();
		if (rc = evaluate_frameset(crit,method,Ai,0,rank_result)) {
			cst_on = cst_global;
			return dtl_error(rc);
			}
		if (belief_level < MIN_SUPPORT_LEVEL) { // interpolate cdf midpoint
			if (rc = dtl_cdf_to_ev(crit,MIN_SUPPORT_LEVEL,&lo_val,&up_val)) {
				cst_on = cst_global;
				return dtl_error(rc);
				}
			// linear interpolation (since pdf derivative is 0 at inflexion point)
			ip = (MIN_SUPPORT_LEVEL-belief_level)/(2.0*MIN_SUPPORT_LEVEL);
			// cannot use rank_result[E_MID][0] due to possible skew (this is
			// the "///" argument: EV mid and 50% cdf do not always coincide)
			lo_value[Ai] = (1.0-ip)*lo_val + ip*up_val;
			up_value[Ai] = (1.0-ip)*up_val + ip*lo_val;
			}
		else if (belief_level > MAX_SUPPORT_LEVEL) { // interpolate EV range
			if (rc = dtl_cdf_to_ev(crit,MAX_SUPPORT_LEVEL,&lo_val,&up_val)) {
				cst_on = cst_global;
				return dtl_error(rc);
				}
			// quadratic interpolation (since cdf is converging steeper than linear)
			ip = sq((belief_level-MAX_SUPPORT_LEVEL)/(1.0-MAX_SUPPORT_LEVEL));
			lo_value[Ai] = (1.0-ip)*lo_val + ip*rank_result[E_MIN][0];
			up_value[Ai] = (1.0-ip)*up_val + ip*rank_result[E_MAX][0];
			}
		else { // within cdf conversion limits
			if (rc = dtl_cdf_to_ev(crit,belief_level,lo_value+Ai,up_value+Ai)) {
				cst_on = cst_global;
				return dtl_error(rc);
				}
			}
		}
	cst_on = cst_global;
	/* Log function result */
	if (cst_ext)
		for (Ai=1; Ai<=df->n_alts; Ai++) {
			sprintf(msg," A%-2d [%.3lf %.3lf]\n",Ai,lo_value[Ai],up_value[Ai]);
			cst_log(msg);
			}
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


/* Matrix with cdf mass of deltas for each pair [Ai,Aj] of alternatives
 *
 * Interpolation modes:
 * mode=0: standard (raw) matrix = no interpolation or rebalancing
 * mode=1: in upper triangle: (i) no row may decrease going to the right &
 *         (ii) no column may increase going down (row prioritisation)
 * mode=2: in upper triangle: (i) no column may increase going down &
 *         (ii) no row may decrease going to the right (column prioritisation)
 * mode=3: in upper triangle: no row may decrease going to the right (row only)
 * mode=-1 in upper triangle: no column may increase going down (undocumented)
 *
 * (in lower triangle: the opposite way around because [Aj,Ai] mirrors [Ai,Aj]) */

rcode DTLAPI DTL_delta_mass(int crit, int mode, ar_matrix delta_value, ar_matrix delta_mass) {
	rcode rc;
	int Ai,Aj,cst_global;
	struct d_frame *df;
	double pos_gamma,neg_gamma,mass_lvl;

	/* Begin single thread semaphore */
	_smx_begin("DMASS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_delta_mass(%d,%d)\n",crit,mode);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(delta_value,1);
	_certify_ptr(delta_mass,2);
	_dtl_assert(delta_value!=delta_mass,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if ((mode < -1) || (mode > 3))
		return dtl_error(DTL_INPUT_ERROR);
	/* Handle partial MC eval */
	if (load_df00(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	df = uf->df;
	/* Evaluate each pair of alternatives */
	cst_global = cst_on;
	cst_on = FALSE;
	dtl_abort_init();
	for (Ai=1; Ai<=df->n_alts; Ai++) {
		delta_value[Ai][Ai] = delta_mass[Ai][Ai] = 0.0;
		for (Aj=Ai+1; Aj<=df->n_alts; Aj++) {
			dtl_abort_check();
			if (rc = evaluate_frameset(crit,E_DELTA,Ai,Aj,rank_result)) {
				cst_on = cst_global;
				return dtl_error(rc);
				}
			delta_value[Ai][Aj] =  rank_result[E_MID][0];
			delta_value[Aj][Ai] = -rank_result[E_MID][0];
			rc = dtl_ev_to_cdf(crit,1.0E-6,&pos_gamma);
			if (!rc)
				rc = dtl_ev_to_cdf(crit,-1.0E-6,&neg_gamma);
			/* Take care of collected rc */
			if (rc) {
				cst_on = cst_global;
				return dtl_error(rc);
				}
			delta_mass[Ai][Aj] = (pos_gamma+neg_gamma)/2.0;
			delta_mass[Aj][Ai] = 1.0-delta_mass[Ai][Aj];
			}
		}
	if (mode) { // interpolate belief
		/* Fetch omega EV */
		for (Ai=1; Ai<=df->n_alts; Ai++) {
			dtl_abort_check();
			if (rc = evaluate_frameset(crit,E_PSI,Ai,0,rank_result)) {
				cst_on = cst_global;
				return dtl_error(rc);
				}
			dummy_value[Ai] = rank_result[E_MID][0];
			omega_order[Ai] = Ai;
			}
		/* Obtain omega ordering */
		sort_b(omega_order,dummy_value,1,df->n_alts,TRUE);
		/* Phase 0: rebalance mass such that higher PSI has at least 50% */
		for (Ai=1; Ai<df->n_alts; Ai++)
			for (Aj=Ai+1; Aj<=df->n_alts; Aj++)
				if (delta_mass[omega_order[Ai]][omega_order[Aj]] < 0.5) {
				/* Two possible causes:
				 * 1. Weak mass distribution from unmodal value base
				 * 2. Omega EV does not adequately rank distributions */
					delta_mass[omega_order[Ai]][omega_order[Aj]] = 0.5;
					delta_mass[omega_order[Aj]][omega_order[Ai]] = 0.5;
					}
		/* Interpolate in phases depending on row/col order (can yield different results).
		 * The term 'upper triangle' refers to values above the diagonal in the matrix
		 * and the term 'lower triangle' refers to values below the diagonal (mirrored). */
		if (mode == 1)
			/* Phase 1: interpolate horizontally in upper triangle (row prio) */
			for (Ai=1; Ai<df->n_alts; Ai++) {
				mass_lvl = 0.0;
				for (Aj=Ai+1; Aj<=df->n_alts; Aj++)
					if (delta_mass[omega_order[Ai]][omega_order[Aj]] < mass_lvl) {
						delta_mass[omega_order[Ai]][omega_order[Aj]] = mass_lvl;
						delta_mass[omega_order[Aj]][omega_order[Ai]] = 1.0-mass_lvl;
						}
					else
						mass_lvl = delta_mass[omega_order[Ai]][omega_order[Aj]];
				}
		if (mode < 3)
			/* Phase 2: interpolate vertically in upper triangle */
			for (Aj=2; Aj<=df->n_alts; Aj++) {
				mass_lvl = 1.0;
				for (Ai=1; Ai<Aj; Ai++)
					if (delta_mass[omega_order[Ai]][omega_order[Aj]] > mass_lvl) {
						delta_mass[omega_order[Ai]][omega_order[Aj]] = mass_lvl;
						delta_mass[omega_order[Aj]][omega_order[Ai]] = 1.0-mass_lvl;
						}
					else
						mass_lvl = delta_mass[omega_order[Ai]][omega_order[Aj]];
				}
		if (mode > 1)
			/* Phase 3: interpolate horizontally in upper triangle (column prio) */
			for (Ai=1; Ai<df->n_alts; Ai++) {
				mass_lvl = 0.0;
				for (Aj=Ai+1; Aj<=df->n_alts; Aj++)
					if (delta_mass[omega_order[Ai]][omega_order[Aj]] < mass_lvl) {
						delta_mass[omega_order[Ai]][omega_order[Aj]] = mass_lvl;
						delta_mass[omega_order[Aj]][omega_order[Ai]] = 1.0-mass_lvl;
						}
					else
						mass_lvl = delta_mass[omega_order[Ai]][omega_order[Aj]];
				}
		}
	cst_on = cst_global;
	if (cst_ext) {
		for (Ai=1; Ai<=df->n_alts; Ai++) {
			sprintf(msg," A%-2d ",Ai);
			cst_log(msg);
			for (Aj=1; Aj<=df->n_alts; Aj++)
				if (Ai==Aj) // self is undefined
					cst_log("       ");
				else {
					sprintf(msg,"%6.3lf ",delta_value[Ai][Aj]);
					cst_log(msg);
					}
			cst_log("\n");
			}
		for (Ai=1; Ai<=df->n_alts; Ai++) {
			sprintf(msg," A%-2d ",Ai);
			cst_log(msg);
			for (Aj=1; Aj<=df->n_alts; Aj++)
				if (Ai==Aj) // self is undefined
					cst_log("     ");
				else {
					sprintf(msg,"%3.0lf%%%% ",100.0*delta_mass[Ai][Aj]);
					cst_log(msg);
					}
			cst_log("\n");
			}
		}
	eval_cache_invalidate();
	_smx_end();
	return DTL_OK;
	}


 /*******************************************************
  *
  *  Ordering/ranking terminology
  *  ----------------------------
  *
  *  x is the content of an alt-vector cell
  *  y is an index position in an alt-vector
  *
  *  Ordering: alt Ax is ranked as y:th best
  *  Ranking:  alt Ay is ranked as x:th best
  *
  *       altvector
  *      +---------+
  *   .  |   ...   |
  *   .  |---------|
  *  [y]:|    x    |
  *   .  |---------|
  *   .  |   ...   |
  *
  *  The order/rank sorting is invertible ->
  *  x == order[rank[x]] == rank[order[x]]
  *
  ***********************************************************/


/* General ranking function with different ranking modes
 * (negative mode numbers are non-traditional rankings)
 *
 * Mode: -3 = delta dominance, strict ranking (not SML)
 *       -2 = psi support level, strict ranking
 *       -1 = gamma cdf, strict ranking
 *        0 = gamma EV, olympic ranking
 *        1 = gamma EV, strict ranking
 *        2 = gamma EV, strict ranking with tiebreaker
 *        3 = gamma EV, group ranking
 *
 * NOTE: gamma cdf (-1) can give unintuitive results when EVs
 * are clustered and spaced far apart (e.g. many alternatives
 * become 0.0 or 1.0 despite being quite different and having
 * rather differing gamma values). Use gamma EV (+1) instead. */

rcode DTLAPI DTL_rank_alternatives(int crit, int mode, double gamma_tolerance, double omega_tolerance, 
		ai_col gamma_rank, ai_col omega_rank, ar_col gamma_value, ar_col omega_value) {
	rcode rc;
	int Ai,Aj,level,cst_global;
#if defined(DOM_2024) && !defined(C_SML)
	int i;
#endif
	double pos_gamma,neg_gamma;
	struct d_frame *df;
	a_row rm1,cm2,cm3;

	/* Begin single thread semaphore */
	_smx_begin("RANK");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_rank_alternatives(%d,%d,%.3lf,%.3lf)\n",
				crit,mode,gamma_tolerance,omega_tolerance);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(gamma_rank,1);
	_certify_ptr(omega_rank,2);
	_certify_ptr(gamma_value,3);
	_certify_ptr(omega_value,4);
	_dtl_assert(gamma_rank!=omega_rank,1);
	_dtl_assert(gamma_value!=omega_value,2);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Phase 0: check input parameters */
	if (abs(mode) > 3)
		return dtl_error(DTL_INPUT_ERROR);
	if (mode == -2) {
		// gamma_tolerance contains support target - not tolerance level
		if ((gamma_tolerance < 0.0) || (gamma_tolerance > 1.0))
			return dtl_error(DTL_INPUT_ERROR);
		if (fabs(gamma_tolerance-0.5) < MIN_SUPPORT_LEVEL) {
			// no cdf available -> convert to omega-only call
			gamma_tolerance = omega_tolerance;
			mode = 4; // internal omega-only ordering
			}
		}
	else
		if ((gamma_tolerance < 0.0) || (gamma_tolerance > 0.1))
			return dtl_error(DTL_WRONG_TOLERANCE);
	if ((omega_tolerance < 0.0) || (omega_tolerance > 0.1))
		return dtl_error(DTL_WRONG_TOLERANCE);
	/* Handle partial MC eval */
	if (load_df00(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	df = uf->df;
	/* Prepare tie-breaker for non-MC df */
	if (crit > 0)
		if (call(TCL_get_moments(df,rm1,cm2,cm3),"TCL_get_moments"))
			return dtl_kernel_error();
	/* Phase 1: generate evaluations */
	cst_global = cst_on;
	cst_on = FALSE;
	dtl_abort_init();
	for (Ai=1; Ai<=df->n_alts; Ai++) {
		dtl_abort_check();
		if (rc = evaluate_frameset(crit,E_PSI,Ai,0,rank_result)) {
			cst_on = cst_global;
			return dtl_error(rc);
			}
		omega_value[Ai] = rank_result[E_MID][0];
		dtl_abort_check();
		if (mode == 4) { // internal mode, copy omega
			gamma_value[Ai] = omega_value[Ai];
			}
		else if (mode == -2) { // rank belief support
			if (gamma_tolerance > MAX_SUPPORT_LEVEL)
				// no cdf available -> convert to upper hull
				gamma_value[Ai] = rank_result[E_MAX][0];
			else if (1.0-gamma_tolerance > MAX_SUPPORT_LEVEL)
				// no cdf available -> convert to lower hull
				gamma_value[Ai] = rank_result[E_MIN][0];
			else if (gamma_tolerance > 0.5) // upper half of cdf
				rc = dtl_cdf_to_ev(crit,2.0*gamma_tolerance-1.0,&pos_gamma,gamma_value+Ai);
			else // lower half of cdf
				rc = dtl_cdf_to_ev(crit,1.0-2.0*gamma_tolerance,gamma_value+Ai,&pos_gamma);
			}
		else { // rank gamma
			if (rc = evaluate_frameset(crit,E_GAMMA,Ai,0,rank_result)) {
				cst_on = cst_global;
				return dtl_error(rc);
				}
			if (mode == -1) { // cdf instead of EV
				rc = dtl_ev_to_cdf(crit,1.0E-6,&pos_gamma);
				if (!rc)
					rc = dtl_ev_to_cdf(crit,-1.0E-6,&neg_gamma);
				gamma_value[Ai] = (pos_gamma+neg_gamma)/2.0;
				}
			else { // EV - prepare for tiebreak
				gamma_value[Ai] = rank_result[E_MID][0];
				// tie-break threshold is set relative to sorting horizon in sort_b
				dummy_value[Ai] = rank_result[E_MID][0] + (crit>0?1.0E-2*cm2[Ai]:0.0);
				}
			}
		if (rc) { // catch conversion rc
			cst_on = cst_global;
			return dtl_error(rc);
			}
		/* Set up order sort vectors */
		gamma_order[Ai] = Ai;
		omega_order[Ai] = Ai;
		}
	/* Phase 2: sort in descending order */
	if (mode == 2) // tiebreak
		sort_b(gamma_order,dummy_value,1,df->n_alts,TRUE);
	else
		sort_b(gamma_order,gamma_value,1,df->n_alts,TRUE);
	sort_b(omega_order,omega_value,1,df->n_alts,TRUE);
	/* Convert order -> ranking */
	if (mode && (mode != 3))
		/* Hard (strict/enforced) ranking, even if same value */
		for (Ai=1; Ai<=df->n_alts; Ai++) {
			gamma_rank[gamma_order[Ai]] = Ai;
			omega_rank[omega_order[Ai]] = Ai;
			}
	else { // 0=olympic or 3=group
		/* Soft ranking = same rank if less than 'tolerance' diff, sorted from
		   the beginning of the ordering, i.e. from the "best" alternative
		   (which constitutes the most important part of the ranking).
		   The internal tolerance in the sorting algorithm is DTL_EPS. */
		Ai = level = 1;
		omega_rank[omega_order[Ai]] = Ai;
		for (Aj=2; Aj<=df->n_alts; Aj++)
			if (omega_value[omega_order[Ai]]-omega_value[omega_order[Aj]] < omega_tolerance+DTL_EPS)
				omega_rank[omega_order[Aj]] = mode?level:Ai;
			else {
				omega_rank[omega_order[Aj]] = mode?++level:Aj;
				Ai = Aj;
				}
		Ai = level = 1;
		gamma_rank[gamma_order[Ai]] = Ai;
		for (Aj=2; Aj<=df->n_alts; Aj++)
			if (gamma_value[gamma_order[Ai]]-gamma_value[gamma_order[Aj]] < gamma_tolerance+DTL_EPS)
				gamma_rank[gamma_order[Aj]] = mode?level:Ai;
			else {
				gamma_rank[gamma_order[Aj]] = mode?++level:Aj;
				Ai = Aj;
				}
		}
#if defined(DOM_2024) && !defined(C_SML)
	if (mode == -3) {
		/* Dominance evaluation -> replace gamma */
		for (i=1; i<df->n_alts; i++) {
			dtl_abort_check();
			Ai = omega_order[i];
			Aj = omega_order[i+1];
			if (rc = dtl_get_dominance(crit,Ai,Aj,gamma_value+Ai,gamma_rank+Ai)) {
				cst_on = cst_global;
				return dtl_error(rc);
				}
			}
		gamma_rank[omega_order[df->n_alts]] = -1; // unused
		gamma_value[omega_order[df->n_alts]] = -1.0; // unused
		}
#endif
	cst_on = cst_global;
	/* Log function result */
	if (cst_ext)
		for (Ai=1; Ai<=df->n_alts; Ai++) {
			sprintf(msg," A%-2d %2d %6.3lf  %2d %6.3lf\n",Ai,gamma_rank[Ai],
					gamma_value[Ai],omega_rank[Ai],omega_value[Ai]);
			cst_log(msg);
			}
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	/* Phase 3: check if same ranking */
#if defined(DOM_2024) && !defined(C_SML)
	if (mode > -3) // dominance does not have rank
#endif
		for (Ai=1; Ai<=df->n_alts; Ai++)
			if (gamma_rank[Ai] != omega_rank[Ai])
				return DTL_DIFFERING_RANKS;
	return DTL_OK;
	}


/* Daisy chain modes:
 *  0 = return absolute omega EV values (default)
 *  1 = return relative omega EV values
 * +2 = mix belief mass and omega EV w/i radius */

/* The outer mix cut point is called the radius. It can range between 0.0
 * and 0.5, where 0.0 means no mixing and 0.5 means mixing deltas being up
 * to half the value scale apart. Larger radii than that are unreasonable. */

rcode DTLAPI DTL_daisy_chain2(int crit, int mode, double radius, 
		ai_col omega_rank, ar_col daisy_value, ar_col omega_value) {
	rcode rc;
	int Ai,i,mixed,cst_global;
	struct d_frame *df;
	double pos_gamma,neg_gamma;

	/* Begin single thread semaphore */
	_smx_begin("DAISY");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_daisy_chain(%d,%d,%.3lf)\n",crit,mode,radius);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(omega_rank,1);
	_certify_ptr(daisy_value,2);
	_certify_ptr(omega_value,3);
	_dtl_assert(daisy_value!=omega_value,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df00(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if ((mode < 0) || (mode > 3))
		return dtl_error(DTL_INPUT_ERROR);
	mixed = (mode&0x02)>>1;
	mode &= 0x01;
	if (mixed) // radius validation
		if ((radius < 0.0) || (radius > 0.5))
			return dtl_error(DTL_INPUT_ERROR);
	df = uf->df;
	/* Evaluate each alternative */
	cst_global = cst_on;
	cst_on = FALSE;
	dtl_abort_init();
	for (Ai=1; Ai<=df->n_alts; Ai++) {
		dtl_abort_check();
		if (rc = evaluate_frameset(crit,E_PSI,Ai,0,rank_result)) {
			cst_on = cst_global;
			return dtl_error(rc);
			}
		omega_value[Ai] = rank_result[E_MID][0];
		/* Set up order sort vectors */
		omega_order[Ai] = Ai;
		}
	/* Sort omega in descending order as basis for daisy chain */
	sort_b(omega_order,omega_value,1,df->n_alts,TRUE);
	/* Evaluate daisy chain */
	for (i=1; i<df->n_alts; i++) {
		dtl_abort_check();
		if (rc = evaluate_frameset(crit,E_DELTA,omega_order[i],omega_order[i+1],rank_result)) {
			cst_on = cst_global;
			return dtl_error(rc);
			}
		rc = dtl_ev_to_cdf(crit,1.0E-6,&pos_gamma);
		if (!rc)
			rc = dtl_ev_to_cdf(crit,-1.0E-6,&neg_gamma);
		if (rc) {
			cst_on = cst_global;
			return dtl_error(rc);
			}
		daisy_value[omega_order[i]] = (pos_gamma+neg_gamma)/2.0;
		daisy_value[omega_order[i]] = max(daisy_value[omega_order[i]],0.5); // trim roundoff error
		}
	daisy_value[omega_order[df->n_alts]] = -1.0; // unused
	if (mixed && (radius > 0.0)) {
		/* Mix belief mass and EV value output */
		for (Ai=1; Ai<df->n_alts; Ai++)
			daisy_value[omega_order[Ai]] -= // reduce cdf mass proportionally within a radius
				max(1.0-(omega_value[omega_order[Ai]]-omega_value[omega_order[Ai+1]])/radius,0.0) *
					(daisy_value[omega_order[Ai]]-0.5);
		}
	if (mode) {
		/* Create relative EV value output */
		for (Ai=1; Ai<df->n_alts; Ai++)
			omega_value[omega_order[Ai]] -= omega_value[omega_order[Ai+1]];
		omega_value[omega_order[df->n_alts]] = -1.0; // unused
		}
	/* Convert order to hard ranking (soft/olympic not useful here).
	 * NOTE: daisy values are belief strengths, cannot rank on them.
	 * Rather, the daisy chain's ranking must be based on omega EV. */
	for (Ai=1; Ai<=df->n_alts; Ai++)
		omega_rank[omega_order[Ai]] = Ai;
	cst_on = cst_global;
	if (cst_ext)
		for (Ai=1; Ai<=df->n_alts; Ai++) {
			sprintf(msg," A%-2d %2d %6.3lf %6.3lf\n",Ai,omega_rank[Ai],daisy_value[Ai],omega_value[Ai]);
			cst_log(msg);
			}
	eval_cache_invalidate();
	_smx_end();
	return DTL_OK;
	}


#define DAISY_RADIUS 0.1 // default radius

/* The original (classic) calls - for simplicity and compatibility
 * DTL_daisy_chain:  no parameters, no radius
 * DTL_daisy_chain1: one mode parameter (see DTL_daisy_chain2) */

rcode DTLAPI DTL_daisy_chain(int crit, ai_col omega_rank, ar_col daisy_value, ar_col omega_value) {

	return DTL_daisy_chain2(crit,0,0.0,omega_rank,daisy_value,omega_value);
	}


rcode DTLAPI DTL_daisy_chain1(int crit, int mode, ai_col omega_rank, ar_col daisy_value, ar_col omega_value) {

	return DTL_daisy_chain2(crit,mode,DAISY_RADIUS,omega_rank,daisy_value,omega_value);
	}


/* Proportional ranking for pie charts (unnormalised ranking is
 * available from DTL_daisy_chain and DTL_rank_alternatives).
 *
 * Based on heuristics: this is experienced as a good pie chart
 * that conveys the belief mass distribution among alternatives.
 * Mode=1 is a faithful pie chart, mode=0 is "cosy"/old compat
 * For frames with only two alternatives, the two modes yield
 * the same result (*).
 *
 * Parameters for mode=1 (ineffectual for mode=0)
 * ----------------------------------------------
 * Moderation1 controls how much of its mass the best alternative
 * distributes along the daisy chain. 0.0 means keep all (default).
 * Moderation2 controls how much of their mass the other alternatives
 * distribute along the daisy chain. 0.0 means keep all (default).
 * Through moderation, the soft effect is obtained more reasonably.
 *
 * mode+2 = mix belief mass and EV w/i radius (see daisy chain call)
 *
 * (*) = For moderation1 > 0, the result for two alternatives in mode
 *       1 is not compatible with mode 0, for moderation2 > 0 it is. */

static ai_col card_order,card_rank;

// DTL layer 0: above DTL proper

rcode DTLAPI DTL_pie_chart2(int crit, int mode, double moderation1, double moderation2, ar_col pie_value) {
	rcode rc;
	int i,Ai;
	struct d_frame *df;
	double sum,pos,cur_val;

	/* Check input parameters */
	if ((mode < 0) || (mode > 3))
		return DTL_INPUT_ERROR;
	if (mode&0x01) {
		// moderations of 0.0 yield the true algorithm, else softer
		if ((moderation1 < 0.0) || (moderation1 > 1.0))
			return DTL_INPUT_ERROR;
		moderation1 += 1.0; // divisor -> 1.0 is neutral
		if ((moderation2 < 0.0) || (moderation2 > 1.0))
			return DTL_INPUT_ERROR;
		moderation2 /= 2.0; // offset -> 0.0 is neutral
		}
	/* Evaluate alternatives and find a ranking by piggybacking */
	if (rc = DTL_daisy_chain1(crit,mode&0x02,card_rank,pie_value,dummy_value))
		return rc;
	df = uf->df;
	/* Convert hard ranking back to order */
	for (Ai=1; Ai<=df->n_alts; Ai++)
		card_order[card_rank[Ai]] = Ai;
	/* Traverse belief chain in order */
	if (mode&0x01) {
		// from largest belief downwards = the most reasonable algorithm
		pos = min(1.0-pie_value[card_order[1]]/moderation1,0.5);
		sum = pie_value[card_order[1]];
		for (i=2; i<=df->n_alts; i++) { // calculate all except the first
			cur_val = pie_value[card_order[i]];
			pie_value[card_order[i]] = pos;
			pos *= moderation2 + 2.0*(1.0-moderation2)*(1.0-cur_val);
			sum += pie_value[card_order[i]];
			}
		}
	else {
		// from smallest belief upwards = less reasonable but "desired" due to "cosiness"
		pie_value[card_order[df->n_alts]] = 1.0-pie_value[card_order[df->n_alts-1]]; // recreate last
		sum = 1.0; // the last two entries sum to one and are already done
		for (i=df->n_alts-2; i; i--) { // calculate all except the last two
			pie_value[card_order[i]] = pie_value[card_order[i+1]] + 2.0*(pie_value[card_order[i]]-0.5);
			sum += pie_value[card_order[i]];
			}
		}
	/* Normalise for pie chart */
	for (Ai=1; Ai<=df->n_alts; Ai++)
		pie_value[Ai] /= sum;
	if (cst_ext) {
		sprintf(msg,"DTL_pie_chart(%.3lf,%.3lf)\n",moderation1-1.0,2.0*moderation2);
		cst_log(msg);
		for (Ai=1; Ai<=df->n_alts; Ai++) {
			sprintf(msg," A%-2d %6.3lf\n",Ai,pie_value[Ai]);
			cst_log(msg);
			}
		}
	return DTL_OK;
	}


/* The original (classic) calls - for compatibility and simplicity.
 * DTL_pie_chart:  no parameters + compat with changed DMC (for ZOR batch)
 * DTL_pie_chart1: one combined moderation parameter
 *  Positive moderation modifies the daisy chain itself as a basis for the chart.
 *  Negative moderation modifies the starting point (anchor) of the pie chart. */

rcode DTLAPI DTL_pie_chart(int crit, ar_col pie_value) {

	return DTL_pie_chart2(crit,1,0.0,0.0,pie_value);
	}


rcode DTLAPI DTL_pie_chart1(int crit, double moderation, ar_col pie_value) {

	return DTL_pie_chart2(crit,1,max(-moderation,0.0),max(moderation,0.0),pie_value);
	}


 /*********************************************************
  *
  *  Security levels
  *
  *********************************************************/

/* Security level cache */
static a_vector strong,marked,weak;

 /*
  * Call semantics: All alternatives are evaluated using sec_level .
  *
  * The result has the form of a matrix {alt} x {min,mid,max} x {contraction},
  * with values from increasing contraction. Currently there are 6-21 values
  * corresponding to contractions of 0-100% in 5-20% steps, same as evaluation.
  *
  * The sec_contract and sec_mass functions deliver one single measure on how large
  * parts of the consistent bases are affected by violation of the thresholds.
  *
  */

static int dtl_sec_level(int crit, double v_min, s_matrix s_result) {
	int Ai,j;
	struct d_frame *df;

	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (dtl_error_count)
		return dtl_error(DTL_OUTPUT_ERROR);
	if ((v_min < 0.0) || (v_min > 1.0))
		return dtl_error(DTL_INPUT_ERROR);
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	df = uf->df;
	/* Get TCL results */
	if (call(TCL_security_level(df,v_min,strong,marked,weak),"TCL_security_level"))
		return dtl_kernel_error();
	/* Transfer result to caller */
	for (Ai=1; Ai<=df->n_alts; Ai++) {
			s_result[Ai][E_MIN] = strong[Ai];
			s_result[Ai][E_MID] = marked[Ai];
			s_result[Ai][E_MAX] = weak[Ai];
			}
	if (cst_ext) {
		for (Ai=1; Ai<=df->n_alts; Ai++) {
			sprintf(msg," A%d:",Ai);
			cst_log(msg);
			for (j=E_MAX; j>=E_MIN; j--) {
				sprintf(msg," %6.3lf",s_result[Ai][j]);
				cst_log(msg);
				}
			cst_log("\n");
			}
		}
	return DTL_OK;
	}


rcode DTLAPI DTL_sec_level(int crit, double v_min, s_matrix s_result) {
	int rc;

	_smx_begin("SEL");
	_certify_ptr(s_result,1);
	if (cst_on) {
		sprintf(msg,"DTL_sec_level(%d,%.3lf)\n",crit,v_min);
		cst_log(msg);
		}
	if (rc = dtl_sec_level(crit,v_min,s_result)) {
		_smx_end();
		return dtl_error(rc);
		}
	_smx_end();
	return DTL_OK;
	}


 /*************************************************************
  *
  *  Three add-on packages are defined at the end of DTLeval.c
  *  in order to have access to internal eval data structures.
  *
  *************************************************************/

#ifdef C_SML

 /*************************************************************
  *
  *  SML = Simplified Multi-criteria/stakeholder Layer
  *
  *************************************************************/

#include "SMLlayer.c"
#endif

#ifdef DOMINANCE

 /*************************************************************
  *
  *  Dominance-based evaluations
  *
  *************************************************************/

#include "DTLdominance.c" // develop and finish in 2024
#endif
