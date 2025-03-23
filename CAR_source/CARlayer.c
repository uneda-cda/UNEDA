/*
 *
 *                 _/_/_/_/       _/        _/_/_/_/
 *               _/             _/ _/      _/     _/
 *              _/            _/    _/    _/      _/
 *             _/           _/      _/   _/      _/
 *            _/           _/_/_/_/_/   _/_/_/_/
 *           _/           _/      _/   _/    _/
 *           _/          _/      _/   _/      _/
 *            _/_/_/    _/      _/   _/       _/
 *
 *
 *        The Cardinal Alternative Ranking Method (CAR)
 *        ---------------------------------------------
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
 *   File: CARlayer.c
 *
 *   Purpose: CAR interface to DTL API
 *
 *   Methods described in M. Danielson: The CAR Method
 *   and M. Danielson: The DURENO Techniques
 *
 */


 /***************************************************************
  *
  *  NOTE: Never mix CAR-generated statements with non-CAR data.
  *  This can lead to inconsistencies from data having different
  *  origins. Remember that CAR is only a syntactic sugar layer.
  *
  *  NOTE2: If CAR data has been generated with CAR_light=OFF,
  *  P/W-tornados with odd mode numbers will be used (no mid),
  *  since a full set of midpoints cannot move around and would
  *  only result in empty tornados.
  *
  ***************************************************************/

#include "DTL.h"
#include "DTLinternal.h"
#include "CAR.h"


 /********************************************************
  *
  *  Layer configuration parameters
  *
  ********************************************************/

/* Number of '>' allowed in statements (increased in 1.23.1).
 * Suggested max are PW=3 V=9, allowed max are now PW=9 V=19.
 *
 * This obfuscates the intended meaning of the terminology
 * since no decision maker can reasonably tell 8 '>' from 9
 * in W or P, or 18 '>' from 19 in V. But V is sorta like a
 * thermometer so many more steps are ok in V than W or P. */

#define noEXTENDED_CAR_STEPS

/* Allowing impure criteria trees is rather muddy waters. It
 * is inconceivable that a decision-maker is able to compare
 * worst-to-best for compound weights. Those comparisons will
 * instead invariably and unfortunately be absolute weights.
 *
 * The immediate and obvious alleviation to this dire dilemma
 * is to allow the UI to draw impure weights in a tree format
 * but keeping it a one-level structure internally. If it is
 * still included, it is mainly for the tidiness of the UI/UX.
 * That is a legitimate and not too bad argument for impurity. */

#define noALLOW_CAR_W_IMPURE_TREES

#ifdef ALLOW_CAR_W_IMPURE_TREES
#pragma message("Note! Impure weight trees allowed in CAR")
#endif

/* Allowing a degenerated V-base with all values equal or not.
 * The case with all values equal yields a meaningles scale
 * with an infinite trade-off -> the weight must be zero.
 * It is also computationally bad since all mids are 0.5
 * and the intervals [0,1] which is an unreasonable slack. */

#define V_DEGEN_SCALE

/* CAR is a layer on top of DTL at API level 0. This means
 * that it normally does not log its calls but rather lets
 * DTL log functions being called. In the same manner, CAR
 * cannot use the error handling within DTL - which is why
 * CAR logging is done after successful parameter checks. */

#define LOG_CAR // logging CAR invocations on top of DTL


 /********************************************************
  *
  *  Layer constants
  *
  ********************************************************/

#define CAR_EPS 1.0E-6

#ifdef EXTENDED_CAR_STEPS // extended cardinal step ranges
#define MAX_STEPS_PW 9
#define MAX_STEPS_V 19
#else  // standard cardinal step ranges
#define MAX_STEPS_PW 3
#define MAX_STEPS_V  9
#endif


 /********************************************************
  *
  *  Internal layer data
  *
  ********************************************************/

static d_row crc;
static int crc_method=0;
static int compat_w_mode=0;
static int compat_v_mode=0;
static int car_light=0;
#define COMPAT_W 0.10 // percent around CAR point
#define COMPAT_V 0.05 // percent of entire scale
static double compat_w=COMPAT_W;
static double compat_v=COMPAT_V;
static int car_activated=FALSE;
static int phull_open=FALSE;


 /**********************************************************
  *
  *  NOTE: The CAR layer runs in user mode, not system mode.
  *  Thus, the layer uses DTL but has no access to its data.
  *
  *  DTL layer 0: all functions are above DTL proper
  *
  **********************************************************/

rcode DTLAPI CAR_init(int method, int mode) {

	/* method=0 -> default weight & prob method, strongly recommended
	 * method>0 -> user selected weight & prob method, research only
	 * mode=0   -> normal modern CAR
	 * mode+1   -> backwards compatible weights with older IIASA Excel models
	 * mode+2   -> backwards compatible values with older IIASA Excel models
	 * mode+4   -> car_light = no midpoints set (no need to decouple
	 *             the midpoint which yields simpler partial hull calls) */

	/* Check if function can start */
	if (car_activated)
		return CAR_STATE_ERROR;
	if (frame_loaded)
		return CAR_STATE_ERROR;
	/* Check input parameters */
	if ((method < 0) || (method > 7))
		return CAR_INPUT_ERROR;
	crc_method = method;
	if ((mode < 0) || (mode > 7))
		return CAR_INPUT_ERROR;
#ifdef LOG_CAR
	/* Log error-free call */
	if (cst_ext) {
		sprintf(msg,"CAR_init(%d,%d)\n",method,mode);
		cst_log(msg);
		}
#endif
	/* Set run parameters */
	compat_w_mode = mode&0x01;
	compat_v_mode = (mode&0x02)>>1;
	car_light = (mode&0x04)>>2;
	car_activated = TRUE;
	return CAR_OK;
	}


rcode DTLAPI CAR_exit() {

	/* Check if function can start */
	if (!car_activated)
		return CAR_STATE_ERROR;
	if (frame_loaded)
		return CAR_STATE_ERROR;
#ifdef LOG_CAR
	/* Log error-free call */
	if (cst_ext)
		cst_log("CAR_exit()\n");
#endif
	/* Reset all run parameters */
	crc_method = 0;
	compat_w_mode = 0;
	compat_v_mode = 0;
	car_light = 0;
	compat_w = COMPAT_W;
	compat_v = COMPAT_V;
	car_activated = FALSE;
	phull_open = FALSE;
	return CAR_OK;
	}


rcode DTLAPI CAR_set_compat(double w_unc, double v_unc) {

	/* Check if function can start */
	if (!car_activated)
		return CAR_STATE_ERROR;
	if (frame_loaded)
		return CAR_STATE_ERROR;
	/* Check input parameters */
	if ((w_unc < 0.02) || (w_unc > 0.20))
		return CAR_INPUT_ERROR;
	if ((v_unc < 0.01) || (v_unc > 0.10))
		return CAR_INPUT_ERROR;
#ifdef LOG_CAR
	/* Log error-free call */
	if (cst_ext) {
		sprintf(msg,"CAR_set_compat(%.3lf,%.3lf)\n",w_unc,v_unc);
		cst_log(msg);
		}
#endif
	/* Set tolerances for Excel compat mode */
	compat_w = w_unc; // percent around CAR point
	compat_v = v_unc; // percent of entire scale
	return CAR_OK;
	}


 /*********************************************************
  *
  *  CAR weight and probability generation procedures
  *
  *********************************************************/

static void gen_roc(int slots, int offset) {
	int i,steps;
	double sum;

	steps = slots+2*offset;
	sum = 0.0;
	for (i=steps; i>offset; i--) {
		sum += 1.0/(double)i;
		crc[i-offset] = sum/(double)steps;
		}
	if (offset) {
		/* Must normalise CRC cutout */
		sum = 0.0;
		for (i=1; i<=slots; i++)
			sum += crc[i];
		for (i=1; i<=slots; i++)
			crc[i] /= sum;
		for (i=slots+1; i<=steps-offset; i++)
			crc[i] = 0.0;
		}
	}


static void gen_rs(int slots, int offset) {
	int i,steps;
	double sum;

	steps = slots+2*offset;
	for (i=slots+offset; i>offset; i--) {
		crc[i-offset] = (double)2.0*(steps+1.0-i)/(double)(steps*(steps+1.0));
		}
	if (offset) {
		/* Must normalise CRC cutout */
		sum = 0.0;
		for (i=1; i<=slots; i++)
			sum += crc[i];
		for (i=1; i<=slots; i++)
			crc[i] /= sum;
		for (i=slots+1; i<=steps-offset; i++)
			crc[i] = 0.0;
		}
	}


static void gen_rx(int slots, int offset, double z) {
	int i,steps;
	double sum;

	steps = slots+2*offset;
	for (i=slots+offset; i>offset; i--) {
		crc[i-offset] = pow(steps+1.0-i,z);
		}
	/* Must normalise */
	sum = 0.0;
	for (i=1; i<=slots; i++)
		sum += crc[i];
	for (i=1; i<=slots; i++)
		crc[i] /= sum;
	for (i=slots+1; i<=steps-offset; i++)
		crc[i] = 0.0;
	}


static void gen_rr(int slots, int offset) {
	int i,steps;
	double sum;

	steps = slots+2*offset;
	sum = 0.0;
	for (i=steps; i>offset; i--) {
		sum += 1.0/(double)i;
		}
	for (i=steps; i>offset; i--) {
		crc[i-offset] = 1.0/sum/(double)i;
		}
	if (offset) {
		/* Must normalise CRC cutout */
		sum = 0.0;
		for (i=1; i<=slots; i++)
			sum += crc[i];
		for (i=1; i<=slots; i++)
			crc[i] /= sum;
		for (i=slots+1; i<=steps-offset; i++)
			crc[i] = 0.0;
		}
	}


static void gen_sr(int slots, int offset) {
	int i,steps;
	double sum;

	/* SR = sum and reciprocal weights */
	steps = slots+2*offset;
	sum = 0.0;
	for (i=slots+offset; i>offset; i--) {
		sum += 1.0/(double)i+(double)(steps+1.0-i)/steps;
		}
	for (i=slots+offset; i>offset; i--) {
		crc[i-offset] = (1.0/(double)i+(double)(steps+1.0-i)/steps)/sum;
		}
	if (offset) {
		/* Must normalise cutout */
		sum = 0.0;
		for (i=1; i<=slots; i++)
			sum += crc[i];
		for (i=1; i<=slots; i++)
			crc[i] /= sum;
		for (i=slots+1; i<=steps-offset; i++)
			crc[i] = 0.0;
		}
	}


static void gen_xr(int slots, int offset, double z) {
	int i,steps;
	double sum;

	/* SRX = RX and reciprocal weights */
	steps = slots+2*offset;
	sum = 0.0;
	for (i=slots+offset; i>offset; i--) {
		sum += 1.0/(double)i+pow(steps+1.0-i,z)/(double)steps;
		}
	for (i=slots+offset; i>offset; i--) {
		crc[i-offset] = (1.0/(double)i+pow(steps+1.0-i,z)/(double)steps)/sum;
		}
	if (offset) {
		/* Must normalise cutout */
		sum = 0.0;
		for (i=1; i<=slots; i++)
			sum += crc[i];
		for (i=1; i<=slots; i++)
			crc[i] /= sum;
		for (i=slots+1; i<=steps-offset; i++)
			crc[i] = 0.0;
		}
	}


 /*********************************************************
  *
  *  Ordinal weight generation
  *
  *********************************************************/


rcode DTLAPI CAR_get_W_ordinal(int n_nodes, cr_col ord_wts) {
	int i;

	if (!car_activated)
		return CAR_NOT_ACTIVATED;
	if (!frame_loaded)
		return CAR_FRAME_NOT_LOADED;
	/* Generate weights from implicit ordinal ranking */
	gen_rx(n_nodes,0,1.0+min((double)n_nodes/60.0,0.25));
	for (i=1; i<=n_nodes; i++)
		ord_wts[i] = crc[i];
	return CAR_OK;
	}
		


 /*********************************************************
  *
  *  CAR load functions: weights, probabilities, values
  *
  *********************************************************/

static h_matrix elobox, eupbox;
static h_vector lobox, upbox;
static h_vector eloboxw, eupboxw;
static h_vector lobosw, midsw, upbosw;
static struct user_stmt_rec ustmt;
static struct user_w_stmt_rec uwstmt;


 /*********************************************************
  *
  *  Weight base
  *
  *********************************************************/

static int W_mark;

static int mark_W_base() {

	W_mark = DTL_nbr_of_W_stmts();
	return W_mark < DTL_OK;
	}


static void rollback_W_base() {
	int i;

	for (i=DTL_nbr_of_W_stmts(); i>W_mark; i--)
		if (DTL_error2(DTL_delete_W_statement(i)))
			return;
	}


rcode DTLAPI CAR_set_W_base(int n_nodes, car_vector ord_crit, car_vector rel) {
	rcode rc;
	int i,k,n_wts,n_act_nodes,tot,inx;
	double rsum;

	/* Check if function can start */
	_certify_ptr(ord_crit,101);
	_certify_ptr(rel,102);
	if (!car_activated)
		return CAR_NOT_ACTIVATED;
	if (!frame_loaded) // protecting PS test
		return CAR_FRAME_NOT_LOADED;
	if (phull_open)
		return CAR_NOT_ALLOWED;
	if (PS)
		return CAR_WRONG_FRAME_TYPE;
	if (mark_W_base())
		return CAR_SYS_CORRUPT;
	if (car_light && dtl_nbr_W_midpoints())
		// created previously with !car_light
		return CAR_NOT_ALLOWED;
	/* Check input parameters */
	if ((n_nodes < 1) || (n_nodes > MAX_NODES))
		return CAR_INPUT_ERROR;
	n_wts = DTL_nbr_of_weights();
	if (n_wts < DTL_OK)
		return n_wts;
	if (n_nodes > n_wts)
		return CAR_INPUT_ERROR;
	for (k=1; k<=n_nodes; k++)
		if ((ord_crit[k] < 1) || (ord_crit[k] > n_wts))
			return CAR_INPUT_ERROR;
#ifndef ALLOW_CAR_W_IMPURE_TREES
	if (!dtl_pure_W_tree())
		return CAR_ILLEGAL_TREE;
#endif
	for (k=2; k<=n_nodes; k++)
		if (dtl_W_node_parents(ord_crit[1],ord_crit[k]))
			return CAR_INPUT_ERROR;
	if (dtl_W_nbr_of_siblings(ord_crit[1]) != n_nodes)
		return CAR_INPUT_ERROR;
	if (n_nodes == 1)
		return CAR_OK;
	/* Initialise box vector */
	for (i=1; i<=n_wts; i++) {
		eloboxw[i] = -2.0;
		eupboxw[i] = -2.0;
		}
	/* Convert weights to [0,1] scale */
	tot = 1;
	for (k=1; k<n_nodes; k++) {
		if (rel[k] == -1) // nullified
			break;
		if ((rel[k] < 0) || (rel[k] > MAX_STEPS_PW))
			return CAR_INPUT_ERROR;
		tot += rel[k];
		}
	n_act_nodes = k; // active nodes
#ifdef LOG_CAR
	/* CAR log starts after input checks */
	if (cst_ext) {
		cst_log("CAR_set_W_base(");
		for (k=1; k<n_nodes; k++) {
			sprintf(msg,"W%d",ord_crit[k]);
			cst_log(msg);
			if (rel[k] > 0) // diff
				for (i=0; i<min(rel[k],MAX_STEPS_PW); i++)
					// min to guard from crap beyond nullifier
					cst_log(">");
			else if (rel[k]) // cut
				cst_log("|");
			else // equal
				cst_log("=");
			}
		sprintf(msg,"W%d) -->\n",ord_crit[n_nodes]);
		cst_log(msg);
		}
#endif
	switch (crc_method) {
		case 0:  gen_rx(tot,0,1.0+min((double)n_act_nodes/60.0,0.25)); break; // adaptive CAR
		case 1:  gen_rs(tot,0);         break;
		case 2:  gen_rr(tot,0);         break;
		case 3:  gen_xr(tot,0,1.35);    break;
		case 4:  gen_sr(tot,0);         break;
		case 5:  gen_roc(tot,0);        break;
		default: return CAR_INPUT_ERROR;
		}
	inx = 1;
	rsum = 0.0;
	for (k=1; k<=n_act_nodes; k++) {
		// active criteria
		eloboxw[ord_crit[k]] = crc[inx];
		eupboxw[ord_crit[k]] = crc[inx];
		lobox[ord_crit[k]] = inx<tot?(crc[inx]+crc[inx+1])/2.0:crc[inx]/2.0;
		upbox[ord_crit[k]] = inx>1?(crc[inx-1]+crc[inx])/2.0:tot>1?0.5+crc[inx]/2.0:1.5;
		rsum += crc[inx];
		inx += rel[k];
		}
	for (; k<=n_nodes; k++) {
		// nullified criteria
		eloboxw[ord_crit[k]] = 0.0;
		eupboxw[ord_crit[k]] = 0.0;
		lobox[ord_crit[k]] = 0.0;
		upbox[ord_crit[k]] = 0.0;
		}
	/* Normalise CRC selection and enter */
	uwstmt.n_terms = 1;
	uwstmt.sign[1] = 1;
	for (k=1; k<=n_nodes; k++) {
		eloboxw[ord_crit[k]] /= rsum;
		eupboxw[ord_crit[k]] /= rsum;
		uwstmt.crit[1] = ord_crit[k];
		if (compat_w_mode) { // Excel compatibility mode
			uwstmt.lobo = (1.0-compat_w)*eloboxw[ord_crit[k]];
			uwstmt.upbo = (1.0+compat_w)*eloboxw[ord_crit[k]];
			}
		else { // default modern mode
			uwstmt.lobo = lobox[ord_crit[k]]/rsum;
			uwstmt.upbo = min(upbox[ord_crit[k]]/rsum,1.0); // because of interpolation upwards
			}
		rc = DTL_add_W_statement(&uwstmt);
		if (rc < DTL_OK) {
			rollback_W_base();
			return rc;
			}
		}
	if (!car_light && (n_act_nodes > 1)) {
		/* For TCL stability */
		for (k=1; k<=n_act_nodes; k++) {
			eloboxw[ord_crit[k]] = max(eloboxw[ord_crit[k]]-CAR_EPS,0.0);
			eupboxw[ord_crit[k]] = min(eloboxw[ord_crit[k]]+CAR_EPS,1.0);
			}
		if (rc = dtl_set_W_mbox_auto(eloboxw,eupboxw)) {
			rollback_W_base();
			return rc;
			}
		}
	/* Return number of statements */
	ord_crit[0] = DTL_nbr_of_W_stmts()-W_mark;
#ifdef LOG_CAR
	/* Log error-free completion */
	if (cst_ext) {
		sprintf(msg,"--> end of CAR_set_W_base(%d)\n",ord_crit[0]);
		cst_log(msg);
		}
#endif
	return DTL_nbr_of_W_stmts()-W_mark;
	}


 /**********************************************************
  *
  *  Partial hull weight verification using DURENO-II
  *  ------------------------------------------------
  *  The partial hull weight verification method works by
  *  pruning the orthogonal hull around the two weights
  *  in the verification statement. It could be seen as a
  *  hull-trimming operator in the same way that standard
  *  normalisation works. And in analogy, since ordinary
  *  normalisation does not guarantee that every pair of
  *  weight combinations taken from the hull are feasible,
  *  in the same way there is no guarantee that weight
  *  combinations from the pruned/cut hull are together
  *  compliant with the verification statement. Here, the
  *  analogy breaks since normalisation can detect which
  *  combinations of weights are feasible while there is
  *  no similar operator to detect which hull combinations
  *  have been verified or not. Such an operator would have
  *  to be non-linear by design since the verifications for
  *  intervals are ratios (w1/w2 >= v2/v1) which are of an
  *  inherently non-linear kind. Thus, while CAR-layer can
  *  verify and prune or cut the hull, it is up to the API
  *  caller to ensure that the weights ultimately are kept
  *  within the specifications supplied by the user.
  *
  *  NOTE: Since this is a user layer on top the caller must
  *  ensure that the phull is closed before unloading frame.
  *
  **********************************************************/

static void pos_first_w_stmt(struct user_w_stmt_rec* swp) {
	int itmp;
	double rtmp;

	/* Ensure positive term first */
	if (swp->sign[1] < 0) {
		itmp = swp->crit[1];
		swp->crit[1] = swp->crit[2];
		swp->crit[2] = itmp;
		itmp = swp->sign[1];
		swp->sign[1] = swp->sign[2];
		swp->sign[2] = itmp;
		rtmp = swp->lobo;
		swp->lobo = swp->upbo;
		swp->upbo = rtmp;
		}
	}


static void swap_w_stmt(struct user_w_stmt_rec* swp) {
	int itmp;

	/* Swap 1st and 2nd crit */
	itmp = swp->sign[1];
	swp->sign[1] = swp->sign[2];
	swp->sign[2] = itmp;
	}


rcode DTLAPI CAR_open_W_phull() {
	rcode rc;
	int i,j,k;

	/* Check if function can start */
	if (!car_activated)
		return CAR_NOT_ACTIVATED;
	if (!frame_loaded) // protecting PS test
		return CAR_FRAME_NOT_LOADED;
	if (PS)
		return CAR_WRONG_FRAME_TYPE;
	if (phull_open)
		return CAR_NOT_ALLOWED;
	if (uf->WP_autogen[0]) {
		for (k=1,i=1; i<=uf->df->n_alts; i++)
			for (j=1; j<=uf->df->tot_cons[i]; j++,k++) {
				eloboxw[k] = -1.0;
				eupboxw[k] = -1.0;
				}
		/* Clear mhull */
		if (rc = dtl_set_W_mbox_auto(eloboxw,eupboxw))
			return rc;
		}
	phull_open = TRUE;
	return CAR_OK;
	}


rcode DTLAPI CAR_close_W_phull() {
	rcode rc;
	int i,j,k;

	/* Check if function can start */
	if (!car_activated)
		return CAR_NOT_ACTIVATED;
	if (!frame_loaded) // protecting PS test
		return CAR_FRAME_NOT_LOADED;
	if (PS)
		return CAR_WRONG_FRAME_TYPE;
	if (!phull_open)
		return CAR_NOT_ALLOWED;
	if (uf->WP_autogen[0]) {
		if (rc = DTL_get_W_hull(0,lobosw,midsw,upbosw))
			return rc;;
		for (k=1,i=1; i<=uf->df->n_alts; i++)
			for (j=1; j<=uf->df->tot_cons[i]; j++,k++) {
				eloboxw[k] = max(midsw[k]-CAR_EPS,0.0);
				eupboxw[k] = min(midsw[k]+CAR_EPS,1.0);
				}
		/* Set new mhull */
		if (rc = dtl_set_W_mbox_auto(eloboxw,eupboxw))
			return rc;
		}
	phull_open = FALSE;
	return CAR_OK;
	}


/* Check the partial hull syntax (swp content) and trade-off */

static rcode car_check_W_phull(struct user_w_stmt_rec* swp, double *tradeoff) {
	rcode rc;
	int n_wts;

	/* Check if function can start */
	_certify_ptr(swp,101);
	_certify_ptr(tradeoff,102);
	if (!car_activated)
		return CAR_NOT_ACTIVATED;
	if (!phull_open)
		return CAR_NOT_ALLOWED;
	if (!frame_loaded) // protecting PS test
		return CAR_FRAME_NOT_LOADED;
	if (PS)
		return CAR_WRONG_FRAME_TYPE;
#ifdef PHULL_LIGHT
	if (!car_light)
		return CAR_NOT_ALLOWED;
	if (dtl_nbr_W_midpoints())
		// created previously with !car_light
		return CAR_NOT_ALLOWED;
#endif
	n_wts = DTL_nbr_of_weights();
	if (n_wts < DTL_OK)
		return n_wts;
#ifndef ALLOW_CAR_W_IMPURE_TREES
	if (!dtl_pure_W_tree())
		return CAR_ILLEGAL_TREE;
#endif
	/* Validate swp */
	if (swp->n_terms != 2)
		return CAR_INPUT_ERROR;
	pos_first_w_stmt(swp);
	if ((swp->crit[1] < 1) || (swp->crit[1] > n_wts))
		return CAR_INPUT_ERROR;
	if (swp->sign[1] != 1)
		return CAR_INPUT_ERROR;
	if ((swp->crit[2] < 1) || (swp->crit[2] > n_wts))
		return CAR_INPUT_ERROR;
	if (swp->sign[2] != -1)
		return CAR_INPUT_ERROR;
	if (swp->crit[1] == swp->crit[2])
		return CAR_INPUT_ERROR;
	if ((swp->lobo < CAR_EPS) || (swp->lobo > 1.0))
		return CAR_INPUT_ERROR;
	if ((swp->upbo < CAR_EPS) || (swp->upbo > 1.0))
		return CAR_INPUT_ERROR;
	if (dtl_W_node_parents(swp->crit[1],swp->crit[2]))
		return CAR_INPUT_ERROR;
	if (!dtl_real_W_crit(swp->crit[1]))
		return CAR_INPUT_ERROR;
	if (!dtl_real_W_crit(swp->crit[2]))
		return CAR_INPUT_ERROR;
#ifdef LOG_CAR
	/* CAR log starts after input checks */
	if (cst_ext) {
		sprintf(msg,"CAR_check_W_phull(%d,%d,%.3lf,%.3lf,%.3lf) -->\n",
				swp->crit[1],swp->crit[2],swp->lobo,swp->upbo,*tradeoff);
		cst_log(msg);
		}
#endif
	/* Get current local hull */
	if (rc = DTL_get_W_hull(0,lobosw,midsw,upbosw))
		return rc;

	/* Return the best possible (max) trade-off ratio on the mutual [0,1] MC scale.
	 * If the ratio is <= 1.0 then there is no trade-off possible in this direction. */

	if (*tradeoff == -1.0) { // also return max trade-off components on the MC scale in swp
		swp->lobo *= upbosw[swp->crit[1]];
		swp->upbo *= lobosw[swp->crit[2]];
		if (swp->upbo > CAR_EPS)
			*tradeoff = swp->lobo/swp->upbo;
		}
	else // leave swp alone (possibly for subsequent use in prune/cut/equal calls)
		if (lobosw[swp->crit[2]] > CAR_EPS)
			*tradeoff = swp->lobo*upbosw[swp->crit[1]]/swp->upbo/lobosw[swp->crit[2]];
		else // not generated by CAR (or nullified)
			*tradeoff = -2.0;
	return CAR_OK;
	}


rcode DTLAPI CAR_check_W_phull(struct user_w_stmt_rec* swp, double *tradeoff) {
	rcode rc;

	if (rc = car_check_W_phull(swp,tradeoff))
		return rc;
#ifdef LOG_CAR
	/* Log error-free completion */
	if (cst_ext) {
		if ((*tradeoff > 0.0) && (upbosw[swp->crit[2]] > 0.0))
			sprintf(msg,"--> end of CAR_check_W_phull(%.3lf,%.3lf)\n",
					(swp->lobo*lobosw[swp->crit[1]]/swp->upbo/upbosw[swp->crit[2]]),*tradeoff);
		else
			sprintf(msg,"--> end of CAR_check_W_phull(INF,%.3lf)\n",*tradeoff);
		cst_log(msg);
		}
#endif
	return CAR_OK;
	}


/* Partial hull verification */

rcode DTLAPI CAR_prune_W_phull(struct user_w_stmt_rec* swp) {
	rcode rc;
	double troff=0.0;

	/* Check input parameters */
	if (rc = car_check_W_phull(swp,&troff))
		return rc;
	if (troff < 0.0) // not generated by CAR
		return CAR_NOT_ALLOWED;

	/* Effectuate W-base changes */
	if (mark_W_base())
		return CAR_SYS_CORRUPT;
	uwstmt.n_terms = 1;
	uwstmt.sign[1] = +1;
	if (swp->lobo*lobosw[swp->crit[1]] - swp->upbo*lobosw[swp->crit[2]] < -CAR_EPS) {
		/* Cr.2's lower bound eats into cr.1's hull on the lower (left) side,
		 * forcing an increase of the lobo of cr.1 to stay consistent. */
		uwstmt.crit[1] = swp->crit[1];
		uwstmt.lobo = swp->upbo*lobosw[swp->crit[2]]/swp->lobo;
		uwstmt.upbo = upbosw[swp->crit[1]];
		if ((uwstmt.upbo-uwstmt.lobo < 2.0*CAR_EPS) || (uwstmt.upbo > 1.0)) {
			// need to catch, otherwise input error
			rollback_W_base();
#ifdef LOG_CAR
			if (cst_ext)
				cst_log("--> end of CAR_prune_W_phull(0,INCONSISTENT)\n");
#endif
			return CAR_INCONSISTENT;
			}
		rc = DTL_add_W_statement(&uwstmt);
		if (rc < DTL_OK) {
			rollback_W_base();
			return rc;
			}
		}
	if (swp->lobo*upbosw[swp->crit[1]] - swp->upbo*upbosw[swp->crit[2]] < -CAR_EPS) {
		/* Cr.1's upper bound eats into cr.2's hull on the upper (right) side,
		 * forcing a decrease of the upbo of cr.2 to stay consistent. */
		uwstmt.crit[1] = swp->crit[2];
		uwstmt.lobo = lobosw[swp->crit[2]];
		uwstmt.upbo = swp->lobo*upbosw[swp->crit[1]]/swp->upbo;
		if ((uwstmt.upbo-uwstmt.lobo < 2.0*CAR_EPS) || (uwstmt.upbo > 1.0)) {
			// need to catch, otherwise input error
			rollback_W_base();
#ifdef LOG_CAR
			if (cst_ext)
				cst_log("--> end of CAR_prune_W_phull(0,INCONSISTENT)\n");
#endif
			return CAR_INCONSISTENT;
			}
		rc = DTL_add_W_statement(&uwstmt);
		if (rc < DTL_OK) {
			rollback_W_base();
			return rc;
			}
		}
#ifdef LOG_CAR
	/* Log error-free completion */
	if (cst_ext) {
		sprintf(msg,"--> end of CAR_prune_W_phull(%d,%.3lf)\n",DTL_nbr_of_W_stmts()-W_mark,troff);
		cst_log(msg);
		}
#endif
	return DTL_nbr_of_W_stmts()-W_mark;
	}


rcode DTLAPI CAR_cut_W_phull(struct user_w_stmt_rec* swp) {
	rcode rc;
	double gap,troff=0.0;

	/* Check input parameters */
	rc = car_check_W_phull(swp,&troff);
	if (rc)
		return rc;
	if (troff < 0.0)
		return CAR_NOT_ALLOWED; // not generated by CAR

	/* Effectuate W-base hull cut changes (more radical) */
	if (mark_W_base())
		return CAR_SYS_CORRUPT;
	uwstmt.n_terms = 1;
	uwstmt.sign[1] = 1;
	gap = swp->upbo*upbosw[swp->crit[2]] - swp->lobo*lobosw[swp->crit[1]];
	if (gap > CAR_EPS) {
		/* Cr.2's upper bound overlaps cr.1's lower bound on the MC scale */
		uwstmt.crit[1] = swp->crit[1];
		uwstmt.lobo = lobosw[swp->crit[1]] + gap/(2.0*swp->lobo);
		uwstmt.upbo = upbosw[swp->crit[1]];
		if ((uwstmt.upbo-uwstmt.lobo) < 3.0*CAR_EPS) {
			// need to catch, otherwise input error
			rollback_W_base();
#ifdef LOG_CAR
			if (cst_ext)
				cst_log("--> end of CAR_cut_W_phull(0,INCONSISTENT)\n");
#endif
			return CAR_INCONSISTENT;
			}
		rc = DTL_add_W_statement(&uwstmt);
		if (rc < DTL_OK) {
			rollback_W_base();
			return rc;
			}
		uwstmt.crit[1] = swp->crit[2];
		uwstmt.lobo = lobosw[swp->crit[2]];
		uwstmt.upbo = upbosw[swp->crit[2]] - gap/(2.0*swp->upbo);
		if ((uwstmt.upbo-uwstmt.lobo) < 3.0*CAR_EPS) {
			// need to catch, otherwise input error
			rollback_W_base();
#ifdef LOG_CAR
			if (cst_ext)
				cst_log("--> end of CAR_cut_W_phull(0,INCONSISTENT)\n");
#endif
			return CAR_INCONSISTENT;
			}
		rc = DTL_add_W_statement(&uwstmt);
		if (rc < DTL_OK) {
			rollback_W_base();
			return rc;
			}
		}
#ifdef LOG_CAR
	/* Log error-free completion */
	if (cst_ext) {
		sprintf(msg,"--> end of CAR_cut_W_phull(%d,%.3lf)\n",DTL_nbr_of_W_stmts()-W_mark,troff);
		cst_log(msg);
		}
#endif
	return DTL_nbr_of_W_stmts()-W_mark;
	}


/* Since _prune is analogous to '>=' (and _cut is analogous to '>'), then as a consequence
 * _equal is implemented as _prune in both directions, i.e. '>=' and '<=' yielding '='. */

rcode DTLAPI CAR_equal_W_phull(struct user_w_stmt_rec* swp) {
	rcode rc,rc2;

	rc = CAR_prune_W_phull(swp);
	if (rc < DTL_OK)
		return rc;
	swap_w_stmt(swp); // GEQ '>=' -> LEQ '<='
	rc2 = CAR_prune_W_phull(swp);
	swap_w_stmt(swp); // reverse statement back ->
	pos_first_w_stmt(swp); // (semi-)non-destructive
	if (rc2 < DTL_OK)
		return rc2;
	return rc+rc2; // total nbr of stmts added
	}


 /*********************************************************
  *
  *  Probability base
  *
  *********************************************************/

static int P_mark;

static int mark_P_base(int crit) {

	P_mark = DTL_nbr_of_P_stmts(crit);
	return P_mark < DTL_OK;
	}


static void rollback_P_base(int crit) {
	int i;

	for (i=DTL_nbr_of_P_stmts(crit); i>P_mark; i--)
		if (DTL_error2(DTL_delete_P_statement(crit,i)))
			return;
	}


rcode DTLAPI CAR_set_P_base(int crit, int alt, int n_nodes, car_vector ord_nodes, car_vector rel) {
	rcode rc;
	int i,j,k,n_act_nodes,t_nodes,tot,inx;
	int n_alts;
	double rsum;

	/* Check if function can start */
	_certify_ptr(ord_nodes,101);
	_certify_ptr(rel,102);
	if (!car_activated)
		return CAR_NOT_ACTIVATED;
	if (!frame_loaded)
		return CAR_FRAME_NOT_LOADED;
	/* Check input parameters */
	if ((n_nodes < 1) || (n_nodes > MAX_NODES))
		return CAR_INPUT_ERROR;
	if (mark_P_base(crit))
		return CAR_CRIT_UNKNOWN;
	if (car_light && dtl_nbr_P_midpoints(crit))
		// created previously with !car_light
		return CAR_NOT_ALLOWED;
	/* Check input parameters */
	n_alts = DTL_nbr_of_alts();
	if (n_alts < DTL_OK)
		return n_alts;
	if ((alt < 1) || (alt > n_alts))
		return CAR_ALT_UNKNOWN;
	t_nodes = DTL_nbr_of_nodes(crit,alt);
	if (t_nodes < DTL_OK)
		return t_nodes;
	if (n_nodes > t_nodes)
		return CAR_INPUT_ERROR;
	for (k=1; k<=n_nodes; k++)
		if ((ord_nodes[k] < 1) || (ord_nodes[k] > t_nodes))
			return CAR_INPUT_ERROR;
	for (k=2; k<=n_nodes; k++)
		if (dtl_P_node_parents(crit,alt,ord_nodes[1],ord_nodes[k]))
			return CAR_INPUT_ERROR;
	if (dtl_P_nbr_of_siblings(crit,alt,ord_nodes[1]) != n_nodes)
		return CAR_INPUT_ERROR;
	if (n_nodes == 1)
		return CAR_OK;
	/* Initialise box matrix */
	for (i=1; i<=n_alts; i++) {
		k = DTL_nbr_of_nodes(crit,i);
		if (k < DTL_OK)
			return k;
		for (j=1; j<=k; j++) {
			elobox[i][j] = -2.0;
			eupbox[i][j] = -2.0;
			}
		}
	/* Convert probabilities to [0,1] scale */
	tot = 1;
	for (k=1; k<n_nodes; k++) {
		if (rel[k] == -1) // nullified
			break;
		if ((rel[k] < 0) || (rel[k] > MAX_STEPS_PW))
			return CAR_INPUT_ERROR;
		tot += rel[k];
		}
	n_act_nodes = k; // active nodes
#ifdef LOG_CAR
	/* CAR log starts after input checks */
	if (cst_ext) {
		sprintf(msg,"CAR_set_P_base(%d,",crit);
		cst_log(msg);
		for (k=1; k<n_nodes; k++) {
			sprintf(msg,"P%d.%d",alt,ord_nodes[k]);
			cst_log(msg);
			if (rel[k] > 0) // diff
				for (i=0; i<min(rel[k],MAX_STEPS_PW); i++)
					cst_log(">");
			else if (rel[k]) // cut
				cst_log("|");
			else // equal
				cst_log("=");
			}
		sprintf(msg,"P%d.%d) -->\n",alt,ord_nodes[n_nodes]);
		cst_log(msg);
		}
#endif
	switch (crc_method) {
		case 0:  gen_rx(tot,0,1.0+min((double)n_act_nodes/60.0,0.25)); break; // adaptive CAR
		case 1:  gen_rs(tot,0);         break;
		case 2:  gen_rr(tot,0);         break;
		case 3:  gen_xr(tot,0,1.35);    break;
		case 4:  gen_sr(tot,0);         break;
		case 5:  gen_roc(tot,0);        break;
		default: return CAR_INPUT_ERROR;
		}
	inx = 1;
	rsum = 0.0;
	for (k=1; k<=n_act_nodes; k++) {
		// active criteria
		elobox[alt][ord_nodes[k]] = crc[inx];
		eupbox[alt][ord_nodes[k]] = crc[inx];
		lobox[ord_nodes[k]] = inx<tot?(crc[inx]+crc[inx+1])/2.0:crc[inx]/2.0;
		upbox[ord_nodes[k]] = inx>1?(crc[inx-1]+crc[inx])/2.0:tot>1?0.5+crc[inx]/2.0:1.5;
		rsum += crc[inx];
		inx += rel[k];
		}
	for (; k<=n_nodes; k++) {
		// nullified criteria
		elobox[alt][ord_nodes[k]] = 0.0;
		eupbox[alt][ord_nodes[k]] = 0.0;
		lobox[ord_nodes[k]] = 0.0;
		upbox[ord_nodes[k]] = 0.0;
		}
	/* Normalise CRC selection and enter */
	ustmt.n_terms = 1;
	ustmt.sign[1] = 1;
	ustmt.alt[1]  = alt;
	for (k=1; k<=n_nodes; k++) {
		elobox[alt][ord_nodes[k]] /= rsum;
		eupbox[alt][ord_nodes[k]] /= rsum;
		ustmt.cons[1] = ord_nodes[k];
		ustmt.lobo = lobox[ord_nodes[k]]/rsum;
		ustmt.upbo = min(upbox[ord_nodes[k]]/rsum,1.0); // limit upwards interpolation
		rc = DTL_add_P_statement(crit,&ustmt);
		if (rc < DTL_OK) {
			rollback_P_base(crit);
			return rc;
			}
		}
	/* For TCL stability */
	if (!car_light && (n_act_nodes > 1)) {
		for (k=1; k<=n_act_nodes; k++) {
			eupbox[alt][ord_nodes[k]] = min(elobox[alt][ord_nodes[k]]+CAR_EPS,1.0);
			elobox[alt][ord_nodes[k]] = max(elobox[alt][ord_nodes[k]]-CAR_EPS,0.0);
			}
		if (rc = dtl_set_P_mbox_auto(crit,elobox,eupbox)) {
			rollback_P_base(crit);
			return rc;
			}
		}
	/* Return number of statements */
	ord_nodes[0] = DTL_nbr_of_P_stmts(crit)-P_mark;
#ifdef LOG_CAR
	/* Log error-free completion */
	if (cst_ext) {
		sprintf(msg,"--> end of CAR_set_P_base(%d,%d)\n",crit,ord_nodes[0]);
		cst_log(msg);
		}
#endif
	return DTL_nbr_of_P_stmts(crit)-P_mark;
	}


 /*********************************************************
  *
  *  Value base
  *
  *********************************************************/

static int V_mark;

static int mark_V_base(int crit) {

	V_mark = DTL_nbr_of_V_stmts(crit);
	return V_mark < DTL_OK;
	}


static void rollback_V_base(int crit) {
	int i;

	for (i=DTL_nbr_of_V_stmts(crit); i>V_mark; i--)
		if (DTL_error2(DTL_delete_V_statement(crit,i)))
			return;
	}


/* CAR_set_V_base takes only full rankings on a criterion scale */

rcode DTLAPI CAR_set_V_base(int crit, car_vector ord_alts, car_vector ord_nodes, car_vector rel) {
	rcode rc;
	int i,j,k,tot,tot_cons,tot_nodes,n_alts,n_nodes,sum;
	double step,compat_v2;

	/* Check if function can start */
	_certify_ptr(ord_alts,101);
	_certify_ptr(ord_nodes,102);
	_certify_ptr(rel,103);
	if (!car_activated)
		return CAR_NOT_ACTIVATED;
	if (!frame_loaded)
		return CAR_FRAME_NOT_LOADED;
	/* Check input parameters */
	if (car_light && dtl_nbr_V_midpoints(crit))
		// created previously with !car_light
		return CAR_NOT_ALLOWED;
	/* CAR_set_V_base takes only full value rankings on a criterion scale.
	 * Thus, only one value base instance can be current at the same time. */
	if (rc = DTL_reset_V_base(crit))
		return rc;
	if (mark_V_base(crit))
		return CAR_CRIT_UNKNOWN;
	/* Initalise box matrix */
	n_alts = DTL_nbr_of_alts();
	if (n_alts < DTL_OK)
		return n_alts;
	for (i=1; i<=n_alts; i++) {
		n_nodes = DTL_nbr_of_nodes(crit,i);
		if (n_nodes < DTL_OK)
			return n_nodes;
		for (j=1; j<=n_nodes; j++) {
			// only full ranking is ok -> init with invalid marker
			elobox[i][j] = -3.0;
			eupbox[i][j] = -3.0;
			}
		}
	tot_nodes = DTL_total_nodes(crit);
	if (tot_nodes < DTL_OK)
		return tot_nodes;
	tot_cons = DTL_total_cons(crit);
	if (tot_cons < DTL_OK)
		return tot_cons;
	for (k=1; k<=tot_cons; k++) {
		/* In theory, this allows for repeated alternatives but that
		 * will be caught at consistency checks (unless rel==0). CAR
		 * is a layer on top of DTL with less firm responsibilities. */
		if ((ord_alts[k] < 1) || (ord_alts[k] > n_alts))
			return CAR_ALT_UNKNOWN;
		if ((ord_nodes[k] < 1) || (ord_nodes[k] > tot_nodes))
			return CAR_INPUT_ERROR;
		}
	/* Map values to [0,1] scale */
	tot = 0;
	for (k=1; k<tot_cons; k++) {
		if ((rel[k] < 0) || (rel[k] > MAX_STEPS_V))
			return CAR_INPUT_ERROR;
		tot += rel[k];
		}
#ifdef LOG_CAR
	/* CAR log starts after input checks */
	if (cst_ext) {
		sprintf(msg,"CAR_set_V_base(%d,",crit);
		cst_log(msg);
		for (k=1; k<tot_cons; k++) {
			sprintf(msg,"V%d.%d",ord_alts[k],ord_nodes[k]);
			cst_log(msg);
			if (rel[k])
				for (i=0; i<rel[k]; i++)
					cst_log(">");
			else
				cst_log("=");
			}
		sprintf(msg,"V%d.%d) -->\n",ord_alts[tot_cons],ord_nodes[tot_cons]);
		cst_log(msg);
		}
#endif
	if (tot) {
		/* Ensure no overlap in Excel compatibility mode */
		compat_v2 = min(compat_v,0.5/(double)(tot+1));
		step = (1.0-2.0*compat_v2)/(double)tot;
		/* Enter statements */
		ustmt.n_terms = 1;
		ustmt.sign[1] = +1;
		sum = 0;
		for (k=1; k<=tot_cons; k++) {
			ustmt.alt[1]  = ord_alts[k];
			ustmt.cons[1] = ord_nodes[k];
			if (compat_v_mode) { // Excel compatibility mode, use old scale
				elobox[ord_alts[k]][ord_nodes[k]] = (tot-sum)*step+compat_v2;
				eupbox[ord_alts[k]][ord_nodes[k]] = elobox[ord_alts[k]][ord_nodes[k]];
				ustmt.lobo = (tot-sum)*step;
				ustmt.upbo = min((tot-sum)*step+2.0*compat_v2,1.0);
				}
			else { // use standard value scale
				elobox[ord_alts[k]][ord_nodes[k]] = 1.0-(double)(2.0*sum+1.0)/(double)(2.0*tot+2.0);
				eupbox[ord_alts[k]][ord_nodes[k]] = elobox[ord_alts[k]][ord_nodes[k]];
				ustmt.lobo = 1.0-(double)(2.0*sum+2.0)/(double)(2.0*tot+2.0);
				ustmt.upbo = 1.0-(double)(2.0*sum)/(double)(2.0*tot+2.0);
				}
			rc = DTL_add_V_statement(crit,&ustmt);
			if (rc < DTL_OK) {
				rollback_V_base(crit);
				return rc;
				}
			sum += rel[k];
			}
		}
	else {
#ifdef V_DEGEN_SCALE
		/* Degenerated scale, no dimension */
		for (k=1; k<=tot_cons; k++) {
			elobox[ord_alts[k]][ord_nodes[k]] = 0.5;
			eupbox[ord_alts[k]][ord_nodes[k]] = 0.5;
			}
#else
		/* Degenerated scale, not allowed */
		return CAR_INPUT_ERROR;
#endif
		}
	/* For TCL stability */
	if (!car_light) {
		for (k=1; k<=tot_cons; k++) {
			if (elobox[ord_alts[k]][ord_nodes[k]] > CAR_EPS)
				elobox[ord_alts[k]][ord_nodes[k]] -= CAR_EPS;
			if (eupbox[ord_alts[k]][ord_nodes[k]] < 1.0-CAR_EPS)
				eupbox[ord_alts[k]][ord_nodes[k]] += CAR_EPS;
			}
		/* Set entire mbox plus number of V-relation steps */
		if (rc = dtl_set_V_mbox_rels(crit,tot,elobox,eupbox)) {
			rollback_V_base(crit);
			return rc;
			}
		}
	/* Return number of statements */
	ord_nodes[0] = DTL_nbr_of_V_stmts(crit)-V_mark;
#ifdef LOG_CAR
	/* Log error-free completion */
	if (cst_ext) {
		sprintf(msg,"--> end of CAR_set_V_base(%d,%d)\n",crit,ord_nodes[0]);
		cst_log(msg);
		}
#endif
	if (tot) // there exist a meaningful ranking
		return DTL_nbr_of_V_stmts(crit)-V_mark;
	else // flat-ranked = meaningless criterion
		return CAR_SAME_RANKINGS;
	}


 /*********************************************************
  *
  *  CAR distance ranking
  *
  *********************************************************/

#include "CARrank.c"
