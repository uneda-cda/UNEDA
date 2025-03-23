/*
 *
 *                   _/_/_/      _/_/_/_/_/   _/
 *                  _/    _/        _/       _/
 *                 _/      _/      _/       _/
 *                _/      _/      _/       _/
 *               _/      _/      _/       _/
 *              _/      _/      _/       _/
 *             _/     _/       _/       _/
 *            _/_/_/_/        _/       _/_/_/_/_/
 *
 *
 *              The Decision Tree Library (DTL)
 *              -------------------------------
 *
 *  +---- o o o --------------------------------------------+
 *  |   o       o            Prof. Mats Danielson           |
 *  |  o  STHLM  o           DECIDE Research Group          |
 *  |  o         o  Dept. of Computer and Systems Sciences  |
 *  |  o   UNI   o           Stockholm University           |
 *  |   o       o    PO Box 7003, SE-164 07 Kista, SWEDEN   |
 *  +---- o o o --------------------------------------------+
 *
 *           Copyright (c) 2022-2024 Mats Danielson
 *                Email: mats.danielson@su.se
 *                    Phone: +46 816 1540
 *
 *
 *   This code is proprietary, NOT free or shared!
 *   It may under NO circumstances be used for any
 *   purpose without a written agreement with the
 *   author.
 *
 */

/*
 *   File: SMLlayer.c
 *
 *   Purpose: SML interface code
 *
 *
 *   Functions exported outside SML
 *   ------------------------------
 *   scall
 *   ucall
 *   SML_init/2
 *   SML_exit
 *   SML_new_DM_flat_frame
 *   SML_new_DM_tree_frame
 *   SML_new_PM_flat_frame
 *   SML_new_PM_tree_frame
 *   SML_new_SM_tree_frame
 *   SML_read_frame
 *   SML_load_frame
 *   SML_unload_frame
 *   SML_dispose_frame
 *   SML_frame_type
 *   SML_set_W_base/2
 *   SML_set_P_base/2
 *   SML_set_V_base/2
 *   SML_set_W_correlations (experimental)
 *   SML_get_V_hull
 *   SML_evaluate_frame
 *   SML_get_mass_range
 *   SML_get_mass_above
 *   SML_get_mass_below
 *   SML_get_support_mass
 *   SML_get_support_lower
 *   SML_get_support_upper
 *   SML_evaluate_cdf
 *   SML_compare_alternatives
 *   SML_evaluate_mid
 *   SML_evaluate_omega/1/2
 *   SML_delta_mass/2
 *   SML_rank_alternatives
 *   SML_rank_gamma
 *   SML_rank_omega
 *   SML_daisy_chain/2
 *   SML_pie_chart/1/2
 *   SML_get_W_tornado/2
 *   SML_get_P_tornado/2
 *   SML_get_MCP_tornado/2
 *   SML_get_V_tornado/2
 *   SML_get_MCV_tornado/2
 *   SML_get_cons_influence
 *   SML_get_errtxt/2
 *   SML_crit_nbr
 *   SML_is_stakeholder
 *
 *   Functions outside module, inside SML
 *   ------------------------------------
 *   sml_nc
 *
 *   Functions internal to module
 *   ----------------------------
 *   sml_error_check
 *   sml_error_code
 *   sml_new_flat_frame
 *   sml_new_tree_frame
 *   sml_get_frame_type
 *   sml_scale_type
 *   sml_get_mass_point
 *   sml_get_support_mass
 *   sml_evaluate_omega
 *   sml_get_PV_tornado
 *
 *
 *   Version history
 *
 *   Ver.   Date   Main reasons
 *   ----  ------  ------------
 *   1.22  220520  SML introduced
 *   1.23  220802  Pointer logging
 *   1.24  221010  Excel VBA support
 *   1.27  230606  Rcode interpreter
 *
 */


 /**********************************************************
  *
  *  SML - Simplified Multi-criteria/stakeholder Layer
  *  -------------------------------------------------
  *
  *  SML has four main purposes:
  *
  *  1. Simplify access to the complete DTL package & API
  *  2. Facilitate easy access to DTL input bases (W,P,V)
  *  3. Integrate autoscale functionality with in/output
  *  4. Maintain the stakeholder interface for the package
  *
  *  This layer facilitates easier implementation of DTL
  *  functionality in different applications and research
  *  projects. A caller accesses max half of the DTL calls.
  *  In fact, there are only 104 published interface calls.
  *  On top of that, the included calls are made simpler.
  *
  *  This is DTL layer 0: all functions are above DTL proper
  *
  **********************************************************/


 /**********************************************************
  *
  *  Configuration parameters
  *  ------------------------
  *
  *  RENORM_W_BASE  renormalise W-base after V-scale change
  *  LOG_SML        logging of SML calls (esp. pointers)
  *  PRE_CHECK      pre-check of inconsistent variables
  *  ERR_TEST       include caller random error testing
  *  PIE_COMPAT     pie chart with 2 alts coincide with eval
  *  EXCEL_VBA      include special support for Excel VBA
  *  SML_STRG       0=null-terminated, 1=length-preceded
  *
  **********************************************************/

#define noRENORM_W_BASE
#define noLOG_SML
#define PRE_CHECK
#define noERR_TEST
#define PIE_COMPAT
#define EXCEL_VBA
#define SML_STRG 0


 /**********************************************************
  *
  *  Layer compiler directives and constants
  *
  **********************************************************/

#define SMC 0 // multi-criteria indicator
#define SML_EPS 1.0E-7
#define SML_FAIL_FREQ 20 // fraction of simulated errors


 /**********************************************************
  *
  *  NOTE: The SML layer runs in user mode, not system mode.
  *  Thus, the layer uses DTL but has no access to its data.
  *  Instead, the layer keeps its own permanent static data.
  *
  **********************************************************/

static int sml_active=FALSE;
static int sml_loaded=0;
static int sml_vsource=0; // v_source (also DMC compatible)
static int sml_ecrit=0;   // evaluation criterion
static int sml_scale=0;   // evaluation scale type
// emethod also works as a guard for belief mass calls
static int sml_emethod;   // evaluation method
#define NO_EMETHOD -1
#define EMETHOD_MASK 0x0C // no risk portfolios in SML

/* Frame types and properties */
#define SM1_FRAME 1 // SM with mirrored criteria
#define SM2_FRAME 2 // SM with duplicated criteria
#define NSM_FRAME 3 // not SM (= DM and PM frames)

#define SM1 (sml_type[0]==SM1_FRAME)
#define SM2 (sml_type[0]==SM2_FRAME)
#define NSM (sml_type[0]==NSM_FRAME)

static int sml_type[MAX_FRAMES+1];
static int sml_n_sh[MAX_FRAMES+1];
static int sml_n_cr[MAX_FRAMES+1];
static int sml_renorm;
#ifdef ERR_TEST
static int sml_err_test=FALSE;
#endif


 /***********************************************************
  *
  *  Generalised error code interpreter and handler
  *
  *  Fetches the appropriate error code and combines it
  *  with the return code into one single error number.
  *
  *  Error number < -2: DTL return code = real error
  *                 -2: no result but not real error
  *                 -1: inconsistent user input (I_USR)
  *                  0: ok
  *                  1: ok + additional state information
  *                > 1: ok + additional numeric information
  *
  *  The parameter 'mode' signals whether an inconsistency
  *  could be caused by the system or a user (numeric input
  *  or sometimes an incomplete user-specified node tree).
  *
  ***********************************************************/

#define I_SYS 0
#define I_USR 1

static rcode sml_error_check(rcode rc, int mode) {
	int err;

	err = DTL_error2(rc) + (mode?DTL_u_error2(rc):DTL_error2(rc));
	switch (err) {
		case 0:
			if (rc<0) return 1; // ok + state info
			else return rc;     // ok + numeric info
		case 2:  return -2;   // no result but no error
		case 3:  return -1;   // inconsistent user input
		default: return rc;   // system error
		}
	}


#ifdef ERR_TEST // exercise callers error handlers

// run call exit through error simulator
#define SML_EC(rc) sml_error_code(rc)

static rcode sml_error_code(rcode rc) {
	time_t now;

	/* Stress test of caller's error handling */
	if (sml_err_test) // pseudo-random failure
		if (!(time(&now) % SML_FAIL_FREQ))
			return SML_OUTPUT_ERROR; // simulate error
	return rc;
	}

#else

// execute call exit directly, bypassing error simulator
#define SML_EC(rc) (rc)

#endif


 /***********************************************************
  *
  *  Inline wrapper error-code interpreters/handlers
  *  Deviate from calling conventions to be less clumsy
  *  The interpreters are run as post-processing calls
  *
  *  scall: standard error code handler
  *  ucall: user/num error code handler
  *
  *  Usage: err = s/ucall(SML_function(par1,par2,...));
  *  err contains the error number as listed above
  *
  ***********************************************************/

rcode DTLAPI scall(rcode func) {

	return sml_error_check(func,I_SYS);
	}


rcode DTLAPI ucall(rcode func) {

	return sml_error_check(func,I_USR);
	}


 /********************************************************
  *
  *  SML system commands
  *
  *  Init modes: 0 = V-base values are human generated
  *                  (default value base behaviour)
  *              1 = V-base values are machine generated
  *                  (also compatible with old DMC)
  *             +2 = activate error exit tests
  *
  ********************************************************/

rcode DTLAPI SML_init2(int mode) {
	rcode rc;
	int i;

	/* Check if function can start */
	if (sml_active)
		return SML_STATE_ERROR;
	/* Check input parameter */
#ifdef ERR_TEST
	if (mode & 0xFFFC) // unknown flag (error test mode)
#else
	if (mode & 0xFFFE) // unknown flag (normal mode)
#endif
		return SML_INPUT_ERROR;
	/* Attempt to open DTL (API level 1) */
	if (rc = DTL_init())
		return rc;
	/* Open CAR (API level 2 = same as SML) */
	if (rc = CAR_init(0,0))
		return rc;
	/* Post processing */
	sml_active = TRUE;
	sml_vsource = mode&0x01;
	for (i=0; i<=MAX_FRAMES; i++) {
		sml_type[i] = 0;
		sml_n_sh[i] = 0;
		sml_n_cr[i] = 0;
		}
#ifdef ERR_TEST
	sml_err_test = (mode&0x02)>>1;
	/* Log function exit only for logging of err_test */
	if (cst_ext && sml_err_test) { // V-base source is machine or human?
		sprintf(msg,"SML_init flags VS=%c ET=T\n",sml_vsource?'M':'H');
		cst_log(msg);
		}
#endif
	return SML_OK;
	}


rcode DTLAPI SML_init(int mode) {

	return SML_init2(mode&0x02); // no machine vsource mode
	}


rcode DTLAPI SML_exit() {
	rcode rc;

	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	if (sml_loaded)
		return SML_FRAME_IN_USE;
	/* Close CAR (API level 2 = same as SML) */
	if (rc = CAR_exit())
		return rc;
	/* Reset internal flags */
	sml_active = FALSE;
	sml_vsource = FALSE;
#ifdef ERR_TEST
	sml_err_test = FALSE;
#endif
	/* Close DTL (API level 1) */
	return DTL_exit(); // returns nbr of trace entries written (or error)
	}


int sml_is_open() {

	// only for internal SML layer use
	return sml_active;
	}


 /********************************************************
  *
  *  SML frame creation
  *  ------------------
  *
  *  DM-frame is the basic type (deterministic MC). Then
  *  there are two extensions: PM adds probabilistic event
  *  trees to the alternatives, while SM instead adds the
  *  concept of multiple stakeholders to a basic DM-frame.
  *
  *  NOTE: SML handles only MC frames. For single-crit
  *  frames, the original DTL interface should be used.
  *
  *  [Trick: create an MC frame with 1 crit to get PS,
  *   but that loophole has been plugged since otherwise
  *   there would be a PS-frame with a modifiable scale.]
  *
  ********************************************************/

static rcode sml_new_flat_frame(rcode (DTLAPI *func)(int,int,int), int ufnbr, int n_crit, int n_alts) {
	rcode rc;

#ifdef LOG_SML
	/* Log function call only for extended logging */
	if (cst_ext) {
		sprintf(msg,"SML_new_XM_flat_frame(%d,%d,%d)\n",ufnbr,n_crit,n_alts);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	if (sml_loaded)
		return SML_FRAME_IN_USE;
	/* Check input parameter */
	if (n_crit < 2) // plug PS loophole
		return SML_INPUT_ERROR;
	/* Create frame */
	if (rc = (*func)(ufnbr,n_crit,n_alts))
		return rc;
	/* Post processing */
	sml_type[ufnbr] = NSM_FRAME;
	sml_n_sh[ufnbr] = 1;
	sml_n_cr[ufnbr] = n_crit;
	return SML_OK;
	}


static rcode sml_new_tree_frame(rcode (DTLAPI *func)(int,int,int,ta_tree), 
		int ufnbr, int n_alts, int n_wtnodes, ta_tree wtree) {
	rcode rc;

#ifdef LOG_SML
	/* Log function call only for extended pointer logging */
	if (cst_ext) {
		sprintf(msg,"SML_new_XM_tree_frame(%d,%d,%d,%p)\n",ufnbr,n_alts,n_wtnodes,wtree);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	if (sml_loaded)
		return SML_FRAME_IN_USE;
	/* Check input parameter */
	if (n_wtnodes < 2) // PS loophole
		return SML_INPUT_ERROR;
	/* Create tree frame */
	if (rc = (*func)(ufnbr,n_alts,n_wtnodes,wtree))
		return rc;
	/* Post processing */
	sml_type[ufnbr] = NSM_FRAME;
	sml_n_sh[ufnbr] = 1;
	sml_n_cr[ufnbr] = uf_list[ufnbr]->n_crit;
	return SML_OK;
	}


/* Create DM-frame */

rcode DTLAPI SML_new_DM_flat_frame(int ufnbr, int n_crit, int n_alts) {

	/* Call SML pattern function */
	return sml_new_flat_frame(&DTL_new_DM_flat_frame,ufnbr,n_crit,n_alts);
	}


rcode DTLAPI SML_new_DM_tree_frame(int ufnbr, int n_alts, int n_wtnodes, ta_tree wtree) {

	/* Call SML pattern function */
	return sml_new_tree_frame(&DTL_new_DM_tree_frame,ufnbr,n_alts,n_wtnodes,wtree);
	}


/* Create PM-frame */

rcode DTLAPI SML_new_PM_flat_frame(int ufnbr, int n_crit, int n_alts) {

	/* Call SML pattern function */
	return sml_new_flat_frame(&DTL_new_PM_flat_frame,ufnbr,n_crit,n_alts);
	}


rcode DTLAPI SML_new_PM_tree_frame(int ufnbr, int n_alts, int n_wtnodes, ta_tree wtree) {

	/* Call SML pattern function */
	return sml_new_tree_frame(&DTL_new_PM_tree_frame,ufnbr,n_alts,n_wtnodes,wtree);
	}


/* Create SM-frame (must be tree since it must have both sh & crit levels) */

#define SM_MODE 3 // create SM + PS frames and copy sh #1 to all other sh

rcode DTLAPI SML_new_SM_tree_frame(int ufnbr, int type, int n_alts, int n_sh,
		int n_wtnodes, ta_tree wt_tree) {
	rcode rc;

#ifdef LOG_SML
	/* Log function call only if extended pointer logging */
	if (cst_ext) {
		sprintf(msg,"SML_new_SM_tree_frame(%d,%d,%d,%d,%d,%p)\n",ufnbr,type,n_alts,n_sh,n_wtnodes,wt_tree);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	if (sml_loaded)
		return SML_FRAME_IN_USE;
	/* Check input parameters */
	if (n_sh < 2) // not SM
		return SML_INPUT_ERROR;
	if (n_wtnodes < 2) // no tree
		return SML_INPUT_ERROR;
	/* Create tree frame */
	if (type == SM1_FRAME) { // mirrored criteria
		if (rc = DTL_new_SM_tree_frame(ufnbr,SM_MODE,n_alts,n_sh,n_wtnodes,wt_tree))
			return rc;
		}
	else if (type == SM2_FRAME) { // duplicated criteria
		if (rc = DTL_new_DM_tree_frame(ufnbr,n_alts,n_wtnodes,wt_tree))
			return rc;
		}
	else // NSM_FRAME type = DM/PM-frames -> should use SML_new_DM/PM instead
		return SML_INPUT_ERROR;
	/* Post processing */
	sml_type[ufnbr] = type;
	sml_n_sh[ufnbr] = n_sh;
	sml_n_cr[ufnbr] = uf_list[ufnbr]->n_crit/n_sh;
	return SML_OK;
	}


/* Read frame from disk */

rcode DTLAPI SML_read_frame(int ufnbr, int type, int n_sh, char *fn, char *folder) {
	rcode rc;

	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	if (sml_loaded)
		return SML_FRAME_IN_USE;
	/* Check input parameters */
	if ((type < SM2_FRAME) || (type > NSM_FRAME))
		// valid types: SM2,DM,PM (SM & DM stored as PM on file)
		return SML_INPUT_ERROR;
	if (n_sh < 1)
		// protect div and rem operators
		return SML_INPUT_ERROR;
	if ((type == NSM_FRAME) && (n_sh != 1))
		// an NSM_FRAME has only one sh (does not model it)
		return SML_INPUT_ERROR;
	/* Call DTL read function */
	if (DTL_error2(rc = DTL_read_frame(ufnbr,fn,folder,FALSE)))
		return rc;
	/* Check that frame is appropriate */
	if (uf_list[ufnbr]->frame_type != PM_FRAME) {
		DTL_dispose_frame(ufnbr);
		return SML_WRONG_FRAME_TYPE;
		}
	if ((type == SM2_FRAME) && (uf_list[ufnbr]->n_crit%n_sh)) {
		DTL_dispose_frame(ufnbr);
		return SML_FRAME_CORRUPT;
		}
	/* Post processing */
	sml_type[ufnbr] = type;
	sml_n_sh[ufnbr] = n_sh;
	sml_n_cr[ufnbr] = uf_list[ufnbr]->n_crit/n_sh;
	return SML_OK;
	}


 /********************************************************
  *
  *  SML frame handling
  *
  ********************************************************/

static rcode sml_get_frame_type(int ufnbr) {
	rcode rc;
	int dtl_type;

	/* Call underlying DTL function */
	if (rc = DTL_frame_type(ufnbr,&dtl_type))
		return rc;
	/* Evaluate DTL type info */
	if (dtl_type != PM_FRAME)
		return SML_WRONG_FRAME_TYPE;
	return sml_type[ufnbr];
	}


rcode DTLAPI SML_load_frame(int ufnbr) {
	rcode rc,rc2;

	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	if ((ufnbr < 0) || (ufnbr > MAX_FRAMES))
		return SML_FRAME_UNKNOWN;
	if (!sml_type[ufnbr])
		return SML_FRAME_UNKNOWN;
	/* Call DTL function */
	if (DTL_error2(rc = DTL_load_frame(ufnbr)))
		return rc;
	/* Check frame type */
	if (DTL_error2(rc2 = sml_get_frame_type(ufnbr))) {
		DTL_unload_frame();
		return rc2;
		}
	/* Post processing */
	sml_loaded  = ufnbr;
	sml_type[0] = sml_type[ufnbr];
	sml_n_sh[0] = sml_n_sh[ufnbr];
	sml_n_cr[0] = sml_n_cr[ufnbr];
	sml_emethod = NO_EMETHOD;
	sml_renorm = FALSE;
	return rc; // nbr of prob trees
	}


rcode DTLAPI SML_unload_frame() {
	rcode rc;

	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	/* Call CAR (same layer level) */
	CAR_close_W_phull();
	/* Call DTL function */
	if (rc = DTL_unload_frame())
		return rc;
	/* Post processing */
	sml_loaded  = 0;
	sml_type[0] = 0;
	sml_n_sh[0] = 0;
	sml_n_cr[0] = 0;
	return SML_OK;
	}


rcode DTLAPI SML_dispose_frame(int ufnbr) {
	rcode rc;

	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	/* Check input parameter */
	if ((ufnbr < 0) || (ufnbr > MAX_FRAMES))
		return SML_FRAME_UNKNOWN;
	if (!sml_type[ufnbr])
		return SML_FRAME_UNKNOWN;
	/* Call underlying DTL function */
	if (rc = DTL_dispose_frame(ufnbr))
		return rc;
	/* Post processing */
	sml_type[ufnbr] = 0;
	sml_n_sh[ufnbr] = 0;
	sml_n_cr[ufnbr] = 0;
	return SML_OK;
	}


rcode DTLAPI SML_frame_type(int ufnbr) {

	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	/* Check input parameter (ufnbr 0 is current frame) */
	if ((ufnbr < 0) || (ufnbr > MAX_FRAMES))
		return SML_FRAME_UNKNOWN;
	/* Return frame type */
	return sml_get_frame_type(ufnbr);
	}


 /********************************************************
  *
  *  SML base input/output
  *
  ********************************************************/

/* The use of -1.0 and -2.0 in mbox is not supported for SML V-bases due to
 * autoscale but still supported for W-bases and P-bases due to normalisation.
 * For V-bases, it is straightforward just to supply the arithmetical midpoint.
 *
 * Since Dirichlet does not have a conceivable modal value, the workable "most
 * likely" definition consists of mean values for W & P and modal values for V.
 * This is because the psychological concept of "most likely" is ambiguous.
 *
 * Two versions of SML_set_X_base: with/without inconsistent variable detection.
 * The number is zero if not inconsistent or if it does not point to single var. */

rcode DTLAPI SML_set_W_base2(h_vector lobox, h_vector mbox, h_vector upbox, int *inc_var) {
	rcode rc;

#ifdef LOG_SML
	/* Log function call only for extended pointer logging */
	if (cst_ext) {
		sprintf(msg,"SML_set_W_base(%p,%p,%p,%p)\n",lobox,mbox,upbox,inc_var);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	_certify_ptr(inc_var,301);
	if (!sml_active)
		return SML_STATE_ERROR;
	*inc_var = 0;
#ifdef PRE_CHECK
	/* Pre-check weight box inconsistency */
	if ((rc = dtl_set_W_check(lobox,mbox,upbox)) < DTL_OK)
		return rc;
	else if (rc) {
		*inc_var = rc;
		return SML_INCONSISTENT;
		}
#endif
	/* Simplified -> clear current W-base data and start all over.
	 * For advanced use with "don't care"/-2.0 markers, use DTL. */
	if (rc = DTL_reset_W_base())
		return rc; // no reset made
	/* Call combined DTL function pair (mid first to catch inco) */
	if (rc = DTL_set_W_mbox1(mbox))
		return rc;
	if (rc = DTL_set_W_box(lobox,upbox))
		return rc;
#ifdef RENORM_W_BASE
	sml_renorm = TRUE; // activate auto renorm after initial weights set
#endif
	return SML_EC(SML_OK);
	}


rcode DTLAPI SML_set_W_base(h_vector lobox, h_vector mbox, h_vector upbox) {
	int inco;

	return SML_set_W_base2(lobox,mbox,upbox,&inco);
	}


rcode DTLAPI SML_set_P_base2(int crit, h_matrix lobox, h_matrix mbox, h_matrix upbox, int *inc_var) {
	rcode rc;

#ifdef LOG_SML
	/* Log function call only for extended pointer logging */
	if (cst_ext) {
		sprintf(msg,"SML_set_P_base(%d,%p,%p,%p,%p)\n",crit,lobox,mbox,upbox,inc_var);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	_certify_ptr(inc_var,301);
	if (!sml_active)
		return SML_STATE_ERROR;
	*inc_var = 0;
#ifdef PRE_CHECK
	/* Pre-check probability box inconsistency */
	if ((rc = dtl_set_P_check(crit,lobox,mbox,upbox)) < DTL_OK)
		return rc;
	else if (rc) {
		*inc_var = rc;
		return SML_INCONSISTENT;
		}
#endif
	/* Simplified -> clear current P-base data and start all over.
	 * For advanced use with "don't care"/-2.0 markers, use DTL. */
	if (rc = DTL_reset_P_base(crit))
		return rc; // no reset made
	/* Call combined DTL function pair (mid first to catch inco) */
	if (rc = DTL_set_P_mbox1(crit,mbox))
		return rc;
	return DTL_set_P_box(crit,lobox,upbox);
	}


rcode DTLAPI SML_set_P_base(int crit, h_matrix lobox, h_matrix mbox, h_matrix upbox) {
	int inco;

	return SML_set_P_base2(crit,lobox,mbox,upbox,&inco);
	}


#define V_MODE 3 // clear old box + add mbox

rcode DTLAPI SML_set_V_base2(int crit, bool rev, int renorm, h_matrix lobox, h_matrix mbox,
			h_matrix upbox,int *inc_var) {
	rcode rc,rc2;
	int i,crit1;

#ifdef LOG_SML
	/* Log function call for extended pointer logging */
	if (cst_ext) {
		sprintf(msg,"SML_set_V_base(%d,%d,%d,%#p,%#p,%#p,%#p)\n",crit,rev,renorm,lobox,mbox,upbox,inc_var);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	_certify_ptr(inc_var,301);
	if (!sml_active)
		return SML_STATE_ERROR;
	if (!sml_loaded) // protect SM1
		return SML_FRAME_NOT_LOADED;
	/* Init parameters */
	if (renorm > 1) // auto
		renorm = sml_renorm;
	*inc_var = 0;
#ifdef PRE_CHECK
	/* Pre-check value box inconsistency, V-base is preserved.
	 * Note that a V-base with all same mid but diff lo/upbox is
	 * accepted -> this must be taken care of in swing proc. */
	if ((rc = dtl_set_V_check(crit,lobox,mbox,upbox)) < DTL_OK)
		return rc;
	else if (rc) { // inc_var > 0 is an inconsistency forecast
		*inc_var = rc;
		return SML_INCONSISTENT;
		}
#endif
	/* Simplified call -> clear current V-base content.
	 * For advanced use with "don't care"/-2.0, see DTL. */
	if (rc = DTL_reset_V_base(crit))
		return rc;
	if (sml_vsource) { // V-base data generated by machine
		if (SM1) // SM1 is inconceivable, has no machine vsource counterpart
			return SML_WRONG_FRAME_TYPE;
		/* Call combined DTL function pair */
		if (DTL_error2(rc = DTL_set_AV_box(crit,rev,renorm,lobox,upbox)))
			return rc; // DTL_SCALE_CHANGE is ok, keep as signal
		if (rc2 = DTL_set_AV_mbox1(crit,mbox))
			return rc2;
		}
	else if (SM1) { // all stakeholders have the same criteria
		/* Find and use first copy of crit */
		crit1 = (crit-1)%sml_n_cr[0]+1;
		/* Call compound DTL function */
		if (DTL_error2(rc = DTL_set_AV_modal(crit1,V_MODE,rev,renorm,lobox,mbox,upbox)))
			return rc;
		/* Map resulting autoscale to mirrored criteria */
		for (i=2; i<=sml_n_sh[0]; i++)
			if (rc = dtl_copy_AV_crit_scale(crit1,(i-1)*sml_n_cr[0]+crit1))
				return rc;
		}
	else // V-base data generated by human + either SM2 or NSM = standard call
		if (DTL_error2(rc = DTL_set_AV_modal(crit,V_MODE,rev,renorm,lobox,mbox,upbox)))
			return rc;
	return SML_EC(rc); // if DTL_SCALE_CHANGE signal and renorm flag -> new weights
	}


rcode DTLAPI SML_set_V_base(int crit, bool rev, int renorm, h_matrix lobox, h_matrix mbox, h_matrix upbox) {
	int inco;

	return SML_set_V_base2(crit,rev,renorm,lobox,mbox,upbox,&inco);
	}


/* Autoscaled V-hull function (caters for modal and reverse scale) */

static h_matrix xlobo,xmid,xupbo; // intermediate

rcode DTLAPI SML_get_V_hull(int crit, h_matrix lobo, h_matrix mid, h_matrix upbo) {
	rcode rc;
	int i,alts,nodes;

#ifdef LOG_SML
	/* Log function call only for extended pointer logging */
	if (cst_ext) {
		sprintf(msg,"SML_get_V_hull(%d,%p,%p,%p)\n",crit,lobo,mid,upbo);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	/* Call DTL functions */
	if (rc = DTL_get_V_hull(crit,xlobo,xmid,xupbo))
		return rc;
	if (!sml_vsource) // if vsource is machine: refrain from modal
		if (rc = DTL_get_V_modal(crit,xmid))
			return rc;
	/* Collect autoscale parameters */
	if ((alts = DTL_nbr_of_alts()) < DTL_OK)
		return alts;
	/* Autoscale output conversion */
	for (i=1; i<=alts; i++) {
		if ((nodes = DTL_nbr_of_nodes(crit,i)) < DTL_OK)
			return nodes;
		if (rc = DTL_get_AV_user_intervals(crit,ABS_SCALE,nodes+1,xlobo[i],xupbo[i],lobo[i],upbo[i]))
			return rc;
		if (rc = DTL_get_AV_user_vector(crit,ABS_SCALE,nodes+1,xmid[i],mid[i]))
			return rc;
		}
	return SML_EC(SML_OK);
	}


 /*******************************************************
  *
  *  Correlation adjustment function (experimental)
  *  ----------------------------------------------
  *
  *  If there is a correlation between two criteria,
  *  their scales are being overemphasised in total.
  *  This should be compensated for by decreasing
  *  their total weight influence on the MC problem.
  *
  *  The pairwise correlations reside in the upper
  *  right triangle of the input matrix corr_mx.
  *
  *  NOTE: There is nothing to prevent an API user
  *  to call this twice, but the function is provided
  *  as a courtesy at the discretion of the user.
  *
  *  NOTE: This will be generalised to more than two
  *  criteria when the need arises.
  *
  *******************************************************/

static h_vector rn_lobo,rn_mid,rn_upbo;

// MS VC erroneously detects "illegal variable use"
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4701) // keep quiet
#endif

rcode DTLAPI SML_set_W_correlations(cr_matrix corr_mx) {
	rcode rc;
	int i,i2,j;
	int n_found,max_nodes,tnode,lnode;
	double com_lo,com_mid,com_up,renorm;

	/* Load MC uf pointer and get current weights */
	if (rc = DTL_get_W_hull(0,rn_lobo,rn_mid,rn_upbo))
		return rc;
	if ((max_nodes = DTL_nbr_of_weights()) < DTL_OK)
		return max_nodes;
	/* Scan correlation matrix rows */
	for (i=1; i<max_nodes; i++) {
		n_found = 0;
		for (j=i+1; j<=max_nodes; j++)
			/* Seek nbr of corrs for node i */
			if (corr_mx[i][j]) {
				/* Correlaton found, store info */
				n_found++;
				i2 = j; // save correlating crit
				}
		if (!n_found)
			/* No correlation found */
			continue;
		if (n_found > 1)
			/* Not a pure correlation pair */
			return DTL_INPUT_ERROR;
		/* Only one correlaton allowed in a row */
		for (j=i2+1; j<=max_nodes; j++)
			/* Look for more corrs for i2 */
			if (corr_mx[i2][j])
				return SML_INPUT_ERROR;
		/* Have found a correlation pair (i,i2) */
		if (dtl_W_node_parents(i,i2))
			return SML_INPUT_ERROR; // different parents
		if (!rn_mid[i] || !rn_mid[i2])
			return SML_INPUT_ERROR; // zero weights not allowed, no ratio possible
		if (corr_mx[i][i2] * max(rn_mid[i],rn_mid[i2])/min(rn_mid[i],rn_mid[i2]) > 2.0)
			return SML_INPUT_ERROR; // too big difference, rogue user, should be closer to 1.0
		/* Calculate common correlation weight limits */
		com_lo  = corr_mx[i][i2] * (rn_lobo[i]+rn_lobo[i2])/2.0;
		com_mid = corr_mx[i][i2] * (rn_mid[i] +rn_mid[i2]) /2.0;
		com_up  = corr_mx[i][i2] * (rn_upbo[i]+rn_upbo[i2])/2.0;
		/* Split common correlation between the pair */
		rn_lobo[i] = (1.0-corr_mx[i][i2])*rn_lobo[i] + com_lo /2.0;
		rn_mid[i]  = (1.0-corr_mx[i][i2])*rn_mid[i]  + com_mid/2.0;
		rn_upbo[i] = (1.0-corr_mx[i][i2])*rn_upbo[i] + com_up /2.0;
		rn_lobo[i2]= (1.0-corr_mx[i][i2])*rn_lobo[i2]+ com_lo /2.0;
		rn_mid[i2] = (1.0-corr_mx[i][i2])*rn_mid[i2] + com_mid/2.0;
		rn_upbo[i2]= (1.0-corr_mx[i][i2])*rn_upbo[i2]+ com_up /2.0;
		/* Renormalise for all nodes with same parent */
		renorm = 1.0-com_mid;
		/* Seek first sibling node */
		for (tnode=i; tnode; tnode=uf->df->prev[1][tnode])
			lnode=tnode; // i>0 -> will do at least one loop backwards
		/* Renormalise all siblings (lnode is the leftmost one) */
		for (tnode=lnode; tnode; tnode=uf->df->next[1][tnode]) {
			rn_lobo[tnode] /= renorm;
			rn_mid[tnode]  /= renorm;
			rn_upbo[tnode]  = min(rn_upbo[tnode]/renorm,1.0); // overflow protection
			}
		}
	/* All ok -> update renormalised weights
	 * NOTE: weights without a declared midpoint
	 * will have one assigned as a side effect */
	if (rc = DTL_set_W_mbox1(rn_mid))
		return rc;
	return DTL_set_W_box(rn_lobo,rn_upbo);
	}

#ifdef _MSC_VER
#pragma warning(pop)
#endif


 /********************************************************
  *
  *  SML basic evaluations (incl. belief mass functions)
  *
  ********************************************************/

static e_matrix x_result,c_result; // intermediate

static int sml_scale_type(int method, int Aj) {
	int e_method;

	/* PSI is a direct evaluation -> absolute scale
	 * while the others are diffs -> relative scale */
	e_method = method&EMETHOD_MASK;
	if ((e_method==E_DIGAMMA) && !Aj) // PSI in disguise
		return ABS_SCALE;
	else
		return e_method==E_PSI?ABS_SCALE:DIFF_SCALE;
	}


rcode DTLAPI SML_evaluate_frame(int crit, int method, int Ai, int Aj, e_matrix e_result) {
	rcode rc;

#ifdef LOG_SML
	/* Log function call only for extended pointer logging */
	if (cst_ext) {
		sprintf(msg,"SML_evaluate_frame(%d,%d,%d,%d,%p)\n",crit,method,Ai,Aj,e_result);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	_init_assert();
	_dtl_assert(e_result!=x_result,301);
	if (!sml_active)
		return SML_STATE_ERROR;
	sml_emethod = NO_EMETHOD;
	/* Call underlying DTL function */
	if (rc = DTL_evaluate_full(crit,method,Ai,Aj,x_result))
		return rc;
	/* Autoscale conversion parameter */
	sml_scale = sml_scale_type(method,Aj); // absolute or relative scale
	/* Autoscale output conversion */
	if (rc = DTL_get_AV_user_intervals(max(crit,0),sml_scale,MAX_RESULTSTEPS,
			x_result[E_MIN],x_result[E_MAX],e_result[E_MIN],e_result[E_MAX]))
		return rc;
	if (rc = DTL_get_AV_user_vector(max(crit,0),sml_scale,MAX_RESULTSTEPS,
			x_result[E_MID],e_result[E_MID]))
		return rc;
	sml_emethod = method&EMETHOD_MASK;
	sml_ecrit = crit;
	return SML_EC(SML_OK);
	}


rcode DTLAPI SML_get_mass_range(double lo_level, double up_level, double *mass) {
	rcode rc;
	int crit;
	double lo_level01,up_level01;

	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	if (sml_emethod == NO_EMETHOD)
		return SML_OUTPUT_ERROR;
	/* Autoscale input conversion */
	crit  = sml_ecrit;
	if (rc = DTL_get_AV_norm_interval(max(crit,0),sml_scale,lo_level,up_level,&lo_level01,&up_level01))
		return rc;
	/* Call underlying DTL function */
	return DTL_get_mass_range(crit,lo_level01,up_level01,mass);
	}


static rcode sml_get_mass_point(rcode (DTLAPI *func)(int,double,double*), double level, double *mass) {
	rcode rc;
	int crit;
	double level01;

	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	if (sml_emethod == NO_EMETHOD)
		return SML_OUTPUT_ERROR;
	/* Autoscale input conversion */
	crit  = sml_ecrit;
	if (rc = DTL_get_AV_norm_value(max(crit,0),sml_scale,level,&level01))
		return rc;
	/* Call supplied DTL function */
	return (*func)(crit,level01,mass);
	}


rcode DTLAPI SML_get_mass_above(double lo_level, double *mass) {

	/* Call SML pattern function */
	return sml_get_mass_point(&DTL_get_mass_above,lo_level,mass);
	}


rcode DTLAPI SML_get_mass_below(double up_level, double *mass) {

	/* Call SML pattern function */
	return sml_get_mass_point(&DTL_get_mass_below,up_level,mass);
	}


static rcode sml_get_support_mass(rcode (DTLAPI *func)(int,double,double*,double*), double belief_level, double *lobo, double *upbo) {
	rcode rc;
	int crit;
	double lobo01,upbo01;

	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	if (sml_emethod == NO_EMETHOD)
		return SML_OUTPUT_ERROR;
	crit = sml_ecrit;
	/* Call supplied DTL function */
	if (DTL_error2(rc = (*func)(crit,belief_level,&lobo01,&upbo01)))
		return rc;
	/* Autoscale output conversion */
	return DTL_get_AV_user_interval(max(crit,0),sml_scale,lobo01,upbo01,lobo,upbo);
	}


rcode DTLAPI SML_get_support_mass(double belief_level, double *lobo, double *upbo) {

	/* Call SML pattern function */
	return sml_get_support_mass(&DTL_get_support_mass,belief_level,lobo,upbo);
	}


rcode DTLAPI SML_get_support_lower(double belief_level, double *lobo, double *upbo) {

	/* Call SML pattern function */
	return sml_get_support_mass(&DTL_get_support_lower,belief_level,lobo,upbo);
	}


rcode DTLAPI SML_get_support_upper(double belief_level, double *lobo, double *upbo) {

	/* Call SML pattern function */
	return sml_get_support_mass(&DTL_get_support_upper,belief_level,lobo,upbo);
	}


rcode DTLAPI SML_evaluate_cdf(int crit, int Ai, c_vector level, c_vector cdf) {
	rcode rc;
	int i;
	double lvl,step,v_min,v_max;

	/* Check if function can start */
	_init_assert();
	_certify_ptr(level,301);
	_certify_ptr(cdf,302);
	_dtl_assert(level!=cdf,301);
	/* Call SML function */
	if (rc = SML_evaluate_frame(crit,E_PSI,Ai,0,c_result))
		return rc;
	/* Loop through MC scale and collect cdf */
	if (rc = DTL_get_AV_MC_scale(&v_min,&v_max))
		return rc;
	step = (v_max-v_min)/(double)MAX_CDF;
	for (i=0, lvl=v_min; i<=MAX_CDF; i++, lvl+=step) {
		if (i==MAX_CDF)
			lvl = v_max; // catch round-off errors
		if (rc = SML_get_mass_below(lvl,cdf+i))
			return rc;
		level[i] = lvl;
		}
	return SML_EC(SML_OK);
	}


 /********************************************************
  *
  *  SML compound evaluations
  *  ------------------------
  *
  *  Most evaluation functions come in two versions
  *
  *  Basic:    simplest possible parameter set
  *  Extended: larger parameter set (name ends in '2')
  *
  ********************************************************/

static ar_col lo_value01,up_value01; // intermediate

/* Method field: use E_PSI or E_GAMMA. Calling with E_DELTA yields "alternative unknown"
 * since there is no alt pair, while E_DIGAMMA will yield the same as E_PSI (no neg term) */

rcode DTLAPI SML_compare_alternatives(int crit, int method, double belief_level, ar_col lo_value, ar_col up_value) {
	rcode rc;
	int alts,scale;

#ifdef LOG_SML
	/* Log function call only for extended pointer logging */
	if (cst_ext) {
		sprintf(msg,"SML_compare_alternatives(%d,%d,%.3lf,%p,%p)\n",crit,method,belief_level,lo_value,up_value);
		cst_log(msg);
		}
#endif
	/* Call SML function */
	if (rc = DTL_compare_alternatives(crit,method,belief_level,lo_value01,up_value01))
		return rc;
	/* Autoscale conversion parameters */
	scale = sml_scale_type(method,0); // abs or rel scale
	if ((alts = DTL_nbr_of_alts()) < DTL_OK)
		return alts;
	/* Autoscale output conversion */
	return DTL_get_AV_user_intervals(max(crit,0),scale,alts+1,lo_value01,up_value01,lo_value,up_value);
	}


/* Mode: ordering/ranking (controls o_pos output)
 *       0 = ordering (default)
 *       1 = ranking  (olympic) [only valid for _mid, not _omega]
 *
 * Mode: output scale type (controls o_result output)
 *       0 = absolute scale (DTL default)
 *      +2 = output in percent of scale (SML default)
 *      +4 = output in percent of EV range (sums to 100%)
 *
 * Note that mode +2 is used differently in SML and DTL
 *
 * Note the name shifts SML:mid       -> DTL:omega
 *                      SML:omega/1/2 -> DTL:omega1 */

static cr_col xo_result; // intermediate

#define T_OMEGA  0 // eval type DTL:omega  (SML:mid)
#define T_OMEGA1 1 // eval type DTL:omega1 (SML:omega)

static rcode sml_evaluate_omega(rcode (DTLAPI *func)(int,int,cr_col,ci_col),
		 int Ai,int mode, int type, cr_col o_result, ci_col o_pos) {
	rcode rc;
	int o_mode,percent;
	double omega;

#ifdef LOG_SML
	/* Log function call only for extended pointer logging */
	if (cst_ext) {
		sprintf(msg,"SML_evaluate_omega(%d,%d,%d,%p,%p)\n",Ai,mode,type,o_result,o_pos);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	if (type) { // omega1 -> check nbr of sh
		if (!sml_loaded) // protect NSM check
			return SML_FRAME_NOT_LOADED;
		if (!NSM) // does not work with multiple sh
			return SML_WRONG_FRAME_TYPE;
		}
	/* Obtain parameters (note that 0x02 order mode is masked away -> scale type) */
	o_mode = mode&(type?0xFC:0xFD); // ranking [not for omega1] + % of EV range
	percent = mode&0x06; // calculate percent (else absolute -> AV)
	/* Call underlying DTL function */
	if (rc = (*func)(Ai,o_mode,percent?o_result:xo_result,o_pos))
		return rc;
	/* Autoscale output conversion */
	if (!percent) // default = non-percent mode
		if (rc = DTL_get_AV_user_vector(SMC,DIST_SCALE,o_pos[0],xo_result+1,o_result+1))
			return rc;
	if (rc = DTL_get_AV_user_value(SMC,ABS_SCALE,percent?o_result[0]:xo_result[0],&omega))
		return rc;
	o_result[0] = omega;
	return SML_EC(SML_OK);
	}


rcode DTLAPI SML_evaluate_mid(int Ai, int mode, cr_col o_result, ci_col o_rank) {

	/* Call SML pattern function */
	return sml_evaluate_omega(&DTL_evaluate_omega,Ai,mode,T_OMEGA,o_result,o_rank);
	}


#define O1_MODE 2 // default SML result scale is percent of MC

static ci_col x_node;

rcode DTLAPI SML_evaluate_omega(int Ai, cr_col o_result) {

	/* Call SML pattern function */
	return sml_evaluate_omega(&DTL_evaluate_omega1,Ai,O1_MODE,T_OMEGA1,o_result,x_node);
	}


rcode DTLAPI SML_evaluate_omega1(int Ai, cr_col o_result, ci_col o_node) {

	/* Call SML pattern function */
	return sml_evaluate_omega(&DTL_evaluate_omega1,Ai,O1_MODE,T_OMEGA1,o_result,o_node);
	}


rcode DTLAPI SML_evaluate_omega2(int Ai, int mode, cr_col o_result, ci_col o_node) {

	/* Call SML pattern function (omega1 with mode parameter) */
	return sml_evaluate_omega(&DTL_evaluate_omega1,Ai,mode,T_OMEGA1,o_result,o_node);
	}


static ar_matrix xdelta_value,xdelta_mass; // placeholders
static ai_col xgamma_rank,xomega_rank;
static ar_col xgamma_value,xomega_value;

#define DMASS_MODE 2 // GAMMA rank with tiebreak

rcode DTLAPI SML_delta_mass2(int crit, int mode, ar_matrix delta_mass, ai_col delta_order) {
	rcode rc;
	int Ai,Aj,alts;

#ifdef LOG_SML
	/* Log function call only for the purpose of extended pointer logging */
	if (cst_ext) {
		sprintf(msg,"SML_delta_mass(%d,%d,%p,%p)\n",crit,mode,delta_mass,delta_order);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	/* Call underlying DTL function (simplified modes: only 0 or 1) */
	if (rc = DTL_delta_mass(crit,mode&0x01,xdelta_value,xdelta_mass))
		return rc;
	/* Collect mass order conversion parameter */
	if ((alts = DTL_nbr_of_alts()) < DTL_OK)
		return alts;
	/* Rank with tiebreak to get a crisp and definitive ranking order */
	if (rc = DTL_rank_alternatives(crit,DMASS_MODE,0.0,0.0,xgamma_rank,xomega_rank,xgamma_value,xomega_value))
		return rc;
	/* Rearrange according to ranking order */
	for (Ai=1; Ai<=alts; Ai++) {
		delta_order[xgamma_rank[Ai]] = Ai;
		for (Aj=1; Aj<=alts; Aj++)
			delta_mass[xgamma_rank[Ai]][xgamma_rank[Aj]] = xdelta_mass[Ai][Aj];
		}
	return SML_EC(SML_OK);
	}


#define R_IPOL 1 // row interpolation is default

rcode DTLAPI SML_delta_mass(int crit, ar_matrix delta_mass, ai_col delta_order) {

	return SML_delta_mass2(crit,R_IPOL,delta_mass,delta_order);
	}


/* Mode field: 0 = soft (olympic) ranking
 *             1 = strict/hard ranking
 *            +2 = tolerances are values (default: percent) */

static ar_col gamma_value01,omega_value01; // intermediate

rcode DTLAPI SML_rank_alternatives(int crit, int mode, double gamma_tolerance, double omega_tolerance, 
		ai_col gamma_rank, ai_col omega_rank, ar_col gamma_value, ar_col omega_value) {
	rcode rc;
	int alts;
	double gamma_tolerance01,omega_tolerance01;

#ifdef LOG_SML
	/* Log function call only for extended pointer logging */
	if (cst_ext) {
		sprintf(msg,"SML_rank_alternatives(%d,%d,%.3lf,%.3lf,%p,%p,%p,%p)\n",crit,mode,gamma_tolerance,omega_tolerance, 
				gamma_rank,omega_rank,gamma_value,omega_value);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	if (mode&0x02) {
		/* Autoscale input conversion (tolerances expressed as values) */
		if (rc = DTL_get_AV_norm_value(max(crit,0),DIST_SCALE,gamma_tolerance,&gamma_tolerance01))
			return rc;
		if (rc = DTL_get_AV_norm_value(max(crit,0),DIST_SCALE,omega_tolerance,&omega_tolerance01))
			return rc;
		}
	else {
		/* Tolerances are percent per default = coincides with [0,1] norm scale */
		gamma_tolerance01 = gamma_tolerance;
		omega_tolerance01 = omega_tolerance;
		}
	/* Call underlying DTL function (simplified mode parameter: only 0 and 1 used) */
	if (DTL_error2(rc = DTL_rank_alternatives(crit,mode&0x01,gamma_tolerance01,omega_tolerance01, 
			gamma_rank,omega_rank,gamma_value01,omega_value01))) // simplified mode
		return rc;
	/* Collect autoscale conversion parameter */
	if ((alts = DTL_nbr_of_alts()) < DTL_OK)
		return alts;
	/* Autoscale output conversion */
	if (rc = DTL_get_AV_user_vector(max(crit,0),DIFF_SCALE,alts+1,gamma_value01,gamma_value))
		return rc;
	if (rc = DTL_get_AV_user_vector(max(crit,0),ABS_SCALE,alts+1,omega_value01,omega_value))
		return rc;
	return SML_EC(SML_OK);
	}


#define RANK_MODE 1 // strict ranking

rcode DTLAPI SML_rank_gamma(int crit, ai_col gamma_rank, ar_col gamma_value) {
	rcode rc;

	rc = SML_rank_alternatives(crit,RANK_MODE,0.0,0.0,gamma_rank,xomega_rank,gamma_value,xomega_value);
	return rc==SML_DIFFERING_RANKS?SML_OK:rc;
	}


rcode DTLAPI SML_rank_omega(int crit, ai_col omega_rank, ar_col omega_value) {
	rcode rc;

	rc =  SML_rank_alternatives(crit,RANK_MODE,0.0,0.0,xgamma_rank,omega_rank,xgamma_value,omega_value);
	return rc==SML_DIFFERING_RANKS?SML_OK:rc;
	}


rcode DTLAPI SML_daisy_chain2(int crit, int mode, ai_col omega_rank, ar_col daisy_value, ar_col omega_value) {
	rcode rc;
	int alts;

#ifdef LOG_SML
	/* Log function call only for extended pointer logging */
	if (cst_ext) {
		sprintf(msg,"SML_daisy_chain(%d,%d,%p,%p,%p)\n",crit,mode,omega_rank,daisy_value,omega_value);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	/* Call underlying DTL function (note: radius not used, else must convert) */
	if (rc = DTL_daisy_chain1(crit,mode&0x01,omega_rank,daisy_value,omega_value01))
		return rc;
	/* Collect autoscale conversion parameter */
	if ((alts = DTL_nbr_of_alts()) < DTL_OK)
		return alts;
	/* Autoscale output conversion */
	if (rc = DTL_get_AV_user_vector(max(crit,0),mode&0x01?DIFF_SCALE:ABS_SCALE,alts+1,omega_value01,omega_value))
		return rc;
	return SML_EC(SML_OK);
	}


#define DAISY_MODE 1

rcode DTLAPI SML_daisy_chain(int crit, ai_col daisy_rank, ar_col daisy_value) {

	return SML_daisy_chain2(crit,DAISY_MODE,daisy_rank,daisy_value,xomega_value);
	}


#define PIE_MODE 1 // modern (not compat with old DMC)

rcode DTLAPI SML_pie_chart2(int crit, double moderation, ar_col pie_value) {
	int n_alts;

#ifdef LOG_SML
	/* Log function call only for extended pointer logging */
	if (cst_ext) {
		sprintf(msg,"SML_pie_chart(%d,%.3lf,%p)\n",crit,moderation,pie_value);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	/* Sanitise input parameter */
	if (fabs(moderation) > 0.5)
		return SML_INPUT_ERROR;
#ifdef PIE_COMPAT
	/* Make pie and delta_mass equal for 2 alts */
	if ((n_alts = DTL_nbr_of_alts()) < DTL_OK)
		return n_alts;
	if (n_alts < 3)
		moderation = 0.0;
#endif
	/* Call underlying DTL function */
	if (moderation < 0.0) // only anchor modified (basic)
		return DTL_pie_chart1(crit,moderation,pie_value);
	else // both anchor and daisy chain modified (recommended)
		return DTL_pie_chart2(crit,PIE_MODE,moderation,moderation,pie_value);
	}


#define PIE_MOD 0.2 // 20% default mix of crisp results for display (benefit of doubt)

rcode DTLAPI SML_pie_chart1(int crit, ar_col pie_value) {

	return SML_pie_chart2(crit,PIE_MOD,pie_value); // moderated mass displayed
	}


rcode DTLAPI SML_pie_chart(int crit, ar_col pie_value) {

	return SML_pie_chart2(crit,0.0,pie_value); // raw crisp mass displayed
	}


 /********************************************************
  *
  *  SML sensitivity analyses
  *
  ********************************************************/

#define T_LOCAL  0 // tornado on local crit scale
#define T_GLOBAL 1 // tornado on global SMC scale
#define T_MODE 0   // default mode is EV + midpoint kept

static h_matrix t_lobo01,t_upbo01; // intermediate

rcode DTLAPI SML_get_W_tornado2(int mode, h_matrix t_lobo, h_matrix t_upbo) {
	rcode rc;
	int i,mass,wts,alts;

#ifdef LOG_SML
	/* Log function call only for extended pointer logging */
	if (cst_ext) {
		sprintf(msg,"SML_get_W_tornado(%d,%p,%p)\n",mode,t_lobo,t_upbo);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	/* Trim parameter */
	mass = mode&0x02;
	/* Call underlying DTL function */
	if (rc = DTL_get_W_tornado(mode,mass?t_lobo:t_lobo01,mass?t_upbo:t_upbo01))
		return rc;
	/* Collect autoscale conversion parameters */
	if ((wts = DTL_nbr_of_weights()) < DTL_OK)
		return wts;
	if ((alts = DTL_nbr_of_alts()) < DTL_OK)
		return alts;
	/* Autoscale output conversion */
	if (!mass) // EV
		for (i=1; i<=alts; i++)
 			if (rc = DTL_get_AV_user_intervals(SMC,DIFF_SCALE,wts+1,t_lobo01[i],t_upbo01[i],t_lobo[i],t_upbo[i]))
				return rc;
	return SML_EC(SML_OK);
	}


rcode DTLAPI SML_get_W_tornado(h_matrix t_lobo, h_matrix t_upbo) {

	return SML_get_W_tornado2(T_MODE,t_lobo,t_upbo);
	}


static rcode sml_get_PV_tornado(rcode (DTLAPI *func)(int,int,h_matrix,h_matrix), 
		int global, int crit, int mode, h_matrix t_lobo, h_matrix t_upbo) {
	rcode rc;
	int i,mass,nodes,alts;

#ifdef LOG_SML
	/* Log function call only for extended pointer logging */
	if (cst_ext) {
		sprintf(msg,"SML_get_PV_tornado(%d,%d,%d,%p,%p)\n",global,crit,mode,t_lobo,t_upbo);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	/* Trim parameter */
	mass = mode&0x02;
	/* Call underlying DTL function */
	if (rc = (*func)(crit,mode,mass?t_lobo:t_lobo01,mass?t_upbo:t_upbo01))
		return rc;
	/* Collect autoscale conversion parameter */
	if ((alts = DTL_nbr_of_alts()) < DTL_OK)
		return alts;
	/* Autoscale output conversion */
	if (!mass) // EV
		for (i=1; i<=alts; i++) {
			/* One more conversion parameter */
			if ((nodes = DTL_nbr_of_nodes(crit,i)) < DTL_OK)
				return nodes;
			if (rc = DTL_get_AV_user_intervals(global?SMC:crit,DIFF_SCALE,nodes+1,t_lobo01[i],t_upbo01[i],t_lobo[i],t_upbo[i]))
				return rc;
			}
	return SML_EC(SML_OK);
	}


rcode DTLAPI SML_get_P_tornado2(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo) {

	/* Call SML pattern function */
	return sml_get_PV_tornado(&DTL_get_P_tornado,T_LOCAL,crit,mode,t_lobo,t_upbo);
	}


rcode DTLAPI SML_get_P_tornado(int crit, h_matrix t_lobo, h_matrix t_upbo) {

	/* Call SML pattern function */
	return sml_get_PV_tornado(&DTL_get_P_tornado,T_LOCAL,crit,T_MODE,t_lobo,t_upbo);
	}


rcode DTLAPI SML_get_MCP_tornado2(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo) {

	/* Call SML pattern function */
	return sml_get_PV_tornado(&DTL_get_MCP_tornado,T_GLOBAL,crit,mode,t_lobo,t_upbo);
	}


rcode DTLAPI SML_get_MCP_tornado(int crit, h_matrix t_lobo, h_matrix t_upbo) {

	/* Call SML pattern function */
	return sml_get_PV_tornado(&DTL_get_MCP_tornado,T_GLOBAL,crit,T_MODE,t_lobo,t_upbo);
	}


rcode DTLAPI SML_get_V_tornado2(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo) {

	/* Call SML pattern function */
	return sml_get_PV_tornado(&DTL_get_V_tornado,T_LOCAL,crit,mode,t_lobo,t_upbo);
	}


rcode DTLAPI SML_get_V_tornado(int crit, h_matrix t_lobo, h_matrix t_upbo) {

	/* Call SML pattern function */
	return sml_get_PV_tornado(&DTL_get_V_tornado,T_LOCAL,crit,T_MODE,t_lobo,t_upbo);
	}


rcode DTLAPI SML_get_MCV_tornado2(int crit, int mode, h_matrix t_lobo, h_matrix t_upbo) {

	/* Call SML pattern function */
	return sml_get_PV_tornado(&DTL_get_MCV_tornado,T_GLOBAL,crit,mode,t_lobo,t_upbo);
	}


rcode DTLAPI SML_get_MCV_tornado(int crit, h_matrix t_lobo, h_matrix t_upbo) {

	/* Call SML pattern function */
	return sml_get_PV_tornado(&DTL_get_MCV_tornado,T_GLOBAL,crit,T_MODE,t_lobo,t_upbo);
	}


/* Cons influence is similar - but not identical - to the alt contribution omega1.
 * This function asks "how important is one crit to various alts' EVs?" while the
 * contribution function asks "which crits contribute to one EV and by how much?"
 *
 * mode=0 is local crit scale, mode=1 is global MC scale */

rcode DTLAPI SML_get_cons_influence(int crit, int mode, h_matrix result) {
	rcode rc;
	int i,nodes,alts;

#ifdef LOG_SML
	/* Log function call only for extended pointer logging */
	if (cst_ext) {
		sprintf(msg,"SML_get_cons_influence(%d,%d,%p)\n",crit,mode,result);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	/* Call underlying DTL function */
	if (rc = DTL_get_cons_influence(crit,mode,t_lobo01))
		return rc;
	/* Collect autoscale conversion parameter */
	if ((alts = DTL_nbr_of_alts()) < DTL_OK)
		return alts;
	/* Autoscale output conversion */
	for (i=1; i<=alts; i++) {
		/* One more conversion parameter */
		if ((nodes = DTL_nbr_of_nodes(crit,i)) < DTL_OK)
			return nodes;
		if (rc = DTL_get_AV_user_vector(mode?SMC:crit,DIFF_SCALE,nodes+1,t_lobo01[i],result[i]))
			return rc;
		}
	return SML_EC(SML_OK);
	}


 /********************************************************
  *
  *  Error text handling (null-terminated/length-preceded)
  *
  ********************************************************/

char* DTLAPI SML_get_errtxt2(rcode drc, int style) {

	if (style)
		return DTL_get_errtxt_p(drc); // Pascal-style text
	else
		return DTL_get_errtxt(drc);   // C-style text
	}


char* DTLAPI SML_get_errtxt(rcode drc) {

	return SML_get_errtxt2(drc,SML_STRG);
	}


 /**********************************************************
  *
  *  Conversion from SML shold & crit nbrs to DTL crit nbr
  *  -----------------------------------------------------
  *
  *  Two versions - with and without type and range checks.
  *  The variant without checks is complemented by a macro
  *  and intended for situations when a subsequent call to
  *  another SML function will validate the type and range.
  *  None of the variants use the standard error code range.
  *
  **********************************************************/

int DTLAPI sml_nc() {

	if (!sml_loaded)
		return 0; // no frame
	return sml_n_cr[0];
	}


int DTLAPI SML_crit_nbr(int sh, int crit) {

	/* Check if function can start */
	if (!sml_active)
		return 0;
	if (!sml_loaded)
		return 0; // no frame
	/* Try to return DTL crit nbr */
	if (SM1 || SM2) { // SM-frame
		if ((sh < 1) || (sh > sml_n_sh[0]))
			return 0; // out of range
		if ((crit < 1) || (crit > sml_n_cr[0]))
			return 0; // out of range
		return (sh-1)*sml_n_cr[0]+crit;
		}
	else if (NSM) // DM/PM-frame
		return crit;
	else  // unknown frame
		return 0;
	}


 /********************************************************
  *
  *  Stakeholder identification (only for test/debug)
  *
  ********************************************************/

rcode DTLAPI SML_is_stakeholder(int node) {

	/* Check if function can start */
	if (!sml_active)
		return SML_STATE_ERROR;
	if (!sml_loaded) // protect NSM
		return SML_FRAME_NOT_LOADED;
	if (NSM) // does not contain sh
		return SML_WRONG_FRAME_TYPE;
	/* Check input parameter */
	if ((node < 1) || (node > DTL_nbr_of_weights()))
		return SML_INPUT_ERROR;
	/* Return re/im-status */
	return !dtl_real_W_crit(node);
	}
