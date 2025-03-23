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
 *   File: DTLvbase.c
 *
 *   Purpose: value base interface to TCL
 *
 *
 *   Functions exported outside DTL
 *   ------------------------------
 *   DTL_add_V_statement
 *   DTL_change_V_statement
 *   DTL_replace_V_statement
 *   DTL_delete_V_statement
 *   DTL_add_V_mid_statement
 *   DTL_delete_V_mid_statement
 *   DTL_set_V_box
 *   DTL_set_V_mbox/1
 *   DTL_remove_V_mbox
 *   DTL_set_V_modal
 *   DTL_get_V_hull
 *   DTL_get_V_modal
 *   DTL_check_V_modality
 *   DTL_get_V_modality_matrix
 *   DTL_reset_V_base
 *
 *   Functions outside module, debug use
 *   -----------------------------------
 *   DTI_real_V_node
 *   DTI_show_V_base
 *   DTI_show_V_box
 *   DTI_show_V_mbox
 *
 *   Functions outside module, within DTL
 *   ------------------------------------
 *   dtl_mid_to_modal
 *   dtl_set_V_check
 *   dtl_set_V_mbox_rels
 *   dtl_nbr_V_midpoints
 *
 *   Functions internal to module
 *   ----------------------------
 *   dtl_modal_to_mid
 *   dtl_check_V_modality
 *   print_V_stmt
 *
 */

#include "DTL.h"
#include "DTLinternal.h"


 /********************************************************
  *
  *  Internal data
  *
  ********************************************************/

/* Hull caches */
static d_row vhlobo,vhupbo,V_mid;


 /*********************************************************
  *
  *  Configuration parameters
  *
  *  V_DEGEN_SCALE - default=OFF (used only by SML)
  *
  *  Setting the flag is to accept a degenerated V-base
  *  with all statements exactly the same = meaningless
  *
  *  V_MODAL_RANGE - default=OFF
  *
  *  Setting it to ON makes DTL incompatible with DMC by
  *  only accepting midpoints with modal representations.
  *
  *  NOTE: The optimal place for such a check routine is
  *  in load_V in TCLvbase.c but since there are not yet
  *  any intra-variable dependencies in the value base,
  *  it works rather well to have the check here instead.
  *  The weak spots are shrinking hulls after mid is set.
  *  So in the long run, it should migrate to TCLvbase.c.
  *  Probably at the next overhaul, if there ever is one.
  *
  *  NOTE2: The safest way is using DTL_check_V_modality.
  *
  *  NOTE3: The best UI solution is to use "most likely"
  *  as the terminology and map that onto V-modal values
  *  here, but onto mean values for W and P where modals
  *  are not even well-defined for a Dirichlet function.
  *
  *  V_MODAL_CHECK - default=OFF (used only by SML)
  *
  *  Similar to V_MODAL_RANGE but for pre-check in SML
  *
  *********************************************************/

#define V_DEGEN_SCALE
#define noV_MODAL_RANGE
#define noV_MODAL_CHECK


 /*********************************************************
  *
  *  Modality functions
  *
  *********************************************************/

/* Conversion functions: midpoint <-> triangular modal value.
   These functions depend on the entire interval [lobo,upbo].

   Return information on error:
   -1.0: intermediate V-base node
   -1.0: box entry is empty
   -2.0: no input info
   -2.0: statement inconsistent
   -3.0: mid out of range wrt triangle base (only dtl_mid_to_modal) */

// not static - also used by W & P bases for modal output (using oor==0)

double dtl_mid_to_modal(double lobo, double mid, double upbo, int oor) {
	double mean;

	/* Check input before conversion */
	if (mid == -2.0)
		// no input info
		return -2.0;
	if ((lobo < -DTL_EPS) || (mid < -DTL_EPS) || (upbo < -DTL_EPS))
		// im-node (for V-base only) or empty
		return -1.0;
	if ((lobo > mid+DTL_EPS) || (mid > upbo+DTL_EPS))
		// inconsistent
		return -2.0;
	if (oor) // out-of-range detection
		if ((3.0*mid < 2.0*lobo+upbo-DTL_EPS) || (3.0*mid > lobo+2.0*upbo+DTL_EPS))
			// modal value ('most likely') for triangle has overhang (lies beyond base)
			return -3.0;
	/* Return triangular modal value (trim out-of-range) */
	mean = max(mid,(2.0*lobo+upbo)/3.0);
	mean = min(mean,(lobo+2.0*upbo)/3.0);
	return min(max(3.0*mean-lobo-upbo,0.0),1.0); // catch round-off errors
	}


double dtl_modal_to_mid(double lobo, double modal, double upbo) {

	/* Check input before conversion */
	if (modal == -2.0)
		// no input info
		return -2.0;
	if ((lobo < -DTL_EPS) || (modal < -DTL_EPS) || (upbo < -DTL_EPS))
		// im-node or empty
		return -1.0;
	if ((lobo > modal+DTL_EPS) || (modal > upbo+DTL_EPS))
		// inconsistent
		return -2.0;
	/* Return triangular mid (no out-of-range in this direction) */
	return (lobo+modal+upbo)/3.0; // no round-off errors
	}


 /*********************************************************
  *
  *  Value base functionality
  *
  *********************************************************/

 /*
  * Call semantics: Add the user value constraint v(crit:alt:cons) =
  * [lobo,upbo] to the value base within the decision frame. The base
  * is checked for consistency with respect to the new interval. In
  * case of inconsistency, nothing is added to the base.
  */

rcode DTLAPI DTL_add_V_statement(int crit, struct user_stmt_rec* ustmtp) {
	struct stmt_rec stmt;

	/* Begin single thread semaphore */
	_smx_begin("AVS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_add_V_statement(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(ustmtp,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if (load_PV_stmt(crit,ustmtp,&stmt,'V'))
		return dtl_error(DTL_STMT_ERROR);
	/* Add statement */
	if (call(TCL_add_V_constraint(uf->df,&stmt),"TCL_add_V_constraint"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return uf->df->V_base->n_stmts;
	}


 /*
  * Call semantics: Change the existing user constraint v(crit:alt:cons) =
  * [old_lobo,old_upbo] to v(crit:alt:cons) = [lobo,upbo] in the value base.
  * The base is checked for consistency with respect to the change. In case
  * of inconsistency, nothing is changed in the base.
  */

rcode DTLAPI DTL_change_V_statement(int crit, int stmt_number, double lobo, double upbo) {

	/* Begin single thread semaphore */
	_smx_begin("CVS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_change_V_statement(%d,%d,%.3lf,%.3lf)\n",crit,stmt_number,lobo,upbo);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	/* Change statement */
	if (call(TCL_change_V_constraint(uf->df,stmt_number,lobo,upbo),"TCL_change_V_constraint"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*
  * Call semantics: Replace the user value constraint (interval)
  * with v(crit:alt:cons) = [lobo,upbo] in the value base. The base
  * is checked for consistency with respect to the new interval. In
  * case of inconsistency, nothing is replaced in the base.
  */

rcode DTLAPI DTL_replace_V_statement(int crit, int stmt_number, struct user_stmt_rec* ustmtp) {
	struct stmt_rec stmt;

	/* Begin single thread semaphore */
	_smx_begin("RVS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_replace_V_statement(%d,%d)\n",crit,stmt_number);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(ustmtp,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if (load_PV_stmt(crit,ustmtp,&stmt,'V'))
		return dtl_error(DTL_STMT_ERROR);
	/* Replace statement */
	if (call(TCL_replace_V_constraint(uf->df,stmt_number,&stmt),"TCL_replace_V_constraint"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*
  * Call semantics: The user value constraint (interval) with
  * position stmt_nbr in the value base is deleted from the base.
  * All constraints with higher positions within the base are
  * shifted one position down.
  */

rcode DTLAPI DTL_delete_V_statement(int crit, int stmt_number) {

	/* Begin single thread semaphore */
	_smx_begin("DVS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_delete_V_statement(%d,%d)\n",crit,stmt_number);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	/* Delete statement */
	if (call(TCL_delete_V_constraint(uf->df,stmt_number),"TCL_delete_V_constraint"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return uf->df->V_base->n_stmts;
	}


 /*
  * Call semantics: Enter the user value core midpoint
  * v(crit:alt:cons) = [lobo,upbo] into the decision frame.
  * The base is checked for consistency with respect to the
  * new core midpoint. In case of inconsistency, nothing is
  * added to the base. Midpoint is mean, not modal (cf. BRS).
  */

rcode DTLAPI DTL_add_V_mid_statement(int crit, struct user_stmt_rec* ustmtp) {
#ifdef V_MODAL_RANGE
	int k;
#endif
	struct stmt_rec stmt;

	/* Begin single thread semaphore */
	_smx_begin("AVM");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_add_V_mid_statement(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(ustmtp,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if (load_PV_stmt(crit,ustmtp,&stmt,'V'))
		return dtl_error(DTL_STMT_ERROR);
#ifdef V_MODAL_RANGE
	/* Get hull and index for modal check */
	if (call(TCL_get_V_hull(uf->df,vhlobo,vhupbo),"TCL_get_V_hull"))
		return dtl_kernel_error();
	if ((k = TCL_get_tot_index(ustmtp->alt[1],ustmtp->cons[1])) < 1)
		return dtl_error(DTL_INPUT_ERROR);
	/* Perform modal check on input */
	if (dtl_mid_to_modal(vhlobo[k],(ustmtp->lobo+ustmtp->upbo)/2.0,vhupbo[k],1) == -3.0)
		return dtl_error(DTL_INCONSISTENT);
#endif
	/* Set midpoint */
	if (call(TCL_add_V_mstatement(uf->df,&stmt),"TCL_add_V_mstatement"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*
  * Call semantics: Remove the user value core midpoint
  * for v(crit:alt:cons) from the decision frame.
  */

rcode DTLAPI DTL_delete_V_mid_statement(int crit, struct user_stmt_rec* ustmtp) {
	struct stmt_rec stmt;

	/* Begin single thread semaphore */
	_smx_begin("DVM");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_delete_V_mid_statement(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(ustmtp,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if (load_PV_stmt(crit,ustmtp,&stmt,'V'))
		return dtl_error(DTL_STMT_ERROR);
	/* Remove midpoint */
	if (call(TCL_delete_V_mstatement(uf->df,&stmt),"TCL_delete_V_mstatement"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*********************************************************
  *
  *  Set box/mbox functions
  *
  *********************************************************/

static d_row box_lobo,box_upbo;


 /*********************************************************
  *
  *  dtl_set_V_check checks for inconsistencies (only SML)
  *
  *  Returns rc = DTL_OK if consistent and spanning a scale
  *             > 0 if inconsistent (offending variable)
  *             = DTL_INPUT_ERROR if degenerated scale
  *             < 0 if other error
  *
  *  Compilation flag DEGEN_SCALE disables DTL_INPUT_ERROR
  *
  *  NOTE: This works since SML does not allow for entering
  *        V-statements one at a time but only all at once.
  *
  *********************************************************/

rcode dtl_set_V_check(int crit, h_matrix lobox, h_matrix mbox, h_matrix upbox) {
	int i,j,k;
	struct d_frame *df;

	/* Check if function can start */
	_init_assert();
	_certify_ptr(lobox,1);
	_certify_ptr(mbox,2);
	_certify_ptr(upbox,3);
	_dtl_assert(lobox!=upbox,1);
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameters */
	if (load_df1(crit))
		return DTL_CRIT_UNKNOWN;
	df = uf->df;
	/* Check box intervals (empty slots not permitted) */
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++)
#ifdef V_MODAL_CHECK // modality check
			if (dtl_mid_to_modal(lobox[i][j],mbox[i][j],upbox[i][j],1) == -3.0)
#else // standard range check
			if ((lobox[i][j]>mbox[i][j]) || (mbox[i][j]>upbox[i][j]))
#endif
				return k;
#ifdef V_DEGEN_SCALE
	return DTL_OK;
#endif
	/* Check non-equal intervals (exactly same not permitted) */
	for (i=1; i<df->n_alts; i++) // check first consequence
		if ((lobox[i][1]!=lobox[i+1][1]) || (mbox[i][1]!=mbox[i+1][1]) || (upbox[i][1]!=upbox[i+1][1]))
				return DTL_OK;
	for (i=1; i<=df->n_alts; i++)
		for (j=1; j<df->tot_cons[i]; j++) // check rest of consequences
			if ((lobox[i][j]!=lobox[i][j+1]) || (mbox[i][j]!=mbox[i][j+1]) || (upbox[i][j]!=upbox[i][j+1]))
				return DTL_OK;
	/* No difference found = scale not spanned */
	return DTL_INPUT_ERROR;
	}


rcode DTLAPI DTL_set_V_box(int crit, h_matrix lobox, h_matrix upbox) {
	int i,j,k;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("SVB");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_set_V_box(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(lobox,1);
	_certify_ptr(upbox,2);
	_dtl_assert(lobox!=upbox,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	df = uf->df;
	/* Transfer and log box */
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			box_lobo[k] = lobox[i][j];
			box_upbo[k] = upbox[i][j];
			if (cst_on && ((lobox[i][j]!=0.0) || (upbox[i][j]!=1.0))) {
				sprintf(msg,"    V%d.%d.%-2d [%.3lf %.3lf] (%le)\n",crit,i,j,
						lobox[i][j],upbox[i][j],upbox[i][j]-lobox[i][j]);
				cst_log(msg);
				}
			}
	if (call(TCL_set_V_box(df,box_lobo,box_upbo),"TCL_set_V_box"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


static d_row mbox_lobo,mbox_upbo;

rcode DTLAPI DTL_set_V_mbox(int crit, h_matrix lobox, h_matrix upbox) {
	int i,j,k;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("SVMB");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_set_V_mbox(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(lobox,1);
	_certify_ptr(upbox,2);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	df = uf->df;
#ifdef V_MODAL_RANGE
	/* Get hull for modal reference */
	if (call(TCL_get_V_hull(df,vhlobo,vhupbo),"TCL_get_V_hull"))
		return dtl_kernel_error();
#endif
	/* Transfer and log midpoint box */
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
#ifdef V_MODAL_RANGE
			if (dtl_mid_to_modal(vhlobo[k],(lobox[i][j]+upbox[i][j])/2.0,vhupbo[k],1) == -3.0)
				return dtl_error(DTL_INCONSISTENT);
#endif
			mbox_lobo[k] = lobox[i][j];
			mbox_upbo[k] = upbox[i][j];
			if (cst_on && ((lobox[i][j]!=-2.0) || (upbox[i][j]!=-2.0))) {
				sprintf(msg,"    V%d.%d.%-2d [%.3lf %.3lf] (%le)\n",crit,i,j,
						lobox[i][j],upbox[i][j],upbox[i][j]-lobox[i][j]);
				cst_log(msg);
				}
			}
	if (call(TCL_set_V_mbox(df,mbox_lobo,mbox_upbo),"TCL_set_V_mbox"))
		return dtl_kernel_error();
	uf->V_n_rels[crit] = 0;
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


/* DTL_set_V_mbox1: pointwise mid, not wide */

rcode DTLAPI DTL_set_V_mbox1(int crit, h_matrix mbox) {

	return DTL_set_V_mbox(crit,mbox,mbox);
	}


/* dtl_set_V_mbox_rels: store nbr of V-relations */

rcode dtl_set_V_mbox_rels(int crit, int V_n_rels, h_matrix lobox, h_matrix upbox) {
	rcode rc;

	rc = DTL_set_V_mbox(crit,lobox,upbox);
	if (!rc)
		uf->V_n_rels[crit] = V_n_rels;
	return rc;
	}


/* Mode: 0 = default
 *      +1 = clear mbox
 *      +2 = set box */

rcode DTLAPI DTL_set_V_modal(int crit, int mode, h_matrix lobox, h_matrix modalx, h_matrix upbox) {
	int i,j,k,mcount=0;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("SVM");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_set_V_modal(%d,%d)\n",crit,mode);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(lobox,1);
	_certify_ptr(modalx,2);
	_certify_ptr(upbox,3);
	_dtl_assert(lobox!=modalx,1);
	_dtl_assert(modalx!=upbox,2);
	_dtl_assert(lobox!=upbox,3);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if (mode&0xFC)
		return dtl_error(DTL_INPUT_ERROR);
	df = uf->df;
	/* First: clear current mbox */
	if (mode&0x01) {
		for (k=1,i=1; i<=df->n_alts; i++)
			for (j=1; j<=df->tot_cons[i]; j++,k++)
				mbox_lobo[k] = -1.0;
		if (call(TCL_set_V_mbox(df,mbox_lobo,mbox_lobo),"TCL_set_V_mbox"))
			return dtl_kernel_error();
		}
	/* Next: create new box & mbox */
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			if (modalx[i][j] >= 0.0) { // mbox entry active
				if ((lobox[i][j] > modalx[i][j]) || (modalx[i][j] > upbox[i][j]))
					return dtl_error(DTL_INCONSISTENT);
				mcount++;
				}
			box_lobo[k] = lobox[i][j];
			mbox_lobo[k] = dtl_modal_to_mid(lobox[i][j],modalx[i][j],upbox[i][j]);
			box_upbo[k] = upbox[i][j];
			if (cst_on && TCL_get_V_index(df,i,j) && 
				 ((lobox[i][j]!=0.0) || (modalx[i][j]!=-2.0) || (upbox[i][j]!=1.0))) {
				if (upbox[i][j]-lobox[i][j] > 0.0)
					sprintf(msg,"    V%d.%d.%d = [%.3lf %.3lf %.3lf] (%le)\n",crit,i,j,
							lobox[i][j],modalx[i][j],upbox[i][j],(upbox[i][j]-lobox[i][j]));
				else
					sprintf(msg,"    V%d.%d.%d = [%.3lf %.3lf %.3lf]\n",crit,i,j,
							lobox[i][j],modalx[i][j],upbox[i][j]);
				cst_log(msg);
				}
			}
	if (!mcount) // no info entered
		return dtl_error(DTL_INPUT_ERROR);
	/* Last: enter new mbox & box */
	if (call(TCL_set_V_mbox(df,mbox_lobo,mbox_lobo),"TCL_set_V_mbox"))
		return dtl_kernel_error();
	if (mode&0x02)
		if (call(TCL_set_V_box(df,box_lobo,box_upbo),"TCL_set_V_box"))
			return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return mcount;
	}


rcode DTLAPI DTL_remove_V_mbox(int crit) {
	int i,j,k;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("RVMB");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_remove_V_mbox(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	df = uf->df;
	/* Remove midpoint box */
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			mbox_lobo[k] = -1.0;
			mbox_upbo[k] = -1.0;
			}
	if (call(TCL_set_V_mbox(df,mbox_lobo,mbox_upbo),"TCL_set_V_mbox"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*********************************************************
  *
  *  Get value base contents
  *
  *********************************************************/

 /*
  * Call semantics: The hull and midpoint are returned in
  * triples [lohull,midpoint,uphull] for each consequence
  * indexed by its alternative and consequence numbers.
  */

rcode DTLAPI DTL_get_V_hull(int crit, h_matrix lobo, h_matrix mid, h_matrix upbo) {
	int i,j,k;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("GVH");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_V_hull(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(lobo,1);
	_certify_ptr(mid,2);
	_certify_ptr(upbo,3);
	_dtl_assert(lobo!=mid,1);
	_dtl_assert(mid!=upbo,2);
	_dtl_assert(lobo!=upbo,3);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	df = uf->df;
	/* Get hull and midpoint */
	if (call(TCL_get_V_hull(df,vhlobo,vhupbo),"TCL_get_V_hull"))
		return dtl_kernel_error();
	if (call(TCL_get_V_masspoint(df,V_mid),"TCL_get_V_masspoint"))
		return dtl_kernel_error();
	/* Transfer and log */
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			lobo[i][j] = vhlobo[k];
			mid[i][j] = V_mid[k];
			upbo[i][j] = vhupbo[k];
			if (cst_ext && (lobo[i][j]!=-1.0)) {
				sprintf(msg,"    V%d.%d.%-2d [%.3lf %.3lf %.3lf] (%le)\n",crit,i,j,
						lobo[i][j],mid[i][j],upbo[i][j],upbo[i][j]-lobo[i][j]);
				cst_log(msg);
				}
			}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


rcode DTLAPI DTL_get_V_modal(int crit, h_matrix modal) {
	int i,j,k;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("GVM");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_V_modal(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(modal,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	df = uf->df;
	/* Get hull and midpoint */
	if (call(TCL_get_V_hull(df,vhlobo,vhupbo),"TCL_get_V_hull"))
		return dtl_kernel_error();
	if (call(TCL_get_V_masspoint(df,V_mid),"TCL_get_V_masspoint"))
		return dtl_kernel_error();
	/* Get modals and transfer */
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			if (vhlobo[k] > -1.0) {
				modal[i][j] = dtl_mid_to_modal(vhlobo[k],V_mid[k],vhupbo[k],0);
				if (cst_ext) {
					sprintf(msg,"    V%d.%d.%d = %.3lf\n",crit,i,j,modal[i][j]);
					cst_log(msg);
					}
				}
			else
				modal[i][j] = -1.0;
			}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


/* DTL_check_V_modality signals whether the value base is modal or unmodal.
 * An unmodal V-base cannot guarantee the correctness of belief functions.
 * While P and W belief is generated from means, V is generated from modals.
 * This is a very important modelling assumption that has to be heeded.
 *
 * If a list of unmodal variables is desired, use DTL_get_V_modal instead
 * and look for -3.0 entries in the output matrix [not yet implemented]
 * or check the output written to the cst_ext log if it has been enabled.
 * Alternatively, get a full mapping from DTL_get_V_modality_matrix below.
 *
 * For MC, use crit=0 to check all the value bases in one function call.
 * The call works since the unmodality lies entirely in the value bases.
 *
 * If Ai=0       -> check all alternatives (for GAMMA and DIGAMMA >2 alts)
 * If Ai>0, Aj=0 -> check one alternative  (for PSI   and DIGAMMA 1 alt)
 * If Ai>0, Aj>0 -> check two alternatives (for DELTA and DIGAMMA 2 alts) */

rcode dtl_check_V_modality(int crit, int Ai, int Aj) {
	int i,j,k,unmodal=0;
	struct d_frame *df;

	/* Check input parameters */
	if (load_df1(crit))
		return DTL_CRIT_UNKNOWN;
	df = uf->df;
	if ((Ai < 0) || (Ai > df->n_alts))
		return DTL_ALT_UNKNOWN;
	if ((Aj < 0) || (Aj > df->n_alts))
		return DTL_ALT_UNKNOWN;
	/* Collect hull & mid for all values */
	if (call(TCL_get_V_hull(df,vhlobo,vhupbo),"TCL_get_V_hull"))
		return dtl_kernel_error();
	if (call(TCL_get_V_masspoint(df,V_mid),"TCL_get_V_masspoint"))
		return dtl_kernel_error();
	/* Check each modal value against its own hull */
	for (i=1, k=1; i<=df->n_alts; i++)
		if (!Ai || (i==Ai) || (i==Aj)) {
			for (j=1; j<=df->tot_cons[i]; j++, k++)
				if (dtl_mid_to_modal(vhlobo[k],V_mid[k],vhupbo[k],1) == -3.0) {
					unmodal++;
					if (cst_ext) {
						sprintf(msg,"    V%d.%d.%d = [%.3lf %.3lf %.3lf]\n",crit,i,j,vhlobo[k],V_mid[k],vhupbo[k]);
						cst_log(msg);
						}
					}
			}
		else // step 1-dim counter
			k += df->tot_cons[i];
	/* Log result */
	if (cst_ext) {
		if (!Ai)
			Aj = 0; // neat printing
		if (unmodal)
			sprintf(msg,"    K%d(%d,%d): %d unmodalit%s found\n",crit,Ai,Aj,unmodal,unmodal>1?"ies":"y");
		else
			sprintf(msg,"    K%d(%d,%d): complete modality found\n",crit,Ai,Aj);
		cst_log(msg);
		}
	return unmodal;
	}


rcode DTLAPI DTL_check_V_modality(int crit, int Ai, int Aj) {
	rcode rc;
	int k,start,stop,unmodal=0;

	/* Begin single thread semaphore */
	_smx_begin("CVMOD");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_check_V_modality(%d,%d,%d)\n",crit,Ai,Aj);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Scan value base(s) - one or all */
	if (crit)
		start = stop = crit;
	else { // crit=0 means all
		start = 1;
		stop = uf->n_crit;
		}
	for (k=start; k<=stop; k++)
		if (rc = dtl_check_V_modality(k,Ai,Aj))
			if (!crit && (rc == DTL_CRIT_UNKNOWN));
				// shadow criterion -> harmless
			else if (rc < DTL_OK)
				return dtl_error(rc);
			else
				unmodal += crit?rc:1;
	/* End single thread semaphore */
	_smx_end();
	return unmodal;
	}


/* DTL_get_V_modality_matrix returns a matrix telling which parts of
 * the value base are modal or unmodal. In the matrix, 0 means unmodal
 * and 1 means modal. An unmodal base cannot guarantee the correctness
 * of belief functions since the distributions extend beyond the hulls.
 *
 * The matrix is indexed with alternative numbers and composed as follows:
 * Row 0 is GAMMA
 * Col 0 is PSI (as are entires [i][i])
 *  (entry [0][0] is the entire V-base)
 * Entries [i][j] are DELTA if i!=j and i,j>0
 *  (entries [i][i] duplicate col 0)
 * Either above is DIGAMMA, depending on 2nd param */

rcode DTLAPI DTL_get_V_modality_matrix(int crit, ai_matrix modal_mx) {
	rcode rc;
	int i,j,k,start,stop;

	/* Begin single thread semaphore */
	_smx_begin("VMODMX");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_V_modality_matrix(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Init matrix - presuming modality */
	for (i=0; i<=uf->n_alts; i++)
		for (j=0; j<=uf->n_alts; j++)
			modal_mx[i][j] = 1;
	/* Scan value base(s) - one or all */
	if (crit)
		start = stop = crit;
	else { // crit=0 means all
		start = 1;
		stop = uf->n_crit;
		}
	for (k=start; k<=stop; k++)
		for (i=1; i<=uf->n_alts; i++) {
			if (rc = dtl_check_V_modality(k,i,0))
				if (!crit && (rc == DTL_CRIT_UNKNOWN))
					continue; // shadow criterion
				else if (rc < DTL_OK)
					return dtl_error(rc);
				else // unmodal
					for (j=0; j<=uf->n_alts; j++) {
						modal_mx[0][j] = 0; // reset base and GAMMA
						modal_mx[i][j] = 0; // reset original DELTA
						modal_mx[j][i] = 0; // reset mirrored DELTA
						}
			}
	/* Log matrix results per alternative if unmodal */
	if (cst_ext && !modal_mx[0][0]) {
		cst_log(" Alt T");
		for (i=1; i<=uf->n_alts; i++) {
			sprintf(msg," %1d",i%10);
			cst_log(msg);
			}
		cst_log("\n ---");
		for (i=0; i<=uf->n_alts; i++)
			cst_log("--");
		cst_log("\n");
		for (i=0; i<=uf->n_alts; i++) {
			if (i)
				sprintf(msg," A%-2d",i);
			else
				sprintf(msg," Tot");
			cst_log(msg);
			for (j=0; j<=uf->n_alts; j++) {
				sprintf(msg," %d",modal_mx[i][j]);
				cst_log(msg);
				}
			cst_log("\n");
			}
		}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


rcode DTLAPI DTL_reset_V_base(int crit) {

	/* Begin single thread semaphore */
	_smx_begin("RSTV");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_reset_V_base(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	/* Clear base */
	if (call(TCL_reset_V_base(uf->df),"TCL_reset_V_base"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*********************************************************
  *
  *  For debug use, limited error handling
  *
  *********************************************************/

/* Check if the tree node is a real/final one or not (otherwise intermediate).
 * Returns >0 if real tree node, =0 if im-node or don't know due to error. */

int DTLAPI DTI_real_V_node(int crit, int alt, int node) {

	/* Check if function can start */
	if (!frame_loaded)
		return 0;
	if (load_df1(crit))
		return 0;
	return TCL_get_V_index(uf->df,alt,node);
	}


 /*********************************************************
  *
  *  For internal use, limited error handling
  *
  *********************************************************/

/* Check if there are any active midpoints in the probability base. Returns >0 if
 * there are midpoints, =0 if no active midpoints, <0 if don't know due to error. */

int dtl_nbr_V_midpoints(int crit) {
	int i,j,k,n_midpoints;
	struct d_frame *df;

	/* Check if function can start */
	if (!frame_loaded)
		return -1;
	if (load_df1(crit))
		return -1;
	df = uf->df;
	if (TCL_get_V_mbox(df,mbox_lobo,mbox_upbo))
		return -1;
	/* Calculate number of midpoints */
	n_midpoints = 0;
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++)
			if (mbox_lobo[k] > -1.0)
				n_midpoints++;
	return n_midpoints;
	}


#ifdef CONSOLE

 /*********************************************************
  *
  *  Console output for test and debug
  *  Internal debug use, no error handling
  *
  *********************************************************/

static void print_V_stmt(int crit, struct base *V, int snbr) {

	if (V->stmt[snbr].n_terms == 1)
		printf("%2d: V%d.%d.%-2d [%.3lf %.3lf]\n",
				snbr,crit,V->stmt[snbr].alt[1],V->stmt[snbr].cons[1],
				V->stmt[snbr].lobo,V->stmt[snbr].upbo);
	else
		printf("%2d: error: %d terms in V statement\n",
					snbr,V->stmt[snbr].n_terms);
	}


void DTLAPI DTI_show_V_base(int crit) {
	int i;
	struct base *V;

	/* Check if function can start */
	if (!frame_loaded)
		return;
	/* Check input parameters */
	if (load_df1(crit))
		return;
	/* Display base */
	V = uf->df->V_base;
	printf("The value base contains %d constraint%s\n",V->n_stmts,(V->n_stmts==1?"":"s"));
	for (i=1; i<=V->n_stmts; i++)
		print_V_stmt(crit,V,i);
	}


void DTLAPI DTI_show_V_box(int crit) {
	int i,j,k,h,n_entries;
	struct d_frame *df;

	/* Check if function can start */
	if (!frame_loaded)
		return;
	/* Check input parameter */
	if (load_df0(crit))
		return;
	df = uf->df;
	if (!df->V_base->box) {
		printf("The value box contains no entries\n");
		return;
		}
	if (TCL_get_V_box(df,vhlobo,vhupbo))
		return;
	/* Calculate number of ranges */
	n_entries = 0;
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++)
			if (TCL_get_V_index(df,i,j)) // re-node
				if ((vhlobo[k] > 0.0) || (vhupbo[k] < 1.0))
					n_entries++;
	printf("The value box contains %d entr%s\n",n_entries,(n_entries==1?"y":"ies"));
	/* Print ranges in order */
	for (k=1,i=1,h=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			if (TCL_get_V_index(df,i,j)) // re-node
				if ((vhlobo[k] > 0.0) || (vhupbo[k] < 1.0)) {
					if (vhupbo[k]-vhlobo[k] > DTL_EPS)
						printf("%2d: V%d.%-2d= [%.3lf %.3lf]\n",h++,i,j,vhlobo[k],vhupbo[k]);
					else
						printf("%2d: V%d.%-2d= %.3lf\n",h++,i,j,vhlobo[k]);
					}
			}
	}


void DTLAPI DTI_show_V_mbox(int crit) {
	int i,j,k,h,n_midpoints;
	struct d_frame *df;

	/* Check if function can start */
	if (!frame_loaded)
		return;
	/* Check input parameter */
	if (load_df1(crit))
		return;
	df = uf->df;
	if (TCL_get_V_mbox(df,mbox_lobo,mbox_upbo))
		return;
	/* Calculate number of midpoints */
	n_midpoints = 0;
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++)
			if (mbox_lobo[k] > -1.0)
				n_midpoints++;
	printf("The value mbox contains %d midpoint%s\n",n_midpoints,(n_midpoints==1?"":"s"));
	/* Print midpoints in order */
	for (k=1,i=1,h=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			if (mbox_lobo[k] > -1.0) {
				if (mbox_upbo[k]-mbox_lobo[k] > DTL_EPS)
					printf("%2d: V%d.%d.%-2d [%.3lf %.3lf]  <- INTERVAL (%.3le)\n",h++,crit,i,j,
							mbox_lobo[k],mbox_upbo[k],mbox_upbo[k]-mbox_lobo[k]);
				else
					printf("%2d: V%d.%d.%-2d %.3lf\n",h++,crit,i,j,mbox_lobo[k]);
				}
			}
	}

#endif

/* Autoscale add-in */
#include "DTLautoscale.c"
