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
 *   File: DTLpbase.c
 *
 *   Purpose: probability base interface to TCL
 *
 *
 *   Functions exported outside DTL
 *   ------------------------------
 *   DTL_add_P_statement
 *   DTL_change_P_statement
 *   DTL_replace_P_statement
 *   DTL_delete_P_statement
 *   DTL_add_P_mid_statement
 *   DTL_delete_P_mid_statement
 *   DTL_set_P_box
 *   DTL_set_P_mbox/1
 *   DTL_remove_P_mbox
 *   DTL_get_P_hull
 *   DTL_reset_P_base
 *
 *   Functions outside module, debug use
 *   -----------------------------------
 *   DTI_P_node_parents
 *   DTI_show_P_base
 *   DTI_show_P_box
 *   DTI_show_P_mbox
 *
 *   Functions outside module, inside DTL
 *   ------------------------------------
 *   dtl_set_P_mbox_auto
 *   dtl_P_node_parents
 *   dtl_P_nbr_of_siblings
 *   dtl_nbr_P_midpoints
 *
 *   Functions internal to module
 *   ----------------------------
 *   print_P_stmt
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
static d_row phlobo,phupbo,pllobo,plupbo;
static d_row P_mid,LP_mid;


 /*********************************************************
  *
  *  Probability base
  *
  *********************************************************/

 /*
  * Call semantics: Add the user constraint p(crit:alt:cons) =
  * [lobo,upbo] to the probability base within the decision frame.
  * The base is checked for consistency with respect to the new
  * interval. In case of inconsistency, nothing is added to the base.
  */

rcode DTLAPI DTL_add_P_statement(int crit, struct user_stmt_rec* ustmtp) {
	struct stmt_rec stmt;

	/* Begin single thread semaphore */
	_smx_begin("APS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_add_P_statement(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(ustmtp,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if (load_PV_stmt(crit,ustmtp,&stmt,'P'))
		return dtl_error(DTL_STMT_ERROR);
	/* Add statement */
	if (call(TCL_add_P_constraint(uf->df,&stmt),"TCL_add_P_constraint"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return uf->df->P_base->n_stmts;
	}


 /*
  * Call semantics: Change the existing user probability constraint
  * p(crit:alt:cons) = [old_lobo,old_upbo]
  * to p(crit:alt:cons) = [lobo,upbo] in the probability base.
  * The base is checked for consistency with respect to the change.
  * In case of inconsistency, nothing is changed in the base.
  */

rcode DTLAPI DTL_change_P_statement(int crit, int stmt_number, double lobo, double upbo) {

	/* Begin single thread semaphore */
	_smx_begin("CPS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_change_P_statement(%d,%d,%.3lf,%.3lf)\n",
				crit,stmt_number,lobo,upbo);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	/* Change statement */
	if (call(TCL_change_P_constraint(uf->df,stmt_number,lobo,upbo),"TCL_change_P_constraint"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*
  * Call semantics: Replace the user probability constraint
  * with p(crit:alt:cons) = [lobo,upbo] in the probability base.
  * The base is checked for consistency with respect to the new
  * interval. In case of inconsistency, nothing is replaced in the
  * base.
  */

rcode DTLAPI DTL_replace_P_statement(int crit, int stmt_number, struct user_stmt_rec* ustmtp) {
	struct stmt_rec stmt;

	/* Begin single thread semaphore */
	_smx_begin("RPS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_replace_P_statement(%d,%d)\n",crit,stmt_number);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(ustmtp,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if (load_PV_stmt(crit,ustmtp,&stmt,'P'))
		return dtl_error(DTL_STMT_ERROR);
	/* Replace statement */
	if (call(TCL_replace_P_constraint(uf->df,stmt_number,&stmt),"TCL_replace_P_constraint"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*
  * Call semantics: The user probability constraint (interval) with
  * position stmt_nbr in the probability base is deleted from the
  * base. All constraints with higher positions within the base are
  * shifted one position down.
  */

rcode DTLAPI DTL_delete_P_statement(int crit, int stmt_number) {

	/* Begin single thread semaphore */
	_smx_begin("DPS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_delete_P_statement(%d,%d)\n",crit,stmt_number);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	/* Delete statement */
	if (call(TCL_delete_P_constraint(uf->df,stmt_number),"TCL_delete_P_constraint"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return uf->df->P_base->n_stmts;
	}


 /*
  * Call semantics: Enter the user probability midpoint
  * p(crit:alt:cons) = [lobo,upbo] into the decision frame.
  * The base is checked for consistency with respect to the
  * new midpoint. In case of inconsistency, nothing is added
  * to the base. NOTE: midpoint is mean, not mode (cf. BRS).
  */

rcode DTLAPI DTL_add_P_mid_statement(int crit, struct user_stmt_rec* ustmtp) {
	struct stmt_rec stmt;

	/* Begin single thread semaphore */
	_smx_begin("APM");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_add_P_mid_statement(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(ustmtp,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if (ustmtp->n_terms != 1)
		return dtl_error(DTL_WRONG_STMT_TYPE);
	if (load_PV_stmt(crit,ustmtp,&stmt,'P'))
		return dtl_error(DTL_STMT_ERROR);
	/* Set midpoint */
	if (call(TCL_add_P_mstatement(uf->df,&stmt),"TCL_add_P_mstatement"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*
  * Call semantics: Remove the user probability core midpoint
  * for p(crit:alt:cons) from the decision frame.
  */

rcode DTLAPI DTL_delete_P_mid_statement(int crit, struct user_stmt_rec* ustmtp) {
	struct stmt_rec stmt;

	/* Begin single thread semaphore */
	_smx_begin("DPM");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_delete_P_mid_statement(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(ustmtp,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if (ustmtp->n_terms != 1)
		return dtl_error(DTL_WRONG_STMT_TYPE);
	if (load_PV_stmt(crit,ustmtp,&stmt,'P'))
		return dtl_error(DTL_STMT_ERROR);
	/* Remove midpoint */
	if (call(TCL_delete_P_mstatement(uf->df,&stmt),"TCL_delete_P_mstatement"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


static d_row box_lobo,box_upbo;

rcode DTLAPI DTL_set_P_box(int crit, h_matrix lobox, h_matrix upbox) {
	int i,j,k;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("SPB");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_set_P_box(%d)\n",crit);
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
	/* Assemble and set box intervals */
	df = uf->df;
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			box_lobo[k] = lobox[i][j];
			box_upbo[k] = upbox[i][j];
			if (cst_on && ((lobox[i][j]!=0.0) || (upbox[i][j]!=1.0))) {
				sprintf(msg,"    P%d.%d.%-2d [%.3lf %.3lf] (%le)\n",crit,i,j,
						lobox[i][j],upbox[i][j],upbox[i][j]-lobox[i][j]);
				cst_log(msg);
				}
			}
	if (call(TCL_set_P_box(df,box_lobo,box_upbo),"TCL_set_P_box"))
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

// NOTE: midpoint is mean, not mode (cf. BRS)
static d_row mbox_lobo,mbox_upbo;

rcode dtl_set_P_check(int crit, h_matrix lobox, h_matrix mbox, h_matrix upbox) {
	int i,j,k;
	struct d_frame *df;

	/* Check if function can start */
	_certify_ptr(lobox,1);
	_certify_ptr(mbox,2);
	_certify_ptr(upbox,3);
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameters */
	if (load_df1(crit))
		return DTL_CRIT_UNKNOWN;
	df = uf->df;
	/* Check box intervals */
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++)
			if (mbox[i][j] < 0.0) { // mbox entry is empty
				if (lobox[i][j] > upbox[i][j])
					return k;
				}
			else { // mbox entry should also be checked
				if ((lobox[i][j] > mbox[i][j]) || (mbox[i][j] > upbox[i][j]))
					return k;
				}
	return DTL_OK;
	}


rcode DTLAPI DTL_set_P_mbox(int crit, h_matrix lobox, h_matrix upbox) {
	int i,j,k;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("SPMB");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_set_P_mbox(%d)\n",crit);
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
	/* Assemble and set mbox midpoints */
	df = uf->df;
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			mbox_lobo[k] = lobox[i][j];
			mbox_upbo[k] = upbox[i][j];
			if (cst_on && ((lobox[i][j]!=-2.0) || (upbox[i][j]!=-2.0))) {
				sprintf(msg,"    P%d.%d.%-2d [%.3lf %.3lf] (%le)\n",crit,i,j,
						lobox[i][j],upbox[i][j],upbox[i][j]-lobox[i][j]);
				cst_log(msg);
				}
			}
	if (call(TCL_set_P_mbox(df,mbox_lobo,mbox_upbo),"TCL_set_P_mbox"))
		return dtl_kernel_error();
	uf->WP_autogen[crit] = FALSE;
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


/* DTL_set_P_mbox1: pointwise mid, not wide */

rcode DTLAPI DTL_set_P_mbox1(int crit, h_matrix mbox) {

	return DTL_set_P_mbox(crit,mbox,mbox);
	}


/* dtl_set_P_mbox_auto: mark as autogenerated */

rcode dtl_set_P_mbox_auto(int crit, h_matrix lobox, h_matrix upbox) {
	rcode rc;

	rc = DTL_set_P_mbox(crit,lobox,upbox);
	if (!rc)
		uf->WP_autogen[crit] = TRUE;
	return rc;
	}


rcode DTLAPI DTL_remove_P_mbox(int crit) {
	int i,j,k;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("RPMB");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_remove_P_mbox(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	/* Remove box midpoint */
	df = uf->df;
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			mbox_lobo[k] = -1.0;
			mbox_upbo[k] = -1.0;
			}
	if (call(TCL_set_P_mbox(df,mbox_lobo,mbox_upbo),"TCL_set_P_mbox"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*********************************************************
  *
  *  Get probability base contents
  *
  *********************************************************/

 /*
  * Call semantics: The hull and mass point are returned
  * for each consequence in three h_matrix structures
  * indexed by alternative and consequence.
  */

rcode DTLAPI DTL_get_P_hull(int crit, int global, h_matrix lobo, h_matrix mid, h_matrix upbo) {
	int i,j,k;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("GPH");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_P_hull(%d,%d)\n",crit,global);
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
	/* Get hull and midpoint */
	df = uf->df;
	if (call(TCL_get_P_hull(df,phlobo,phupbo,pllobo,plupbo),"TCL_get_P_hull"))
		return dtl_kernel_error();
	if (call(TCL_get_P_masspoint(df,P_mid,LP_mid),"TCL_get_P_masspoint"))
		return dtl_kernel_error();
	/* Transfer hull and midpoint */
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			lobo[i][j] = global?phlobo[k]:pllobo[k];
			mid[i][j]  = global?P_mid[k]: LP_mid[k];
			upbo[i][j] = global?phupbo[k]:plupbo[k];
			if (cst_ext) {
				sprintf(msg,"    P%d.%d.%-2d [%.3lf %.3lf %.3lf] (%le)\n",crit,i,j,
						lobo[i][j],mid[i][j],upbo[i][j],upbo[i][j]-lobo[i][j]);
				cst_log(msg);
				}
			}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


rcode DTLAPI DTL_reset_P_base(int crit) {

	/* Begin single thread semaphore */
	_smx_begin("RSTP");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_reset_P_base(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df1(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	/* Reset information */
	if (call(TCL_reset_P_base(uf->df),"TCL_reset_P_base"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*********************************************************
  *
  *  For internal use, limited error handling
  *
  *********************************************************/

int dtl_P_node_parents(int crit, int alt, int node1, int node2) {
	struct d_frame *df;

	if (!frame_loaded)
		return -1;
	if (load_df1(crit))
		return -1;
	df = uf->df;
	if ((alt<1) || (alt>df->n_alts))
		return -1;
	return -1*TCL_different_parents(df,alt,node1,node2);
	}


int DTLAPI DTI_P_node_parents(int crit, int alt, int node1, int node2) {

	/* Call stack refurbishing */
	return dtl_P_node_parents(crit,alt,node1,node2);
	}


int dtl_P_nbr_of_siblings(int crit, int alt, int node) {
	struct d_frame *df;

	if (!frame_loaded)
		return -1;
	if (load_df1(crit))
		return -1;
	df = uf->df;
	if ((alt<1) || (alt>df->n_alts))
		return -1;
	return -1*TCL_nbr_of_siblings(df,alt,node);
	}


/* Check if there are any active midpoints in the probability base. Returns >0 if
 * there are midpoints, =0 if no active midpoints, <0 if don't know due to error. */

int dtl_nbr_P_midpoints(int crit) {
	int i,j,k,n_midpoints;
	struct d_frame *df;

	/* Check if function can start */
	if (!frame_loaded)
		return -1;
	if (load_df1(crit))
		return -1;
	df = uf->df;
	if (TCL_get_P_mbox(df,mbox_lobo,mbox_upbo))
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
  *
  *********************************************************/

static void print_P_stmt(int crit, struct base *P, int snbr) {

	if (P->stmt[snbr].n_terms == 1)
		printf("%2d: P%d.%d.%-2d [%.3lf %.3lf]\n",
				snbr,crit,P->stmt[snbr].alt[1],P->stmt[snbr].cons[1],
				P->stmt[snbr].lobo,P->stmt[snbr].upbo);
	else
		printf("%2d: error: %d terms in P statement\n",
				snbr,P->stmt[snbr].n_terms);
	}


void DTLAPI DTI_show_P_base(int crit) {
	int i;
	struct base *P;

	/* Check if function can start */
	if (!frame_loaded)
		return;
	/* Check input parameter */
	if (load_df1(crit))
		return;
	/* Display base */
	P = uf->df->P_base;
	printf("The probability base contains %d constraint%s\n",P->n_stmts,(P->n_stmts==1?"":"s"));
	for (i=1; i<=P->n_stmts; i++)
		print_P_stmt(crit,P,i);
	}


void DTLAPI DTI_show_P_box(int crit) {
	int i,j,k,h,n_entries;
	struct d_frame *df;

	/* Check if function can start */
	if (!frame_loaded)
		return;
	/* Check input parameter */
	if (load_df1(crit))
		return;
	df = uf->df;
	if (!df->P_base->box) {
		printf("The probability box contains no entries\n");
		return;
		}
	if (TCL_get_P_box(df,pllobo,plupbo))
		return;
	/* Calculate number of estimates */
	n_entries = 0;
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++)
			if ((pllobo[k] > 0.0) || (plupbo[k] < 1.0))
				n_entries++;
	printf("The probability box contains %d entr%s\n",n_entries,(n_entries==1?"y":"ies"));
	/* Print estimates in order */
	for (k=1,i=1,h=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			if ((pllobo[k] > 0.0) || (plupbo[k] < 1.0)) {
				if (plupbo[k]-pllobo[k] > DTL_EPS)
					printf("%2d: P%d.%-2d= [%.3lf %.3lf]\n",h++,i,j,pllobo[k],plupbo[k]);
				else
					printf("%2d: P%d.%-2d= %.3lf\n",h++,i,j,pllobo[k]);
				}
			}
	}


void DTLAPI DTI_show_P_mbox(int crit) {
	int i,j,k,h,n_midpoints;
	struct d_frame *df;

	/* Check if function can start */
	if (!frame_loaded)
		return;
	/* Check input parameter */
	if (load_df1(crit))
		return;
	df = uf->df;
	if (TCL_get_P_mbox(df,mbox_lobo,mbox_upbo))
		return;
	/* Calculate number of midpoints */
	n_midpoints = 0;
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++)
			if (mbox_lobo[k] > -1.0)
				n_midpoints++;
	printf("The probability mbox contains %d midpoint%s\n",n_midpoints,(n_midpoints==1?"":"s"));
	/* Print midpoints in order */
	for (k=1,i=1,h=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			if (mbox_lobo[k] > -1.0) {
				if (mbox_upbo[k]-mbox_lobo[k] > DTL_EPS)
					printf("%2d: P%d.%d.%-2d [%.3lf %.3lf]  <- INTERVAL (%.3le)\n",h++,crit,i,j,
							mbox_lobo[k],mbox_upbo[k],mbox_upbo[k]-mbox_lobo[k]);
				else
					printf("%2d: P%d.%d.%-2d %.3lf\n",h++,crit,i,j,mbox_lobo[k]);
				}
			}
	}

#endif
