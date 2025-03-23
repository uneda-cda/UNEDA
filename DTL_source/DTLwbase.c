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
 *   File: DTLwbase.c
 *
 *   Purpose: weight base interface to TCL
 *
 *
 *   Functions exported outside DTL
 *   ------------------------------
 *   DTL_add_W_statement
 *   DTL_change_W_statement
 *   DTL_replace_W_statement
 *   DTL_delete_W_statement
 *   DTL_add_W_mid_statement
 *   DTL_delete_W_mid_statement
 *   DTL_set_W_box
 *   DTL_set_W_mbox/1
 *   DTL_remove_W_mbox
 *   DTL_get_W_hull
 *   DTL_reset_W_base
 *
 *   Functions outside module, debug use
 *   -----------------------------------
 *   DTI_pure_W_tree
 *   DTI_W_node_parents
 *   DTI_real_W_crit
 *   DTI_nbr_W_midpoints
 *   DTI_show_W_base
 *   DTI_show_W_box
 *   DTI_show_W_mbox
 *
 *   Functions outside module, inside DTL
 *   ------------------------------------
 *   dtl_set_W_mbox_auto
 *   dtl_pure_W_tree
 *   dtl_W_node_parents
 *   dtl_W_nbr_of_siblings
 *   dtl_real_W_crit
 *   dtl_nbr_W_midpoints
 *
 *   Functions internal to module
 *   ----------------------------
 *   print_W_stmt
 *
 */

#include "DTL.h"
#include "DTLinternal.h"


 /********************************************************
  *
  *  Internal data
  *
  ********************************************************/

static d_row box_lobo,box_upbo;
static d_row mbox_lobo,mbox_upbo;
static d_row hlobo,hupbo,llobo,lupbo,W_mid,LW_mid;


 /*********************************************************
  *
  *  Weight base
  *
  *********************************************************/

 /*
  * Call semantics: Add the user weight statement w(crit) =
  * [lobo,upbo] to the weight base within the decision frame. The
  * base is checked for consistency with respect to the new interval.
  * In case of inconsistency, nothing is added to the base.
  */

rcode DTLAPI DTL_add_W_statement(struct user_w_stmt_rec* uwstmtp) {
	struct stmt_rec stmt;

	/* Begin single thread semaphore */
	_smx_begin("AWS");
	/* Log function call */
	if (cst_on)
		cst_log("DTL_add_W_statement()\n");
	/* Check if function can start */
	_certify_ptr(uwstmtp,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (PS)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if (load_df0(0))
		return dtl_error(DTL_SYS_CORRUPT);
	if (load_W_stmt(uwstmtp,&stmt))
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
  * Call semantics: Change the existing user weight statement w(crit) =
  * [old_lobo,old_upbo] to w(crit) = [lobo,upbo] in the weight base.
  * The base is checked for consistency with respect to the change.
  * In case of inconsistency, nothing is changed in the base.
  */

rcode DTLAPI DTL_change_W_statement(int stmt_number, double lobo, double upbo) {

	/* Begin single thread semaphore */
	_smx_begin("CWS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_change_W_statement(%d,%.3lf,%.3lf)\n",stmt_number,lobo,upbo);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (PS)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if (load_df0(0))
		return dtl_error(DTL_SYS_CORRUPT);
	/* Change statement */
	if (call(TCL_change_P_constraint(uf->df,stmt_number,lobo,upbo),"TCL_change_P_constraint"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*
  * Call semantics: Replace the user weight statement with
  * w(crit) = [lobo,upbo] in the weight base. The base is
  * checked for consistency with respect to the new interval.
  * In case of inconsistency, nothing is replaced in the base.
  */

rcode DTLAPI DTL_replace_W_statement(int stmt_number, struct user_w_stmt_rec* uwstmtp) {
	struct stmt_rec stmt;

	/* Begin single thread semaphore */
	_smx_begin("RWS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_replace_W_statement(%d)\n",stmt_number);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(uwstmtp,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (PS)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if (load_df0(0))
		return dtl_error(DTL_SYS_CORRUPT);
	if (load_W_stmt(uwstmtp,&stmt))
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
  * Call semantics: The user weight constraint with position
  * wstmt_nbr in the weight base is deleted from the base.
  * All constraints with higher positions within the base are
  * shifted one position down.
  */

rcode DTLAPI DTL_delete_W_statement(int stmt_number) {

	/* Begin single thread semaphore */
	_smx_begin("DWS");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_delete_W_statement(%d)\n",stmt_number);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (PS)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if (load_df0(0))
		return dtl_error(DTL_SYS_CORRUPT);
	/* Delete statement */
	if (call(TCL_delete_P_constraint(uf->df,stmt_number),"TCL_delete_P_constraint"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return uf->df->P_base->n_stmts;
	}


 /*
  * Call semantics: Enter the user weight midpoint
  * w(crit) = [lobo,upbo] into the decision frame.
  * The base is checked for consistency with respect
  * to the midpoint. In case of inconsistency,
  * nothing is added to the base.
  */

rcode DTLAPI DTL_add_W_mid_statement(struct user_w_stmt_rec* uwstmtp) {
	struct stmt_rec stmt;

	/* Begin single thread semaphore */
	_smx_begin("AWM");
	/* Log function call */
	if (cst_on)
		cst_log("DTL_add_W_mid_statement()\n");
	/* Check if function can start */
	_certify_ptr(uwstmtp,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (PS)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if (load_df0(0))
		return dtl_error(DTL_SYS_CORRUPT);
	if (uwstmtp->n_terms != 1)
		return dtl_error(DTL_WRONG_STMT_TYPE);
	if (load_W_stmt(uwstmtp,&stmt))
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
  * Call semantics: Remove the user weight core midpoint for
  * w(crit) from the decision frame.
  */

rcode DTLAPI DTL_delete_W_mid_statement(struct user_w_stmt_rec* uwstmtp) {
	struct stmt_rec stmt;

	/* Begin single thread semaphore */
	_smx_begin("DWM");
	/* Log function call */
	if (cst_on)
		cst_log("DTL_delete_W_mid_statement()\n");
	/* Check if function can start */
	_certify_ptr(uwstmtp,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (PS)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if (load_df0(0))
		return dtl_error(DTL_SYS_CORRUPT);
	if (uwstmtp->n_terms != 1)
		return dtl_error(DTL_WRONG_STMT_TYPE);
	if (load_W_stmt(uwstmtp,&stmt))
		return dtl_error(DTL_STMT_ERROR);
	/* Remove midpoint */
	if (call(TCL_delete_P_mstatement(uf->df,&stmt),"TCL_delete_P_mstatement"))
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

rcode dtl_set_W_check(h_vector lobox, h_vector mbox, h_vector upbox) {
	int j;
	struct d_frame *df;

	/* Check if function can start */
	_certify_ptr(lobox,1);
	_certify_ptr(mbox,2);
	_certify_ptr(upbox,3);
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	if (PS)
		return DTL_WRONG_FRAME_TYPE;
	/* Check input parameters */
	if (load_df0(0))
		return DTL_SYS_CORRUPT;
	df = uf->df;
	/* Alternative 1: the weight container */
	for (j=1; j<=df->tot_cons[1]; j++)
		if (mbox[j] < 0.0) { // mbox entry is empty
			if (lobox[j] > upbox[j])
				return j;
			}
		else { // mbox entry should also be checked
			if ((lobox[j] > mbox[j]) || (mbox[j] > upbox[j]))
				return j;
			}
	return DTL_OK;
	}


rcode DTLAPI DTL_set_W_box(h_vector lobox, h_vector upbox) {
	int j,j2;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("SWB");
	/* Log function call */
	if (cst_on)
		cst_log("DTL_set_W_box()\n");
	/* Check if function can start */
	_certify_ptr(lobox,1);
	_certify_ptr(upbox,2);
	_dtl_assert(lobox!=upbox,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (PS)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if (load_df0(0))
		return dtl_error(DTL_SYS_CORRUPT);
	df = uf->df;
	/* Alternative 1: the weight container */
	for (j=1; j<=df->tot_cons[1]; j++) {
		box_lobo[j] = lobox[j];
		box_upbo[j] = upbox[j];
		if (cst_on && ((lobox[j]!=0.0) || (upbox[j]!=1.0))) {
			sprintf(msg,"    W%-2d [%.3lf %.3lf] (%le)\n",j,lobox[j],upbox[j],upbox[j]-lobox[j]);
			cst_log(msg);
			}
		}
	/* Alternatives 2..n: dummies */
	for (j2=2; j2<=df->n_alts; j++,j2++) {
		box_lobo[j] = 0.0;
		box_upbo[j] = 1.0;
		}
	if (call(TCL_set_P_box(df,box_lobo,box_upbo),"TCL_set_P_box"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


rcode DTLAPI DTL_set_W_mbox(h_vector lobox, h_vector upbox) {
	int j,j2;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("SWMB");
	/* Log function call */
	if (cst_on)
		cst_log("DTL_set_W_mbox()\n");
	/* Check if function can start */
	_certify_ptr(lobox,1);
	_certify_ptr(upbox,2);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (PS)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if (load_df0(0))
		return dtl_error(DTL_SYS_CORRUPT);
	df = uf->df;
	/* Assemble and set mbox midpoints */
	for (j=1; j<=df->tot_cons[1]; j++) {
		mbox_lobo[j] = lobox[j];
		mbox_upbo[j] = upbox[j];
		if (cst_on && ((lobox[j]!=-2.0) || (upbox[j]!=-2.0))) {
			sprintf(msg,"    W%-2d [%.3lf %.3lf] (%le)\n",j,lobox[j],upbox[j],upbox[j]-lobox[j]);
			cst_log(msg);
			}
		}
	/* Alternatives 2..n: dummies */
	for (j2=2; j2<=df->n_alts; j++,j2++) {
		mbox_lobo[j] = -1.0;
		mbox_upbo[j] = -1.0;
		}
	if (call(TCL_set_P_mbox(df,mbox_lobo,mbox_upbo),"TCL_set_P_mbox"))
		return dtl_kernel_error();
	uf->WP_autogen[0] = FALSE;
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


/* DTL_set_W_mbox1: pointwise mid, not wide */

rcode DTLAPI DTL_set_W_mbox1(h_vector mbox) {

	return DTL_set_W_mbox(mbox,mbox);
	}


/* dtl_set_W_mbox_auto: mark as autogenerated */

rcode dtl_set_W_mbox_auto(h_vector lobox, h_vector upbox) {
	rcode rc;

	rc = DTL_set_W_mbox(lobox,upbox);
	if (!rc)
		uf->WP_autogen[0] = TRUE;
	return rc;
	}


rcode DTLAPI DTL_remove_W_mbox() {
	int j,j2;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("RWMB");
	/* Log function call */
	if (cst_on)
		cst_log("DTL_remove_W_mbox()\n");
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (PS)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if (load_df0(0))
		return dtl_error(DTL_SYS_CORRUPT);
	df = uf->df;
	/* Remove box midpoint */
	for (j=1; j<=df->tot_cons[1]; j++) {
		mbox_lobo[j] = -1.0;
		mbox_upbo[j] = -1.0;
		}
	/* Alternatives 2..n: dummies */
	for (j2=2; j2<=df->n_alts; j++,j2++) {
		mbox_lobo[j] = -1.0;
		mbox_upbo[j] = -1.0;
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
  *  Get weight base contents
  *
  *********************************************************/

 /*
  * Call semantics: The hull and mass point are returned in
  * triples [lobo,mid,upbo] for each weight in the hull.
  */

rcode DTLAPI DTL_get_W_hull(int global, h_vector lobo, h_vector mid, h_vector upbo) {
	int j;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("GWH");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_W_hull(%d)\n",global);
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
	if (PS)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if (load_df0(0))
		return dtl_error(DTL_SYS_CORRUPT);
	df = uf->df;
	/* Get hull and midpoint */
	if (call(TCL_get_P_hull(df,hlobo,hupbo,llobo,lupbo),"TCL_get_P_hull"))
		return dtl_kernel_error();
	if (call(TCL_get_P_masspoint(df,W_mid,LW_mid),"TCL_get_P_masspoint"))
		return dtl_kernel_error();
	/* Transfer hull and midpoint */
	for (j=1; j<=df->tot_cons[1]; j++) {
		lobo[j] = global?hlobo[j]:llobo[j];
		mid[j]  = global?W_mid[j]:LW_mid[j];
		upbo[j] = global?hupbo[j]:lupbo[j];
		if (cst_ext) {
			sprintf(msg,"    W%-2d [%.3lf %.3lf %.3lf] (%le)\n",j,lobo[j],mid[j],upbo[j],upbo[j]-lobo[j]);
			cst_log(msg);
			}
		}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


rcode DTLAPI DTL_reset_W_base() {

	/* Begin single thread semaphore */
	_smx_begin("RSTW");
	/* Log function call */
	if (cst_on)
		cst_log("DTL_reset_W_base()\n");
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (PS)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if (load_df0(0))
		return dtl_error(DTL_SYS_CORRUPT);
	/* Reset information */
	if (call(TCL_reset_P_base(uf->df),"TCL_reset_P_base"))
		return dtl_kernel_error();
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


static ci_col xorder,xnode;
static cr_col xval;

/* Obtain a ranking of the children nodes to a flat weight structure, or to
 * the top node or an intermediate weight node in case of a weight tree */

DTL_get_W_node_ranking(int snode, ci_col Wnode, cr_col Wval) {
	struct d_frame *df;
	int i,tnode,rc;

	_smx_begin("GWNR");
	_certify_ptr(Wnode,1);
	_certify_ptr(Wval,2);
	if (cst_on) {
		sprintf(msg,"DTL_get_W_node_ranking(%d)\n",snode);
		cst_log(msg);
		}
	/* Check input parameters */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (PS)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	if (rc = load_df0(0))
		return dtl_error(rc);
	df = uf->df;
	/* Verify the starting node is legit */
	if ((snode < 0) || (snode > df->tot_cons[1]))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if (!df->down[1][snode])
		return dtl_error(DTL_CRIT_UNKNOWN);
	/* Obtain focal point values for the children of snode */
	if (call(TCL_get_P_masspoint(df,W_mid,LW_mid),"TCL_get_P_masspoint"))
		return dtl_kernel_error();
	for (i=1, tnode=df->down[1][snode]; tnode; i++, tnode=df->next[1][tnode]) {
		xorder[i] = i;
		xnode[i] = tnode;
		xval[i] = W_mid[tnode];
		}
	/* Obtain a ranking order for the children's focal point values */
	Wnode[0] = i-1; // nbr of node entries
	sort_b(xorder,xval,1,Wnode[0],TRUE);
	for (i=1; i<=Wnode[0]; i++) {
		Wnode[i] = xnode[xorder[i]];
		Wval[i] = xval[xorder[i]];
		if (cst_ext) {
			sprintf(msg,"    W%d = %.3lf\n",Wnode[i],Wval[i]);
			cst_log(msg);
			}
		}
	_smx_end();
	return DTL_OK;
	}


 /*********************************************************
  *
  *  For internal DTL use, limited error handling
  *
  *********************************************************/

/* Check if the weight tree is pure or not. Pure is defined as all nodes
 * with the same parent being either re-nodes or im-nodes, not a mixture.
 * Returns 1 if pure tree, 0 if impure tree, <0 if don't know due to error. */

int dtl_pure_W_tree() {

	if (!frame_loaded)
		return -1;
	if (PS) // no weights
		return -1;
	if (load_df0(0))
		return -1;
	return -1*TCL_pure_tree(uf->df,1);
	}


int DTLAPI DTI_pure_W_tree() {

	/* Call stack refurbishing */
	return dtl_pure_W_tree();
	}


/* Check if two weight nodes have the same parent or not. Returns 0 if
 * same parent, 1 if different parents, <0 if don't know due to error. */

int dtl_W_node_parents(int node1, int node2) {

	if (!frame_loaded)
		return -1;
	if (!PM)
		return -1;
	if (load_df0(0))
		return -1;
	return -1*TCL_different_parents(uf->df,1,node1,node2);
	}


int DTLAPI DTI_W_node_parents(int node1, int node2) {

	/* Call stack refurbishing */
	return dtl_W_node_parents(node1,node2);
	}


/* Check the number of siblings a weight node has (incl. itself). Returns the
 * number of siblings (1 if no other sibling), <0 if don't know due to error. */

int dtl_W_nbr_of_siblings(int node) {

	if (!frame_loaded)
		return -1;
	if (load_df0(0))
		return -1;
	return -1*TCL_nbr_of_siblings(uf->df,1,node);
	}


/* Check if the crit node is a real/final one or not (otherwise intermediate).
 * Returns >0 if real crit node, =0 if im-node or don't know due to error. */

int dtl_real_W_crit(int node) {

	/* Check if function can start */
	if (!frame_loaded)
		return 0;
	if (PS) // no weights
		return 0;
	if (load_df0(0))
		return 0;
	return TCL_get_V_index(uf->df,1,node);
	}


int DTLAPI DTI_real_W_crit(int node) {

	/* Call stack refurbishing */
	return dtl_real_W_crit(node);
	}


/* Check if there are any active midpoints in the weight base. Returns >0 if
 * there are midpoints, =0 if no active midpoints, <0 if don't know due to error. */

int dtl_nbr_W_midpoints() {
	int k,n_midpoints;
	struct d_frame *df;

	/* Check if function can start */
	if (!frame_loaded)
		return -1;
	if (PS)
		return -1;
	if (load_df0(0))
		return -1;
	df = uf->df;
	if (TCL_get_P_mbox(df,mbox_lobo,mbox_upbo))
		return -1;
	/* Calculate number of midpoints */
	n_midpoints = 0;
	for (k=1; k<=df->tot_cons[1]; k++)
		if (mbox_lobo[k] > -1.0)
			n_midpoints++;
	return n_midpoints;
	}


int DTLAPI DTI_nbr_W_midpoints() {

	/* Call stack refurbishing */
	return dtl_nbr_W_midpoints();
	}


#ifdef CONSOLE

 /*********************************************************
  *
  *  Console output for test and debug
  *
  *********************************************************/

static void print_W_stmt(struct base *P, int snbr) {
	int cr;

	if (P->stmt[snbr].n_terms == 1)
		if (cr=TCL_get_V_index(uf->df,1,P->stmt[snbr].cons[1]))
			printf("%2d: W%-2d= [%.3lf %.3lf]  K%d%s\n",
					snbr,P->stmt[snbr].cons[1],
					P->stmt[snbr].lobo,P->stmt[snbr].upbo,cr,
					P->stmt[snbr].alt[1]==1?"":" #alt[1]!=1");
		else
			printf("%2d: W%-2d= [%.3lf %.3lf]  *%s\n",
					snbr,P->stmt[snbr].cons[1],
					P->stmt[snbr].lobo,P->stmt[snbr].upbo,
					P->stmt[snbr].alt[1]==1?"":" #alt[1]!=1");
	else
		printf("%2d: error: %d terms in W statement\n",
				snbr,P->stmt[snbr].n_terms);
	}


void DTLAPI DTI_show_W_base() {
	int i;
	struct base *W;

	/* Check if function can start */
	if (!frame_loaded)
		return;
	if (PS)
		return;
	/* Check input parameters */
	if (load_df0(0))
		return;
	/* Display base */
	W = uf->df->P_base;
	printf("The weight base contains %d constraint%s\n",W->n_stmts,(W->n_stmts==1?"":"s"));
	for (i=1; i<=W->n_stmts; i++)
		print_W_stmt(W,i);
	}


void DTLAPI DTI_show_W_base_crit(int crit) {
	int i;
	struct base *W;

	if (!frame_loaded)
		return;
	if (PS)
		return;
	if (load_df0(0))
		return;
	W = uf->df->P_base;
	for (i=1; i<=W->n_stmts; i++)
		if (W->stmt[i].n_terms == 1)
			if (W->stmt[i].cons[1] == crit)
				print_W_stmt(W,i);
	}


void DTLAPI DTI_show_W_box() {
	int i,j,k,h,cr,n_entries;
	struct d_frame *df;

	/* Check if function can start */
	if (!frame_loaded)
		return;
	if (PS)
		return;
	/* Check input parameter */
	if (load_df0(0))
		return;
	df = uf->df;
	if (!df->P_base->box) {
		printf("The weight box contains no entries\n");
		return;
		}
	if (TCL_get_P_box(df,llobo,lupbo))
		return;
	/* Calculate number of entries */
	n_entries = 0;
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++)
			if ((llobo[k] > 0.0) || (lupbo[k] < 1.0))
				n_entries++;
	printf("The weight box contains %d entr%s\n",n_entries,(n_entries==1?"y":"ies"));
	/* Print entries in order */
	for (k=1,i=1,h=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			if ((llobo[k] > 0.0) || (lupbo[k] < 1.0)) {
				if (cr = TCL_get_V_index(df,i,j)) // for crawler
					if (lupbo[k]-llobo[k] > DTL_EPS)
						printf("%2d: W%-2d= [%.3lf %.3lf]  K%d %s\n",h++,j,llobo[k],lupbo[k],cr,i==1?"":"#ineffectual");
					else
						printf("%2d: W%-2d= %.3lf  K%d %s\n",h++,j,llobo[k],cr,i==1?"":"#ineffectual");
				else
					printf("%2d: W%-2d= %.3lf  *  %s\n",h++,j,llobo[k],i==1?"":"#ineffectual");
				}
			}
	}


void DTLAPI DTI_show_W_mbox() {
	int k,h,cr,n_midpoints;
	struct d_frame *df;

	/* Check if function can start */
	if (!frame_loaded)
		return;
	if (PS)
		return;
	/* Check input parameter */
	if (load_df0(0))
		return;
	df = uf->df;
	if (TCL_get_P_mbox(df,mbox_lobo,mbox_upbo))
		return;
	/* Calculate number of midpoints */
	n_midpoints = 0;
	for (k=1; k<=df->tot_cons[1]; k++)
		if (mbox_lobo[k] > -1.0)
			n_midpoints++;
	printf("The weight mbox contains %d midpoint%s\n",n_midpoints,(n_midpoints==1?"":"s"));
	/* Print midpoints in order */
	for (k=1,h=1; k<=df->tot_cons[1]; k++) {
		if (mbox_lobo[k] > -1.0) {
			if (cr=TCL_get_P_index(df,1,k))
				if (mbox_upbo[k]-mbox_lobo[k] > DTL_EPS)
					printf("%2d: W%-2d [%.3lf %.3lf]  K%d\n",h++,k,mbox_lobo[k],mbox_upbo[k],cr);
				else
					printf("%2d: W%-2d %.3lf  K%d\n",h++,k,mbox_lobo[k],cr);
			else
				printf("%2d: W%-2d %.3lf  *\n",h++,k,mbox_lobo[k]);
			}
		}
	}

#endif