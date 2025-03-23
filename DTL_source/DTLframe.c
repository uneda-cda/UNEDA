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
 *   File: DTLframe.c
 *
 *   Purpose: frame handling in DTL
 *
 *
 *   Functions exported outside DTL
 *   ------------------------------
 *   DTL_new_PS_flat_frame
 *   DTL_new_PS_tree_frame
 *   DTL_new_PM_flat_frame
 *   DTL_new_PM_tree_frame
 *   DTL_new_PM_crit_tree
 *   DTL_load_PM_crit
 *   DTL_unload_PM_crit
 *   DTL_delete_PM_crit
 *   DTL_new_DM_flat_frame
 *   DTL_new_DM_tree_frame
 *   DTL_new_SM_tree_frame
 *   DTL_load_frame
 *   DTL_unload_frame
 *   DTL_dispose_frame
 *   DTL_frame_name
 *   DTL_frame_type
 *   DTL_load_status
 *   DTL_PM_crit_exists
 *
 *   Functions outside module, debug use
 *   -----------------------------------
 *   DTI_split_DM_frame
 *   DTI_crit_exists
 *
 *   Functions outside module, inside DTL
 *   ------------------------------------
 *   dtl_dispose_frame
 *
 *   Functions internal to module
 *   ----------------------------
 *   dtl2tcl_tree
 *   dtl2tcl_dlevel
 *
 */

#include "DTL.h"
#include "DTLinternal.h"


 /********************************************************
  *
  *  Internal data
  *
  ********************************************************/

static t_matrix tnext,tdown;


 /*****************************************************
  *
  *  Create frames (PS,PM)
  *
  *****************************************************/

/* dtl2tcl: DTL to TCL conversion of general tree to vectors */

static int dtl2tcl_dlevel(int alt, int snode, ta_tree tree, 
		t_row next, t_row down, int *re_cons);

static int dtl2tcl_tree(int alt, int snode, ta_tree tree, 
		t_row next, t_row down, int *re_cons) {
	int tnode,last_node;

	last_node = DTL_TREE_ERROR; // last node is also error code if < 0
	for (tnode=tree[snode].down; tnode; tnode=tree[tnode].next) {
		if (tnode < 0)
			return DTL_TREE_ERROR;
		next[tnode] = tree[tnode].next;
		down[tnode] = tree[tnode].down;
		if (tree[tnode].down) {
			if (toupper(tree[tnode].type) == 'E') {
				/* Event node */
				if ((last_node = dtl2tcl_tree(alt,tnode,tree,next,down,re_cons)) < 0)
					return last_node;
				}
			else if ((toupper(tree[tnode].type) == 'D') || (toupper(tree[tnode].type) == 'F')) {
				/* Decison node */
				if ((last_node = dtl2tcl_dlevel(alt,tnode,tree,next,down,re_cons)) < 0)
					return last_node;
				}
			else {
				if (cst_on) {
					sprintf(msg," dtl2tcl_tree failed at alt=%d node=%d type=%c next=%d down=%d\n",
							alt,tnode,tree[tnode].type,tree[tnode].next,tree[tnode].down);
					cst_log(msg);
					}
				return DTL_TREE_ERROR;
				}
			}
		else {
			last_node = tnode;
			(*re_cons)++;
			}
		}
	return last_node;
	}


static int dtl2tcl_dlevel(int alt, int snode, ta_tree tree, 
		t_row next, t_row down, int *re_cons) {
	int tnode,last_node;

	last_node = DTL_TREE_ERROR; // last node is also error code if < 0
	for (tnode=tree[snode].down; tnode; tnode=tree[tnode].next) {
		if (tnode < 0)
			return DTL_TREE_ERROR;
		next[tnode] = tree[tnode].next;
		down[tnode] = tree[tnode].down;
		if (tree[tnode].down) {
			if (toupper(tree[tnode].type) == 'E') {
				if ((last_node = dtl2tcl_tree(alt,tnode,tree,next,down,re_cons)) < 0)
					return last_node;
				}
			else /* D or F */ {
				if (cst_on) {
					sprintf(msg," dtl2tcl_dlevel failed at alt=%d node=%d type=%c next=%d down=%d\n",
							alt,tnode,tree[tnode].type,tree[tnode].next,tree[tnode].down);
					cst_log(msg);
					}
				return DTL_TREE_ERROR;
				}
			}
		else {
			last_node = tnode;
			(*re_cons)++;
			}
		}
	return last_node;
	}


rcode DTLAPI DTL_new_PS_flat_frame(int ufnbr, int n_alts, int n_cons[]) {
	int j;
	struct user_frame *tmp_uf = NULL;
	char name[20];

	/* Begin single thread semaphore */
	_smx_begin("PSF");
	/* Protect log function */
	_certify_ptr(n_cons,1);
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_new_PS_flat_frame(%d,%d",ufnbr,n_alts);
		cst_log(msg);
		if (n_alts <= MAX_ALTS)
			for (j=1; j<=n_alts; j++) {
				sprintf(msg,",%d",n_cons[j]);
				cst_log(msg);
				}
		cst_log(")\n");
		}
	/* Check if function can start */
	if (!dtl_is_init())
		return dtl_error(DTL_STATE_ERROR);
	/* Check input parameters */
	if ((ufnbr < 1) || (ufnbr > MAX_FRAMES))
		return dtl_error(DTL_FRAME_UNKNOWN);
	if (n_alts > MAX_ALTS)
		return dtl_error(DTL_ALT_OVERFLOW);
	for (j=1; j<=n_alts; j++)
		if (n_cons[j] > MAX_COPA)
			return dtl_error(DTL_CONS_OVERFLOW);
	/* Allocate new user frame */
	if (!(tmp_uf = new_uf(ufnbr))) {
		return dtl_error(DTL_FRAME_EXISTS);
		}
	sprintf(name,"PS-%03d",ufnbr);
	strcpy(tmp_uf->frame_name,name);
	if (call(TCL_create_flat_frame(&(tmp_uf->df),n_alts,n_cons),"TCL_create_flat_frame")) {
		dispose_uf(ufnbr);
		return dtl_kernel_error();
		}
	/* Set df name */
	sprintf(tmp_uf->df->name,"%s-01F",tmp_uf->frame_name);
	/* Set up uf */
	tmp_uf->frame_type = PS_FRAME;
	tmp_uf->frame_nbr = ufnbr;
	tmp_uf->n_alts = n_alts;
	tmp_uf->n_crit = 1;
	tmp_uf->n_sh = 1;
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


rcode DTLAPI DTL_new_PS_tree_frame(int ufnbr, int n_alts, int n_nodes[], tt_tree xtree) {
	int i,j,k,errc,re_cons;
	struct user_frame *tmp_uf = NULL;
	char name[20];

	/* Begin single thread semaphore */
	_smx_begin("PST");
	/* Protect log function */
	_certify_ptr(n_nodes,1);
	_certify_ptr(xtree,2);
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"%s(%d,%d)","DTL_new_PS_tree_frame",ufnbr,n_alts);
		cst_log(msg);
		if (n_alts <= MAX_ALTS)
			for (j=1; j<=n_alts; j++) {
				sprintf(msg,"\n    A%d:%d",j,n_nodes[j]);
				cst_log(msg);
				for (k=1; k<=n_nodes[j]; k++) {
					if (xtree[j][k].type && (n_nodes[j] <= MAX_NOPA))
						sprintf(msg,"(%d,%c,%d,%d)",k,xtree[j][k].type,xtree[j][k].next,xtree[j][k].down);
					else
						sprintf(msg,"(%d,*,*,*)",k);
					cst_log(msg);
					}
				}
		cst_log("\n");
		}
	/* Check if function can start */
	if (!dtl_is_init())
		return dtl_error(DTL_STATE_ERROR);
	/* Check input parameters */
	if ((ufnbr < 1) || (ufnbr > MAX_FRAMES))
		return dtl_error(DTL_FRAME_UNKNOWN);
	if (n_alts > MAX_ALTS)
		return dtl_error(DTL_ALT_OVERFLOW);
	/* Check nbr of nodes per alt */
	for (j=1; j<=n_alts; j++) {
		/* Protect t-data structures in dtl2tcl */
		if (n_nodes[j] > MAX_NOPA)
			return dtl_error(DTL_NODE_OVERFLOW);
		}
	/* Allocate new user frame */
	if (!(tmp_uf = new_uf(ufnbr))) {
		return dtl_error(DTL_FRAME_EXISTS);
		}
	sprintf(name,"PS-%03d",ufnbr);
	strcpy(tmp_uf->frame_name,name);
	for (i=1; i<=n_alts; i++) {
		/* Set up root node if missing */
		xtree[i][0].type ='D';
		xtree[i][0].next = 0;
		xtree[i][0].down = 1;
		/* Fill vectors for next & down */
		re_cons = 0;
		if ((errc = dtl2tcl_tree(i,0,xtree[i],tnext[i],tdown[i],&re_cons)) < 0) {
			if (cst_on) {
				sprintf(msg," %s alt=%d rc=%d\n","DTL_new_PS_tree_frame",i,errc);
				cst_log(msg);
				}
			dispose_uf(ufnbr);
			return dtl_error(DTL_TREE_ERROR);
			}
		if (errc != n_nodes[i]) {
			dispose_uf(ufnbr);
			return dtl_error(DTL_TREE_ERROR);
			}
		if (re_cons > MAX_COPA) {
			dispose_uf(ufnbr);
			return dtl_error(DTL_CONS_OVERFLOW);
			}
		}
	if (call(TCL_create_tree_frame(&(tmp_uf->df),n_alts,n_nodes,tnext,tdown),"TCL_create_tree_frame")) { /* tmp_uf.df !! */
		dispose_uf(ufnbr);
		return dtl_kernel_error();
		}
	/* Set df name */
	sprintf(tmp_uf->df->name,"%s-01T",tmp_uf->frame_name);
	/* Finish uf initialisation */
	tmp_uf->frame_type = PS_FRAME;
	tmp_uf->frame_nbr = ufnbr;
	tmp_uf->n_alts = n_alts;
	tmp_uf->n_crit = 1;
	tmp_uf->n_sh = 1;
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


rcode DTLAPI DTL_new_PM_flat_frame(int ufnbr, int n_crit, int n_alts) {
	int j,n_cons[MAX_ALTS+1];
	struct user_frame *tmp_uf = NULL;
	char name[20];

	/* Begin single thread semaphore */
	_smx_begin("PMF");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_new_PM_flat_frame(%d,%d,%d)\n",ufnbr,n_crit,n_alts);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!dtl_is_init())
		return dtl_error(DTL_STATE_ERROR);
	/* Check input parameters */
	if ((ufnbr < 1) || (ufnbr > MAX_FRAMES))
		return dtl_error(DTL_FRAME_UNKNOWN);
	if (n_alts > MAX_ALTS)
		return dtl_error(DTL_ALT_OVERFLOW);
	if (n_crit > MAX_CRIT)
		return dtl_error(DTL_CRIT_OVERFLOW);
	/* Allocate new user frame */
	if (!(tmp_uf = new_uf(ufnbr))) {
		return dtl_error(DTL_FRAME_EXISTS);
		}
	sprintf(name,"PM-%03d",ufnbr);
	strcpy(tmp_uf->frame_name,name);
	n_cons[1] = n_crit;
	for (j=2; j<=n_alts; j++)
		n_cons[j] = 1;
	if (call(TCL_create_flat_frame(&(tmp_uf->df),n_alts,n_cons),"TCL_create_flat_frame")) {
		dispose_uf(ufnbr);
		return dtl_kernel_error();
		}
	sprintf(tmp_uf->df->name,"%s-MCF",tmp_uf->frame_name);
	tmp_uf->frame_type = PM_FRAME;
	tmp_uf->frame_nbr = ufnbr;
	tmp_uf->n_alts = n_alts;
	tmp_uf->n_crit = n_crit;
	tmp_uf->n_sh = 1;
	tmp_uf->df_list[0] = tmp_uf->df;
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


rcode DTLAPI DTL_new_PM_tree_frame(int ufnbr, int n_alts, int n_wtnodes, ta_tree wtree) {
	int n_nodes[MAX_ALTS+1];
	int i,k,errc,re_cons;
	struct user_frame *tmp_uf = NULL;
	char name[20];

	/* Begin single thread semaphore */
	_smx_begin("PMT");
	/* Protect log function */
	_certify_ptr(wtree,1);
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"%s(%d,%d,%d)","DTL_new_PM_tree_frame",ufnbr,n_alts,n_wtnodes);
		cst_log(msg);
		sprintf(msg,"\n    W:%d",n_wtnodes);
		cst_log(msg);
		if (n_wtnodes <= MAX_NOPA)
			for (k=1; k<=n_wtnodes; k++) {
				if (wtree[k].type)
					sprintf(msg,"(%d,%c,%d,%d)",k,wtree[k].type,wtree[k].next,wtree[k].down);
				else
					sprintf(msg,"(%d,*,*,*)",k);
				cst_log(msg);
				}
		cst_log("\n");
		}
	/* Check if function can start */
	if (!dtl_is_init())
		return dtl_error(DTL_STATE_ERROR);
	/* Check input parameters */
	if ((ufnbr < 1) || (ufnbr > MAX_FRAMES))
		return dtl_error(DTL_FRAME_UNKNOWN);
	if (n_alts < 2)
		return dtl_error(DTL_TOO_FEW_ALTS);
	if (n_alts > MAX_ALTS)
		return dtl_error(DTL_ALT_OVERFLOW);
	/* Check nbr of weight nodes */
	if (n_wtnodes > MAX_NOPA)
		return dtl_error(DTL_NODE_OVERFLOW);
	/* Allocate new user frame */
	if (!(tmp_uf = new_uf(ufnbr))) {
		return dtl_error(DTL_FRAME_EXISTS);
		}
	sprintf(name,"PM-%03d",ufnbr);
	strcpy(tmp_uf->frame_name,name);
	/* Fill vectors for next & down */
	re_cons = 0;
	if ((errc = dtl2tcl_tree(1,0,wtree,tnext[1],tdown[1],&re_cons)) < 0) {
		if (cst_on) {
			sprintf(msg," %s wt rc=%d\n","DTL_new_PM_tree_frame",errc);
			cst_log(msg);
			}
		dispose_uf(ufnbr);
		return dtl_error(DTL_TREE_ERROR);
		}
	if (errc != n_wtnodes) {
		dispose_uf(ufnbr);
		return dtl_error(DTL_TREE_ERROR);
		}
	if (re_cons > MAX_CRIT) {
		dispose_uf(ufnbr);
		return dtl_error(DTL_CRIT_OVERFLOW);
		}
	n_nodes[1] = n_wtnodes;
	/* Pseudo-alternatives */
	for (i=2; i<=n_alts; i++) {
		tnext[i][0] = 0;
		tdown[i][0] = 1;
		tnext[i][1] = 0;
		tdown[i][1] = 0;
		n_nodes[i] = 1;
		}
	if (call(TCL_create_tree_frame(&(tmp_uf->df),n_alts,n_nodes,tnext,tdown),"TCL_create_tree_frame")) { /* tmp_uf.df !! */
		dispose_uf(ufnbr);
		return dtl_kernel_error();
		}
	/* Set df name */
	sprintf(tmp_uf->df->name,"%s-MCT",tmp_uf->frame_name);
	/* Finish uf initialisation */
	tmp_uf->frame_type = PM_FRAME;
	tmp_uf->frame_nbr = ufnbr;
	tmp_uf->n_alts = n_alts;
	tmp_uf->n_crit = tmp_uf->df->n_cons[1];
	tmp_uf->n_sh = 1;
	tmp_uf->df_list[0] = tmp_uf->df;
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


rcode DTLAPI DTL_new_PM_crit_tree(int crit, int n_nodes[], tt_tree xtree) {
	int i,j,k,errc,re_cons;

	/* Begin single thread semaphore */
	_smx_begin("PMCT");
	/* Protect log function */
	_certify_ptr(n_nodes,1);
	_certify_ptr(xtree,2);
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_new_PM_crit_tree(%d)",crit);
		cst_log(msg);
		if (frame_loaded && n_nodes) {
			for (j=1; j<=uf->n_alts; j++) {
				sprintf(msg,"\n    A%d:%d",j,n_nodes[j]);
				cst_log(msg);
				for (k=1; k<=n_nodes[j]; k++) {
					if (xtree[j][k].type && (n_nodes[j] <= MAX_NOPA))
						sprintf(msg,"(%d,%c,%d,%d)",k,xtree[j][k].type,xtree[j][k].next,xtree[j][k].down);
					else
						sprintf(msg,"(%d,*,*,*)",k);
					cst_log(msg);
					}
				}
			}
		cst_log("\n");
		}
	/* Check if function can start */
	if (!dtl_is_init())
		return dtl_error(DTL_STATE_ERROR);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (!PM)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if ((crit < 1) || (crit > uf->n_crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if (uf->df_list[crit])
		return dtl_error(DTL_CRIT_EXISTS);
	/* Check nbr of nodes per alt */
	for (j=1; j<=uf->n_alts; j++) {
		/* Protect t-data structures in dtl2tcl */
		if (n_nodes[j] > MAX_NOPA)
			return dtl_error(DTL_NODE_OVERFLOW);
		}
	for (i=1; i<=uf->n_alts; i++) {
		/* Set up root node if missing */
		xtree[i][0].type ='D';
		xtree[i][0].next = 0;
		xtree[i][0].down = 1;
		/* Fill vectors for next, down */
		re_cons = 0;
		if ((errc = dtl2tcl_tree(i,0,xtree[i],tnext[i],tdown[i],&re_cons)) < 0) {
			if (cst_on) {
				sprintf(msg," DTL_new_PM_crit_tree alt=%d rc=%d\n",i,errc);
				cst_log(msg);
				}
			return dtl_error(DTL_TREE_ERROR);
			}
		if (errc != n_nodes[i])
			return dtl_error(DTL_TREE_ERROR);
		if (re_cons > MAX_COPA)
			return dtl_error(DTL_CONS_OVERFLOW);
		}
	if (call(TCL_create_tree_frame(&(uf->df_list[crit]),uf->n_alts,n_nodes,tnext,tdown),"TCL_create_tree_frame")) {
		uf->df_list[crit] = NULL;
		return dtl_kernel_error();
		}
	if (load_df1(crit)) {
		/* Error in frame */
		call(TCL_dispose_frame(uf->df_list[crit]),"TCL_dispose_frame");
		uf->df_list[crit] = NULL;
		return dtl_error(DTL_FRAME_CORRUPT);
		}
	/* Set df name */
	sprintf(uf->df_list[crit]->name,"%s-%03dT",uf->frame_name,crit%1000); // exactly 2 digits
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


rcode DTLAPI DTL_load_PM_crit(int crit, int ufnbr) {
	struct user_frame *ps_uf;

	/* Begin single thread semaphore */
	_smx_begin("LPMC");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_load_PM_crit(%d,%d)\n",crit,ufnbr);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (!PM)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if ((crit < 1) || (crit > uf->n_crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if (uf->df_list[crit])
		return dtl_error(DTL_CRIT_EXISTS);
	/* Find PS-tree user frame */
	if ((ps_uf = get_uf(ufnbr)) == NULL) {
		return dtl_error(DTL_FRAME_UNKNOWN);
		}
	if (ps_uf->frame_type != PS_FRAME)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check same number of alts */
	if (uf->n_alts != ps_uf->n_alts)
		return dtl_error(DTL_ALT_MISMATCH);
	/* Incorporate into PM-frame */
	uf->df_list[crit] = ps_uf->df;
	if (load_df1(crit)) {
		/* Bogus frame */
		uf->df_list[crit] = NULL;
		return dtl_error(DTL_FRAME_CORRUPT);
		}
	if (!dispose_uf(ufnbr)) {
		/* Could not dispose, don't want multiple copies */
		uf->df_list[crit] = NULL;
		return dtl_error(DTL_FRAME_CORRUPT);
		}
	/* Set df name */
	if (uf->df_list[crit]->tree)
		sprintf(uf->df_list[crit]->name,"%s-%03dT",uf->frame_name,crit%1000);
	else
		sprintf(uf->df_list[crit]->name,"%s-%03dF",uf->frame_name,crit%1000);
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


rcode DTLAPI DTL_unload_PM_crit(int crit, int new_ufnbr) {
	struct user_frame *ps_uf;
	char name[20];

	/* Begin single thread semaphore */
	_smx_begin("UPMC");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_unload_PM_crit(%d,%d)\n",crit,new_ufnbr);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (!PM)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if ((new_ufnbr < 1) || (new_ufnbr > MAX_FRAMES))
		return dtl_error(DTL_FRAME_UNKNOWN);
	if ((crit < 1) || (crit > uf->n_crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if (!(uf->df_list[crit]))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if (load_df0(0))
		return dtl_error(DTL_SYS_CORRUPT);
	/* Allocate new user frame */
	if (!(ps_uf = new_uf(new_ufnbr))) {
		return dtl_error(DTL_FRAME_EXISTS);
		}
	ps_uf->df = uf->df_list[crit];
	sprintf(name,"PS-%03d",new_ufnbr);
	strcpy(ps_uf->frame_name,name);
	if (ps_uf->df->tree)
		sprintf(ps_uf->df->name,"%s-01T",ps_uf->frame_name);
	else
		sprintf(ps_uf->df->name,"%s-01F",ps_uf->frame_name);
	/* Finish uf initialisation */
	ps_uf->frame_type = PS_FRAME;
	ps_uf->frame_nbr = new_ufnbr;
	ps_uf->n_alts = uf->n_alts;
	ps_uf->n_crit = 1;
	ps_uf->n_sh = 1;
	uf->df_list[crit] = NULL;
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


rcode DTLAPI DTL_delete_PM_crit(int crit) {

	/* Begin single thread semaphore */
	_smx_begin("DPMC");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_delete_PM_crit(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (!PM)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if ((crit < 1) || (crit > uf->n_crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if (!(uf->df_list[crit]))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if (load_df0(0))
		return dtl_error(DTL_SYS_CORRUPT);
	/* Remove user frame */
	if (call(TCL_dispose_frame(uf->df_list[crit]),"TCL_dispose_frame")) {
		return dtl_kernel_error();
		}
	uf->df_list[crit] = NULL;
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


 /*****************************************************
  *
  *  Create pseudo frames (DM,SM)
  *
  *****************************************************/

/* Note that DM create/split functions (unlike PM or PS) require that
 * no frame is loaded and will leave the new frame unloaded as well.
 * Further, the frame with the highest number (MAX_FRAMES) is used. */

// DTL layer 0: above DTL proper

rcode DTLAPI DTL_new_DM_flat_frame(int ufnbr, int n_crit, int n_alts) {
	rcode rc;
	int i,n_cons[MAX_ALTS+1];

	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_new_DM_flat_frame(%d,%d,%d) -->\n",ufnbr,n_crit,n_alts);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!dtl_is_init())
		return dtl_error(DTL_STATE_ERROR);
	if (frame_loaded)
		return dtl_error(DTL_FRAME_IN_USE);
	/* Create and load mother frame */
	if (rc = DTL_new_PM_flat_frame(ufnbr,n_crit,n_alts))
		return rc;
	if (rc = DTL_load_frame(ufnbr)) {
		DTL_dispose_frame(ufnbr);
		return rc;
		}
	uf->frame_name[0] = uf->df->name[0] = 'D'; // "PM"->"DM"
	/* Prepare singleton template */
	for (i=1; i<=n_alts; i++)
		n_cons[i] = 1;
	for (i=1; i<=n_crit; i++) {
		/* Create and attach singleton frame */
		if (rc = DTL_new_PS_flat_frame(MAX_FRAMES,n_alts,n_cons)) {
			DTL_unload_frame();
			DTL_dispose_frame(ufnbr);
			return rc;
			}
		if (rc = DTL_load_PM_crit(i,MAX_FRAMES)) {
			DTL_unload_frame();
			DTL_dispose_frame(ufnbr);
			return rc;
			}
		}
	/* Finish up */
	rc = DTL_unload_frame();
	if (cst_on) {
		sprintf(msg,"--> DTL_new_DM_flat_frame(%d,%d,%d) END\n",ufnbr,n_crit,n_alts);
		cst_log(msg);
		}
	return rc;
	}


// DTL layer 0: above DTL proper

rcode DTLAPI DTL_new_DM_tree_frame(int ufnbr, int n_alts, int n_wtnodes, ta_tree wtree) {
	rcode rc;
	int i,c,n_crit=0,n_cons[MAX_ALTS+1];

	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_new_DM_tree_frame(%d,%d,%d) -->\n",ufnbr,n_wtnodes,n_alts);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(wtree,1);
	if (!dtl_is_init())
		return dtl_error(DTL_STATE_ERROR);
	if (frame_loaded)
		return dtl_error(DTL_FRAME_IN_USE);
	/* Create and load mother tree frame */
	if (rc = DTL_new_PM_tree_frame(ufnbr,n_alts,n_wtnodes,wtree))
		return rc;
	if (rc = DTL_load_frame(ufnbr)) {
		DTL_dispose_frame(ufnbr);
		return rc;
		}
	uf->frame_name[0] = uf->df->name[0] = 'D'; // "PM"->"DM"
	/* Prepare singleton template */
	for (i=1; i<=n_alts; i++)
		n_cons[i] = 1;
	for (i=1; i<=n_wtnodes; i++)
		if (c = dtl_real_W_crit(i)) {
			/* Create and attach singleton frame */
			if (rc = DTL_new_PS_flat_frame(MAX_FRAMES,n_alts,n_cons)) {
				DTL_unload_frame();
				DTL_dispose_frame(ufnbr);
				return rc;
				}
			if (rc = DTL_load_PM_crit(c,MAX_FRAMES)) {
				DTL_unload_frame();
				DTL_dispose_frame(ufnbr);
				return rc;
				}
			n_crit++;
			}
	/* Finish up */
	rc = DTL_unload_frame();
	if (cst_on) {
		sprintf(msg,"--> DTL_new_DM_tree_frame(%d,%d,%d) END\n",ufnbr,n_crit,n_alts);
		cst_log(msg);
		}
	return rc;
	}


/* Split a DM-frame read from file into DTL internal format.
 * This is only used for testing. For prod, use create.
 *
 * NOTE: MC pseudo-alts 2..n must be manually set to 1 node
 *       if the splitted DM-frame is written to disk as PM. */

// DTL layer 0: above DTL proper

rcode DTLAPI DTI_split_DM_frame(int ufnbr) {
	rcode rc;
	int i,crit,n_alts,n_crit,n_cons[MAX_ALTS+1];
	struct base *V0;
	struct user_stmt_rec *u_stmt;
	struct user_frame *ufx;

	/* Log function call */
	if (cst_ext) {
		sprintf(msg,"DTI_split_DM_frame(%d) -->\n",ufnbr);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!dtl_is_init())
		return dtl_error(DTL_STATE_ERROR);
	if (frame_loaded)
		return DTL_FRAME_IN_USE;
	/* Check input parameter */
	if ((ufx = get_uf(ufnbr)) == NULL) {
		return DTL_FRAME_UNKNOWN;
		}
	if (ufx->frame_type != DM_FRAME)
		return DTL_WRONG_FRAME_TYPE;
	/* Semi-proper PM-frame (alts 2..n too big) */
	ufx->frame_type = PM_FRAME;
	strcat(ufx->df->name,ufx->df->tree?"-MCT":"-MCF");
	ufx->df_list[0] = ufx->df;
#ifdef ADJUST_DM
	for (i=2; i<=ufx->df->n_alts; i++) {
		/* Adjust MC pseudo-alts 2..n (has side effects) */
		ufx->df->tot_cons[i] = ufx->df->n_cons[i] = 1;
		ufx->df->im_cons[i] = 0;
		}
#endif
	/* Load as PM-frame */
	if (rc = DTL_load_frame(ufnbr))
		return rc;
	if (load_df0(0))
		return DTL_SYS_CORRUPT;
	/* Prepare for DM -> PM split */
	V0 = ufx->df->V_base;
	n_alts = ufx->n_alts;
	n_crit = ufx->n_crit;
	for (i=1; i<=n_alts; i++)
		n_cons[i] = 1;
	for (i=1; i<=n_crit; i++) {
		/* Create and connect singleton frame */
		if (rc = DTL_new_PS_flat_frame(MAX_FRAMES,n_alts,n_cons)) {
			DTL_unload_frame();
			return rc;
			}
		if (rc = DTL_load_PM_crit(i,MAX_FRAMES)) {
			DTL_unload_frame();
			return rc;
			}
		}
	/* Transfer statements DM -> PM */
	for (i=1; i<=V0->n_stmts; i++) {
		u_stmt = (struct user_stmt_rec*)&V0->stmt[i];
		crit = u_stmt->cons[1]; // select criterion
		u_stmt->cons[1] = 1; // from DM = one cons in each alt
		if ((rc = DTL_add_V_statement(crit,u_stmt)) < DTL_OK) {
			DTL_unload_frame();
			return rc;
			}
		}
	/* Finish up */
	V0->n_stmts = 0; // reset source
	rc = DTL_unload_frame();
	if (cst_ext) {
		sprintf(msg,"--> DTI_split_DM_frame(%d) END\n",ufnbr);
		cst_log(msg);
		}
	return rc;
	}


/* Create an SM-frame (stakeholders) and prepare according to mode flag:
 * Mode: 0 = only copy sh #1 to all other sh
 *       1 = create SM mother frame + copy
 *       2 = create PS crit frames + copy
 *       3 = create SM + PS frames + copy
 * NOTE: It is up to the caller to supply a symmetric stakeholder hierarchy in
 *       which the lowest level are the criteria (i.e. many stakeholder levels
 *       but only one criteria level) since this is mostly stakeholder-focused
 */

// DTL layer 0: above DTL proper

rcode DTLAPI DTL_new_SM_tree_frame(int ufnbr, int mode, int n_alts, int n_sh, int n_wtnodes, ta_tree wtree) {
	rcode rc;
	int i,c,n_crit1,n_cons[MAX_ALTS+1];

	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_new_SM_tree_frame(%d,%d,%d,%d) -->\n",ufnbr,n_alts,n_sh,n_wtnodes);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(wtree,1);
	if (!dtl_is_init())
		return dtl_error(DTL_STATE_ERROR);
	if (frame_loaded)
		return DTL_FRAME_IN_USE;
	/* Check input parameter */
	if (n_sh < 2)
		return DTL_INPUT_ERROR;
	/* Modes 1 & 3 = Create mother tree frame */
	if (mode & 0x01)
		if (rc = DTL_new_PM_tree_frame(ufnbr,n_alts,n_wtnodes,wtree))
			return rc;
	/* Load mother tree frame */
	if (rc = DTL_load_frame(ufnbr)) {
		DTL_dispose_frame(ufnbr);
		return rc;
		}
	uf->frame_name[0] = uf->df->name[0] = 'S'; // "PM"->"SM"
	/* Investigate mother tree frame */
	if (uf->n_crit%n_sh) { // unbalanced
		DTL_unload_frame();
		DTL_dispose_frame(ufnbr);
		return DTL_TREE_ERROR;
		}
	n_crit1 = uf->n_crit/n_sh;
	/* Prepare singleton template */
	for (i=1; i<=n_alts; i++)
		n_cons[i] = 1;
	/* Modes 2 & 3 = inject PS frames */
	if (mode & 0x02) {
		/* First sh: create and attach singleton frames */
		for (i=1,c=0; c<n_crit1; i++)
			if (c = dtl_real_W_crit(i)) {
				if (rc = DTL_new_PS_flat_frame(MAX_FRAMES,n_alts,n_cons)) {
					DTL_unload_frame();
					DTL_dispose_frame(ufnbr);
					return rc;
					}
				if (rc = DTL_load_PM_crit(c,MAX_FRAMES)) {
					DTL_unload_frame();
					DTL_dispose_frame(ufnbr);
					return rc;
					}
				}
		/* Other sh: copy singleton frame pointer from 1st sh */
		for (; i<=n_wtnodes; i++)
			if (c = dtl_real_W_crit(i))
				uf->df_list[c] = uf->df_list[(c-1)%n_crit1+1];
		}
	uf->n_sh = n_sh;
	/* Finish up */
	rc = DTL_unload_frame();
	if (cst_on) { // note that 4th parameter differs
		sprintf(msg,"--> DTL_new_SM_tree_frame(%d,%d,%d,%d) END\n",ufnbr,n_alts,n_sh,n_crit1);
		cst_log(msg);
		}
	return rc;
	}


 /*****************************************************
  *
  *  Load/unload/dispose of frames
  *
  *****************************************************/

 /*
  * DTL_load_frame call semantics: Attach
  * the frame to TCL and check for consistency
  */

rcode DTLAPI DTL_load_frame(int ufnbr) {
	int i,nbr_p_frames=0;

	/* Begin single thread semaphore */
	_smx_begin("LOAD");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_load_frame(%d)\n",ufnbr);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!dtl_is_init())
		return dtl_error(DTL_STATE_ERROR);
	if (frame_loaded)
		return dtl_error(DTL_FRAME_IN_USE);
	/* Check input parameters */
	if ((uf = get_uf(ufnbr)) == NULL) {
		return dtl_error(DTL_FRAME_UNKNOWN);
		}
	if (DM)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	else if (PM) {
		/* Verify consistency for each criterion */
		for (i=1; i<=uf->n_crit; i++)
			if (uf->df_list[i]) {
				/* Exercise load-unload cycle */
				if (call(TCL_attach_frame(uf->df_list[i]),"TCL_attach_frame"))
					return dtl_kernel_error();
				if (uf->df_list[i]->tot_cons[0] > uf->df_list[i]->n_alts)
					nbr_p_frames++;
				if (call(TCL_detach_frame(uf->df_list[i]),"TCL_detach_frame"))
					return dtl_kernel_error();
				}
		/* Load weight tree */
		if (uf->n_alts != uf->df_list[0]->n_alts)
			return dtl_error(DTL_FRAME_CORRUPT);
		if (call(TCL_attach_frame(uf->df_list[0]),"TCL_attach_frame"))
			return dtl_kernel_error();
		uf->df = uf->df_list[0];
		uf->load_crit = 0;
		}
	else { // (PS)
		/* Load full problem */
		if (uf->n_alts != uf->df->n_alts)
			return dtl_error(DTL_FRAME_CORRUPT);
		if (call(TCL_attach_frame(uf->df),"TCL_attach_frame"))
			return dtl_kernel_error();
		// could do: nbr_p_frames = uf->df->tot_cons[0]>uf->df->n_alts;
		}
	dtl_error_count = 0;
	frame_loaded = ufnbr;
	eval_cache_invalidate();
	/* End single thread semaphore */
	_smx_end();
	return nbr_p_frames;
	}


 /*
  * Call semantics: Detach the frame and free the interface for new frames
  *                 DTL_unload_frame2 returns the frame that was unloaded
  */

rcode DTLAPI DTL_unload_frame() {

	/* Begin single thread semaphore */
	_smx_begin("UNL");
	/* Log function call */
	if (cst_on)
		cst_log("DTL_unload_frame()\n");
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Detach/unload */
	if (PM) {
		if (uf->load_crit>=0) {
			/* Unload current criterion */
			if (call(TCL_detach_frame(uf->df),"TCL_detach_frame"))
				return dtl_kernel_error();
			uf->load_crit = -1;
			}
		uf->df = NULL;
		}
	else
		/* Unload full problem */
		if (call(TCL_detach_frame(uf->df),"TCL_detach_frame"))
			return dtl_kernel_error();
	uf = NULL;
	frame_loaded = 0;
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


rcode DTLAPI DTL_unload_frame2() {
	rcode rc;
	int fnr;

	fnr = frame_loaded;
	if (rc = DTL_unload_frame())
		return rc;
	return fnr;
	}


 /*
  * Call semantics: Dispose resources
  */

// DTL layer 1: DTL API level

rcode dtl_dispose_frame(int ufnbr) {
	rcode rc;
	int i,n_crit;
	struct user_frame *tmp_uf;

	rc = DTL_OK;
	/* Check if function can start */
	if (ufnbr && (ufnbr == frame_loaded))
		return dtl_error(DTL_FRAME_IN_USE);
	/* Check input parameters */
	if ((tmp_uf=get_uf(ufnbr)) == NULL) {
		return dtl_error(DTL_FRAME_UNKNOWN);
		}
	/* Release resources */
	if (tmp_uf->frame_type == PM_FRAME) {
		n_crit = tmp_uf->n_crit/tmp_uf->n_sh;
		for (i=0; i<=n_crit; i++)
			if (tmp_uf->df_list[i])
				if (call(TCL_dispose_frame(tmp_uf->df_list[i]),"TCL_dispose_frame")) {
					if (!rc)
						rc = dtl_kernel_error();
					}
		}
	else {
		if (call(TCL_dispose_frame(tmp_uf->df),"TCL_dispose_frame")) {
			rc = dtl_kernel_error();
			}
		}
	/* Dispose DTL resource even if TCL failed */
	if (!dispose_uf(ufnbr)) {
		return dtl_error(DTL_FRAME_CORRUPT);
		}
	return rc;
	}


rcode DTLAPI DTL_dispose_frame(int ufnbr) {
	rcode rc;

	/* Begin single thread semaphore */
	_smx_begin("DISP");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_dispose_frame(%d)\n",ufnbr);
		cst_log(msg);
		}
	/* Release resources */
	rc = dtl_dispose_frame(ufnbr);
	/* End single thread semaphore */
	if (!rc)
		_smx_end();
	return rc;
	}


 /*****************************************************
  *
  *  Frame information services
  *
  *****************************************************/

rcode DTLAPI DTL_frame_name(char *fname, unsigned c_size, int *ftype) {
	unsigned u,stop;

	/* Log function call */
	if (cst_ext) {
		sprintf(msg,"DTL_frame_name(strg[%u])\n",c_size);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(fname,1);
	_certify_ptr(ftype,2);
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check that file name fits */
	stop = (unsigned)strlen(uf->frame_name);
	if (stop >= c_size)
		return DTL_BUFFER_OVERRUN;
	strcpy(fname,uf->frame_name);
	/* Replace underscore with space */
	for (u=0U; u<=stop; u++)
		if (fname[u]=='_')
			fname[u]=' ';
	*ftype = uf->frame_type;
	if (cst_ext) {
		sprintf(msg," frame name=\"%.40s\" type=%2d\n",fname,*ftype);
		cst_log(msg);
		}
	return DTL_OK;
	}


rcode DTLAPI DTL_load_status(int *f_loaded) {

	/* Log function call */
	if (cst_ext) {
		sprintf(msg,"DTL_load_status(%d)\n",frame_loaded);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(f_loaded,1);
	*f_loaded = frame_loaded;
	return DTL_OK;
	}


rcode DTLAPI DTL_frame_type(int ufnr, int *type) {

	/* Log function call */
	if (cst_ext) {
		sprintf(msg,"DTL_frame_type(%d)\n",ufnr);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(type,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check index bounds */
	if ((ufnr < 0) || (ufnr > MAX_FRAMES))
		return dtl_error(DTL_FRAME_UNKNOWN);
	else if (!ufnr) // 0 = currently loaded
		*type = uf->frame_type;
	else if (uf_list[ufnr])
		*type = uf_list[ufnr]->frame_type;
	else
		*type = 0;
	return DTL_OK;
	}


rcode DTLAPI DTL_PM_crit_exists(int crit, int *exists) {

	/* Log function call */
	if (cst_ext) {
		sprintf(msg,"DTL_PM_crit_exists(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(exists,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	if (!PM)
		return dtl_error(DTL_WRONG_FRAME_TYPE);
	/* Check input parameters */
	if ((crit < 1) || (crit > uf->n_crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	/* Deliver result */
	if (uf->df_list[crit])
		*exists = TRUE;
	else
		*exists = FALSE;
	return DTL_OK;
	}


/* For debug only. More general than DTL_PM_crit_exists above since
 * it accepts PS-frames as well, plus the MC-crit (crit 0) of PM. */

rcode DTLAPI DTI_crit_exists(int crit) {

	/* Log function call even though debug only - else confusing log */
	if (cst_ext) {
		sprintf(msg,"DTI_crit_exists(%d)\n",crit);
		cst_log(msg);
		}
	/* Check if function can start */
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Crude test, could look in df_list instead for PM */
	return !load_df0(crit);
	}
