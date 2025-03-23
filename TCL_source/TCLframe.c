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
 *   File: TCLframe.c
 *
 *   Purpose: TCL open layer: manage modifications of the frame
 *
 *
 *   Functions exported outside TCL
 *   ------------------------------
 *   TCL_create_flat_frame
 *   TCL_create_tree_frame
 *   TCL_dispose_frame
 *   TCL_attach_frame
 *   TCL_detach_frame
 *   TCL_get_real_index
 *   TCL_get_tot_index
 *   TCL_pure_tree
 *   TCL_different_parents
 *   TCL_nbr_of_siblings
 *
 *   Functions outside module, inside TCL
 *   ------------------------------------
 *   NONE
 *
 *   Functions internal to module
 *   ----------------------------
 *   tree_end
 *   lonely_im_child
 *   init_global_tree
 *   pure_node
 *
 */

#include "TCLinternal.h"


 /*********************************************************
  *
  *  Tree data structures in permanent storage
  *
  *********************************************************/

t_matrix t2f,t2r,t2i,r2t,i2t;
tn_row f2r,f2i,r2f,i2f,i2end;
int n_alts;
int alt_inx[MAX_ALTS+1];
int n_vars;
int im_alt_inx[MAX_ALTS+1];
int im_vars;
int tot_alt_inx[MAX_ALTS+1];
int tot_vars;


 /*********************************************************
  *
  *  Check tree structure and collect end node
  *
  *********************************************************/

static int tree_end(struct d_frame *df, int alt, int snode) {
	int tnode,end,inx;

	end = snode;
	for (tnode=df->down[alt][snode]; tnode; tnode=df->next[alt][tnode]) {
		if (tnode != end+1)
			/* Error in user tree */
			return 0;
		if (df->down[alt][tnode])	{
			/* Im-node */
			end = tree_end(df,alt,tnode);
			inx = at2i(alt,tnode);
			i2end[inx] = tot_alt_inx[alt-1] + end;
			}
		else
			end = tnode;
		}
	return end;
	}


/* Check tree to spot lonely im (intermediate) nodes
 * (children is defined here as "at least one child") */

static bool lonely_im_child(struct d_frame *df, int alt, int snode) {
	int tnode,t1;

	/* Check own first child node */
	if (t1=df->down[alt][snode]) // is this an im-node?
		/* A node is im and lonely if it has children but no siblings */
		if ((df->down[alt][t1]) && !df->next[alt][t1])
			return TRUE;
	/* Ask all children to check their children */
	for (tnode=t1; tnode; tnode=df->next[alt][tnode])
		if (df->down[alt][tnode])
			if (lonely_im_child(df,alt,tnode))
				return TRUE;
	return FALSE;
	}


 /*********************************************************
  *
  *  Initialise global tree structure
  *
  *********************************************************/

static rcode init_global_tree(struct d_frame *df) {
	int h,i,j,k1,k2,f1,f2;

	n_alts = df->n_alts;
	alt_inx[0] = 0;
	im_alt_inx[0] = 0;
	tot_alt_inx[0] = 0;
	for (i=1; i<=n_alts; i++)	{
		alt_inx[i] = alt_inx[i-1] + df->n_cons[i];
		im_alt_inx[i] = im_alt_inx[i-1] + df->im_cons[i];
		tot_alt_inx[i] = alt_inx[i] + im_alt_inx[i];
		}
	n_vars = alt_inx[n_alts];
	im_vars = im_alt_inx[n_alts];
	tot_vars = n_vars + im_vars;
	/* Map tree up */
	f1 = f2 = 1;
	for (h=1,i=1; i<=df->n_alts; i++)	{
		k1 = k2 = 1;
		for (j=1; j<=df->tot_cons[i]; j++,h++) {
			t2f[i][j] = h;
			if (df->down[i][j]) {
				/* Im-node */
				t2r[i][j] = 0;
				t2i[i][j] = k2;
				i2t[i][k2++] = j;
				f2r[h] = 0;
				f2i[h] = f2;
				i2f[f2++] = h;
				}
			else {
				/* Re-node */
				t2r[i][j] = k1;
				t2i[i][j] = 0;
				r2t[i][k1++] = j;
				f2r[h] = f1;
				f2i[h] = 0;
				r2f[f1++] = h;
				}
			}
		if ((df->n_cons[i] != k1-1) || (df->im_cons[i] != k2-1))
			return TCL_TREE_ERROR;
		if (tree_end(df,i,0) != df->tot_cons[i])
			return TCL_TREE_ERROR;
		}
	return TCL_OK;
	}


 /*********************************************************
  *
  *  Create or destruct data frame
  *
  *********************************************************/

/* Flat frame */

rcode TCL_create_flat_frame(struct d_frame **dfp, int n_alts, int n_cons[]) {
	rcode rc;
	int i,j;

	/* Check input parameters */
	if (n_alts < 2)
		return TCL_TOO_FEW_ALTS;
	if (n_alts > MAX_ALTS)
		return TCL_TOO_MANY_ALTS;
	n_cons[0] = 0;
	for (i=1; i<=n_alts; i++) {
		if (n_cons[i] < 1)
			return TCL_INPUT_ERROR;
		if (n_cons[i] > MAX_COPA)
			return TCL_TOO_MANY_CONS;
		n_cons[0] += n_cons[i];
		}
	if (n_cons[0] > MAX_CONS)
		return TCL_TOO_MANY_CONS;

	/* Allocate a decision frame */
	*dfp = (struct d_frame *)mem_alloc(sizeof(struct d_frame),"struct d_frame","TCL_create_flat_frame");
	if (!*dfp)
		return TCL_OUT_OF_MEMORY;
	(*dfp)->watermark = D_MARK;
	(*dfp)->name[0] = '\0';
	(*dfp)->tree = FALSE;
	(*dfp)->n_alts = n_alts;
	for (i=0; i<=n_alts; i++)	{
		(*dfp)->n_cons[i] = (*dfp)->tot_cons[i] = n_cons[i];
		(*dfp)->im_cons[i] = 0;
		}
	for (i=n_alts+1; i<=MAX_ALTS; i++)
		(*dfp)->n_cons[i] = (*dfp)->tot_cons[i] = (*dfp)->im_cons[i] = 0;
	/* Create next and down pointers plus reverse pointers. */
	for (i=1; i<=n_alts; i++)	{
		(*dfp)->down[i][0] = 1;
		(*dfp)->up[i][0] = 0;
		(*dfp)->next[i][0] = 0;
		(*dfp)->prev[i][0] = 0;
		for (j=1; j<=n_cons[i]; j++)	{
			(*dfp)->down[i][j] = 0;
			(*dfp)->up[i][j] = 0;
			(*dfp)->next[i][j] = ((j<n_cons[i])?j+1:0);
			(*dfp)->prev[i][j] = j-1;
			}
		}
	(*dfp)->attached = FALSE;

	/* Allocate bases in memory */
	if (rc = create_P(*dfp))
		return rc;
	if (rc = create_V(*dfp))
		return rc;
	return TCL_OK;
	}


/* Create tree frame */

rcode TCL_create_tree_frame(struct d_frame **dfp, int n_alts, int tot_cons[], 
				t_matrix next, t_matrix down) {
	rcode rc;
	int i,j,max_next,re_cons[MAX_ALTS+1],im_cons[MAX_ALTS+1];

	/* Check input parameters */
	if (n_alts < 2)
		return TCL_INPUT_ERROR;
	if (n_alts > MAX_ALTS)
		return TCL_TOO_MANY_ALTS;
	/* Find nbr of real and intermediate cons */
	for (i=1; i<=n_alts; i++) {
		if (tot_cons[i] > MAX_NOPA)
			return TCL_INPUT_ERROR;
		im_cons[i] = 0;
		max_next = 1; /* = down[i][0] implicit */
		for (j=1; j<=tot_cons[i]; j++)	{
			if (down[i][j])
				im_cons[i]++;
			if (next[i][j] > max_next)
				max_next = next[i][j];
			}
		if (tot_cons[i] != max_next)
			return TCL_TREE_ERROR;
		re_cons[i] = tot_cons[i] - im_cons[i];
		}
	re_cons[0] = 0;
	im_cons[0] = 0;
	for (i=1; i<=n_alts; i++) {
		if (re_cons[i] < 1)
			return TCL_INPUT_ERROR;
		if (re_cons[i] > MAX_COPA)
			return TCL_TOO_MANY_CONS;
		re_cons[0] += re_cons[i];
		im_cons[0] += im_cons[i];
		}
	if (re_cons[0] > MAX_CONS)
		return TCL_TOO_MANY_CONS;
	tot_cons[0] = re_cons[0] + im_cons[0];

	/* Allocate a data frame */
	*dfp = (struct d_frame *)mem_alloc(sizeof(struct d_frame),"struct d_frame","TCL_create_tree_frame");
	if (!*dfp)
		return TCL_OUT_OF_MEMORY;
	(*dfp)->watermark = D_MARK;
	(*dfp)->name[0] = '\0';
	(*dfp)->tree = (im_cons[0] > 0);
	(*dfp)->n_alts = n_alts;
	for (i=0; i<=n_alts; i++)	{
		(*dfp)->n_cons[i] = re_cons[i];
		(*dfp)->im_cons[i] = im_cons[i];
		(*dfp)->tot_cons[i] = tot_cons[i];
		}
	for (i=n_alts+1; i<=MAX_ALTS; i++) {
		(*dfp)->n_cons[i] = 0;
		(*dfp)->im_cons[i] = 0;
		}

	/* Copy next and down pointers. Create reverse pointers. */
	/* If prev==0 then up>0, except for top node where up==0 */
	for (i=1; i<=n_alts; i++)	{
		(*dfp)->down[i][0] = 1;
		(*dfp)->up[i][0] = 0;
		(*dfp)->next[i][0] = 0;
		(*dfp)->prev[i][0] = 0;
		for (j=1; j<=tot_cons[i]; j++)	{
			(*dfp)->prev[i][j] = 0;
			(*dfp)->up[i][j] = 0;
			}
		for (j=1; j<=tot_cons[i]; j++)	{
			(*dfp)->next[i][j] = next[i][j];
			if (next[i][j])
				(*dfp)->prev[i][next[i][j]] = j;
			(*dfp)->down[i][j] = down[i][j];
			if (down[i][j])
				(*dfp)->up[i][down[i][j]] = j;
			}
		/* Im-node must have >1 children (except implicit node 0) */
		if (lonely_im_child(*dfp,i,0)) {
			(*dfp)->watermark = 0;
			mem_free((void *)(*dfp));
			return TCL_TREE_ERROR;
			}
		}
	(*dfp)->attached = FALSE;

	/* Create bases */
	if (rc = create_P(*dfp))
		return rc;
	if (rc = create_V(*dfp))
		return rc;
	return TCL_OK;
	}


/* Delete frame */

rcode TCL_dispose_frame(struct d_frame *df) {
	rcode rc,rc2;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	/* Release base memory */
	rc  = dispose_P(df);
	rc2 = dispose_V(df);
	if (rc+rc2)
		return max(rc,rc2);
	df->watermark = 0;
	/* Release own memory */
	return mem_free((void *)df);
	}


 /*********************************************************
  *
  *  Attach/detach data frame
  *
  *********************************************************/

/* Load frame and bases */

rcode TCL_attach_frame(struct d_frame *df) {
	rcode rc;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (df->attached)
		return TCL_ATTACHED;
	/* Load all internal data structures */
	if (rc = init_global_tree(df))
		return rc;
	if (rc = load_P(df))
		return rc;
	if (rc = load_V(df))
		return rc;
	df->attached = TRUE;
	return TCL_OK;
	}

/* Unload frame */

rcode TCL_detach_frame(struct d_frame *df) {

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	/* Mark data structures invalid */
	df->attached = FALSE;
	return TCL_OK;
	}


 /*********************************************************
  *
  *  Index conversion functions between indexing types
  *
  *********************************************************/

/* Get real index (B2) from alt and cons (A2) */

int TCL_get_real_index(int alt, int cons) {

	if ((alt<1) || (alt>n_alts))
		return 0;
	if ((cons<1) || (cons>tot_alt_inx[alt]-tot_alt_inx[alt-1]))
		return 0;
	return tot_alt_inx[alt-1]+cons;
	}


/* Get total index (B1) from alt and cons (A2) */

int TCL_get_tot_index(int alt, int cons) {

	/* Check input parameters */
	if ((alt<1) || (alt>n_alts))
		return 0;
	if ((cons<1) || (cons>tot_alt_inx[alt]-tot_alt_inx[alt-1]))
		return 0;
	/* Return B1 index */
	return tot_alt_inx[alt-1]+r2t[alt][cons];
	}


 /*********************************************************
  *
  *  Property checks
  *
  *********************************************************/

/* Check if tree is pure or not. Pure is defined as all nodes with
 * the same parent being either re-nodes or im-nodes, not a mixture. */

static int pure_node(struct d_frame *df, int alt, int snode) {
	int tnode,re=0,im=0;

	/* Check all children of this node */
	for (tnode=df->down[alt][snode]; tnode; tnode=df->next[alt][tnode])
		if (df->down[alt][tnode]) {
			/* Intermediate node */
			if (!pure_node(df,alt,tnode))
				return FALSE;
			im++;
			}
		else
			/* Real/final node */
			re++;
	return !re || !im;
	}


/* Returns -1 if pure, 0 if impure, >0 if error */

int TCL_pure_tree(struct d_frame *df, int alt) {

	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;
	if ((alt < 1) || (alt > df->n_alts))
		return TCL_INPUT_ERROR;
	return -1*pure_node(df,alt,0);
	}


/* Check if two tree nodes have the same parent. Returns
 * 0 if same parent, -1 if different parents, >0 if error */

int TCL_different_parents(struct d_frame *df, int alt, int node1, int node2) {
	int tc,tc1,tc2;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;
	if ((alt < 1) || (alt > df->n_alts))
		return TCL_INPUT_ERROR;
	if ((node1 < 1) || (node1 > df->tot_cons[alt]))
		return TCL_INPUT_ERROR;
	if ((node2 < 1) || (node2 > df->tot_cons[alt]))
		return TCL_INPUT_ERROR;
	if (node1 == node2)
		return 0;
	/* Same level if higher is on lower's tail */
	tc1 = min(node1,node2);
	tc2 = max(node1,node2);
	for (tc=df->next[alt][tc1]; tc<=tc2; tc=df->next[alt][tc]) {
		if (!tc)
			return -1;
		if (tc==tc2)
			return 0;
		}
	return -1;
	}


/* Count nbr of siblings incl. self: <0 (-#) if ok, >0 if error.
 * Note: never returns 0, do not use boolean conditional. */

int TCL_nbr_of_siblings(struct d_frame *df, int alt, int node) {
	int tc,n_nodes=1;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;
	if ((alt < 1) || (alt > df->n_alts))
		return TCL_INPUT_ERROR;
	if ((node < 1) || (node > df->tot_cons[alt]))
		return TCL_INPUT_ERROR;
	/* Look at head and tail */
	for (tc=df->prev[alt][node]; tc; tc=df->prev[alt][tc])
		n_nodes++;
	for (tc=df->next[alt][node]; tc; tc=df->next[alt][tc])
		n_nodes++;
	return -n_nodes;
	}
