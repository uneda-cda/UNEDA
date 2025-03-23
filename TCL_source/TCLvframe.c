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
 *   File: TCLvframe.c
 *
 *   Purpose: Modifications to V-base
 *
 *
 *   Functions exported outside TCL
 *   ------------------------------
 *   TCL_reset_V_base
 *   TCL_add_V_constraint
 *   TCL_replace_V_constraint
 *   TCL_change_V_constraint
 *   TCL_delete_V_constraint
 *   TCL_add_V_mstatement
 *   TCL_delete_V_mstatement
 *   TCL_set_V_box
 *   TCL_unset_V_box
 *   TCL_set_V_mbox
 *   TCL_get_V_box
 *   TCL_get_V_mbox
 *   TCL_get_V_index
 *   TCL_get_V_hull
 *   TCL_get_V_masspoint
 *
 *   Functions outside module, inside TCL
 *   ------------------------------------
 *   NONE
 *
 *   Functions internal to module
 *   ----------------------------
 *   NONE
 *
 */

#include "TCLinternal.h"

rcode TCL_reset_V_base(struct d_frame *df) {
	rcode rc;
	int i;
	struct base *V;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	V = df->V_base;

	/* Remove constraint set from base */
	V->n_stmts = 0;
	/* Remove range box intervals from base */
	V->box = FALSE;
	/* Remove midpoint box intervals from base */
	for (i=1; i<=tot_vars; i++)
		if (f2r[i])	{
			V->lo_midbox[f2r[i]] = -1.0;
			V->up_midbox[f2r[i]] = -1.0;
			}

	rc = TCL_OK;
	if (df->attached)
		rc = load_V(df);
	return rc;
	}


 /*********************************************************
  *
  *  Data input procedures, single constraint
  *
  *********************************************************/

rcode TCL_add_V_constraint(struct d_frame *df, struct stmt_rec *V_stmt) {
	rcode rc,rc2;
	struct base *V;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	V = df->V_base;
	if (V->n_stmts >= MAX_STMTS)
		return TCL_TOO_MANY_STMTS;

	/* Add the constraint last to the base */
	V->n_stmts++;
	memcpy(&(V->stmt[V->n_stmts]),V_stmt,sizeof(struct stmt_rec));

	rc = TCL_OK;
	if (df->attached) {
		/* Try to load new value base */
		rc = load_V(df);
		if (rc) {
			/* Failed to load, inconsistent */
			V->n_stmts--;
			rc2 = load_V(df);
			if (rc2)
				df->attached = FALSE;
			}
		}
	return rc;
	}


rcode TCL_replace_V_constraint(struct d_frame *df, int stmt_nbr, struct stmt_rec *V_stmt) {
	rcode rc,rc2;
	struct stmt_rec tmp;
	struct base *V;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;
	V = df->V_base;
	if ((stmt_nbr < 1) || (stmt_nbr > V->n_stmts))
		return TCL_INPUT_ERROR;

	/* Copy the constraint on top of the old one in the base */
	memcpy(&tmp,&(V->stmt[stmt_nbr]),sizeof(struct stmt_rec));
	memcpy(&(V->stmt[stmt_nbr]),V_stmt,sizeof(struct stmt_rec));

	/* Try to load new base */
	rc = load_V(df);
	if (rc) {
		/* Failed to load, inconsistent */
		memcpy(&(V->stmt[stmt_nbr]),&tmp,sizeof(struct stmt_rec));
		rc2 = load_V(df);
		if (rc2)
			df->attached = FALSE;
		}
	return rc;
	}


rcode TCL_change_V_constraint(struct d_frame *df, int stmt_nbr, double lobo, double upbo) {
	rcode rc,rc2;
	struct stmt_rec tmp;
	struct base *V;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;
	V = df->V_base;
	if ((stmt_nbr < 1) || (stmt_nbr > V->n_stmts))
		return TCL_INPUT_ERROR;
	rc = TCL_OK;

	/* Copy the bounds to the constraint */
	memcpy(&tmp,&(V->stmt[stmt_nbr]),sizeof(struct stmt_rec));
	V->stmt[stmt_nbr].lobo = lobo;
	V->stmt[stmt_nbr].upbo = upbo;

	/* Try to load new base */
	rc = load_V(df);
	if (rc) {
		/* Failed to load, inconsistent */
		memcpy(&(V->stmt[stmt_nbr]),&tmp,sizeof(struct stmt_rec));
		rc2 = load_V(df);
		if (rc2)
			df->attached = FALSE;
		}
	return rc;
	}


rcode TCL_delete_V_constraint(struct d_frame *df, int stmt_nbr) {
	rcode rc,rc2;
	int i;
	struct stmt_rec tmp;
	struct base *V;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	V = df->V_base;
	if ((stmt_nbr < 1) || (stmt_nbr > V->n_stmts))
		return TCL_INPUT_ERROR;

	/* Move constraints up */
	memcpy(&tmp,&(V->stmt[stmt_nbr]),sizeof(struct stmt_rec));
	for (i=stmt_nbr; i<V->n_stmts; i++)
		memcpy(&(V->stmt[i]),&(V->stmt[i+1]),sizeof(struct stmt_rec));
	V->n_stmts--;

	rc = TCL_OK;
	if (df->attached) {
		/* Try to attach new base */
		rc = load_V(df);
		if (rc) {
			/* Failed to load, inconsistent -> restore */
			for (i=V->n_stmts; i>=stmt_nbr; i--)
				memcpy(&(V->stmt[i+1]),&(V->stmt[i]),sizeof(struct stmt_rec));
			memcpy(&(V->stmt[stmt_nbr]),&tmp,sizeof(struct stmt_rec));
			V->n_stmts++;
			rc2 = load_V(df);
			if (rc2)
				df->attached = FALSE;
			}
		}
	return rc;
	}


rcode TCL_add_V_mstatement(struct d_frame *df, struct stmt_rec *V_stmt) {
	rcode rc,rc2;
	int index;
	struct base *V;
	double save_lo, save_up;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	V = df->V_base;

	if ((index = get_V_index(V_stmt->alt[1],V_stmt->cons[1])) == 0)
		return TCL_INPUT_ERROR; // im-node
	if ((V_stmt->lobo < 0.0) || (V_stmt->lobo > 1.0) || 
			(V_stmt->upbo < 0.0) || (V_stmt->upbo > 1.0))
		return TCL_INPUT_ERROR;
	save_lo = V->lo_midbox[index];
	save_up = V->up_midbox[index];
	V->lo_midbox[index] = V_stmt->lobo;
	V->up_midbox[index] = V_stmt->upbo;

	rc = TCL_OK;
	if (df->attached) {
		/* Try to load new base */
		rc = load_V(df);
		if (rc) {
			/* Failed to load, inconsistent */
			V->lo_midbox[index] = save_lo;
			V->up_midbox[index] = save_up;
			rc2 = load_V(df);
			if (rc2)
				df->attached = FALSE;
			}
		}
	return rc;
	}


rcode TCL_delete_V_mstatement(struct d_frame *df, struct stmt_rec *V_stmt) {
	rcode rc,rc2;
	int index;
	struct base *V;
	double save_lo, save_up;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	V = df->V_base;

	if ((index = get_V_index(V_stmt->alt[1],V_stmt->cons[1])) == 0)
		return TCL_INPUT_ERROR;
	if ((V->lo_midbox[index] == -1.0) && (V->up_midbox[index] == -1.0))
		return TCL_OK;
	save_lo = V->lo_midbox[index];
	save_up = V->up_midbox[index];
	V->lo_midbox[index] = -1.0;
	V->up_midbox[index] = -1.0;

	rc = TCL_OK;
	if (df->attached) {
		/* Try to load new base */
		rc = load_V(df);
		if (rc) {
			/* Failed to load, inconsistent */
			V->lo_midbox[index] = save_lo;
			V->up_midbox[index] = save_up;
			rc2 = load_V(df);
			if (rc2)
				df->attached = FALSE;
			}
		}
	return rc;
	}


 /*********************************************************
  *
  *  Data input procedures, entire box
  *
  *********************************************************/

rcode TCL_set_V_box(struct d_frame *df, d_row tbox_lobo, d_row tbox_upbo) {
	rcode rc,rc2;
	int i;
	struct base *V;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	for (i=1; i<=tot_vars; i++)
		if ((tbox_lobo[i] < 0.0) || (tbox_lobo[i] > 1.0) || 
				(tbox_upbo[i] < 0.0) || (tbox_upbo[i] > 1.0))
			return TCL_INPUT_ERROR;

	/* Add the box to the base */
	V = df->V_base;
	for (i=1; i<=tot_vars; i++)
		if (f2r[i])	{
			V->box_lobo[f2r[i]] = tbox_lobo[i];
			V->box_upbo[f2r[i]] = tbox_upbo[i];
			}
	V->box = TRUE;

	rc = TCL_OK;
	if (df->attached) {
		/* Try to load the new base */
		rc = load_V(df);
		if (rc) {
			/* Failed to load, inconsistent */
			V->box = FALSE;
			rc2 = load_V(df);
			if (rc2)
				df->attached = FALSE;
			}
		}
	return rc;
	}


rcode TCL_unset_V_box(struct d_frame *df) {
	rcode rc;
	struct base *V;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	V = df->V_base;
	if (!V->box)
		return TCL_OK;

	/* Remove the intervals from the base */
	V->box = FALSE;
	rc = TCL_OK;
	if (df->attached) {
		/* Try to load original base */
		rc = load_V(df);
		if (rc) {
			df->attached = FALSE;
			}
		}
	return rc;
	}


static d_row save_lobo,save_upbo;

rcode TCL_set_V_mbox(struct d_frame *df, d_row tmbox_lobo, d_row tmbox_upbo) {
	rcode rc,rc2;
	int i;
	struct base *V;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	for (i=1; i<=tot_vars; i++)
		if (f2r[i])
			/* Check within range [0,1], empty (-1) or unoccupied (-2) */
			if (((tmbox_lobo[i] < 0.0) || (tmbox_lobo[i] > 1.0) || 
					(tmbox_upbo[i] < 0.0) || (tmbox_upbo[i] > 1.0)) && 
					(tmbox_lobo[i] != -1.0) && (tmbox_lobo[i] != -2.0))
				return TCL_INPUT_ERROR;

	V = df->V_base;

	/* Add the box to the base */
	for (i=1; i<=tot_vars; i++)
		if (f2r[i])	{
			save_lobo[i] = V->lo_midbox[f2r[i]];
			save_upbo[i] = V->up_midbox[f2r[i]];
			if (tmbox_lobo[i] > -2.0) {
				V->lo_midbox[f2r[i]] = tmbox_lobo[i];
				V->up_midbox[f2r[i]] = tmbox_upbo[i];
				}
			}

	rc = TCL_OK;
	if (df->attached) {
		/* Try to load new base */
		rc = load_V(df);
		if (rc) {
			/* Failed to load <- inconsistent */
			for (i=1; i<=tot_vars; i++)
				if (f2r[i])	{
					V->lo_midbox[f2r[i]] = save_lobo[i];
					V->up_midbox[f2r[i]] = save_upbo[i];
					}
			rc2 = load_V(df);
			if (rc2)
				df->attached = FALSE;
			}
		}
	return rc;
	}


 /*********************************************************
  *
  *  Data output procedures
  *
  *********************************************************/

rcode TCL_get_V_box(struct d_frame *df, d_row lo_box, d_row up_box) {
	int i;
	struct base *V;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;

	/* If box=TRUE  then return combined stmt + box */
	/* If box=FALSE then return [0,1] i.e. w/o stmts */
	V = df->V_base;
	for (i=1; i<=df->tot_cons[0]; i++)
		if (f2r[i])	{
			if (V->box) {
				lo_box[i] = V->box_lobo[f2r[i]];
				up_box[i] = V->box_upbo[f2r[i]];
				}
			else {
				lo_box[i] = 0.0;
				up_box[i] = 1.0;
				}
			}
		else {
			lo_box[i] = -1.0;
			up_box[i] = -1.0;
			}
	return TCL_OK;
	}


rcode TCL_get_V_mbox(struct d_frame *df, d_row lo_mid, d_row up_mid) {
	int i;
	struct base *V;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;

	V = df->V_base;
	for (i=1; i<=df->tot_cons[0]; i++)
		if (f2r[i])	{
			lo_mid[i] = V->lo_midbox[f2r[i]];
			up_mid[i] = V->up_midbox[f2r[i]];
			}
		else	{
			lo_mid[i] = -1.0;
			up_mid[i] = -1.0;
			}
	return TCL_OK;
	}


int TCL_get_V_index(struct d_frame *df, int alt, int node) {

	/* Check input parameters */
	if (df == NULL)
		return 0;
	if (df->watermark != D_MARK)
		return 0;
	/* Return B2 index */
	return get_V_index(alt,node);
	}


/* Get hull for user */

rcode TCL_get_V_hull(struct d_frame *df, d_row lobo, d_row upbo) {

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;
	/* Return hull boundaries */
	hull_V(lobo,upbo);
	return TCL_OK;
	}


/* Mass point */

rcode TCL_get_V_masspoint(struct d_frame *df, d_row V_mid) {

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;

	cpoint_V(V_mid); // B1 indexing
	return TCL_OK;
	}
