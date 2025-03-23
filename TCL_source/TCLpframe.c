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
 *   File: TCLpframe.c
 *
 *   Purpose: manage modifications to the P-base
 *
 *
 *   Functions exported outside TCL
 *   ------------------------------
 *   TCL_reset_P_base
 *   TCL_add_P_constraint
 *   TCL_replace_P_constraint
 *   TCL_change_P_constraint
 *   TCL_delete_P_constraint
 *   TCL_add_P_mstatement
 *   TCL_delete_P_mstatement
 *   TCL_set_P_box
 *   TCL_unset_P_box
 *   TCL_set_P_mbox
 *   TCL_get_P_box
 *   TCL_get_P_mbox
 *   TCL_get_P_index
 *   TCL_get_P_hull
 *   TCL_get_P_masspoint
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

rcode TCL_reset_P_base(struct d_frame *df) {
	rcode rc;
	int i;
	struct base *P;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	P = df->P_base;

	/* Remove constraint set from base */
	P->n_stmts = 0;
	/* Remove intervals from base */
	P->box = FALSE;
	/* Remove midpoints from base */
	for (i=1; i<=tot_vars; i++)
		if (f2r[i])	{
			P->lo_midbox[f2r[i]] = -1.0;
			P->up_midbox[f2r[i]] = -1.0;
			}
		else {
			P->lo_im_midbox[f2i[i]] = -1.0;
			P->up_im_midbox[f2i[i]] = -1.0;
			}

	rc = TCL_OK;
	if (df->attached) {
		rc = load_P(df);
		}
	return rc;
	}


 /*********************************************************
  *
  *  Data input procedures, single constraint
  *
  *********************************************************/

rcode TCL_add_P_constraint(struct d_frame *df, struct stmt_rec *P_stmt) {
	rcode rc,rc2;
	struct base *P;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	P = df->P_base;
	if (P->n_stmts >= MAX_STMTS)
		return TCL_TOO_MANY_STMTS;

	/* Add the constraint last to the base */
	P->n_stmts++;
	memcpy(&(P->stmt[P->n_stmts]),P_stmt,sizeof(struct stmt_rec));

	rc = TCL_OK;
	if (df->attached) {
		/* Try to load new base */
		rc = load_P(df);
		if (rc) {
			/* Failed to load, inconsistent */
			P->n_stmts--;
			rc2 = load_P(df);
			if (rc2)
				df->attached = FALSE;
			}
		}
	return rc;
	}


rcode TCL_replace_P_constraint(struct d_frame *df, int stmt_nbr, struct stmt_rec *P_stmt) {
	rcode rc,rc2;
	struct stmt_rec tmp;
	struct base *P;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;
	P = df->P_base;
	if ((stmt_nbr < 1) || (stmt_nbr > P->n_stmts))
		return TCL_INPUT_ERROR;

	/* Not destructive, can be undone */
	memcpy(&tmp,&(P->stmt[stmt_nbr]),sizeof(struct stmt_rec));

	/* Copy the constraint on top of the old one in the base */
	memcpy(&(P->stmt[stmt_nbr]),P_stmt,sizeof(struct stmt_rec));

	/* Try to load new base */
	rc = load_P(df);
	if (rc) {
		/* Failed to load, inconsistent */
		memcpy(&(P->stmt[stmt_nbr]),&tmp,sizeof(struct stmt_rec));
		rc2 = load_P(df);
		if (rc2)
			df->attached = FALSE;
		}
	return rc;
	}


rcode TCL_change_P_constraint(struct d_frame *df, int stmt_nbr, double lobo, double upbo) {
	rcode rc,rc2;
	struct stmt_rec tmp;
	struct base *P;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;
	P = df->P_base;
	if ((stmt_nbr < 1) || (stmt_nbr > P->n_stmts))
		return TCL_INPUT_ERROR;
	rc = TCL_OK;

	/* Copy the bounds to the constraint */
	memcpy(&tmp,&(P->stmt[stmt_nbr]),sizeof(struct stmt_rec));
	P->stmt[stmt_nbr].lobo = lobo;
	P->stmt[stmt_nbr].upbo = upbo;

	/* Try to load new base */
	rc = load_P(df);
	if (rc) {
		/* Failed to load, inconsistent */
		memcpy(&(P->stmt[stmt_nbr]),&tmp,sizeof(struct stmt_rec));
		rc2 = load_P(df);
		if (rc2)
			df->attached = FALSE;
		}
	return rc;
	}


rcode TCL_delete_P_constraint(struct d_frame *df, int stmt_nbr) {
	rcode rc,rc2;
	int i;
	struct stmt_rec tmp;
	struct base *P;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;

	P = df->P_base;
	if ((stmt_nbr < 1) || (stmt_nbr > P->n_stmts))
		return TCL_INPUT_ERROR;

	/* Move constraints up */
	memcpy(&tmp,&(P->stmt[stmt_nbr]),sizeof(struct stmt_rec));
	for (i=stmt_nbr; i<P->n_stmts; i++)
		memcpy(&(P->stmt[i]),&(P->stmt[i+1]),sizeof(struct stmt_rec));
	P->n_stmts--;

	rc = TCL_OK;
	if (df->attached) {
		/* Try to attach new base */
		rc = load_P(df);
		if (rc) {
			/* Failed to load, inconsistent -> restore */
			for (i=P->n_stmts; i>=stmt_nbr; i--)
				memcpy(&(P->stmt[i+1]),&(P->stmt[i]),sizeof(struct stmt_rec));
			memcpy(&(P->stmt[stmt_nbr]),&tmp,sizeof(struct stmt_rec));
			P->n_stmts++;
			rc2 = load_P(df);
			if (rc2)
				df->attached = FALSE;
			}
		}

	return rc;
	}


rcode TCL_add_P_mstatement(struct d_frame *df, struct stmt_rec *P_stmt) {
	rcode rc,rc2;
	int index;
	struct base *P;
	double save_lo, save_up;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	P = df->P_base;

	if ((P_stmt->lobo < 0.0) || (P_stmt->lobo > 1.0) || 
			(P_stmt->upbo < 0.0) || (P_stmt->upbo > 1.0))
		return TCL_INPUT_ERROR;
	if (df->down[P_stmt->alt[1]][P_stmt->cons[1]]) {
		/* Intermediate node */
		if ((index = get_P_im_index(P_stmt->alt[1],P_stmt->cons[1])) == 0)
			return TCL_INPUT_ERROR;
		save_lo = P->lo_im_midbox[index];
		save_up = P->up_im_midbox[index];
		P->lo_im_midbox[index] = P_stmt->lobo;
		P->up_im_midbox[index] = P_stmt->upbo;
		}
	else	{
		/* Real node */
		if ((index = get_P_index(P_stmt->alt[1],P_stmt->cons[1])) == 0)
			return TCL_INPUT_ERROR;
		save_lo = P->lo_midbox[index];
		save_up = P->up_midbox[index];
		P->lo_midbox[index] = P_stmt->lobo;
		P->up_midbox[index] = P_stmt->upbo;
		}

	rc = TCL_OK;
	if (df->attached) {
		/* Try to load new base */
		rc = load_P(df);
		if (rc) {
			/* Failed to load <- inconsistent */
			if (df->down[P_stmt->alt[1]][P_stmt->cons[1]]) {
				/* Intermediate */
				P->lo_im_midbox[index] = save_lo;
				P->up_im_midbox[index] = save_up;
				}
			else	{
				/* Real */
				P->lo_midbox[index] = save_lo;
				P->up_midbox[index] = save_up;
				}
			rc2 = load_P(df);
			if (rc2)
				df->attached = FALSE;
			}
		}
	return rc;
	}


rcode TCL_delete_P_mstatement(struct d_frame *df, struct stmt_rec *P_stmt) {
	rcode rc,rc2;
	int index;
	struct base *P;
	double save_lo, save_up;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	P = df->P_base;

	if (df->down[P_stmt->alt[1]][P_stmt->cons[1]])	{
		if ((index = get_P_im_index(P_stmt->alt[1],P_stmt->cons[1])) == 0)
			return TCL_INPUT_ERROR;
		if ((P->lo_im_midbox[index] == -1.0) || (P->up_im_midbox[index] == -1.0))
			return TCL_INPUT_ERROR;
		save_lo = P->lo_im_midbox[index];
		save_up = P->up_im_midbox[index];
		P->lo_im_midbox[index] = -1.0;
		P->up_im_midbox[index] = -1.0;
		}
	else	{
		if ((index = get_P_index(P_stmt->alt[1],P_stmt->cons[1])) == 0)
			return TCL_INPUT_ERROR;
		if ((P->lo_midbox[index] == -1.0) && (P->up_midbox[index] == -1.0))
			return TCL_OK;
		save_lo = P->lo_midbox[index];
		save_up = P->up_midbox[index];
		P->lo_midbox[index] = -1.0;
		P->up_midbox[index] = -1.0;
		}

	rc = TCL_OK;
	if (df->attached) {
		/* Try to load new base */
		rc = load_P(df);
		if (rc) {
			/* Failed to load, inconsistent */
			if (df->down[P_stmt->alt[1]][P_stmt->cons[1]]) {
				/* Intermediate */
				P->lo_im_midbox[index] = save_lo;
				P->up_im_midbox[index] = save_up;
				}
			else	{
				/* Real */
				P->lo_midbox[index] = save_lo;
				P->up_midbox[index] = save_up;
				}
			rc2 = load_P(df);
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

int TCL_set_P_box(struct d_frame *df, d_row tbox_lobo, d_row tbox_upbo) {
	rcode rc,rc2;
	int i;
	struct base *P;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	for (i=1; i<=tot_vars; i++)
		if ((tbox_lobo[i] < 0.0) || (tbox_lobo[i] > 1.0) || 
				(tbox_upbo[i] < 0.0) || (tbox_upbo[i] > 1.0))
			return TCL_INPUT_ERROR;

	/* Split and add the box to the base */
	P = df->P_base;
	for (i=1; i<=tot_vars; i++)
		if (f2r[i])	{
			P->box_lobo[f2r[i]] = tbox_lobo[i];
			P->box_upbo[f2r[i]] = tbox_upbo[i];
			}
		else {
			P->im_box_lobo[f2i[i]] = tbox_lobo[i];
			P->im_box_upbo[f2i[i]] = tbox_upbo[i];
			}
	P->box = TRUE;

	rc = TCL_OK;
	if (df->attached) {
		/* Try to load new base */
		rc = load_P(df);
		if (rc) {
			/* Failed to load, inconsistent */
			P->box = FALSE;
			rc2 = load_P(df);
			if (rc2)
				df->attached = FALSE;
			}
		}
	return rc;
	}


rcode TCL_unset_P_box(struct d_frame *df) {
	rcode rc;
	struct base *P;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	P = df->P_base;
	if (!P->box)
		return TCL_OK;

	/* Remove the intervals from the base */
	P->box = FALSE;
	rc = TCL_OK;
	if (df->attached) {
		/* Try to load original base */
		rc = load_P(df);
		if (rc) {
			df->attached = FALSE;
			}
		}
	return rc;
	}


static d_row save_lobo,save_upbo;

rcode TCL_set_P_mbox(struct d_frame *df, d_row tmbox_lobo, d_row tmbox_upbo) {
	rcode rc,rc2;
	int i;
	struct base *P;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	for (i=1; i<=tot_vars; i++)
			/* Check within range [0,1], empty (-1) or unoccupied (-2) */
		if (((tmbox_lobo[i] < 0.0) || (tmbox_lobo[i] > 1.0) || 
				(tmbox_upbo[i] < 0.0) || (tmbox_upbo[i] > 1.0)) && 
				(tmbox_lobo[i] != -1.0) && (tmbox_lobo[i] != -2.0))
			return TCL_INPUT_ERROR;

	P = df->P_base;

	/* Split and add the box to the base */
	for (i=1; i<=tot_vars; i++)
		if (f2r[i])	{
			save_lobo[i] = P->lo_midbox[f2r[i]];
			save_upbo[i] = P->up_midbox[f2r[i]];
			if (tmbox_lobo[i] > -2.0) {
				P->lo_midbox[f2r[i]] = tmbox_lobo[i];
				P->up_midbox[f2r[i]] = tmbox_upbo[i];
				}
			}
		else {
			save_lobo[i] = P->lo_im_midbox[f2i[i]];
			save_upbo[i] = P->up_im_midbox[f2i[i]];
			if (tmbox_lobo[i] > -2.0) {
				P->lo_im_midbox[f2i[i]] = tmbox_lobo[i];
				P->up_im_midbox[f2i[i]] = tmbox_upbo[i];
				}
			}

	rc = TCL_OK;
	if (df->attached) {
		/* Try to load new base */
		rc = load_P(df);
		if (rc) {
			/* Failed to load, inconsistent */
			for (i=1; i<=tot_vars; i++)
				if (f2r[i])	{
					P->lo_midbox[f2r[i]] = save_lobo[i];
					P->up_midbox[f2r[i]] = save_upbo[i];
					}
				else {
					P->lo_im_midbox[f2i[i]] = save_lobo[i];
					P->up_im_midbox[f2i[i]] = save_upbo[i];
					}
			rc2 = load_P(df);
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

rcode TCL_get_P_box(struct d_frame *df, d_row lo_box, d_row up_box) {
	int i;
	struct base *P;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;

	/* If box=TRUE  then return combined stmt + box */
	/* If box=FALSE then return [0,1] i.e. w/o stmts */
	P = df->P_base;
	if (P->box) {
		for (i=1; i<=df->tot_cons[0]; i++)
			if (f2r[i])	{
				lo_box[i] = P->box_lobo[f2r[i]];
				up_box[i] = P->box_upbo[f2r[i]];
				}
			else	{
				lo_box[i] = P->im_box_lobo[f2i[i]];
				up_box[i] = P->im_box_upbo[f2i[i]];
				}
			}
	else {
		for (i=1; i<=df->tot_cons[0]; i++) {
			lo_box[i] = 0.0;
			up_box[i] = 1.0;
			}
		}
	return TCL_OK;
	}


rcode TCL_get_P_mbox(struct d_frame *df, d_row lo_mid, d_row up_mid) {
	int i;
	struct base *P;

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;

	P = df->P_base;
	for (i=1; i<=df->tot_cons[0]; i++)
		if (f2r[i])	{
			lo_mid[i] = P->lo_midbox[f2r[i]];
			up_mid[i] = P->up_midbox[f2r[i]];
			}
		else	{
			lo_mid[i] = P->lo_im_midbox[f2i[i]];
			up_mid[i] = P->up_im_midbox[f2i[i]];
			}
	return TCL_OK;
	}


int TCL_get_P_index(struct d_frame *df, int alt, int node) {

	/* Check input parameters */
	if (df == NULL)
		return 0;
	if (df->watermark != D_MARK)
		return 0;
	/* Return B2 index */
	return get_P_index(alt,node);
	}


/* Global and local hull */

rcode TCL_get_P_hull(struct d_frame *df, d_row hlobo, d_row hupbo, d_row llobo, d_row lupbo) {

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;
	/* Return orthogonal hull */
	hull_P(hlobo,hupbo);
	l_hull_P(llobo,lupbo);
	return TCL_OK;
	}


/* Mass point */

rcode TCL_get_P_masspoint(struct d_frame *df, d_row P_mid, d_row P_lmid) {

	/* Check input parameters */
	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;

	cpoint_P(P_mid);
	l_cpoint_P(P_lmid);
	return TCL_OK;
	}
