/*
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
 *   Copyright (c) 2012-2025  Prof. Mats Danielson, Stockholm University
 *
 *   Website: https://people.dsv.su.se/~mad/UNEDA
 *   GitHub:  https://github.com/uneda-cda/UNEDA
 *
 *   Licensed under CC BY 4.0: https://creativecommons.org/licenses/by/4.0/.
 *   Provided "as is", without warranty of any kind, express or implied.
 *   Reuse and modifications are encouraged, with proper attribution.
 *
 */

/*
 *   File: file2.c
 *
 *   Purpose: reading files for user tester
 *
 */

#include "uct.h"

#define WARN_EST

extern int user_mode;
static t_matrix tnext,tdown;
static d_row xlobo,xupbo;

double vv1x(double val, double min, double max) {
	/* convert v -> x for range */
	return val*(max-min)+min;
	}


/* Read user frame from disk */

int read_ufile(char *fn, char *folder, struct user_frame *uf) {
	int i,j,rc,n_alts,n_stmts,frame_type;
	int n_cons[MAX_ALTS+1];
	FILE *fp;
	struct stmt_rec stmt;
	char fn_ddt[64];

	rc = TCL_OK;
	strcpy(fn_ddt,folder);
	strcat(fn_ddt,fn);
	strcat(fn_ddt,".ddt");
	fp = fopen(fn_ddt,"r");
	if (!fp) {
		return TCL_NO_FILE;
		}

	/* User frame header */
	fscanf(fp,"%d.%d.%d ",&uf_ver_main,&uf_ver_func,&uf_ver_tech);
	fscanf(fp,"%120s ",uf->frame_name);
	if (strcmp(fn,uf->frame_name))
		printf("Warning: name mismatch - file '%s' frame '%s'\n",fn,uf->frame_name);

	/* Read frame type if modern frame (assumed!) */
	fscanf(fp,"%d ",&frame_type);
	/* Old flag can be 0 */
	if (frame_type > 1) {
		fclose(fp);
		return TCL_CORRUPTED;
		}

	/* Read level type if multilevel release */
	if (uf_ver_main >= 4)
		fscanf(fp,"%d ",&uf->multilevel);
	else
		uf->multilevel = FALSE;

	/* UCT frame header */
	fscanf(fp,"%d ",&n_alts);
	if (n_alts < 2) {
		fclose(fp);
		return TCL_INCONSISTENT;
		}
	n_cons[0] = 0;
	for (i=1; i<=n_alts; i++)
		fscanf(fp,"%d ",n_cons+i);
	if (uf->multilevel) {
		for (i=1; i<=n_alts; i++) {
			for (j=1; j<=n_cons[i]; j++)
				fscanf(fp,"%d ",tnext[i]+j);
			for (j=1; j<=n_cons[i]; j++)
				fscanf(fp,"%d ",tdown[i]+j);
			}
		if (rc = TCL_create_tree_frame(&(uf->df),n_alts,n_cons,tnext,tdown)) {
			fclose(fp);
			return rc;
			}
		}
	else {
		if (rc = TCL_create_flat_frame(&(uf->df),n_alts,n_cons)) {
			fclose(fp);
			return rc;
			}
		}
	strcpy(uf->df->name,uf->frame_name);
	uf->n_crit = 1;

	/* Activate frame to load structure */
	if (rc = TCL_attach_frame(uf->df)) {
		fclose(fp);
		return rc;
		}

	/* Read names */
	if ((uf_ver_main > 3) || ((uf_ver_main == 3) && (uf_ver_func > 2))) {
		for (i=1; i<=n_alts; i++)
			fscanf(fp,"%29s ",uf->alt_name[i]);
		}
	else {
		for (i=1; i<=n_alts; i++)
			strcpy(uf->alt_name[i],"no_name");
		}
			
	/* Probability base */
	P_links[0].head = 0;
	fscanf(fp,"%d ",&n_stmts);
	if (n_stmts < 0) {
		fclose(fp);
		return TCL_INCONSISTENT;
		}
	for (i=1; i<=n_stmts; i++) {
		fscanf(fp,"%d ",&(stmt.n_terms));
		for (j=1; j<=stmt.n_terms; j++)
			fscanf(fp,"%d %d %d ",stmt.alt+j,stmt.cons+j,stmt.sign+j);
		fscanf(fp,"%lf ",&(stmt.lobo));
		fscanf(fp,"%lf ",&(stmt.upbo));
		if (stmt.n_terms == 1) {
			if (rc = TCL_add_P_constraint(uf->df,&stmt)) {
				fclose(fp);
				return rc;
				}
			}
		else { // unknown type
			fclose(fp);
			return TCL_INPUT_ERROR;
			}
		}

	/* Read value scale */
	if ((uf_ver_main > 3) || ((uf_ver_main == 3) && (uf_ver_func > 2))) {
		fscanf(fp,"%lf ",&(uf->v_lo));
		fscanf(fp,"%lf ",&(uf->v_up));
		}
	else {
		uf->v_lo = 0.0;
		uf->v_up = 1.0;
		}
	if ((uf->v_lo!=0.0) || (uf->v_up!=1.0))
		printf("Warning: value scale [%.3lf %.3lf]\n",uf->v_lo,uf->v_up);

	/* Value base */
	V_links[0].head = 0;
	fscanf(fp,"%d ",&n_stmts);
	if (n_stmts < 0) {
		fclose(fp);
		return TCL_INCONSISTENT;
		}
	for (i=1; i<=n_stmts; i++) {
		fscanf(fp,"%d ",&(stmt.n_terms));
		for (j=1; j<=stmt.n_terms; j++)
			fscanf(fp,"%d %d %d ",stmt.alt+j,stmt.cons+j,stmt.sign+j);
		fscanf(fp,"%lf ",&(stmt.lobo));
		fscanf(fp,"%lf ",&(stmt.upbo));
		if (stmt.n_terms == 1) {
			if (rc = TCL_add_V_constraint(uf->df,&stmt)) {
				fclose(fp);
				return rc;
				}
			}
		else { // unknown type
			fclose(fp);
			return TCL_INPUT_ERROR;
			}
		}

	if (uf_ver_main >= 3) {
		/* Probability midpoints */
		fscanf(fp,"%d ",&n_stmts);
		if (n_stmts < 0) {
			fclose(fp);
			return TCL_INCONSISTENT;
			}
		for (i=1; i<=n_stmts; i++) {
			fscanf(fp,"%d %d ",stmt.alt+1,stmt.cons+1);
			fscanf(fp,"%lf ",&(stmt.lobo));
			fscanf(fp,"%lf ",&(stmt.upbo));
			stmt.n_terms = 1;
			stmt.sign[1] = 1;
			if (rc = TCL_add_P_mstatement(uf->df,&stmt)) {
				fclose(fp);
				return rc;
				}
#ifdef WARN_EST
			if (stmt.upbo-stmt.lobo > EPS)
				printf("Interval core P%d.%d = [%.3lf %.3lf]\n",stmt.alt[1],stmt.cons[1],
								stmt.lobo,stmt.upbo);
#endif
			}
	
		/* Value midpoints */
		fscanf(fp,"%d ",&n_stmts);
		if (n_stmts < 0) {
			fclose(fp);
			return TCL_INCONSISTENT;
			}
		for (i=1; i<=n_stmts; i++) {
			fscanf(fp,"%d %d ",stmt.alt+1,stmt.cons+1);
			fscanf(fp,"%lf ",&(stmt.lobo));
			fscanf(fp,"%lf ",&(stmt.upbo));
			stmt.n_terms = 1;
			stmt.sign[1] = 1;
			if (rc = TCL_add_V_mstatement(uf->df,&stmt)) {
				fclose(fp);
				return rc;
				}
#ifdef WARN_EST
			if (stmt.upbo-stmt.lobo > EPS)
				printf("Interval core V%d.%d = [%.3lf %.3lf]\n",stmt.alt[1],stmt.cons[1],
								vv1x(stmt.lobo,uf->v_lo,uf->v_up),vv1x(stmt.upbo,uf->v_lo,uf->v_up));
#endif
			}
		}

	fclose(fp);
	return rc;
	}


/* Backup a user file if it does exist */

int backup_ufile(char *fn, char *folder) {
	char fn_ddt[64],fn_bkp[64];

	strcpy(fn_ddt,folder);
	strcat(fn_ddt,fn);
	strcat(fn_ddt,".ddt");
	strcpy(fn_bkp,folder);
	strcat(fn_bkp,fn);
	strcat(fn_bkp,".bkp");
	remove(fn_bkp);
	if (rename(fn_ddt,fn_bkp))
		return TCL_NO_FILE;
	return TCL_OK;
	}


/* Write user frame to disk */
static d_row lo_midpoint,up_midpoint;

int write_ufile(char *fn, char *folder, struct user_frame *uf) {
	int i,j,k,rc,n_midpoints;
	FILE *fp;
	struct d_frame *df;
	struct base *P,*V;
	char fn_ddt[64];

	rc = TCL_OK;
	backup_ufile(fn,folder);
	strcpy(fn_ddt,folder);
	strcat(fn_ddt,fn);
	strcat(fn_ddt,".ddt");
	fp = fopen(fn_ddt,"w"); // will overwrite if exists
	if (!fp)
		return TCL_NO_FILE;

	/* User frame header */
	fprintf(fp,"%d.%d.%d\n",DTL_MAIN,DTL_FUNC,DTL_TECH);
	fprintf(fp,"%s\n",uf->frame_name);

	/* Frame parameters */
	df = uf->df;
	fprintf(fp,"1\n"); // PS-type
	fprintf(fp,"%d \n",uf->multilevel);

	/* Decision frame alternatives */
	fprintf(fp,"%d ",df->n_alts);
	for (i=1; i<=df->n_alts; i++)
		fprintf(fp,"%d ",df->tot_cons[i]);
	fprintf(fp,"\n");

	/* Tree description */
	if (uf->multilevel) {
		for (i=1; i<=df->n_alts; i++) {
			for (j=1; j<=df->tot_cons[i]; j++)
				fprintf(fp,"%d ",uf->df->next[i][j]);
			fprintf(fp,"\n");
			for (j=1; j<=df->tot_cons[i]; j++)
				fprintf(fp,"%d ",uf->df->down[i][j]);
			fprintf(fp,"\n");
			}
		}

	/* Alt and crit names */
	for (i=1; i<=df->n_alts; i++)
		fprintf(fp,"%s\n",uf->alt_name[i]);

	/* Probability base */
	P = df->P_base;
	fprintf(fp,"%d\n",P->n_stmts);
	for (i=1; i<=P->n_stmts; i++) {
		fprintf(fp,"%d ",P->stmt[i].n_terms);
		for (j=1; j<=P->stmt[i].n_terms; j++)
			fprintf(fp,"%d %d %d ",P->stmt[i].alt[j],P->stmt[i].cons[j],P->stmt[i].sign[j]);
		fprintf(fp,"%.10le ",P->stmt[i].lobo);
		fprintf(fp,"%.10le\n",P->stmt[i].upbo);
		}

	/* Value base */
	V = df->V_base;
	fprintf(fp,"%.10le ",uf->v_lo);
	fprintf(fp,"%.10le\n",uf->v_up);
	fprintf(fp,"%d\n",V->n_stmts);
	for (i=1; i<=V->n_stmts; i++) {
		fprintf(fp,"%d ",V->stmt[i].n_terms);
		for (j=1; j<=V->stmt[i].n_terms; j++)
			fprintf(fp,"%d %d %d ",V->stmt[i].alt[j],V->stmt[i].cons[j],V->stmt[i].sign[j]);
		fprintf(fp,"%.10le ",V->stmt[i].lobo);
		fprintf(fp,"%.10le\n",V->stmt[i].upbo);
		}

	/* Probability midpoints */
	if (rc = TCL_get_P_mbox(uf->df,lo_midpoint,up_midpoint)) {
		fclose(fp);
		return rc;
		}

	/* Calculate number of midpoints */
	n_midpoints = 0;
	for (k=1,i=1; i<=uf->df->n_alts; i++)
		for (j=1; j<=uf->df->tot_cons[i]; j++,k++)
			if (lo_midpoint[k] >= 0.0)
				n_midpoints++;
	fprintf(fp,"%d\n",n_midpoints);

	/* Write midpoints in order */
	for (k=1,i=1; i<=uf->df->n_alts; i++)
		for (j=1; j<=uf->df->tot_cons[i]; j++,k++)
			if (lo_midpoint[k] >= 0.0) {
				fprintf(fp,"%d %d ",i,j);
				fprintf(fp,"%.10le ",lo_midpoint[k]);
				fprintf(fp,"%.10le\n",up_midpoint[k]);
				}

	/* Value midpoints */
	if (rc = TCL_get_V_mbox(uf->df,lo_midpoint,up_midpoint)) {
		fclose(fp);
		return rc;
		}

	/* Calculate number of midpoints */
	n_midpoints = 0;
	for (k=1,i=1; i<=uf->df->n_alts; i++)
		for (j=1; j<=uf->df->tot_cons[i]; j++,k++)
			if (lo_midpoint[k] >= 0.0)
				n_midpoints++;
	fprintf(fp,"%d\n",n_midpoints);

	/* Write midpoints in order */
	for (k=1,i=1; i<=uf->df->n_alts; i++)
		for (j=1; j<=uf->df->tot_cons[i]; j++,k++)
			if (lo_midpoint[k] >= 0.0) {
				fprintf(fp,"%d %d ",i,j);
				fprintf(fp,"%.10le ",lo_midpoint[k]);
				fprintf(fp,"%.10le\n",up_midpoint[k]);
				}
	fclose(fp);
	return rc;
	}
