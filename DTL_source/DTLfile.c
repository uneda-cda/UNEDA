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
 *   File: DTLfile.c
 *
 *   Purpose: reading and writing files from the user interface
 *
 *
 *   Functions exported outside DTL
 *   ------------------------------
 *   DTL_read_frame
 *   DTL_write_frame
 *
 *   Functions outside of module, inside DTL
 *   ---------------------------------------
 *   NONE
 *
 *   Functions internal to module
 *   ----------------------------
 *   read_ufile
 *   read_dfile
 *   dtl_read_file
 *   backup_ufile
 *   nospace
 *   write_ufile
 *   write_dfile
 *   dtl_write_file
 *
 */

#include "DTL.h"
#include "DTLinternal.h"


 /********************************************************
  *
  *  Configuration parameters
  *
  ********************************************************/

#ifdef CONSOLE
// Warning if criterion 0 (MC) contains a value base with entries
#define V_CRIT0
// Warning if a criterion has probability or value midpoints in file
#define WARN_MIDPT
#define noWARN_MIDPT_EXT
// Warning if MC weight statement use alt>1 as carrier - ineffectual
#define WARN_MC
// Warning if MC pseudo-alternative is too big
#define WARN_ERR
// Warning if links encountered (not supported in basic UNEDA)
#define WARN_LINK
#endif

#define FOSIZE 256 // folder name size


 /********************************************************
  *
  *  Internal data
  *
  ********************************************************/

static t_matrix tnext,tdown;
static int uf_dtl_main,uf_dtl_func;
static int links_skipped;


 /*********************************************************
  *
  *  File format
  *
  *********************************************************/

 /*
  * .dmc file format for PM, PS, DM file types
  *
  * ufile structure
  * ---------------
  * main.func
  * frame_name
  * type (1=PS, 2=DM, 3=PM)
  * #alts #crit
  * {critmap} (one per crit incl. crit0/wt: 1=crit present 0=not present)
  * {dfile structure} (one per crit present incl. crit0/wt)
  *
  * dfile structure
  * ---------------
  * crit_name
  * #alts {nodes} (one per alt)
  * tree (0=flat 1=multilevel)
  * if (multilevel)
  *   tree structure [see below]
  * #p_base-stmts
  * {p_base-stmt} (one row per stmt)
  * #v_base-stmts
  * {v_base-stmt} (one row per stmt)
  * if (main.func > 6.04)
  *   #p_box-stmts
  *   {p_box-stmt} (one row per stmt)
  *   #v_box-stmts
  *   {v_box-stmt} (one row per stmt)
  * #p_est-stmts
  * {p_est-stmt} (one row per stmt)
  * #v_est-stmts
  * {v_est-stmt} (one row per stmt)
  *
  * tree structure
  * --------------
  * {tnext
  *  tdown} (one pair per alt, one number per node in each row)
  *
  * NOTE: the error handling is much more limited since this is an
  * internal function in the sense that it only reads its own data.
  *
  */


 /*********************************************************
  *
  *  Read user frame from file
  *
  *********************************************************/

static struct d_frame *read_dfile(FILE *fp, int crit);

static rcode read_ufile(char *fn, char *folder, struct user_frame *uframe) {
	int i,i2;
	int crit_map[MAX_CRIT+1];
	FILE *fp;
	char fn_dmc[FOSIZE+FNSIZE+6];

	/* Check input parameters */
	strcpy(fn_dmc,folder);
	strcat(fn_dmc,fn);
	strcat(fn_dmc,".dmc");
	fp = fopen(fn_dmc,"r");
	if (!fp)
		return DTL_FILE_UNKNOWN;
	/* User frame header */
	fscanf(fp,"%d.%d ",&uf_dtl_main,&uf_dtl_func);
	fscanf(fp,"%120s ",uframe->frame_name);
	/* Read frame type */
	fscanf(fp,"%d ",&uframe->frame_type);
	/* Frame header */
	fscanf(fp,"%d ",&uframe->n_alts);
	if ((uframe->n_alts < 2) || (uframe->n_alts > MAX_ALTS)) {
		fclose(fp);
		return DTL_ALT_OVERFLOW;
		}
	fscanf(fp,"%d ",&uframe->n_crit);
	if ((uframe->n_crit < 1) || (uframe->n_crit > MAX_CRIT)) {
		fclose(fp);
		return DTL_CRIT_OVERFLOW;
		}
	if (uframe->frame_type == PS_FRAME)
		if (uframe->n_crit != 1) {
			fclose(fp);
			return DTL_FRAME_CORRUPT;
			}
	uframe->n_sh = 1;
	if (uframe->frame_type == PM_FRAME) {
		for (i=0; i<=uframe->n_crit; i++)
			fscanf(fp,"%d ",crit_map+i);
		if (!crit_map[0]) { // must have MC frame
			fclose(fp);
			return DTL_FRAME_CORRUPT;
			}
		for (i=0; i<=uframe->n_crit; i++)
			if (crit_map[i]) {
				if (!(uframe->df_list[i] = read_dfile(fp,i))) {
					// Release what's been collected
					for (i2=0; i2<i; i2++)
						if (crit_map[i2])
							TCL_dispose_frame(uframe->df_list[i2]);
					fclose(fp);
					return DTL_FRAME_CORRUPT;
					}
				}
			else
				uframe->df_list[i] = NULL;
		}
	else { // PS,DM
		if (!(uframe->df = read_dfile(fp,1))) {
			fclose(fp);
			return DTL_FRAME_CORRUPT;
			}
		}
	fclose(fp);
	return DTL_OK;
	}


static struct d_frame *read_dfile(FILE *fp, int crit) {
	int i,j,n_alts,n_stmts;
	int n_nodes[MAX_ALTS+1];
	int multilevel;
	struct stmt_rec stmt;
	struct d_frame *df;
	char df_name[FNSIZE+4];

	/* Frame header */
	if (fscanf(fp,"%120s ",df_name) < 1)
		return NULL;
	if (fscanf(fp,"%d ",&n_alts) < 1)
		return NULL;
	n_nodes[0] = 0;
	for (i=1; i<=n_alts; i++)
		fscanf(fp,"%d ",n_nodes+i);
	/* Check correct node count for MC */
	if (!crit)
		for (i=2; i<=n_alts; i++)
			if (n_nodes[i] != 1) {
#ifdef WARN_ERR
				printf("Error! MC pseudo-alt %d has node count %d (should be 1)\n",i,n_nodes[i]);
#endif
				return NULL;
				}
	/* Read multilevel */
	fscanf(fp,"%d ",&multilevel);
	if (multilevel) {
		for (i=1; i<=n_alts; i++) {
			for (j=1; j<=n_nodes[i]; j++)
				fscanf(fp,"%d ",tnext[i]+j);
			for (j=1; j<=n_nodes[i]; j++)
				fscanf(fp,"%d ",tdown[i]+j);
			}
		if (TCL_create_tree_frame(&df,n_alts,n_nodes,tnext,tdown)) {
			return NULL;
			}
		}
	else {
		if (TCL_create_flat_frame(&df,n_alts,n_nodes)) {
			return NULL;
			}
		}
	strcpy(df->name,df_name);
	/* Start using frame */
	if (TCL_attach_frame(df)) {
		TCL_dispose_frame(df);
		return NULL;
		}
	/* Probability base */
	fscanf(fp,"%d ",&n_stmts);
	if (n_stmts < 0) {
		TCL_dispose_frame(df);
		return NULL;
		}
	for (i=1; i<=n_stmts; i++) {
		fscanf(fp,"%d ",&(stmt.n_terms));
		for (j=1; j<=stmt.n_terms; j++)
			fscanf(fp,"%d %d %d ",stmt.alt+j,stmt.cons+j,stmt.sign+j);
		fscanf(fp,"%lf ",&(stmt.lobo));
		fscanf(fp,"%lf ",&(stmt.upbo));
#ifdef WARN_MC
		if (!crit)
			for (j=1; j<=stmt.n_terms; j++)
				if (stmt.alt[j] != 1)
					printf("Warning: MC weight statement using alt %d as carrier - ineffectual\n",stmt.alt[j]);
#endif
		if (stmt.n_terms == 1) {
			if (TCL_add_P_constraint(df,&stmt)) {
				TCL_dispose_frame(df);
				return NULL;
				}
			}
		else {
#ifdef WARN_LINK
			if (crit)
				printf("Skipping P-link P%d.%d-P%d.%d = [%lg,%lg] in criterion K%d%s\n",
						stmt.alt[1],stmt.cons[1],stmt.alt[2],stmt.cons[2],stmt.lobo,stmt.upbo,crit,
						stmt.alt[1]!=stmt.alt[2]?"  <- TOXIC":"");
			else
				printf("Skipping W-link W%d-W%d = [%lg,%lg] in weight base%s\n",
						stmt.cons[1],stmt.cons[2],stmt.lobo,stmt.upbo,
						stmt.alt[1]+stmt.alt[2]!=2?"  <- TOXIC":"");
#endif
			links_skipped++;
			}
		}
	/* Value base */
	fscanf(fp,"%d ",&n_stmts);
	if (n_stmts < 0) {
		TCL_dispose_frame(df);
		return NULL;
		}
	if (!crit && n_stmts) {
#ifdef V_CRIT0
		printf("Error! Criterion 0 contains a value base with %d entries\n",n_stmts);
#else
		TCL_dispose_frame(df);
		return NULL;
#endif
		}
	for (i=1; i<=n_stmts; i++) {
		fscanf(fp,"%d ",&(stmt.n_terms));
		for (j=1; j<=stmt.n_terms; j++)
			fscanf(fp,"%d %d %d ",stmt.alt+j,stmt.cons+j,stmt.sign+j);
		fscanf(fp,"%lf ",&(stmt.lobo));
		fscanf(fp,"%lf ",&(stmt.upbo));
		if (stmt.n_terms == 1) {
			if (TCL_add_V_constraint(df,&stmt)) {
				TCL_dispose_frame(df);
				return NULL;
				}
			}
		else {
#ifdef WARN_LINK
			printf("Skipping V-link V%d.%d-V%d.%d = [%lg,%lg] in criterion K%d\n",
					stmt.alt[1],stmt.cons[1],stmt.alt[2],stmt.cons[2],stmt.lobo,stmt.upbo,crit);
#endif
			links_skipped++;
			}
		}
	/* Handle different file formats */
	if ((uf_dtl_main == 5) || ((uf_dtl_main == 6) && (uf_dtl_func < 5)))
		goto no_box;

	/* Note that the box is entered as 1-dim constraints.
	 * This could overload the statement pool compared
	 * to having an actual box. */

	/* Probability box */
	fscanf(fp,"%d ",&n_stmts);
	if (n_stmts < 0) {
		TCL_dispose_frame(df);
		return NULL;
		}
	for (i=1; i<=n_stmts; i++) {
		fscanf(fp,"%d %d ",stmt.alt+1,stmt.cons+1);
		fscanf(fp,"%lf ",&(stmt.lobo));
		fscanf(fp,"%lf ",&(stmt.upbo));
		stmt.n_terms = 1;
		stmt.sign[1] = 1;
#ifdef WARN_MC
		if (!crit)
			if (stmt.alt[1] != 1)
				printf("Warning: MC weight box using alt %d as carrier - ineffectual\n",stmt.alt[1]);
#endif
		if (TCL_add_P_constraint(df,&stmt)) {
			TCL_dispose_frame(df);
			return NULL;
			}
		}
	/* Value box */
	fscanf(fp,"%d ",&n_stmts);
	if (n_stmts < 0) {
		TCL_dispose_frame(df);
		return NULL;
		}
	if (!crit && n_stmts) {
#ifdef V_CRIT0
		printf("Error! Criterion 0 contains a value box with %d entries\n",n_stmts);
#else
		TCL_dispose_frame(df);
		return NULL;
#endif
		}
	for (i=1; i<=n_stmts; i++) {
		fscanf(fp,"%d %d ",stmt.alt+1,stmt.cons+1);
		fscanf(fp,"%lf ",&(stmt.lobo));
		fscanf(fp,"%lf ",&(stmt.upbo));
		stmt.n_terms = 1;
		stmt.sign[1] = 1;
		if (TCL_add_V_constraint(df,&stmt)) {
			TCL_dispose_frame(df);
			return NULL;
			}
		}
no_box:
	/* Probability midpoints */
	fscanf(fp,"%d ",&n_stmts);
	if (n_stmts < 0) {
		TCL_dispose_frame(df);
		return NULL;
		}
#ifdef WARN_MIDPT
	if (n_stmts) {
		if (crit)
			printf("Criterion %d has %d probability midpoint%s\n",crit,n_stmts,n_stmts==1?"":"s");
		else
			printf("Frame has %d weight midpoint%s\n",n_stmts,n_stmts==1?"":"s");
		}
#endif
	for (i=1; i<=n_stmts; i++) {
		fscanf(fp,"%d %d ",stmt.alt+1,stmt.cons+1);
		fscanf(fp,"%lf ",&(stmt.lobo));
		if (fscanf(fp,"%lf ",&(stmt.upbo)) < 1) {
			TCL_dispose_frame(df);
			return NULL;
			}
		stmt.n_terms = 1;
		stmt.sign[1] = 1;
		if (TCL_add_P_mstatement(df,&stmt)) {
			TCL_dispose_frame(df);
			return NULL;
			}
#ifdef WARN_MIDPT
		if (stmt.upbo-stmt.lobo > DTL_EPS)
			if (crit)
				printf("P%d.%d.%d = [%.3lf %.3lf]  <- INTERVAL\n",crit,stmt.alt[1],stmt.cons[1],stmt.lobo,stmt.upbo);
			else
				printf("W%-2d= [%.3lf %.3lf]  <- INTERVAL\n",stmt.cons[1],stmt.lobo,stmt.upbo);
#ifdef WARN_MIDPT_EXT
		else
			if (crit)
				printf("P%d.%d.%d = %.3lf\n",crit,stmt.alt[1],stmt.cons[1],stmt.lobo);
			else
				printf("W%-2d= %.3lf\n",stmt.cons[1],stmt.lobo);
#endif
#endif
		}
	/* Value midpoints */
	fscanf(fp,"%d ",&n_stmts);
	if (n_stmts < 0) {
		TCL_dispose_frame(df);
		return NULL;
		}
	if (!crit && n_stmts) {
#ifdef V_CRIT0
		printf("Error! Criterion 0 contains %d value midpoints\n",n_stmts);
#else
		TCL_dispose_frame(df);
		return NULL;
#endif
		}
#ifdef WARN_MIDPT
	if (crit && n_stmts)
			printf("Criterion %d has %d value midpoint%s\n",crit,n_stmts,n_stmts==1?"":"s");
#endif
	for (i=1; i<=n_stmts; i++) {
		fscanf(fp,"%d %d ",stmt.alt+1,stmt.cons+1);
		fscanf(fp,"%lf ",&(stmt.lobo));
		if (fscanf(fp,"%lf ",&(stmt.upbo)) < 1) {
			TCL_dispose_frame(df);
			return NULL;
			}
		stmt.n_terms = 1;
		stmt.sign[1] = 1;
		if (TCL_add_V_mstatement(df,&stmt)) {
			TCL_dispose_frame(df);
			return NULL;
			}
#ifdef WARN_MIDPT
		if (stmt.upbo-stmt.lobo > DTL_EPS)
			printf("V%d.%d.%d = [%.3lf %.3lf]  <- INTERVAL\n",crit,stmt.alt[1],stmt.cons[1],stmt.lobo,stmt.upbo);
#ifdef WARN_MIDPT_EXT
		else
			printf("V%d.%d.%d = %.3lf\n",crit,stmt.alt[1],stmt.cons[1],stmt.lobo);
#endif
#endif
		}
	if (TCL_detach_frame(df)) {
		TCL_dispose_frame(df);
		return NULL;
		}
	return df;
	}


// DTL layer 1: DTL API level

static rcode dtl_read_file(int ufnbr, char *fn, char *folder, int mode) {
	rcode rc;
	struct user_frame *tmp_uf = NULL;

	/* Check input parameters */
	if ((ufnbr < 1) || (ufnbr > MAX_FRAMES))
		return dtl_error(DTL_FRAME_UNKNOWN);
	if (!fn || !fn[0])
		return dtl_error(DTL_NAME_MISSING);
	if (!folder) // empty folder allowed
		return dtl_error(DTL_NAME_MISSING);
	if (strlen(fn) > FNSIZE)
		return dtl_error(DTL_NAME_TOO_LONG);
	if (strlen(folder) > FOSIZE)
		return dtl_error(DTL_NAME_TOO_LONG);
	/* Allocate new user frame */
	if (!(tmp_uf = new_uf(ufnbr))) {
		return dtl_error(DTL_FRAME_EXISTS);
		}
	if (rc = read_ufile(fn,folder,tmp_uf)) {
		dispose_uf(ufnbr);
		return dtl_error(rc);
		}
	if (mode)
		strcpy(tmp_uf->frame_name,fn);
	tmp_uf->frame_nbr = ufnbr;
	return DTL_OK;
	}


rcode DTLAPI DTL_read_frame(int ufnbr, char *fn, char *folder, int mode) {
	rcode rc;

	/* Begin single thread semaphore */
	_smx_begin("FREAD");
	/* Log function call */
	if (cst_on) {
		if (fn && fn[0]) // does not log folder name
			sprintf(msg,"DTL_read_frame(%d,%.40s.dmc,%d)\n",ufnbr,fn,mode);
		else
			sprintf(msg,"DTL_read_frame(%d,_,%d)\n",ufnbr,mode);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(fn,1);
	_certify_ptr(folder,2);
	if (!dtl_is_init())
		return dtl_error(DTL_STATE_ERROR);
	if (frame_loaded)
		return dtl_error(DTL_FRAME_IN_USE);
	/* Read file */
	links_skipped = 0; // start link count
	rc = dtl_read_file(ufnbr,fn,folder,mode);
	/* End single thread semaphore */
	_smx_end();
	if (rc < DTL_OK)
		return rc;
	else
		return links_skipped;
	}


 /*********************************************************
  *
  *  Write user frame to file
  *
  *********************************************************/

static d_row lo_midbox,up_midbox;


/* Backup a user file if it does exist */

static rcode backup_ufile(char *fn, char *folder) {
	char fn_dmc[FOSIZE+FNSIZE+6],fn_bkp[FOSIZE+FNSIZE+6];

	strcpy(fn_dmc,folder);
	strcat(fn_dmc,fn);
	strcat(fn_dmc,".dmc");
	strcpy(fn_bkp,folder);
	strcat(fn_bkp,fn);
	strcat(fn_bkp,".dbk");
	remove(fn_bkp);
	if (rename(fn_dmc,fn_bkp))
		return DTL_FILE_UNKNOWN;
	return DTL_OK;
	}


static void rollback_ufile(char *fn, char *folder) {
	char fn_dmc[FOSIZE+FNSIZE+6];

	strcpy(fn_dmc,folder);
	strcat(fn_dmc,fn);
	strcat(fn_dmc,".dmc");
	remove(fn_dmc);
	}


static void nospace(char *strg) {

	for (; *strg; strg++)
		if (*strg==' ')
			*strg='_';
	}


static rcode write_dfile(FILE *fp, struct d_frame *df);

static rcode write_ufile(char *fn, char *folder) {
	rcode rc;
	int i;
	FILE *fp;
	char fn_dmc[FOSIZE+FNSIZE+6],frn[FNSIZE+6];

	rc = DTL_OK;
	backup_ufile(fn,folder);
	/* Check input parameters */
	strcpy(fn_dmc,folder);
	strcat(fn_dmc,fn);
	strcat(fn_dmc,".dmc");
	fp = fopen(fn_dmc,"w");
	if (!fp)
		return DTL_FILE_UNKNOWN;

	/* User frame header (note: DTL_MAIN bumped by 7 to be compat with DMC tests) */
	fprintf(fp,"%d.%02d\n",DTL_MAIN+7,DTL_FUNC);
	strcpy(frn,uf->frame_name);
	nospace(frn);
	fprintf(fp,"%s\n",frn);

	/* Frame parameters */
	fprintf(fp,"%d\n",uf->frame_type);
	fprintf(fp,"%d ",uf->n_alts);
	fprintf(fp,"%d\n",uf->n_crit);

	/* Tree description */
	if (PM) {
		// criteria map
		for (i=0; i<=uf->n_crit; i++)
			fprintf(fp,"%d ",uf->df_list[i]?1:0);
		fprintf(fp,"\n");
		// write weights to mirror map
		for (i=0; i<=uf->n_crit; i++)
			if (uf->df_list[i]) {
				if (rc = load_df0(i)) {
					fclose(fp);
					return rc;
					}
				/* Log function call */
				if (cst_ext) {
					sprintf(msg," writing PM-crit %d...\n",i);
					cst_log(msg);
					}
				rc = write_dfile(fp,uf->df);
				/* Log function result */
				if (cst_ext) {
					if (rc)
						sprintf(msg," ...crit %d write error %d\n",i,rc);
					else
						sprintf(msg," ...crit %d written\n",i);
					cst_log(msg);
					}
				if (rc) {
					fclose(fp);
					return rc;
					}
				}
		}
	else {
		/* Log function call */
		if (cst_ext)
			cst_log(" writing frame...\n");
		rc = write_dfile(fp,uf->df);
		/* Log function result */
		if (cst_ext) {
			if (rc)
				sprintf(msg," ...frame write error %d\n",rc);
			else
				sprintf(msg," ...frame written\n");
			cst_log(msg);
			}
		}
	fclose(fp);
	return rc;
	}


static rcode write_dfile(FILE *fp, struct d_frame *df) {
	rcode rc;
	int i,j,k,n_midpoints;
	struct base *P,*V;
	char frn[FNSIZE+6];

	rc = DTL_OK;
	strcpy(frn,df->name);
	nospace(frn);
	fprintf(fp,"%s\n",frn);

	/* Decision frame alternatives */
	fprintf(fp,"%d ",df->n_alts);
	for (i=1; i<=df->n_alts; i++)
		fprintf(fp,"%d ",df->tot_cons[i]);
	fprintf(fp,"\n");

	/* Multilevel flag */
	fprintf(fp,"%d \n",df->tree);

	if (df->tree) {
		/* Tree description */
		for (i=1; i<=df->n_alts; i++) {
			for (j=1; j<=df->tot_cons[i]; j++)
				fprintf(fp,"%d ",df->next[i][j]);
			fprintf(fp,"\n");
			for (j=1; j<=df->tot_cons[i]; j++)
				fprintf(fp,"%d ",df->down[i][j]);
			fprintf(fp,"\n");
			}
		}

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
	fprintf(fp,"%d\n",V->n_stmts);
	for (i=1; i<=V->n_stmts; i++) {
		fprintf(fp,"%d ",V->stmt[i].n_terms);
		for (j=1; j<=V->stmt[i].n_terms; j++)
			fprintf(fp,"%d %d %d ",V->stmt[i].alt[j],V->stmt[i].cons[j],V->stmt[i].sign[j]);
		fprintf(fp,"%.10le ",V->stmt[i].lobo);
		fprintf(fp,"%.10le\n",V->stmt[i].upbo);
		}

	/* Probability box */
	if (rc = TCL_get_P_box(df,lo_midbox,up_midbox))
		return rc;
	/* Calculate number of box entries */
	n_midpoints = 0;
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++)
			if ((df->tot_cons[i] > 1) && ((lo_midbox[k] > DTL_EPS) || (up_midbox[k] < 1.0-DTL_EPS)))
				n_midpoints++;
	fprintf(fp,"%d\n",n_midpoints);
	/* Write box in order */
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			if ((df->tot_cons[i] > 1) && ((lo_midbox[k] > DTL_EPS) || (up_midbox[k] < 1.0-DTL_EPS))) {
				fprintf(fp,"%d %d ",i,j);
				fprintf(fp,"%.10le ",lo_midbox[k]);
				fprintf(fp,"%.10le\n",up_midbox[k]);
				}
			}

	/* Value box */
	if (rc = TCL_get_V_box(df,lo_midbox,up_midbox))
		return rc;
	/* Calculate number of box entries */
	n_midpoints = 0;
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++)
			if ((lo_midbox[k] > DTL_EPS) || ((up_midbox[k] < 1.0-DTL_EPS) && (up_midbox[k] != -1.0)))
				n_midpoints++;
	fprintf(fp,"%d\n",n_midpoints);
	/* Write box in order */
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++)
			if ((lo_midbox[k] > DTL_EPS) || ((up_midbox[k] < 1.0-DTL_EPS) && (up_midbox[k] != -1.0))) {
				fprintf(fp,"%d %d ",i,j);
				fprintf(fp,"%.10le ",lo_midbox[k]);
				fprintf(fp,"%.10le\n",up_midbox[k]);
				}

	/* Probability midpoints */
	if (rc = TCL_get_P_mbox(df,lo_midbox,up_midbox))
		return rc;
	/* Calculate number of midpoints */
	n_midpoints = 0;
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++)
			if ((df->tot_cons[i] > 1) && (lo_midbox[k] > -1.0))
				n_midpoints++;
	fprintf(fp,"%d\n",n_midpoints);
	/* Write midpoints in order */
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++) {
			if ((df->tot_cons[i] > 1) && (lo_midbox[k] > -1.0)) {
				fprintf(fp,"%d %d ",i,j);
				fprintf(fp,"%.10le ",lo_midbox[k]);
				fprintf(fp,"%.10le\n",up_midbox[k]);
				}
			}

	/* Value midpoints */
	if (rc = TCL_get_V_mbox(df,lo_midbox,up_midbox))
		return rc;
	/* Calculate number of midpoints */
	n_midpoints = 0;
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++)
			if (lo_midbox[k] > -1.0)
				n_midpoints++;
	fprintf(fp,"%d\n",n_midpoints);
	/* Write midpoints in order */
	for (k=1,i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++,k++)
			if (lo_midbox[k] > -1.0) {
				fprintf(fp,"%d %d ",i,j);
				fprintf(fp,"%.10le ",lo_midbox[k]);
				fprintf(fp,"%.10le\n",up_midbox[k]);
				}
	return rc;
	}


static rcode dtl_write_file(char *fn, char *folder) {
	rcode rc;

	/* Check input parameters */
	if (!fn || !fn[0])
		return DTL_NAME_MISSING;
	if (!folder || !folder[0])
		return DTL_NAME_MISSING;
	if (strlen(fn) > FNSIZE)
		return DTL_NAME_TOO_LONG;
	if (strlen(folder) > FOSIZE)
		return DTL_NAME_TOO_LONG;
	/* Write file */
	if (rc = write_ufile(fn,folder)) {
		/* Not a proper file created */
		rollback_ufile(fn,folder);
		return DTL_KERNEL_ERROR+rc; // context switch from TCL to DTL
		}
	return DTL_OK;
	}


/* DM-frames are stored verbatim here as PM-frames, while for
 * SM-frames there is no .dmc file format counterpart. Thus,
 * for SM2-frames (duplicated criteria) all is captured except
 * nbr of stakeholders, but for SM1-frames (mirrored criteria)
 * the content of each SM criterion is copied per stakeholder. */

rcode DTLAPI DTL_write_frame(char *fn, char *folder) {
	rcode rc;

	/* Begin single thread semaphore */
	_smx_begin("FWRT");
	/* Log function call */
	if (cst_on) {
		if (fn && fn[0] && folder && folder[0])
			sprintf(msg,"DTL_write_frame(%.40s,%.20s)\n",fn,folder);
		else
			sprintf(msg,"DTL_write_frame(_,_)\n");
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(fn,1);
	_certify_ptr(folder,2);
	if (!dtl_is_init())
		return dtl_error(DTL_STATE_ERROR);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Write file */
	rc = dtl_write_file(fn,folder);
	/* End single thread semaphore */
	_smx_end();
	return rc;
	}
