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
 *   File: DTLfile2.c
 *
 *   Purpose: reading .ddt files from the UCT test tool
 *
 *
 *   Functions exported outside DTL
 *   ------------------------------
 *   DTL_read_ddt_frame
 *
 *   Functions exported outside DTL, debug use
 *   -----------------------------------------
 *   NONE
 *
 *   Functions outside of module, inside DTL
 *   ---------------------------------------
 *   NONE
 *
 *   Functions internal to module
 *   ----------------------------
 *   read_ddt_ufile
 *   dtl_read_ddt_file
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
// warning if a file has probability or value midpoints
#define WARN_MIDPT
#define noWARN_MIDPT_EXT
// warning if .ddt file value scale differs from internal [0,1]
#define WARN_VSCALE
// Warning if alternative is too big
#define WARN_ERR
// Warning if links encountered (not supported in basic UNEDA)
#define WARN_LINK
#endif

#define FOSIZE 128 // folder name size
#define FOSIZE_DDT 192 // folder path size
#define MAX_STMTS_LOCAL 50


 /********************************************************
  *
  *  Internal data
  *
  ********************************************************/

static t_matrix tnext,tdown;
static int uf_ver_main,uf_ver_func,uf_ver_tech,multilevel;
static char alt_name[MAX_ALTS+1][32];
static char crit_name[MAX_CRIT+1][32];
static double v_lo,v_up;
static int links_skipped,P_skipped;
static h_matrix lobo,midpt,upbo;
static struct stmt_rec skip[MAX_STMTS_LOCAL+1];


 /*********************************************************
  *
  *  File format
  *
  *********************************************************/

 /*
  * .ddt file format for PS and DM file types
  *
  * ufile structure
  * ---------------
  * main.func.tech
  * fname
  * type (0=old_PS, 1=PS, 2=DM)
  * if (main > 3)
  *   tree (0=flat 1=multilevel)
  * #alts {nodes} (one per alt)
  * if (multilevel)
  *   tree structure [see below]
  * if (main.func > 3.2)
  *   {alt_name} (one per alt)
  *   if (#crit > 1)
  *     {crit_name} (one per crit)
  * #p_base-stmts
  * {p_base-stmt} (one row per stmt)
  * v_min v_max   (heeded by DTL -> non-compliant vs DMC)
  * #v_base-stmts
  * {v_base-stmt} (one row per stmt)
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
  * NOTE: The error handling is much more limited since this is an
  * internal function in the sense that it only reads its own data.
  *
  */

 /********************************************************
  *
  *  File access functions
  *
  ********************************************************/

/* Read user frame from file */

static rcode read_ddt_ufile(char *fn, char *folder, struct user_frame *uf) {
	rcode rc;
	int i,j,n_stmts;
	int n_cons[MAX_ALTS+1];
	FILE *fp;
	struct stmt_rec stmt;
	char fn_ddt[64];

	rc = DTL_OK;
	/* Check input parameters */
	strcpy(fn_ddt,folder);
	strcat(fn_ddt,fn);
	strcat(fn_ddt,".ddt");
	fp = fopen(fn_ddt,"r");
	if (!fp)
		return DTL_FILE_UNKNOWN;
	/* User frame header */
	fscanf(fp,"%d.%d.%d ",&uf_ver_main,&uf_ver_func,&uf_ver_tech);
	fscanf(fp,"%120s ",uf->frame_name);
	/* Read frame type */
	fscanf(fp,"%d ",&uf->frame_type);
	/* Old flag can be 0 */
	if (!uf->frame_type)
		uf->frame_type = PS_FRAME;
	/* Read level type if multilevel release */
	if (uf_ver_main > 3)
		fscanf(fp,"%d ",&multilevel);
	else
		multilevel = FALSE;
	/* .ddt frame header */
	fscanf(fp,"%d ",&(uf->n_alts));
	if ((uf->n_alts < 2) || (uf->n_alts > MAX_ALTS)) {
		fclose(fp);
		return DTL_FRAME_CORRUPT;
		}
	n_cons[0] = 0;
	for (i=1; i<=uf->n_alts; i++) {
		fscanf(fp,"%d ",n_cons+i);
		if (DM && (n_cons[i]!=n_cons[1])) {
			fclose(fp);
			return DTL_FRAME_CORRUPT;
			}
		}
	if (multilevel) {
		for (i=1; i<=uf->n_alts; i++) {
			for (j=1; j<=n_cons[i]; j++)
				fscanf(fp,"%d ",tnext[i]+j);
			for (j=1; j<=n_cons[i]; j++)
				fscanf(fp,"%d ",tdown[i]+j);
			}
		if (rc = TCL_create_tree_frame(&(uf->df),uf->n_alts,n_cons,tnext,tdown)) {
			fclose(fp);
			return DTL_KERNEL_ERROR+rc;
			}
		}
	else {
		if (rc = TCL_create_flat_frame(&(uf->df),uf->n_alts,n_cons)) {
			fclose(fp);
			return DTL_KERNEL_ERROR+rc;
			}
		}
	strcpy(uf->df->name,uf->frame_name);
	if (DM)
		uf->n_crit = n_cons[1];
	else
		uf->n_crit = 1;
	if ((uf->n_crit < 1) || (uf->n_crit > MAX_CRIT)) {
		fclose(fp);
		return DTL_FRAME_CORRUPT;
		}
	uf->n_sh = 1;
	/* Attach frame to load structure */
	if (rc = TCL_attach_frame(uf->df)) {
		TCL_dispose_frame(uf->df);
		fclose(fp);
		return DTL_KERNEL_ERROR+rc;
		}
	/* Read names */
	if ((uf_ver_main > 3) || ((uf_ver_main == 3) && (uf_ver_func > 2))) {
		for (i=1; i<=uf->n_alts; i++)
			fscanf(fp,"%29s ",alt_name[i]);
		if (uf->n_crit > 1)
			for (i=1; i<=uf->n_crit; i++)
				fscanf(fp,"%29s ",crit_name[i]);
		}
	else {
		for (i=1; i<=uf->n_alts; i++)
			strcpy(alt_name[i],"no_name");
		if (uf->n_crit > 1)
			for (i=1; i<=uf->n_crit; i++)
				strcpy(crit_name[i],"no_name");
		}
	/* Probability base */
	fscanf(fp,"%d ",&n_stmts);
	if (n_stmts < 0) {
		TCL_dispose_frame(uf->df);
		fclose(fp);
		return DTL_FRAME_CORRUPT;
		}
	for (i=1; i<=n_stmts; i++) {
		fscanf(fp,"%d ",&(stmt.n_terms));
		for (j=1; j<=stmt.n_terms; j++) {
			fscanf(fp,"%d %d %d ",stmt.alt+j,stmt.cons+j,stmt.sign+j);
			if (DM && (stmt.alt[j]!=1)) {
				TCL_dispose_frame(uf->df);
				fclose(fp);
				return DTL_FRAME_CORRUPT;
				}
			}
		fscanf(fp,"%lf ",&(stmt.lobo));
		fscanf(fp,"%lf ",&(stmt.upbo));
		if (stmt.n_terms == 1) {
			if (rc = TCL_add_P_constraint(uf->df,&stmt)) {
				TCL_dispose_frame(uf->df);
				fclose(fp);
				return DTL_KERNEL_ERROR+rc;
				}
			}
		else {
#ifdef WARN_LINK
			printf("Skipping P-link P%d.%d-P%d.%d = [%lg %lg]%s\n",
					stmt.alt[1],stmt.cons[1],stmt.alt[2],stmt.cons[2],stmt.lobo,stmt.upbo,
					stmt.alt[1]!=stmt.alt[2]?"  <- TOXIC":"");
#endif
			links_skipped++;
			if (links_skipped <= MAX_STMTS_LOCAL) {
				memcpy(&skip[links_skipped],&stmt,sizeof(struct stmt_rec));
				}
			P_skipped++;
			}
		}
	/* Value base */
	if ((uf_ver_main > 3) || ((uf_ver_main == 3) && (uf_ver_func > 2))) {
		/* Read and install value range */
		fscanf(fp,"%lf ",&v_lo);
		fscanf(fp,"%lf ",&v_up);
		uf->av_min[0] = v_lo;
		uf->av_max[0] = v_up;
#ifdef WARN_VSCALE
		if ((v_lo != 0.0) || (v_up != 1.0))
			printf("Warning: Frame has value range [%.3lf %.3lf]\n",v_lo,v_up);
#endif
		}
	fscanf(fp,"%d ",&n_stmts);
	if (n_stmts < 0) {
		TCL_dispose_frame(uf->df);
		fclose(fp);
		return DTL_FRAME_CORRUPT;
		}
	for (i=1; i<=n_stmts; i++) {
		fscanf(fp,"%d ",&(stmt.n_terms));
		for (j=1; j<=stmt.n_terms; j++)
			fscanf(fp,"%d %d %d ",stmt.alt+j,stmt.cons+j,stmt.sign+j);
		fscanf(fp,"%lf ",&(stmt.lobo));
		fscanf(fp,"%lf ",&(stmt.upbo));
		if (stmt.n_terms == 1) {
			if (rc = TCL_add_V_constraint(uf->df,&stmt)) {
				TCL_dispose_frame(uf->df);
				fclose(fp);
				return DTL_KERNEL_ERROR+rc;
				}
			}
		else {
#ifdef WARN_LINK
			printf("Skipping V-link V%d.%d-V%d.%d = [%lg %lg]\n",
					stmt.alt[1],stmt.cons[1],stmt.alt[2],stmt.cons[2],stmt.lobo,stmt.upbo);
#endif
			links_skipped++;
			if (links_skipped <= MAX_STMTS_LOCAL)
				memcpy(&skip[links_skipped],&stmt,sizeof(struct stmt_rec));
			}
		}
	if (uf_ver_main > 2) {
		/* Probability midpoints */
		fscanf(fp,"%d ",&n_stmts);
		if (n_stmts < 0) {
			TCL_dispose_frame(uf->df);
			fclose(fp);
			return DTL_FRAME_CORRUPT;
			}
#ifdef WARN_MIDPT
		if (n_stmts)
			printf("Frame has %d %s midpoint%s\n",n_stmts,DM?"weight":"probability",n_stmts==1?"":"s");
#endif
		for (i=1; i<=n_stmts; i++) {
			fscanf(fp,"%d %d ",stmt.alt+1,stmt.cons+1);
			if (DM && (stmt.alt[1]!=1)) {
#ifdef WARN_ERR
				printf("Error! DM weight statement using alt %d as carrier - illegal\n",stmt.alt[1]);
#endif
				TCL_dispose_frame(uf->df);
				fclose(fp);
				return DTL_FRAME_CORRUPT;
				}
			fscanf(fp,"%lf ",&(stmt.lobo));
			fscanf(fp,"%lf ",&(stmt.upbo));
			stmt.n_terms = 1;
			stmt.sign[1] = 1;
			if (rc = TCL_add_P_mstatement(uf->df,&stmt)) {
				TCL_dispose_frame(uf->df);
				fclose(fp);
				return DTL_KERNEL_ERROR+rc;
				}
#ifdef WARN_MIDPT
			if (stmt.upbo-stmt.lobo > DTL_EPS)
				printf("P%d.%d = [%.3lf %.3lf]  <- INTERVAL\n",stmt.alt[1],stmt.cons[1],stmt.lobo,stmt.upbo);
#ifdef WARN_MIDPT_EXT
			else
				printf("P%d.%d = %.3lf\n",stmt.alt[1],stmt.cons[1],stmt.lobo);
#endif
#endif
			}
		/* Value midpoints */
		fscanf(fp,"%d ",&n_stmts);
		if (n_stmts < 0) {
			TCL_dispose_frame(uf->df);
			fclose(fp);
			return DTL_FRAME_CORRUPT;
			}
#ifdef WARN_MIDPT
		if (n_stmts)
			printf("Frame has %d value midpoint%s\n",n_stmts,n_stmts==1?"":"s");
#endif
		for (i=1; i<=n_stmts; i++) {
			fscanf(fp,"%d %d ",stmt.alt+1,stmt.cons+1);
			fscanf(fp,"%lf ",&(stmt.lobo));
			fscanf(fp,"%lf ",&(stmt.upbo));
			stmt.n_terms = 1;
			stmt.sign[1] = 1;
			if (rc = TCL_add_V_mstatement(uf->df,&stmt)) {
				TCL_dispose_frame(uf->df);
				fclose(fp);
				return DTL_KERNEL_ERROR+rc;
				}
#ifdef WARN_MIDPT
			if (stmt.upbo-stmt.lobo > DTL_EPS)
				printf("V%d.%d = [%.3lf %.3lf]  <- INTERVAL\n",stmt.alt[1],stmt.cons[1],stmt.lobo,stmt.upbo);
#ifdef WARN_MIDPT_EXT
			else
				printf("V%d.%d = %.3lf\n",stmt.alt[1],stmt.cons[1],stmt.lobo);
#endif
#endif
			}
		}
	/* Unload frame */
	if (rc = TCL_detach_frame(uf->df)) {
		TCL_dispose_frame(uf->df);
		fclose(fp);
		return DTL_KERNEL_ERROR+rc;
		}
	fclose(fp);
	return rc;
	}


// DTL layer 1: DTL API level

static rcode dtl_read_ddt_file(int ufnbr, char *fn, char *folder, int mode) {
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
	/* Read file */
	if (rc = read_ddt_ufile(fn,folder,tmp_uf)) {
		dispose_uf(ufnbr);
		return dtl_error(rc);
		}
	if (mode)
		strcpy(tmp_uf->frame_name,fn);
	tmp_uf->frame_nbr = ufnbr;
	return DTL_OK;
	}


rcode DTLAPI DTL_read_ddt_frame(int ufnbr, char *fn, char *folder, int mode) {
	rcode rc;

	/* Begin single thread semaphore */
	_smx_begin("FRDDT");
	/* Log function call */
	if (cst_on) {
		if (fn && fn[0]) // does not log folder name
			sprintf(msg,"DTL_read_ddt_frame(%d,%.40s.ddt,%d)\n",ufnbr,fn,mode);
		else
			sprintf(msg,"DTL_read_ddt_frame(%d,_,%d)\n",ufnbr,mode);
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
	links_skipped = P_skipped = 0; // start link count
	rc = dtl_read_ddt_file(ufnbr,fn,folder,mode);
	/* End single thread semaphore */
	_smx_end();
	if (rc < DTL_OK)
		return rc;
	else
		return links_skipped;
	}
