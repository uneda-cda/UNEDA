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
 *   File: DTLautoscale.c
 *
 *   Purpose: Automatic scale add-in
 *
 *
 *   Functions exported outside DTL
 *   ------------------------------
 *   DTL_set_AV_box
 *   DTL_set_AV_mbox/1
 *   DTL_set_AV_modal
 *   DTL_get_AV_crit_scale
 *   DTL_set_AV_MC_scale
 *   DTL_copy_AV_MC_scale
 *   DTL_reset_AV_MC_scale
 *   DTL_get_AV_MC_scale
 *   DTL_get_AV_user_vector
 *   DTL_get_AV_user_value
 *   DTL_get_AV_user_intervals
 *   DTL_get_AV_user_interval
 *   DTL_get_AV_norm_vector
 *   DTL_get_AV_norm_value
 *   DTL_get_AV_norm_intervals
 *   DTL_get_AV_norm_interval
 *   DTL_check_AV_user_values
 *   DTL_check_AV_norm_values
 *
 *   Functions outside module, debug use
 *   -----------------------------------
 *   DTI_set_AV_crit_scale
 *   DTI_reset_AV_crit_scale
 *   DTI_AV_scale_ratio
 *   DTI_check_AV_values
 *   DTI_is_AV_default_scale
 *
 *   Functions outside module, inside DTL
 *   ------------------------------------
 *   dtl_copy_AV_crit_scale
 *
 *   Functions internal to module
 *   ----------------------------
 *   dtl_set_AV_scale
 *   dtl_get_AV_scale
 *   dtl_trf_AV_input
 *   dtl_AV_renorm_W_base
 *   dti_AV_rev_scale
 *   dti_AV_opposite_scales
 *   dti_AV_scale_ratio
 *
 */


 /********************************************************
  *
  *  Compiler directives
  *
  ********************************************************/

#include <stdarg.h> // variadic function header file

#ifdef _MSC_VER
/* VC++ compiler mistake: warns for non-existent problem (varargs allowed in C89)
 * "warning: nonstandard extension used: function declaration used ellipsis". The
 * var-argument types must be self-promoting so that the default promotions do not
 * change their types. That rules out array and function types, as well as float,
 * char (signed or not), and short int (signed or not). This is C89/C99 standard. */
#pragma warning(disable:4212)
#pragma warning(disable:4701)
#endif


 /********************************************************
  *
  *  Configuration parameters
  *
  ********************************************************/

#define noLOG_AS // logging autoscale invocations on top of DTL
#define noSEMI_RATIO // inflating Dirac scales to semi-ratio
#define PART_TREE // conversions can be called with partial tree


 /********************************************************
  *
  *  Autoscale for value bases (and thus for DTL results).
  *  This is in practice a code layer on top of DTLvbase,
  *  even though it is included at the end of that file.
  *
  *  It is a companion to the SML package, with the same
  *  purpose of simplifying access to DTL functionality.
  *
  *  NOTE: Autoscale uses only box + mbox, not statements.
  *
  *  DTL layer 0: all functions here are above DTL proper
  *
  ********************************************************/

#define VSCALE_MAXVAL  1.0E+12
#define VSCALE_MINSPAN 1.0E-03


/* Scale setting and retrieval calls. Only MC (output) scale is allowed
 * to be set manually, otherwise the meaning of statements would change.
 * Further, this is only allowed for PM-frames (PS has only one scale).
 *
 * [For test and debug, any scale can be changed through a DTI call.]
 *
 * Autoscale usage sequence
 * ------------------------
 * 1. Set each crit scale by entering box+mbox or modal values
 * 2. For PM: set MC output scale (for PS: nothing to set)
 * 3. Convert to/from DTL user/norm values and vectors */


 /*********************************************************
  *
  *  Internal autoscale functions (set,get)
  *
  *********************************************************/

static rcode dtl_set_AV_scale(int crit, double v_min, double v_max) {

#ifdef LOG_AS
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_set_AV_scale(%d,%.3lf,%.3lf)\n",crit,v_min,v_max);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameters */
	if (check_df0(crit)) // validate crit
		return DTL_CRIT_UNKNOWN;
	if ((_fabs(v_min) > VSCALE_MAXVAL) || (_fabs(v_max) > VSCALE_MAXVAL))
		return DTL_INPUT_ERROR;
	if (_fabs(v_max-v_min) < VSCALE_MINSPAN)
		return DTL_INPUT_ERROR; // must span an interval
	/* Set autoscale endpoints */
	uf->av_min[crit] = v_min;
	uf->av_max[crit] = v_max;
	return DTL_OK;
	}


static rcode dtl_get_AV_scale(int crit, double *v_min, double *v_max) {

#ifdef LOG_AS
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_AV_scale(%d)\n",crit);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	_init_assert();
	_certify_ptr(v_min,201);
	_certify_ptr(v_max,202);
	_dtl_assert(v_min!=v_max,201);
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameter */
	if (check_df0(crit)) // validate crit
		return DTL_CRIT_UNKNOWN;
	/* Return autoscale endpoints */
	*v_min = uf->av_min[crit];
	*v_max = uf->av_max[crit];
	return DTL_OK;
	}


 /***********************************************************
  *
  *  Criteria autoscale functions (set,copy,reset,get)
  *
  *  If the input spans a scale, the endpoints of the input
  *  become the endpoints of the automatic scale. But when
  *  the input does not span a scale (all are the same = a),
  *  if the module configuration flag SEMI_RATIO is defined,
  *  a semi-ratio scale is synthetically spanned with [0,2a]
  *  as the spanned scale if a>0 and [2a,0] if a<0. For a=0,
  *  or if the configuration flag SEMI_RATIO is not defined,
  *  a scale is spanned with range [a-MIN_SPAN,a+MIN_SPAN].
  *
  ***********************************************************/

#define TNO_MODAL 0 // do not handle modal values
#define TRF_MODAL 1 // transform also modal values

static h_matrix av_lobox,av_modalx,av_upbox; // autoscale matrices

static rcode dtl_trf_AV_input(int crit, int mode, bool rev, h_matrix lobox, 
		h_matrix modalx, h_matrix upbox, double *t_min, double *t_max) {
	int i,j,as_chg;
	struct d_frame *df;
	// input scale, set up at opposite extremes
	double vi_min = +VSCALE_MAXVAL;
	double vi_max = -VSCALE_MAXVAL;
	// criterion scale, empty so far
	double v_min,v_max;

	/* Check if subroutine can start */
	_init_assert();
	_certify_ptr(lobox,201);
	_certify_ptr(modalx,202);
	_certify_ptr(upbox,203);
	_certify_ptr(t_min,204);
	_certify_ptr(t_max,205);
	_dtl_assert(lobox!=upbox,201);
	_dtl_assert(lobox!=modalx,202);
	_dtl_assert(modalx!=upbox,203);
	/* NOTE: requires df to be loaded, not only checked, upon calling */
	df = uf->df;
	/* Seek input scale endpoints */
	for (i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++) {
			if (lobox[i][j] < -VSCALE_MAXVAL)
				return DTL_INPUT_ERROR;
			vi_min = min(vi_min,lobox[i][j]);
			if (upbox[i][j] > +VSCALE_MAXVAL)
				return DTL_INPUT_ERROR;
			vi_max = max(vi_max,upbox[i][j]);
			}
	if (vi_min > vi_max)
		return DTL_INCONSISTENT;
	if (vi_max-vi_min < VSCALE_MINSPAN)
		/* Scale not spanned -> should be adjusted in such a way so that the
		 * Dirac point (if it exists) winds up in the middle of the scale */
#ifdef SEMI_RATIO
		if ((vi_max < 0.0) || (vi_min > 0.0)) {
			/* Dirac non-zero scale -> inflate to semi-ratio */
			vi_min = min(vi_min+vi_max,0.0);
			vi_max = max(vi_min+vi_max,0.0);
			}
		else
#endif
			{
			/* Dirac (zero/non-zero) scale -> expand symmetrically */
			vi_min -= VSCALE_MINSPAN;
			vi_max += VSCALE_MINSPAN;
			}
	/* Get scale mapping - normal or reverse */
	v_min = rev?vi_max:vi_min;
	v_max = rev?vi_min:vi_max;
	/* Transform all entries from user (external)
	 * [a,b] scale to DTL internal [0,1] scale */
	for (i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++) {
			av_lobox[i][j] = ((rev?upbox[i][j]:lobox[i][j])-v_min)/(v_max-v_min);
			if (mode) { // modal exists
				if ((modalx[i][j]<vi_min) || (modalx[i][j]>vi_max))
					return DTL_INCONSISTENT;
				av_modalx[i][j] = (modalx[i][j]-v_min)/(v_max-v_min);
				}
			av_upbox[i][j] = ((rev?lobox[i][j]:upbox[i][j])-v_min)/(v_max-v_min);
			}
	/* Check if change in autoscale */
	as_chg = (v_min != uf->av_min[crit]) || (v_max != uf->av_max[crit]);
#ifdef LOG_AS
	/* Log new scale */
	if (cst_ext && as_chg) {
		sprintf(msg,"  scale change [%lg %lg] -> [%lg %lg]\n",
				uf->av_min[crit],uf->av_max[crit],v_min,v_max);
		cst_log(msg);
		}
#endif
	/* Return new scale */
	*t_min = v_min;
	*t_max = v_max;
	/* Return scale change flag */
	return as_chg?DTL_SCALE_CHANGE:DTL_OK;
	}


/* Duality renormalisation using the DURENO-I technique */

static h_vector rn_lobo,rn_mid,rn_upbo;

static rcode dtl_AV_renorm_W_base(int crit, double sfact) {
	rcode rc;
	int snode,tnode,lnode;
	double norm;

	/* Load MC uf pointer and get current weights */
	if (rc = DTL_get_W_hull(0,rn_lobo,rn_mid,rn_upbo))
		return rc;
	/* Find node nbr for modified criterion */
	if ((snode = dtl_crit2node(crit)) < DTL_OK)
		return snode;
	if (!snode) // no crit node found
		return DTL_CRIT_UNKNOWN;
	/* Calculate new weight norm */
	norm = (sfact-1.0)*rn_mid[snode] + 1.0;
	/* Rescale affected node by scaling factor */
	rn_lobo[snode] *= sfact;
	rn_mid[snode]  *= sfact;
	rn_upbo[snode] *= sfact;
	/* Seek first sibling node (Lisp car/cdr-style pointers) */
	for (tnode=snode; tnode; tnode=uf->df->prev[1][tnode])
		lnode=tnode; // snode > 0 -> will always do at least one round
	/* Renormalise all siblings (lnode is the leftmost one) */
	for (tnode=lnode; tnode; tnode=uf->df->next[1][tnode]) {
		rn_lobo[tnode] /= norm;
		rn_mid[tnode]  /= norm;
		rn_upbo[tnode]  = min(rn_upbo[tnode]/norm,1.0); // overflow protection
		}
	/* Set renormalised weights
	 * NOTE: weights without a declared midpoint
	 * will have one assigned as a side effect */
	if (rc = DTL_reset_W_base())
		return rc;
	if (rc = DTL_set_W_mbox1(rn_mid))
		return rc;
	return DTL_set_W_box(rn_lobo,rn_upbo);
	}


/* Box loading calls. NOTE1: If DTL_set_AV_box signals a scale
 * change, the midbox must be reloaded using DTL_set_AV_mbox.
 *
 * NOTE2: To have lower values being preferred, set 'rev'=TRUE.
 *
 * NOTE3: Since this is a layer on top of the value base: using
 * the standard API calls DTL_set_V_box & DTL_set_V_mbox, care
 * must be taken to reset the scales if the standard calls are
 * being used without this layer (i.e. bypassing the layer). A
 * user layer on top must manage its own integrity in the base.
 *
 * NOTE4: Never signal remorm before the weights have been set! */

rcode DTLAPI DTL_set_AV_box(int crit, bool rev, bool renorm, h_matrix lobox, h_matrix upbox) {
	rcode rc,rc2;
	double v_min,v_max,scaling;

#ifdef LOG_AS
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_set_AV_box(%d,%d,%d)\n",crit,rev,renorm);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameter */
	if (load_df1(crit))
		return DTL_CRIT_UNKNOWN;
	/* Transform input to [0,1] scale */
	if (DTL_error2(rc = dtl_trf_AV_input(crit,TNO_MODAL,rev,lobox,av_modalx,upbox,&v_min,&v_max)))
		return rc;
	/* If scale changed, midbox must be cleared now and later
	 * reinstated - or else there might be inconsistencies. */
	if (rc == DTL_SCALE_CHANGE)
		if (rc2 = DTL_remove_V_mbox(crit))
			return rc2;
	/* Enter new [0,1]-scaled box */
	if (rc2 = DTL_set_V_box(crit,av_lobox,av_upbox))
		return rc2;
	if (rc == DTL_SCALE_CHANGE) {
		/* Calculate scaling factor (span > 0 -> no div-by-zero risk) */
		scaling = (v_max-v_min) / (uf->av_max[crit]-uf->av_min[crit]);
		/* Set new autoscale endpoints */
		uf->av_min[crit] = v_min;
		uf->av_max[crit] = v_max;
		if (renorm) { // renormalisation requested & required
			/* Renormalise criteria weights */
			if (DTL_error2(rc2 = dtl_AV_renorm_W_base(crit,scaling)))
				return rc2;
			}
		}
	/* Signal whether autoscale has changed */
	return rc;
	}


/* NOTE: -1.0 and -2.0 are legit mbox inputs in standard DTL, but not
 * in autoscale or simplified SML. If necessary in the future, empty
 * markers could be signalled by a marker in lobox and a lower value
 * in upbox, thereby defying the ordinary DTL V-base syntax rulebook.
 * But that would unnecessarily complicate the user <-> norm functions.
 * The issue is easily circumvented by entering the numerical midpoint. */

rcode DTLAPI DTL_set_AV_mbox(int crit, h_matrix lobox, h_matrix upbox) {
	int i,j,rev;
	struct d_frame *df;
	double v_min,v_max,vi_min,vi_max;

	/* Standard or reverse scale (if higher or lower values are better) */
	rev = frame_loaded?uf->av_min[crit]>uf->av_max[crit]:0;
#ifdef LOG_AS
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_set_AV_mbox(%d,%d)\n",crit,rev);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	_certify_ptr(lobox,201);
	_certify_ptr(upbox,202);
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameter */
	if (load_df1(crit))
		return DTL_CRIT_UNKNOWN;
	df = uf->df;
	/* Autoscale endpoints */
	v_min = uf->av_min[crit];
	v_max = uf->av_max[crit];
	vi_min = min(v_min,v_max);
	vi_max = max(v_min,v_max);
	/* Transform all entries from user [a,b] to norm [0,1] scale */
	for (i=1; i<=df->n_alts; i++)
		for (j=1; j<=df->tot_cons[i]; j++) {
			/* Check endpoint compliance */
			if ((lobox[i][j] < vi_min) || (upbox[i][j] > vi_max))
				return DTL_INPUT_ERROR;
			av_lobox[i][j] = ((rev?upbox[i][j]:lobox[i][j])-v_min)/(v_max-v_min);
			av_upbox[i][j] = ((rev?lobox[i][j]:upbox[i][j])-v_min)/(v_max-v_min);
			}
	return DTL_set_V_mbox(crit,av_lobox,av_upbox);
	}


rcode DTLAPI DTL_set_AV_mbox1(int crit, h_matrix mbox) {

	return DTL_set_AV_mbox(crit,mbox,mbox);
	}


/* Since DTL_set_AV_modal sets both box and modal in the same call,
 * there is no need for scale change signalling or sync. Thus, the
 * change signal should be seen as purely informative in this case.
 * Mode: 0 = set mbox
 *      +1 = clear mbox before setting
 *      +2 = also set box */

rcode DTLAPI DTL_set_AV_modal(int crit, int mode, bool rev, bool renorm, h_matrix lobox, h_matrix modalx, h_matrix upbox) {
	rcode rc,rc2;
	double v_min,v_max,scaling;

#ifdef LOG_AS
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_set_AV_modal(%d,%d,%d,%d)\n",crit,mode,rev,renorm);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameter */
	if (load_df1(crit))
		return DTL_CRIT_UNKNOWN;
	/* Transform input to [0,1] scale */
	if (DTL_error2(rc = dtl_trf_AV_input(crit,TRF_MODAL,rev,lobox,modalx,upbox,&v_min,&v_max)))
		return rc;
	/* If scale changed, mbox must be cleared */
	if (rc == DTL_SCALE_CHANGE)
		mode |= 0x01;
	/* Enter new [0,1]-scaled box and modal values */
	if (DTL_error2(rc2 = DTL_set_V_modal(crit,mode,av_lobox,av_modalx,av_upbox)))
		return rc2;
	if (rc == DTL_SCALE_CHANGE) {
		/* Calculate scaling factor (span > 0 -> no div-by-zero risk) */
		scaling = (v_max-v_min) / (uf->av_max[crit]-uf->av_min[crit]);
		/* Set new autoscale endpoints */
		uf->av_min[crit] = v_min;
		uf->av_max[crit] = v_max;
		/* Renormalise criteria weights if requested */
		if (renorm)
			if (DTL_error2(rc2 = dtl_AV_renorm_W_base(crit,scaling)))
				return rc2;
		}
	/* Signal whether autoscale has changed or not */
	return rc;
	}


/* Set criterion scale - only for test/debug use */

rcode DTLAPI DTI_set_AV_crit_scale(int crit, double v_min, double v_max) {

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameter */
	if (check_df1(crit))
		return DTL_CRIT_UNKNOWN;
	return dtl_set_AV_scale(crit,v_min,v_max);
	}


/* Copy one criterion scale to another - only for internal use */

rcode dtl_copy_AV_crit_scale(int cr_from, int cr_to) {
	double avf_min,avf_max;

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check and get input parameters */
	if (check_df1(cr_from))
		return DTL_CRIT_UNKNOWN;
	avf_min = uf->av_min[cr_from];
	avf_max = uf->av_max[cr_from];
	return dtl_set_AV_scale(cr_to,avf_min,avf_max);
	}


/* While DTL_set_AV_box and DTL_set_AV_modal turns the autoscale
 * on for a particular criterion, this call instead turns it off.
 * Only allowed for testing, changes meaning of value statements. */

rcode DTLAPI DTI_reset_AV_crit_scale(int crit) {

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameter */
	if (check_df1(crit))
		return DTL_CRIT_UNKNOWN;
	return dtl_set_AV_scale(crit,0.0,1.0);
	}


rcode DTLAPI DTL_get_AV_crit_scale(int crit, double *v_min, double *v_max) {

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameter */
	if (check_df1(crit))
		return DTL_CRIT_UNKNOWN;
	/* Return autoscale input limits */
	return dtl_get_AV_scale(crit,v_min,v_max);
	}


 /*********************************************************
  *
  *  MC scale functions (set,copy,reset,get)
  *
  *********************************************************/

#define AMC 0 // autoscale multi-criteria

rcode DTLAPI DTL_set_AV_MC_scale(double v_min, double v_max) {

	return dtl_set_AV_scale(AMC,v_min,v_max);
	}


rcode DTLAPI DTL_copy_AV_MC_scale(int crit) {

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameter */
	if (check_df1(crit))
		return DTL_CRIT_UNKNOWN;
	return dtl_set_AV_scale(AMC,uf->av_min[crit],uf->av_max[crit]);
	}


rcode DTLAPI DTL_reset_AV_MC_scale() {

	/* Reset autoscale MC limits */
	return dtl_set_AV_scale(AMC,0.0,1.0);
	}


rcode DTLAPI DTL_get_AV_MC_scale(double *v_min, double *v_max) {

	/* Return autoscale MC limits */
	return dtl_get_AV_scale(AMC,v_min,v_max);
	}


 /*********************************************************
  *
  *  Scale ratio calculator (a.k.a. trade-off measure)
  *
  *  Mode: 0 = crit level, unweighted scale span ratio
  *        1 = MC level, ratio weighted with midpoints
  *
  *  NOTE: Ratio is negative if one scale is reversed ->
  *        use |ratio| for a leverage measure btw crits
  *
  *********************************************************/

static h_vector Wlobo,Wmid,Wupbo;

static bool dti_AV_rev_scale(int crit) {

	/* Return scale direction */
	return uf->av_min[crit] > uf->av_max[crit];
	}


static bool dti_AV_opposite_scales(int crit1, int crit2) {
	bool rev1,rev2;

	/* Calc scale directions */
	rev1 = dti_AV_rev_scale(crit1);
	rev2 = dti_AV_rev_scale(crit2);
	/* Return scale direction diff */
	return rev1 != rev2;
	}


static rcode dti_AV_scale_ratio(int c_from, int c_to, int mode, double *ratio) {
	rcode rc;
	double span_from,span_to;

	/* Get criteria weights (MC is 1.0) */
	if (rc = DTL_get_W_hull(1,Wlobo,Wmid,Wupbo))
		return rc;
	if (!c_to) { // MC = full weight tree
		r2t[1][0] = 0;
		Wmid[0] = 1.0;
		}
	/* Calculate ratio */
	span_from = _fabs(uf->av_max[c_from]-uf->av_min[c_from]);
	span_to = _fabs(uf->av_max[c_to]-uf->av_min[c_to]);
	if (!span_from)
		*ratio = -1.0; // infinity/impossible (Dirac)
	else {
		*ratio = span_to/span_from;
		if (mode) // global/MC level comparison
			if (!Wmid[r2t[1][c_to]])
				*ratio = -1.0; // infinity/impossible (Dirac)
			else
				*ratio *= Wmid[r2t[1][c_from]]/Wmid[r2t[1][c_to]];
		}
	return DTL_OK;
	}


/* mode 0 = local scales compared, mode 1 = global MC scale via weights */

rcode DTLAPI DTI_AV_scale_ratio(int c_from, int c_to, int mode, double *ratio) {
	rcode rc;

	/* Check if function can start */
	_certify_ptr(ratio,201);
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	if (!PM)
		return DTL_NOT_ALLOWED;
	/* Check input parameters */
	if (!c_from) // MC can only be output
		return DTL_INPUT_ERROR;
	if (c_from == c_to)
		return DTL_INPUT_ERROR;
	if (check_df1(c_from))
		return DTL_CRIT_UNKNOWN;
	if (check_df0(c_to)) // MC allowed as output
		return DTL_CRIT_UNKNOWN;
	/* Get ratio (or infinity) */
	if (rc = dti_AV_scale_ratio(c_from,c_to,mode,ratio))
		return rc;
	if (*ratio == -1.0) { // infinity due to Dirac scale
#ifdef INFINITY
		*ratio = INFINITY; // infinite Dirac delta density
#elif defined(DBL_MAX)
		*ratio = DBL_MAX;
#else
		*ratio = 1.79E308; // "almost infinite" (max is around 1.7977E308 in IEEE754 binary64)
#endif
		return DTL_INFINITE_MASS;
		}
	else if (dti_AV_opposite_scales(c_from,c_to))
		*ratio *= -1.0;
	return DTL_OK;
	}


 /*********************************************************
  *
  *  Scale conversion functions (vector,value,interval)
  *
  *                 rel neg |x|  scale type   acronym
  *                 --- --- ---  ----------   -------
  *  Type field: 1   N   N   N   absolute     ABS_SCALE
  *              2   Y   Y   N   difference   DIFF_SCALE
  *              3   Y   N   Y   distance     DIST_SCALE
  *          (*) 4   Y   Y   Y   reverse diff REVD_SCALE
  *
  *  Legend: rel = relative scale (absolute is default)
  *          neg = allow negative norm input [-1,0]
  *          |x| = trim to non-negative output [0,a]
  *
  *  (*) = Reverse difference scale: it treats a reverse
  *        scale [b,a] (b>a) as if it were a scale [a,b]
  *
  *  NOTE: Interval calls handle reverse scales, use
  *        them instead of vector calls for intervals
  *
  *  NOTE: The parameter 'crit' is a criterion number,
  *        not an evaluation (partial = crit<0) marker,
  *        unless the config param PART_TREE is active.
  *
  *********************************************************/

/* From normalised [0,1] to user [a,b] scale. For type 1 (absolute) and type 3
 * (relative) the normalised scale is [0,1] while for types 2 & 4 (relative) it
 * is [-1,1]. This intentionally leaves im-V-node content undefined for relative
 * scales since it should not be accessed or used (contains only empty markers).
 * For type 1 (absolute), -1.0 indicates an empty slot for compatibility reasons. */

rcode DTLAPI DTL_get_AV_user_vector(int crit, int type, int size, double v_val[], double av_val[]) {
	int i;
	double v_min,v_max;

#ifdef LOG_AS
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_AV_user_vector(%d,%d,%d)\n",crit,type,size);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	_init_assert();
	_certify_ptr(v_val,201);
	_certify_ptr(av_val,202);
	_dtl_assert(v_val!=av_val,201);
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameters */
#ifdef PART_TREE
	crit = max(crit,0); // partial eval crit nbr
#endif
	if (check_df0(crit))
		return DTL_CRIT_UNKNOWN;
	if ((type < 1) || (type > 4))
		return DTL_INPUT_ERROR;
	if ((size < 1) || (size > MAX_NODES))
		return DTL_INPUT_ERROR;
	/* Get autoscale limits */
	v_min = uf->av_min[crit];
	v_max = uf->av_max[crit];
	/* Convert vector from internal [0,1] to external [a,b] scale */
	for (i=0; i<size; i++) {
		if ((type<2) && (v_val[i]==-1.0)) // empty slot maps onto norm=0.0
			av_val[i] = v_min; // to facilitate compatible inverse calling
		else { // occupied slot
			if ((v_val[i] < (type&1?0.0:-1.0)) || (v_val[i] > 1.0))
				return DTL_INPUT_ERROR;
			av_val[i] = v_val[i]*(type<3?v_max-v_min:_fabs(v_max-v_min))+(type<2?v_min:0.0);
#ifdef LOG_AS
			if (cst_ext) {
				sprintf(msg,"  v[%d] %lg -> %lg\n",i,v_val[i],av_val[i]);
				cst_log(msg);
				}
#endif
			}
		}
	return DTL_OK;
	}


rcode DTLAPI DTL_get_AV_user_value(int crit, int type, double v_val, double *av_val) {

	/* Convert single value */
	return DTL_get_AV_user_vector(crit,type,1,&v_val,av_val);
	}


rcode DTLAPI DTL_get_AV_user_intervals(int crit, int type, int size, double v_lobo[], 
		double v_upbo[], double av_lobo[], double av_upbo[]) {
	rcode rc;
	bool rev;

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameter */
#ifdef PART_TREE
	crit = max(crit,0); // partial eval crit nbr
#endif
	if (check_df0(crit))
		return DTL_CRIT_UNKNOWN;
	/* Check if reverse scale */
	rev = uf->av_min[crit] > uf->av_max[crit];
	/* Sort upper and lower bounds (types 3 & 4 already sorted due to fabs above) */
	if (rc = DTL_get_AV_user_vector(crit,type,size,v_lobo,rev&&type<3?av_upbo:av_lobo))
		return rc;
	if (rc = DTL_get_AV_user_vector(crit,type,size,v_upbo,rev&&type<3?av_lobo:av_upbo))
		return rc;
	return DTL_OK;
	}


rcode DTLAPI DTL_get_AV_user_interval(int crit, int type, double v_lobo, double v_upbo, 
		double *av_lobo, double *av_upbo) {

	/* Convert single interval */
	return DTL_get_AV_user_intervals(crit,type,1,&v_lobo,&v_upbo,av_lobo,av_upbo);
	}


/* From user [a,b] to normalised [0,1] scale (range [-1,1] for types 2 & 4) */

rcode DTLAPI DTL_get_AV_norm_vector(int crit, int type, int size, double av_val[], double v_val[]) {
	int i;
	double v_min,v_max;

#ifdef LOG_AS
	/* Log function call */
	if (cst_ext) {
		sprintf(msg,"DTL_get_AV_norm_vector(%d,%d,%d)\n",crit,type,size);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	_init_assert();
	_certify_ptr(av_val,201);
	_certify_ptr(v_val,202);
	_dtl_assert(av_val!=v_val,201);
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameters */
#ifdef PART_TREE
	crit = max(crit,0); // partial eval crit nbr
#endif
	if (check_df0(crit))
		return DTL_CRIT_UNKNOWN;
	if ((type < 1) || (type > 4))
		return DTL_INPUT_ERROR;
	if ((size < 1) || (size > MAX_NODES))
		return DTL_INPUT_ERROR;
	/* Get autoscale limits */
	v_min = uf->av_min[crit];
	v_max = uf->av_max[crit];
	/* Convert vector from external [a,b] to internal [0,1] scale */
	for (i=0; i<size; i++) {
		v_val[i] = (av_val[i]-(type<2?v_min:0.0))/(type<3?v_max-v_min:_fabs(v_max-v_min));
		if ((v_val[i] < (type&1?0.0:-1.0)) || (v_val[i] > 1.0))
			return DTL_INPUT_ERROR;
#ifdef LOG_AS
		if (cst_ext) {
			sprintf(msg,"  v[%d] %lg -> %lg\n",i,av_val[i],v_val[i]);
			cst_log(msg);
			}
#endif
		}
	return DTL_OK;
	}


rcode DTLAPI DTL_get_AV_norm_value(int crit, int type, double av_val, double *v_val) {

	/* Convert single value */
	return DTL_get_AV_norm_vector(crit,type,1,&av_val,v_val);
	}


rcode DTLAPI DTL_get_AV_norm_intervals(int crit, int type, int size, double av_lobo[], double av_upbo[], 
		double v_lobo[], double v_upbo[]) {
	rcode rc;
	bool rev;

	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameter */
#ifdef PART_TREE
	crit = max(crit,0); // partial eval crit nbr
#endif
	if (check_df0(crit))
		return DTL_CRIT_UNKNOWN;
	/* Check if reverse scale */
	rev = uf->av_min[crit] > uf->av_max[crit];
	/* Sort upper and lower bounds (types 3 & 4 already sorted due to fabs above) */
	if (rc = DTL_get_AV_norm_vector(crit,type,size,av_lobo,rev&&type<3?v_upbo:v_lobo))
		return rc;
	if (rc = DTL_get_AV_norm_vector(crit,type,size,av_upbo,rev&&type<3?v_lobo:v_upbo))
		return rc;
	return DTL_OK;
	}


rcode DTLAPI DTL_get_AV_norm_interval(int crit, int type, double av_lobo, double av_upbo, 
		double *v_lobo, double *v_upbo) {

	/* Convert single interval */
	return DTL_get_AV_norm_intervals(crit,type,1,&av_lobo,&av_upbo,v_lobo,v_upbo);
	}


 /*******************************************************
  *
  *  Scale checking functions (user and norm scales)
  *
  *  Accept a list of up to ten values to range check.
  *  Duplicated code since varargs cannot call varargs.
  *
  *  Return code: <0 = error
  *                0 = at least one value out of range
  *                1 = all values are within the range
  *
  *  NOTE1: Use with care since type check is minimal.
  *         The stack is interpreted at stack pointer.
  *
  *  NOTE2: See type specification above for type field.
  *         But differently executed here since this is
  *         applied to output from the functions above.
  *
  *******************************************************/

#define MAX_AV_CHECK 10
#define MAX_AV_CHECK_SHORT 4

rcode DTLAPI DTL_check_AV_user_values(int crit, int type, int count, ...) {
	int i;
	double user,v_min,v_max,vi_min,vi_max;
	va_list ap;

#ifdef LOG_AS
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_check_AV_user_values(%d,%d,%d)\n",crit,type,count);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameters */
#ifdef PART_TREE
	crit = max(crit,0); // partial eval crit nbr
#endif
	if (check_df0(crit))
		return DTL_CRIT_UNKNOWN;
	if ((type < 1) || (type > 4))
		return DTL_INPUT_ERROR;
	if ((count < 1) || (count > MAX_AV_CHECK))
		return DTL_INPUT_ERROR;
	/* Autoscale endpoints */
	v_min = uf->av_min[crit];
	v_max = uf->av_max[crit];
	if (type > 1) {
		/* Relative scale */
		vi_max = _fabs(v_max-v_min);
		vi_min = type&1?0.0:-vi_max;
		}
	else {
		/* Absolute scale */
		vi_min = min(v_min,v_max);
		vi_max = max(v_min,v_max);
		}
	/* Parse list of user values */
	va_start(ap,count);
	for (i=0; i<count; i++) {
		user = va_arg(ap,double);
		if ((user < vi_min) || (user > vi_max))
			return FALSE;
		}
	va_end(ap);
	return TRUE;
	}


rcode DTLAPI DTL_check_AV_norm_values(int type, int count, ...) {
	int i;
	double norm;
	va_list ap;

#ifdef LOG_AS
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_check_AV_norm_values(%d,%d)\n",type,count);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameters */
	if ((type < 1) || (type > 4))
		return DTL_INPUT_ERROR;
	if ((count < 1) || (count > MAX_AV_CHECK))
		return DTL_INPUT_ERROR;
	/* Parse list of norm values */
	va_start(ap,count);
	for (i=0; i<count; i++) {
		norm = va_arg(ap,double);
		if ((type<2) && (norm==-1.0)) // empty slot
			continue;
		if ((norm<(type&1?0.0:-1.0)) || (norm>1.0))
			return FALSE;
		}
	va_end(ap);
	return TRUE;
	}


/* DTI functions for test and debug use */

/* Compound value range check call
 *
 * Handles two interval pairs or four freestanding values.
 * Varargs could be replaced by vectors if types are same.
 * Use type += 0x08 for norm scale instead of user scale.
 *
 * Return code: <0 = error
 *               0 = out of range
 *               1 = within range */

rcode DTLAPI DTI_check_AV_values(int crit, int type, int count, ...) {
	int i,res;
	double val[MAX_AV_CHECK_SHORT];
	va_list ap;

	/* Check input parameter */
	if ((count < 1) || (count > MAX_AV_CHECK_SHORT))
		return DTL_INPUT_ERROR;
	/* Parse list of numbers */
	va_start(ap,count);
	for (i=0; i<count; i++)
		val[i] = va_arg(ap,double);
	for (; i<MAX_AV_CHECK_SHORT; i++)
		val[i] = 0.0;
	va_end(ap);
#ifdef PART_TREE
	crit = max(crit,0); // partial eval crit nbr
#endif
	if (type&0x08) // norm scale indicator
		res = DTL_check_AV_norm_values(type&0x07,count,val[0],val[1],val[2],val[3]);
	else
		res = DTL_check_AV_user_values(crit,type,count,val[0],val[1],val[2],val[3]);
	return res;
	}


/* Check if MC or crit has retained default [0,1] scale or changed it
 *
 * Return code: <0 = error
 *               0 = scale changed
 *               1 = scale preserved */

rcode DTLAPI DTI_is_AV_default_scale(int crit) {
	rcode rc;
	double v_min,v_max;

#ifdef LOG_AS
	/* Log function call */
	if (cst_ext) {
		sprintf(msg,"DTI_is_AV_default_scale(%d)\n",crit);
		cst_log(msg);
		}
#endif
	/* Check if function can start */
	if (!frame_loaded)
		return DTL_FRAME_NOT_LOADED;
	/* Check input parameter */
#ifdef PART_TREE
	crit = max(crit,0); // partial eval crit nbr
#endif
	if (check_df0(crit)) // MC allowed
		return DTL_CRIT_UNKNOWN;
	/* Get and compare scale endpoints */
	if (rc = dtl_get_AV_scale(crit,&v_min,&v_max))
		return rc;
	return v_min==0.0 && v_max==1.0;
	}
