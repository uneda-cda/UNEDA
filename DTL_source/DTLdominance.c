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
 *   File: DTLdominance.c
 *
 *   Purpose: dominance evaluation of alternatives
 *
 *
 *   Functions exported outside DTL
 *   ------------------------------
 *   DTL_get_dominance
 *   DTL_get_dominance_matrix
 *   DTL_get_dominance_nt_matrix
 *   DTL_get_dominance_rank
 *   DTL_get_cardinal_dominance_matrix
 *   DTL_get_abs_dominance_matrix
 *
 *   Functions outside of module, inside DTL
 *   ---------------------------------------
 *   dtl_get_dominance
 *
 *   Functions internal to module
 *   ----------------------------
 *   add_dom
 *   abs_dom
 *
 */


 /*********************************************************
  *
  *  Configuration flags
  *
  *********************************************************/

#define noSTRICT_DOM // strict dominance for non-transitive map
#define noABSDOM_ONLY_PM // compile flag to stop PS frames


 /*********************************************************
  *
  *  Belief dominance evaluations
  *
  *********************************************************/

#define DOMINANCE_LIMIT 1.0E-3 // significant diff (visible rounded to 3 decimals)

static e_matrix e_result1,e_result2;

static double add_dom(double diff, double sum, int *dom) {

	sum += diff;
	if (diff > DOMINANCE_LIMIT)
		*dom |= 0x01; // alt 1
	else if (diff < -DOMINANCE_LIMIT)
		*dom |= 0x02; // alt 2
	return sum;
	}


// DTL layer 1: at DTL API level

/* Comparing cdfs is done in a reverse way in this function.
 * Usually: fixed steps of EV, compare cdfs, lower = better.
 * Here: fixed steps of cdf, compare EVs, higher = better.
 * Same result - they cross (2nd order) or not (1st order).
 *
 * NOTE: if ever porting to DMC, change to fixed steps of EV.
 * Or, at least, cater for both to be able to compare them. */

rcode dtl_get_dominance(int crit, int Ai, int Aj, double *cd_value, int *d_order) {
	rcode rc;
	int i,dom,cst_global;
	double sum;

	/* Check input parameters */
	if (load_df00(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if (Ai == Aj)
		return dtl_error(DTL_INPUT_ERROR);
	/* Collect both cdfs */
	cst_global = cst_on;
	cst_on = FALSE;
	if (rc = evaluate_frameset(crit,E_PSI,Ai,0,e_result1)) {
		cst_on = cst_global;
		return dtl_error(rc); // already at level 1?
		}
	// note: expansion type 1 = towards cdf 50%, not mass point
	expand_eval_result1(crit,0,e_result1);
	if (rc = evaluate_frameset(crit,E_PSI,Aj,0,e_result2)) {
		cst_on = cst_global;
		return dtl_error(rc); // already at level 1?
		}
	expand_eval_result1(crit,0,e_result2);
	/* dom=0: not determined
		 dom=1: Ai found superior
		 dom=2: Aj found superior
		 dom=3: neither dominates fully (not 1st order) */
	dom = 0;
	sum = add_dom(e_result1[E_MID][MAX_RESULTSTEPS-1]-e_result2[E_MID][MAX_RESULTSTEPS-1],0.0,&dom);
	for (i=0; i<MAX_RESULTSTEPS-1; i++) {
		sum = add_dom(e_result1[E_MIN][i]-e_result2[E_MIN][i],sum,&dom);
		sum = add_dom(e_result1[E_MAX][i]-e_result2[E_MAX][i],sum,&dom);
		}
	/* Should roughly yield the EV difference between the alternatives (not exactly, we
	 * sample in only 21 steps + use expand_1=cdf steps, not EV + add max&min, nor mid).
	 * NOTE: we use the diff btw psi instead of the delta to get at the diff in cdf, not EV.*/
	*cd_value = sum/(double)(2*MAX_RESULTSTEPS-1);
	if (fabs(*cd_value) >= DOMINANCE_LIMIT)
		*d_order = dom>1?dom-1:dom; // 1=1st order dom, 2=2nd order dom
	else
		*d_order = 0; // 0=no dom
	eval_cache_invalidate();
	cst_on = cst_global;
	return DTL_OK;
	}


rcode DTLAPI DTL_get_dominance(int crit, int Ai, int Aj, double *cd_value, int *d_order) {
	rcode rc;

	/* Begin single thread semaphore */
	_smx_begin("GDOM");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_dominance(%d,%d,%d)\n",crit,Ai,Aj);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(cd_value,1);
	_certify_ptr(d_order,2);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Fetch the belief dominance */
	if (rc = dtl_get_dominance(crit,Ai,Aj,cd_value,d_order))
		return rc;
	/* Log function result */
	if (cst_ext) {
		sprintf(msg," %+.3lf %d-order dominance\n",*cd_value,*d_order);
		cst_log(msg);
		}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


rcode DTLAPI DTL_get_dominance_matrix(int crit, double threshold, ai_matrix dominance_mx) {
	rcode rc;
	int Ai,Aj,d_order;
	double cd_value;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("GDOMX");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_dominance_matrix(%d,%.3lf)\n",crit,threshold);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(dominance_mx,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df00(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if ((threshold < 0.0) || (threshold > 0.1))
		return dtl_error(DTL_INPUT_ERROR);
	df = uf->df;
	/* Collect the belief dominances */
	dtl_abort_init();
	for (Ai=1; Ai<=df->n_alts; Ai++)
		dominance_mx[Ai][Ai] = 0; // cannot dominate itself
	for (Ai=1; Ai<df->n_alts; Ai++)
		for (Aj=Ai+1; Aj<=df->n_alts; Aj++) {
			dtl_abort_check();
			if (rc = dtl_get_dominance(crit,Ai,Aj,&cd_value,&d_order))
				return rc;
			if (cd_value > threshold) {
				dominance_mx[Ai][Aj] = d_order; // Ai dominates
				dominance_mx[Aj][Ai] = 0;       // Aj dominated
				}
			else if (cd_value < -threshold) {
				dominance_mx[Ai][Aj] = 0;       // Ai dominated
				dominance_mx[Aj][Ai] = d_order; // Aj dominates
				}
			else {
				dominance_mx[Ai][Aj] = 0; // Ai not dominated
				dominance_mx[Aj][Ai] = 0; // Aj not dominated
				}
			}
	/* Log function result */
	if (cst_ext) {
		for (Ai=1; Ai<=df->n_alts; Ai++) {
			for (Aj=1; Aj<=df->n_alts; Aj++)
				if (Ai==Aj)
					cst_log("  ");
				else {
					sprintf(msg," %1d",dominance_mx[Ai][Aj]);
					cst_log(msg);
					}
			cst_log("\n");
			}
		}
	/* End single thread semaphore */
	if (!rc)
		_smx_end();
	return rc;
	}


static ai_matrix t_mx,nt_mx;

// DTL layer 0: above DTL proper

rcode DTLAPI DTL_get_dominance_nt_matrix(int crit, double threshold, ai_matrix dominance_mx) {
	rcode rc;
	int Ai,Aj,Ak;
	struct d_frame *df;

	if (dominance_mx == t_mx) // will not map correctly
		return DTL_INTERNAL_ERROR;
	if (rc = DTL_get_dominance_matrix(crit,threshold,t_mx))
		return rc;
	df = uf->df;
	/* Transform the belief dominances */
	for (Ai=1; Ai<=df->n_alts; Ai++)
		for (Aj=1; Aj<=df->n_alts; Aj++)
			dominance_mx[Ai][Aj] = t_mx[Ai][Aj];
	for (Ai=1; Ai<=df->n_alts; Ai++)
		for (Aj=1; Aj<=df->n_alts; Aj++)
#ifdef STRICT_DOM
			if (dominance_mx[Ai][Aj]==1)
				for (Ak=1; Ak<=df->n_alts; Ak++)
					if ((t_mx[Ai][Ak]==1) && (t_mx[Ak][Aj]==1)) {
#else
			if (dominance_mx[Ai][Aj])
				for (Ak=1; Ak<=df->n_alts; Ak++)
					if (((dominance_mx[Ai][Aj]==1) && (t_mx[Ai][Ak]==1) && (t_mx[Ak][Aj]==1)) ||
							((dominance_mx[Ai][Aj]==2) && t_mx[Ai][Ak] && t_mx[Ak][Aj])) {
#endif
						/* Redundant by transitivity */
						dominance_mx[Ai][Aj] = 0;
						break;
						}
	return DTL_OK;
	}


static ai_vector active;

/* Dominance ranking
 * -----------------
 *
 * Mode parameter:
 * mode: 0 = group ranking
 *       1 = olympic ranking
 *       2 = hard/sharp ranking - use when unique indices are required
 *
 * Dominance parameter:
 * dmode: 0 = weak dominance   -> both 1-order & 2-order counts (recommended)
 *        1 = strong dominance -> only 1-order counts (sometimes unintuitive) */

// DTL layer 0: above DTL proper

rcode DTLAPI DTL_get_dominance_rank(int crit, int mode, int dmode, double threshold, ai_vector dom_rank) {
	rcode rc;
	int Ai,Aj,loop,level,remaining;
	struct d_frame *df;

	if ((mode < 0) || (mode > 2))
		return DTL_INPUT_ERROR;
	if ((dmode < 0) || (dmode > 1))
		return DTL_INPUT_ERROR;
	if (rc = DTL_get_dominance_nt_matrix(crit,threshold,nt_mx))
		return rc;
	df = uf->df;
	/* Reset all levels */
	for (Ai=1; Ai<=df->n_alts; Ai++) {
		active[Ai] = TRUE;
		dom_rank[Ai] = 0;
		}
	remaining = df->n_alts;
	loop = level = 1; // loop & level count
	while (remaining) {
		for (Aj=1; Aj<=df->n_alts; Aj++)
			if (active[Aj]) {
				// scan all active on this level by column
				for (Ai=1; Ai<=df->n_alts; Ai++)
					// if dominated by someone who is active -> not at this level
					if (active[Ai] && nt_mx[Ai][Aj])
						if (!dmode || (nt_mx[Ai][Aj]==1))
							break;
				if (Ai > df->n_alts) {
					// not dominated by anyone active -> belongs to this level
					dom_rank[Aj] = mode>1?df->n_alts+1-remaining:(mode?level:loop);
					remaining--;
					}
				}
		// update active rows and level at end of loop, else mixup
		for (Ai=1; Ai<=df->n_alts; Ai++)
			active[Ai] = (dom_rank[Ai]==0);
		level = df->n_alts+1-remaining;
		if (++loop > MAX_ALTS)
			return DTL_INTERNAL_ERROR;
		}
	return DTL_OK;
	}


rcode DTLAPI DTL_get_cardinal_dominance_matrix(int crit, int dmode, double threshold, ar_matrix cardinal_mx) {
	rcode rc;
	int Ai,Aj,total;
	double dominance;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("GCDOMX");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_cardinal_dominance_matrix(%d,%d,%.3lf)\n",crit,dmode,threshold);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(cardinal_mx,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
	/* Check input parameters */
	if (load_df00(crit))
		return dtl_error(DTL_CRIT_UNKNOWN);
	if ((dmode < 0) || (dmode > 1))
		return dtl_error(DTL_INPUT_ERROR);
	if ((threshold < 0.0) || (threshold > 0.1))
		return dtl_error(DTL_INPUT_ERROR);
	df = uf->df;
	/* Collect the belief dominances */
	for (Ai=1; Ai<=df->n_alts; Ai++)
		cardinal_mx[Ai][Ai] = 0.0; // cannot dominate itself
	for (Ai=1; Ai<df->n_alts; Ai++)
		for (Aj=Ai+1; Aj<=df->n_alts; Aj++) {
			if (rc = dtl_get_dominance(crit,Ai,Aj,&dominance,&total))
				return rc;
			if ((total==1) || (total && !dmode)) {
				if (dominance > threshold) {
					cardinal_mx[Ai][Aj] = dominance;  // Ai dominates
					cardinal_mx[Aj][Ai] = 0.0;        // Aj dominated
					}
				else if (dominance < -threshold) {
					cardinal_mx[Ai][Aj] = 0.0;        // Ai dominated
					cardinal_mx[Aj][Ai] = -dominance; // Aj dominates
					}
				else {
					cardinal_mx[Ai][Aj] = 0.0; // Ai not dominated
					cardinal_mx[Aj][Ai] = 0.0; // Aj not dominated
					}
				}
			else {
				cardinal_mx[Ai][Aj] = 0.0; // Ai not dominated
				cardinal_mx[Aj][Ai] = 0.0; // Aj not dominated
				}
			}
	/* Log function result */
	if (cst_ext) {
		for (Ai=1; Ai<=df->n_alts; Ai++) {
			for (Aj=1; Aj<=df->n_alts; Aj++)
				if (Ai==Aj)
					cst_log("      ");
				else {
					sprintf(msg," %.3lf",cardinal_mx[Ai][Aj]);
					cst_log(msg);
					}
			cst_log("\n");
			}
		}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}


/* Abs dom: 2-dom overshadows 1-dom
 *          no-dom overshadows both 1 and 2
 *          the state chain is 1 -> 2 -> 0

            cur
           0 1 2
          -------
       0 | 0 0 0
   old 1 | 0 1 2
       2 | 0 2 2

 * State machine: old x cur -> new */

static int abs_dom(int old, int cur) {

	switch (old) {
		case 2:  if (cur==1) return old;
		default: return old?cur:old;
		}
	}

/* Abs sum: 2-dom overshadows no-dom
 *          1-dom overshadows both no and 2
 *          the state chain is 0 -> 2 -> 1

            cur
           0 1 2
          -------
       0 | 0 1 2
   old 1 | 1 1 1
       2 | 2 1 2

 * State machine: old x cur -> new */

static int abs_sum(int old, int cur) {

	switch (old) {
		case 2:  if (!cur) return old;
		default: return old==1?old:cur;
		}
	}


/* Absolute dominance
 * ------------------
 *
 * Absolute dominance means regardless of criteria weights.
 * Used to filter out alternatives before an analysis begins.
 * The procedure is asymmetric in that those to be excluded
 * are truly dominated while the contenders often might but
 * need not contain the best alternative. It is only meaningful
 * for MC frames but could let PS frames in for compatibility.
 * The latter is controlled by the ABSDOM_ONLY_PM parameter.
 *
 * Dmode: 0 = 2nd order dominance -> 2nd order also counts (default)
 *        1 = 1st order dominance -> only 1st order counts
 *
 * Output matrix dmx:
 *   dmx[0][0] >0 if any dominance found (1=yes)
 *   dmx[i][0] >0 if alt Ai dominates the other alts
 *   dmx[0][j] >0 if alt Aj is dominated by the other alts
 *   dmx[i][j] >0 if alt Ai dominates alt Aj (1=1st order, 2=2nd order)
 *
 * NOTE: This does not measure stochastic dominance in P-trees but
 *       is instead intended for weight-independent EV dominance.
 */

rcode DTLAPI DTL_get_abs_dominance_matrix(int dmode, double threshold, ai_matrix dominance_mx) {
	rcode rc;
	int Ai,Aj,k,d_order,doms,n_crit;
	double cd_value;
	struct d_frame *df;

	/* Begin single thread semaphore */
	_smx_begin("GADOMX");
	/* Log function call */
	if (cst_on) {
		sprintf(msg,"DTL_get_abs_dominance_matrix(%d,%.3lf)\n",dmode,threshold);
		cst_log(msg);
		}
	/* Check if function can start */
	_certify_ptr(dominance_mx,1);
	if (!frame_loaded)
		return dtl_error(DTL_FRAME_NOT_LOADED);
#ifdef ABSDOM_ONLY_PM
	if (!PM)
		return dtl_error(DTL_NOT_ALLOWED);
#endif
	/* Check input parameters */
	if ((threshold < 0.0) || (threshold > 0.1))
		return dtl_error(DTL_INPUT_ERROR);
	if (threshold < DTL_EPS)
		threshold = DTL_EPS;
	if ((n_crit = DTL_nbr_of_crit()) < DTL_OK)
		return dtl_error(n_crit);
	df = uf->df;
	/* Collect the absolute belief dominances */
	dtl_abort_init();
	/* Init matrix to full dominance - disprove in belief loop */
	for (Ai=1; Ai<=df->n_alts; Ai++)
		for (Aj=1; Aj<=df->n_alts; Aj++)
			if (Ai==Aj)
				dominance_mx[Ai][Ai] = 0; // cannot dominate itself
			else
				dominance_mx[Ai][Aj] = 1; // posit 1-order dominance
	/* Check every pair (Ai,Aj) for belief dominance */
	for (Ai=1; Ai<=df->n_alts; Ai++) {
		dtl_abort_check();
		for (Aj=1; Aj<=df->n_alts; Aj++)
			if (Ai!=Aj) // not self
				for (k=1; k<=n_crit; k++) { // for each criterion
#ifdef ABSDOM_ONLY_PM
					if (dtl_is_shadow_crit(k))
#else
					if ((n_crit>1) && dtl_is_shadow_crit(k))
#endif
						cd_value = 0.0;
					else if (rc = dtl_get_dominance(k,Ai,Aj,&cd_value,&d_order))
						return dtl_error(rc);
					if (cd_value > threshold) // Ai dominates
						dominance_mx[Ai][Aj] = abs_dom(dominance_mx[Ai][Aj],d_order);
					else { // Ai not dominating
						dominance_mx[Ai][Aj] = 0;
						break; // end value of state machine -> final verdict
						}
					}
		}
	/* Scan every alt col Aj for being dominated (= can be excluded) */
	doms = 0; // reset dominated
	for (Aj=1; Aj<=df->n_alts; Aj++) {
		dominance_mx[0][Aj] = 0; // posit Aj not dominated (start of state machine)
		for (Ai=1; Ai<=df->n_alts; Ai++)
			if (Ai!=Aj) { // not self
				dominance_mx[0][Aj] = abs_sum(dominance_mx[0][Aj],dominance_mx[Ai][Aj]);
				if (dominance_mx[0][Aj]==1) // found to be dominated
					break; // end of state machine -> final verdict
				}
		if (dmode && (dominance_mx[0][Aj]==2)) // strong: 2nd order is not dom
			dominance_mx[0][Aj] = 0;
		doms += dominance_mx[0][Aj]?1:0; // increase dominated count
		}
	dominance_mx[0][0] = doms > 0; // dominance occurence flag (>0 = TRUE)
	/* Scan every alt row Ai for exerting dominance (= dominators) */
	for (Ai=1; Ai<=df->n_alts; Ai++) {
		if (!dominance_mx[0][Ai]) { // not itself dominated -> scan for dominance
			dominance_mx[Ai][0] = 0; // posit Ai does not dominate (start of state machine)
			for (Aj=1; Aj<=df->n_alts; Aj++)
				if (Ai!=Aj) { // not self
					dominance_mx[Ai][0] = abs_sum(dominance_mx[Ai][0],dominance_mx[Ai][Aj]);
					if (dominance_mx[Ai][0]==1) // found to be dominating
						break; // end of state machine -> final verdict
					}
			}
		else
			dominance_mx[Ai][0] = 0; // no dominance
		if (dmode && (dominance_mx[Ai][0]==2)) // strong: 2nd order is not dom
			dominance_mx[Ai][0] = 0;
		}
	/* Log function result */
#ifdef ABSDOM_ONLY_PM
	if (cst_ext) {
#else
	if ((n_crit>1) && cst_ext) {
#endif
		for (Ai=0; Ai<=df->n_alts; Ai++) {
			for (Aj=0; Aj<=df->n_alts; Aj++)
				if (Ai==Aj)
					cst_log(" -");
				else {
					sprintf(msg," %1d",dominance_mx[Ai][Aj]);
					cst_log(msg);
					}
			cst_log("\n");
			}
		}
	/* End single thread semaphore */
	_smx_end();
	return DTL_OK;
	}
