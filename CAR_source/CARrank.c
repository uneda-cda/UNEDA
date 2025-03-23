/*
 *
 *                 _/_/_/_/       _/        _/_/_/_/
 *               _/             _/ _/      _/     _/
 *              _/            _/    _/    _/      _/
 *             _/           _/      _/   _/      _/
 *            _/           _/_/_/_/_/   _/_/_/_/
 *           _/           _/      _/   _/    _/
 *           _/          _/      _/   _/      _/
 *            _/_/_/    _/      _/   _/       _/
 *
 *
 *        The Cardinal Alternative Ranking Method (CAR)
 *        ---------------------------------------------
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
 */

/*
 *   File: CARrank.c
 *
 *   Purpose: Distance ranking functions
 *
 */

 /*********************************************************
  *
  *  Weight base distance ranking
  *
  *********************************************************/

rcode DTLAPI CAR_rank_W_base(int n_nodes, car_vector ord_crit, double dist) {
	rcode rc;
	int k,n_wts;
	double dfact;

	/* Check if function can start */
	_certify_ptr(ord_crit,101);
#ifdef LOG_CAR
	/* Log starts after pointer check */
	if (cst_ext) {
		cst_log("CAR_rank_W_base(");
		for (k=1; k<n_nodes; k++) {
			sprintf(msg,"W%d>",ord_crit[k]);
			cst_log(msg);
			}
		sprintf(msg,"W%d) -->\n",ord_crit[n_nodes]);
		cst_log(msg);
		}
#endif
	if (!frame_loaded) // protecting PS test
		return CAR_FRAME_NOT_LOADED;
	if (PS)
		return CAR_WRONG_FRAME_TYPE;
	if (fabs(dist) > 1.0)
		return CAR_INPUT_ERROR;
	if (mark_W_base())
		return CAR_SYS_CORRUPT;
	/* Check input parameters */
	if ((n_nodes < 1) || (n_nodes > MAX_NODES))
		return CAR_INPUT_ERROR;
	n_wts = DTL_nbr_of_weights();
	if (n_wts < CAR_OK)
		return n_wts;
	if (n_nodes > n_wts)
		return CAR_INPUT_ERROR;
	for (k=1; k<=n_nodes; k++)
		if ((ord_crit[k] < 1) || (ord_crit[k] > n_wts))
			return CAR_INPUT_ERROR;
#ifndef ALLOW_CAR_W_IMPURE_TREES
	if (!dtl_pure_W_tree())
		return CAR_ILLEGAL_TREE;
#endif
	for (k=2; k<=n_nodes; k++)
		if (dtl_W_node_parents(ord_crit[1],ord_crit[k]))
			return CAR_INPUT_ERROR;
	if (dtl_W_nbr_of_siblings(ord_crit[1]) != n_nodes)
		return CAR_INPUT_ERROR;
	if (n_nodes == 1)
		return CAR_OK;
	/* Build distance ranking statements and enter */
	switch (crc_method) {
		case 0:  gen_rx(n_nodes,0,1.0+min((double)n_nodes/60.0,0.25)); break; // adaptive CAR
		case 1:  gen_rs(n_nodes,0);      break;
		case 2:  gen_rr(n_nodes,0);      break;
		case 3:  gen_xr(n_nodes,0,1.35); break;
		case 4:  gen_sr(n_nodes,0);      break;
		case 5:  gen_roc(n_nodes,0);     break;
		default: return CAR_INPUT_ERROR;
		}
	for (k=1; k<=n_nodes; k++) {
		/* Positive dist = gap, negative dist = overlap */
		dfact = (dist+1.0)/2.0;
		lobox[ord_crit[k]] = k<n_nodes?dfact*crc[k]+(1.0-dfact)*crc[k+1]:dfact*crc[k];
		upbox[ord_crit[k]] = k>1?(1.0-dfact)*crc[k-1]+dfact*crc[k]:n_nodes>1?(1.0-dfact)+dfact*crc[k]:1.0;
		}
	uwstmt.n_terms = 1;
	uwstmt.sign[1] = 1;
	for (k=1; k<=n_nodes; k++) {
		uwstmt.crit[1] = ord_crit[k];
		uwstmt.lobo = lobox[ord_crit[k]];
		uwstmt.upbo = min(upbox[ord_crit[k]],1.0); // due to interpolation upwards
		rc = DTL_add_W_statement(&uwstmt);
		if (rc < CAR_OK) {
			rollback_W_base();
			return rc;
			}
		}
	/* Return number of statements */
	ord_crit[0] = DTL_nbr_of_W_stmts()-W_mark;
#ifdef LOG_CAR
	/* Log error-free completion */
	if (cst_ext) {
		sprintf(msg,"--> end of CAR_rank_W_base(%d)\n",ord_crit[0]);
		cst_log(msg);
		}
#endif
	return DTL_nbr_of_W_stmts()-W_mark;
	}


 /*********************************************************
  *
  *  Probability base ranking
  *
  *********************************************************/

rcode DTLAPI CAR_rank_P_base(int crit, int alt, int n_nodes, car_vector ord_nodes, double dist) {
	rcode rc;
	int k,t_nodes;
	int n_alts;
	double dfact;

	/* Check if function can start */
	_certify_ptr(ord_nodes,101);
#ifdef LOG_CAR
	/* Log starts after pointer check */
	if (cst_ext) {
		sprintf(msg,"CAR_rank_P_base(%d,%d,",crit,alt);
		cst_log(msg);
		for (k=1; k<n_nodes; k++) {
			sprintf(msg,"P%d.%d>",alt,ord_nodes[k]);
			cst_log(msg);
			}
		sprintf(msg,"P%d.%d) -->\n",alt,ord_nodes[n_nodes]);
		cst_log(msg);
		}
#endif
	if (!frame_loaded)
		return CAR_FRAME_NOT_LOADED;
	/* Check input parameters */
	if ((n_nodes < 1) || (n_nodes > MAX_NODES))
		return CAR_INPUT_ERROR;
	if (fabs(dist) > 1.0)
		return CAR_INPUT_ERROR;
	if (mark_P_base(crit))
		return CAR_CRIT_UNKNOWN;
	n_alts = DTL_nbr_of_alts();
	if (n_alts < DTL_OK)
		return n_alts;
	if ((alt < 1) || (alt > n_alts))
		return CAR_ALT_UNKNOWN;
	t_nodes = DTL_nbr_of_nodes(crit,alt);
	if (t_nodes < DTL_OK)
		return t_nodes;
	if (n_nodes > t_nodes)
		return CAR_INPUT_ERROR;
	for (k=1; k<=n_nodes; k++)
		if ((ord_nodes[k] < 1) || (ord_nodes[k] > t_nodes))
			return CAR_INPUT_ERROR;
	for (k=2; k<=n_nodes; k++)
		if (dtl_P_node_parents(crit,alt,ord_nodes[1],ord_nodes[k]))
			return CAR_INPUT_ERROR;
	if (dtl_P_nbr_of_siblings(crit,alt,ord_nodes[1]) != n_nodes)
		return CAR_INPUT_ERROR;
	if (n_nodes == 1)
		return CAR_OK;
	/* Build distance ranking statements and enter */
	switch (crc_method) {
		case 0:  gen_rx(n_nodes,0,1.0+min((double)n_nodes/60.0,0.25)); break; // adaptive CAR
		case 1:  gen_rs(n_nodes,0);      break;
		case 2:  gen_rr(n_nodes,0);      break;
		case 3:  gen_xr(n_nodes,0,1.35); break;
		case 4:  gen_sr(n_nodes,0);      break;
		case 5:  gen_roc(n_nodes,0);     break;
		default: return CAR_INPUT_ERROR;
		}
	/* Generate statements */
	for (k=1; k<=n_nodes; k++) {
		/* Positive dist = gap, negative dist = overlap */
		dfact = (dist+1.0)/2.0;
		lobox[ord_nodes[k]] = k<n_nodes?dfact*crc[k]+(1.0-dfact)*crc[k+1]:dfact*crc[k];
		upbox[ord_nodes[k]] = k>1?(1.0-dfact)*crc[k-1]+dfact*crc[k]:n_nodes>1?(1.0-dfact)+dfact*crc[k]:1.0;
		}
	/* Enter statements */
	ustmt.n_terms = 1;
	ustmt.sign[1] = 1;
	for (k=1; k<=n_nodes; k++) {
		ustmt.alt[1] = alt;
		ustmt.cons[1] = ord_nodes[k];
		ustmt.lobo = lobox[ord_nodes[k]];
		ustmt.upbo = min(upbox[ord_nodes[k]],1.0); // due to interpolation upwards
		rc = DTL_add_P_statement(crit,&ustmt);
		if (rc < CAR_OK) {
			rollback_P_base(crit);
			return rc;
			}
		}
	/* Return number of statements */
	ord_nodes[0] = DTL_nbr_of_P_stmts(crit)-P_mark;
#ifdef LOG_CAR
	/* Log error-free completion */
	if (cst_ext) {
		sprintf(msg,"--> end of CAR_rank_P_base(%d,%d,%d)\n",crit,alt,ord_nodes[0]);
		cst_log(msg);
		}
#endif
	return DTL_nbr_of_P_stmts(crit)-P_mark;
	}


 /***********************************************************
  *
  *  Value base ranking
  *
  *  No request for a value base rank function has been made.
  *  It will be implemented once there is a demand for it.
  *
  ***********************************************************/
