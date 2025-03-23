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
 *   File: TCLseclevel.c
 *
 *   Purpose: set security constraints/thresholds for the alternatives
 *
 *
 *   Functions exported outside TCL
 *   ------------------------------
 *   TCL_security_level
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

static d_row P_lobo,P_upbo,V_lobo,V_upbo;
static i_row strong_ixset,marked_ixset,weak_ixset;

rcode TCL_security_level(struct d_frame *df, double sec_level, 
					a_vector strong, a_vector marked, a_vector weak) {
	int Ai,Ai_begin,Ai_end,i;
	double P_strong,P_marked,P_weak;

	if (df == NULL)
		return TCL_CORRUPTED;
	if (df->watermark != D_MARK)
		return TCL_CORRUPTED;
	if (!df->attached)
		return TCL_DETACHED;
	if ((sec_level < 0.0) || (sec_level > 1.0))
		return TCL_INPUT_ERROR;

	hull_V(V_lobo,V_upbo);
	for (Ai=1; Ai<=df->n_alts; Ai++) {
		Ai_begin = get_V_start(Ai);
		Ai_end = get_V_end(Ai);
		/* Each cons can take on min and max independently */
		for (i=Ai_begin; i<=Ai_end; i++) {
			strong_ixset[i] = (V_upbo[i] < sec_level-EPS);
			marked_ixset[i] = ((V_upbo[i]+V_lobo[i])/2.0 < sec_level-EPS);
			weak_ixset[i] = (V_lobo[i] < sec_level-EPS);
			}
		/* Check if dangerous set has high probability:
		 * pick the upper hull for the dangerous set
		 * and the lower hull for the complement set. */
		/* Strong S-dominance */
		P_strong = ixset_P_min(Ai,0,strong_ixset);
		if (P_strong < -EPS)
			return TCL_INCONSISTENT;
		/* Marked S-dominance */
		P_marked  = ixset_P_min(Ai,0,marked_ixset);
		P_marked += ixset_P_max(Ai,0,marked_ixset);
		P_marked /= 2.0;
		/* Weak S-dominance */
		P_weak = ixset_P_max(Ai,0,weak_ixset);
		if (P_weak < -EPS)
			return TCL_INCONSISTENT;
		strong[Ai] = P_strong;
		marked[Ai] = P_marked;
		weak[Ai] = P_weak;
		}
	return TCL_OK;
	}
