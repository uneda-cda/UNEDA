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
 *   File: TCLwarp.c
 *
 *   Purpose: calculate mass point in P-base
 *
 *   Algorithm described in Sundgren, D., Danielson, M. and Ekenberg, L.
 *   "Warp effects on calculating interval probabilities", International
 *   Journal of Approximate Reasoning 50 (2009), pp. 1360-1368.
 *
 *
 *   Functions exported outside TCL
 *   ------------------------------
 *   NONE
 *
 *   Functions outside module, inside TCL
 *   ------------------------------------
 *   adjust_vx
 *
 *   Functions internal to module
 *   ----------------------------
 *   f1_T
 *   ijar_warp_T
 *   ijar_warp_sum2
 *
 */

#ifdef TDL_COMPAT
#define VX_MAXDIM 12   // max nbr of prob vars
#define VX_CUTOFFDIM 8 // where cutoff starts
#define VX_MAXVER 2048 // max nbr of vertices inside hyper-pyramide (=2^(VX_MAXDIM-1))
#else
#define VX_MAXDIM 30    // max nbr of prob vars
#define VX_CUTOFFDIM 26 // where cutoff starts
#define VX_MAXVER 65536 // max nbr of vertices inside hyper-pyramide (=2^(VX_MAXDIM-1))
#endif
double sigma[VX_MAXVER+1];
double delta[VX_MAXVER+1];
double s_pow[VX_MAXVER+1];
int upnodes[VX_MAXVER+1];
int s_path[VX_MAXDIM+1][VX_MAXVER+1];
int s_count;
double mp_lobo[VX_MAXDIM+1];
double mp_upbo[VX_MAXDIM+1];


static int f1_T(double value, double target, int cur, int stop, int path[], int active[]) {
	int i,upn,ofl;

	/* Check if still inside hyper-pyramid */
	if (value > target-EPS)
		return 1;
	/* Check if end of path */
	if (cur == stop) {
		if (s_count >= VX_MAXVER)
			return 1;
		sigma[++s_count] = value;
		s_pow[s_count] = pow(target-value,(double)(stop-1));
		upn = 0;
		for (i=1; i<=stop; i++) {
			if (active[i]) {
				/* Record the up/down path... */
				s_path[i][s_count] = path[i];
				/* ...and the nbr of up-nodes */
				if (path[i])
					upn++;
				}
			}
		upnodes[s_count] = upn;
		return 0;
		}
	/* Next active vertex */
	while (!active[cur+1])
		cur++;
	/* Split into two paths and try both if necessary */
	path[cur+1] = 0;
	ofl = f1_T(value+mp_lobo[cur+1],target,cur+1,stop,path,active);
	path[cur+1] = 1;
	if (!ofl) f1_T(value+mp_upbo[cur+1],target,cur+1,stop,path,active);
	return ofl;
	}


/* Vertex based mass point algorithm */

static double ijar_warp_T(int x, int dim, double target, double sum2) {
	int i;
	double cij,sum1,term,tp;

	if (dim <= 1)
		/* Only one dim -> no vertices */
		return target;
	/* Sum weights from all vertices for this dimension x */
	sum1 = 0.0;
	for (i=1; i<=s_count; i++) {
		if (s_path[x][i])
			cij = mp_upbo[x];
		else
			cij = mp_lobo[x];
		term = s_pow[i] * (dim*cij + target - sigma[i]);
		if (upnodes[i] % 2)
			sum1 -= term;
		else
			sum1 += term;
		}
	tp = sum1 / (dim*sum2);
	return tp;
	}


static double ijar_warp_sum2() {
	int i;
	double sum2;

	/* Sum normalisation weights from all vertices */
	sum2 = 0.0;
	for (i=1; i<=s_count; i++) {
		if (upnodes[i] % 2)
			sum2 -= s_pow[i];
		else
			sum2 += s_pow[i];
		}
	return sum2;
	}


/* Adjust mp according to vertex algorithm if dimension permits.
   The algorithm explodes in dimensionality and can only run in small cases. */

static void adjust_vx(int alt, int snode, double lofrac) {
	int j,last,act_dim,sd_path[VX_MAXDIM+1],active[VX_MAXDIM+1];
	double target,mp_factor,sum2;
	int tnode,mp_dim;

	for (tnode=tdown[alt][snode],j=1; tnode; tnode=tnext[alt][tnode],j++) {
		if (j > VX_MAXDIM)
			/* mp_factor = 0 */
			return;
		mp_lobo[j] = (tdown[alt][tnode]?im_L_mhull_lobo[at2i(alt,tnode)]:L_mhull_lobo[at2r(alt,tnode)]);
		mp_upbo[j] = (tdown[alt][tnode]?im_L_mhull_upbo[at2i(alt,tnode)]:L_mhull_upbo[at2r(alt,tnode)]);
		}
	mp_dim = j-1;
	if (mp_dim <= VX_CUTOFFDIM)
		mp_factor = 1.0;
	else
		mp_factor = (double)(VX_MAXDIM+1-mp_dim)/(double)(VX_MAXDIM+1-VX_CUTOFFDIM);
#ifndef TDL_COMPAT
	mp_factor /= 2.0; // max half of mp is from this algorithm
#endif
	/* Collect prob mass and nbr of cons (=dim) to distribute */
	target = 1.0;
	last = 0;
	act_dim = mp_dim;
	for (j=1; j<=mp_dim; j++) {
		/* Skip collapsed dimensions */
		active[j] = (mp_upbo[j]-mp_lobo[j] > EPS100);
		if (active[j])
			last = j;
		else {
			target -= lofrac*mp_lobo[j]+(1.0-lofrac)*mp_upbo[j];
			act_dim--;
			}
		}
	if (act_dim < 2)
		/* Cannot change anything */
		return;
	/* Initialise invariants */
	s_count = 0;
	f1_T(0.0,target,0,last,sd_path,active);
	if (!s_count)
		/* No path */
		return;
	sum2 = ijar_warp_sum2();
#ifdef TDL_COMPAT
	if (sum2 < EPS)
#else // depends on depth
	if (sum2 < EPS/100.0)
#endif
		/* No normaliser */
		return;
	for (tnode=tdown[alt][snode],j=1; tnode; tnode=tnext[alt][tnode],j++)
		if (active[j]) {
			/* Non-collapsed dimension */
			if (tdown[alt][tnode])
				im_L_mass_point[at2i(alt,tnode)] = (1.0-mp_factor) * im_L_mass_point[at2i(alt,tnode)] +
					mp_factor * ijar_warp_T(j,act_dim,target,sum2);
			else
				L_mass_point[at2r(alt,tnode)] = (1.0-mp_factor) * L_mass_point[at2r(alt,tnode)] +
					mp_factor * ijar_warp_T(j,act_dim,target,sum2);
			}
	}
