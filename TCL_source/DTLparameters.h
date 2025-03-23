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
 *                   UNEDA Decision Tree Layer (TCL)
 *                   -------------------------------
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
 *   File: DTLparameters.h
 *
 *   Purpose: public DTL parameters (also used by TCL)
 *
 */

 /*********************************************************
  *
  *  Configuration
  *
  *********************************************************/

/* Allowing TCL and DTL to work with large problem sizes:
 * 14400 alts on 8 GB PM Linux server, 53800 on 128 GB PM */

#define noBIGFOOT


 /*********************************************************
  *
  *  Parameters, fields and directives
  *
  *********************************************************/

#ifdef BIGFOOT

/* Max number of alternatives */
#define MAX_ALTS 3000
/* Max number of re-nodes (consequences) and re+im-nodes in total */
#define MAX_CONS MAX_ALTS+25
#define MAX_NODES MAX_ALTS+30
/* Max number of re-nodes (consequences) and re+im-nodes per alternative */
#define MAX_COPA 22
#define MAX_NOPA 24

#else

/* Max number of alternatives */
#define MAX_ALTS 100
/* Max number of re-nodes (consequences) and re+im-nodes in total */
#define MAX_CONS 920
#define MAX_NODES MAX_CONS+MAX_CONS-2
/* Max number of re-nodes (consequences) and re+im-nodes per alternative */
#define MAX_COPA 512
#define MAX_NOPA MAX_COPA+MAX_COPA-2

#endif

/* Max size of frame name */
#define FNSIZE 128

/* Output fields */
#define E_MIN 0
#define E_MID 1
#define E_MAX 2
#define MAX_ERESULT E_MAX

#ifdef _MSC_VER
/* Assignment within conditional expression is a coding style in DTL/TCL */
#pragma warning(disable:4706) // (this corresponds to /wd4706 in newer VC++)
#endif
