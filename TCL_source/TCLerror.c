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
 *   File: TCLerror.c
 *
 *   Purpose: provide error texts
 *
 *
 *   Functions exported outside TCL
 *   ------------------------------
 *   TCL_get_errtxt
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


 /*********************************************************
  *
  *  Lower case error texts for low level kernel errors
  *
  *********************************************************/

static char *tcl_errtxt[] = {
	"tcl ok",                   //  0
	"inconsistent",             //  1
	"input error",              //  2
	"tree error",               //  3
	"illegal node",             //  4
	"too many alternatives",    //  5
	"too many consequences",    //  6
	"too many statements",      //  7
	"too narrow statement",     //  8
	"too few alternatives",     //  9
	"frame corrupted",          // 10
	"frame attached",           // 11
	"frame detached",           // 12
	"no such file or directory",// 13
	"infinite solution",        // 14
	"out of memory",            // 15
	"memory leak",              // 16
	"- rc out of range -",      // 17+
	};


char *TCL_get_errtxt(rcode rc) {

	if ((rc<0) || (rc>MAX_RCODE))
		return tcl_errtxt[MAX_RCODE+1];
	else
		return tcl_errtxt[rc];
	}
