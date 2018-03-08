/*;DISCLAIMER:
;
;This software was written at the Cosmology Data Analysis Center in
;support of the Cosmic Background Explorer (COBE) Project under NASA
;contract number NAS5-30750.
;
;This software may be used, copied, modified or redistributed so long
;as it is not sold and this disclaimer is distributed along with the
;software.  If you modify the software please indicate your
;modifications in a prominent place in the source code.  
;
;All routines are provided "as is" without any express or implied
;warranties whatsoever.  All routines are distributed without guarantee
;of support.  If errors are found in this code it is requested that you
;contact us by sending email to the address below to report the errors
;but we make no claims regarding timely fixes.  This software has been 
;used for analysis of COBE data but has not been validated and has not 
;been used to create validated data sets of any type.
;
;Please send bug reports to CGIS@COBECL.DNET.NASA.GOV.  */



/*
**++
**  FACILITY:  cpusec
**
**  MODULE DESCRIPTION:
**
**      CPUSec returns the number of cpu seconds consumed
**	by the current process as a floating point value.   
**
**  AUTHORS:
**
**      Pete Kryszak-Servin, General Sciences Corp.
**
**  CREATION DATE:  28-Aug-1992
**
**  DESIGN ISSUES:
**
**      This routine uses the C run-time library routine
**	clock() and the value CLK_TCK.  
**
**  INCLUDED FILES:
**   
**  	stdio.h	    standard c header file
**	time.h	    standard c header file
**	export.h    supplied by RSI with IDL, needed for definitions
**   
**  PORTABILITY ISSUES:
**
**  	VAX C supports the macro for CLK_TCK (clock cycles per second)
**	while GNU C does not.  There is a check for the definition and
**	if CLK_TCK is undefined, then it is defined to 100.  Other 
**	hardware platforms may require another value.
**   
**	This routine may work in plain IDL or COBE's UIDL, but
**      IDL linkimage is supported on only a few hardware platforms.
**
**   
**  SUBSYSTEM:
**   
**  	This routine is part of the CDAC Guest Investigator System
**	subsystem.
**   
**  USAGE:
**  
**	Example 1.    
**
**	UIDL> setlog, 'cpusec_exe', 'cgis$idl:cpusec.exe'
**	UIDL> linkimage, 'cpusec', 'cpusec_exe', 1
**	UIDL> print, cpusec()
**	494.34
**
**	This means that the current process has consumed 8 minutes
**	and 14.34 seconds of cpu time.  Control-T would have given
**	output like the following:
**	NODE::USER 11:42:35 COBE_IDL   CPU=00:08:14.34 PF=367 IO=951 MEM=512
**
**	Example 2. 
**  
**	UIDL> t0 = cpusec()
**	UIDL> ...routine or steps to be timed...
**	UIDL> t1 = cpusec()
**	UIDL> print, 'Routine consumed ', t1-t0, ' CPU seconds.'
**
**	This example shows how to time a routine or several steps
**	of an IDL session, and would normally be used as part of 
**	a procedure rather than run interactively as shown.
**
**
**  MODIFICATION HISTORY:
**       09/95                    Newmark         Upgrade to IDL 4.0
**
**      {@tbs@}...
**--
*/

#include <stdio.h>
#include <time.h>
#include "export_wrap.h"

#ifndef CLK_TCK
/* valid for VAX VMS run time library clock() */
#define CLK_TCK (100)
#endif

 
IDL_VPTR cpusec( argc, argv, argk ) 
int argc;
IDL_VPTR argv[];
char *argk;
{

/*  
**  this IDL_VARIABLE is declared as static so as to avoid
**  returning a pointer to a stack VARIABLE or using
**  IDL's temporary VARIABLEs which need to be cleaned up.
*/
static IDL_VARIABLE cpuseconds;


/* it's a IDL_LONG, and it has a value */

cpuseconds.type    = IDL_TYP_FLOAT;
cpuseconds.value.f = ((float) clock()) / (float)CLK_TCK;

return( &cpuseconds );

}
