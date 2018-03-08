/*
   HEADER: export_wrap.h

   This file is provided for the COBE C software. Please read
   it carefully and make the appropriate changed for the version
   of IDL you are using and your pathways.
*/
#include <stdio.h>	/* export.h needs this */

	/* define problem symbols before including export.h */
#if defined(__GNUC__) && defined(VMS)
#define noshare
#define globalref extern
#endif
	/* IDL's export.h file, at /include=idl_dir:[source] */
#define IDL_UCHAR UCHAR
#include "export.h"	

/* the following line take care of obsolete typedefs */
/* uncomment this section if running versions of IDL */
/* OLDER than 4.0 (e.g. 3.5,3.6)                      */
/* #include "cobe_obs.h"                 */

/* the following 2 lines are only for IDL v4.0  */
/* uncomment if necessary (fixed in 4.0.1)     */
/* #define ULONG IDL_ULONG */
/* #define IDL_DELTMP IDL_Deltmp */
/*  end idl 4.0 section                                */

/* The export.h file was changed in IDL version 3.5 to add   */
/* definitions of LONG and ULONG. If you are running older   */
/* versions, e.g. IDL v3.0 uncomment (remove the leading and */
/* trailing symbols) the next two lines.                     */
/* typedef long int LONG;            */
/* typedef unsigned long int ULONG;  */


/* ;DISCLAIMER:
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
;Please send bug reports to CGIS@COBECL.DNET.NASA.GOV. */

