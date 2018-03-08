;*******************************************************************************
;
;	This is the COBE-specific IDL environment startup procedure.
;
;	AUTHOR:	Song Yom (HUGHES STX) (derived from STARTUP1.PRO)
;	DATE:	03/92 (SER #9525)
;
;	MODIFICATION HISTORY:
;	DATE	SPR/SER		PROGRAMMER	DESCRIPTION
;	----	-------		----------	-----------
;	03/92	08588		Song Yom	Reconfiguration of PROJECTION
;		09616				and UHELP installation.
;       01/93   10513           K. Turpie       Added user defined startup
;       01/93   10479           K. Turpie       Added system variables for FITS
;       05/93   --              K. Turpie       Set mail and edit global
;                                               variables; added @wexpath.inc
;       06/93   11003           J. Newmark      Add system variable for
;                                               multiple operating systems.
;       06/93   --              K. Turpie       Removed mail and edit global
;                                               variables.
;       07/93   11127           J. Newmark      IDL for Windows compatability
;
;*******************************************************************************
;
;Create special system variables for IDL Astronomy User's Library routines:
;
;          Variable        Used in Routine
;          --------        --------------------------------------------------
 DefSysV, '!DEBUG',    0 ; ReadCol,SkyPro,ZParCheck
 DefSysV, '!PRIV',     0 ; ?
 DefSysV, '!TEXTUNIT', 0 ; FITS_Info,TextClose,TextOpen
 DefSysV, '!TEXTOUT',  1 ; FITS_Info,TextClose,TextOpen
;
; Create a local system variable to condense the number of different
; UNIX-like operating systems and hold the current !version.os.
;
 DefSysV,'!cgis_os',!version.os
 IF (!VERSION.OS EQ 'Win32') THEN defsysv,'!cgis_os','windows'
;
; Given below are three versions of IDL that for the purposes of
; IDL are both considered UNIX. If you are running on a different
; UNIX-compatible operating system , e.g. hpux(?) instead of
; sunos, make the appropriate change in the call.
;
 IF (!VERSION.OS EQ 'ultrix') THEN defsysv,'!cgis_os','unix'
 IF (!VERSION.OS EQ 'sunos') THEN defsysv,'!cgis_os','unix'
 IF (!VERSION.OS EQ 'OSF') THEN defsysv,'!cgis_os','unix'
;
 IF ((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN DEVICE, RETAIN=2
;
 !PROMPT = 'UIDL> '
;
 PRINT,'Linkimaging applications...'
;un-comment appropriate line for operating system
;modify path for UNIX systems

@$CSRC/link_unix.pro        ;unix
;@CGIS$IDL:LINK_VMS.PRO      ;vms
;
 PRINT,'Enter ? followed by a carriage return for IDL HELP'
 PRINT,'Enter UHELP followed by a carriage return for UIDL HELP'
;
;un-comment appropriate line for operating system
@$USER_START    ;unix
;@USER_START      ;vms

;the following line may be necessary for compilation of some very large
;routines. If you get an error message trying to compile a routine
;e.g. DCONVERT, try changing 33500 -> 42000, or whatever is necessary
;.size 33500 7500

;DISCLAIMER:
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
;Please send bug reports to CGIS@ZWICKY.GSFC.NASA.GOV.

