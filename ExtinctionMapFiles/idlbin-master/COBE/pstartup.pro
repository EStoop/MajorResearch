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
; Given below are two versions of IDL that for the purposes of
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
; PRINT,'Linkimaging applications...'
;un-comment appropriate line for operating system
;modify path for UNIX systems

;@$CSRC/link_unix.pro        ;unix
;@CGIS$IDL:LINK_VMS.PRO      ;vms
;
 PRINT,'Enter ? followed by a carriage return for IDL HELP'
 PRINT,'Enter UHELP followed by a carriage return for UIDL HELP'
;
;@$USER_START    ;unix
;@USER_START      ;vms

;the following line is necessary for compilation of some very large
;routines. If you get an error message trying to compile a routine
;e.g. DCONVERT, try changing 33500 -> 42000, or whatever is necessary
;.size 33500 7500

;modify path as appropriate
setenv,'cgis_data=c:\rsi\idl50\cgis\qlut\'
setenv,'cgis_fits=c:\rsi\idl50\cgis\data'
setenv,'cgis_ciss=c:\rsi\idl50\cgis\data'
setenv,'home=c:\rsi\idl50'
setenv,'uimage_help=c:\rsi\idl50\cgis\uimage_h.dat'
.run \rsi\idl50\cgis\quimage\xdispla0
.run \rsi\idl50\cgis\quimage\change_0
.run \rsi\idl50\cgis\quimage\change_1
.run \rsi\idl50\cgis\quimage\display0
.run \rsi\idl50\cgis\quimage\select_0
.run \rsi\idl50\cgis\quidl\draw_gri
.run \rsi\idl50\cgis\quimage\draw_gri
.run \rsi\idl50\cgis\quimage\xdisplay
.run \rsi\idl50\cgis\quimage\project_image
