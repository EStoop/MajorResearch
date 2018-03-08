; Function cpusec
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    CPUSec returns the amount of cpu time consumed as a FLOAT.
;
;DESCRIPTION:
;    This routine calls the C run-time library function clock()
;    and returns the number of cpu seconds consumed by the current
;    process as a FLOAT value.
;
;CALLING SEQUENCE:
;    cpu = cpusec()
;
;ARGUMENTS:
;	none
;
;WARNINGS:
;	This is the cumulative total for the current process
;	and does not include subprocesses, but does include 
;	all previous processing including login.  The value
;	matches that given by Control-T, except in format.
;
;EXAMPLE: 
;
;	Example 1.    
;
;	UIDL> setlog, 'cpusec_exe', 'cgis$idl:cpusec.exe'
;	UIDL> linkimage, 'cpusec', 'cpusec_exe', 1
;	UIDL> print, cpusec()
;	494.34
;
;	This means that the current process has consumed 8 minutes
;	and 14.34 seconds of cpu time.  Control-T would have given
;	output like the following:
;	NODE::USER 11:42:35 COBE_IDL   CPU=00:08:14.34 PF=367 IO=951 MEM=512
;
;	Example 2. 
;  
;	UIDL> t0 = cpusec()
;	UIDL> ...routine or steps to be timed...
;	UIDL> t1 = cpusec()
;	UIDL> print, 'Routine consumed ', t1-t0, ' CPU seconds.'
;
;	This example shows how to time a routine or several steps
;	of an IDL session, and would normally be used as part of 
;	a procedure rather than run interactively as shown.
;
;
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;	It is a C program which calls the C run-time library routine clock().
;	The number of clock ticks is then converted to a time in seconds by
;	dividing by CLK_TCK.  The routine is prototyped in and the constant
;	defined in TIME.H which is part of ANSI standard C.
;
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    none 
;
;MODIFICATION HISTORY
;    Written by Pete Kryszak-Servin, General Sciences Corp.   28-August-1992 
;	SER 9926
;
;.TITLE
;Routine CPUSec
;-
;
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


