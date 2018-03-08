;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     GOIRAF invokes the iraf startup command approriate to the OS
;
;DESCRIPTION:  
;     GOIRAF checks the IDL system variable !version.os and
;     branches to the appropriate section of code.  It then invokes
;     a command procedure or a shell script to start up iraf
;
;CALLING SEQUENCE:  
;     GOIRAF
;
;ARGUMENTS (I = input, O = output, [] = optional):
;     NONE
;
;WARNINGS:
;
;EXAMPLE:
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES): 
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     None
;  
;MODIFICATION HISTORY:
;     Written by Dave Bazell,  General Sciences Corp. Mar 1993 spr 10687
;     SPR 11003  Jun 02 93   Change !version.os to !cgis_os.  J Newmark
;
;.TITLE
; Routine GOIRAF
;-
;
; Check on input parameters
;
pro goiraf

on_error, 2

case !cgis_os of
	'vms'	: spawn, '@CGIS$IDL:StartIRAF'
	'unix': spawn, 'source $csrc/qultrix/goiraf.csh'
        'win': spawn, 'source csrc\win\goiraf.csh'
endcase

end
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


