PRO FULLMAP, Filename
;+
;   NAME:
;      FULLMAP
;   PURPOSE:
;      Creates a full DSK7 skymap file from a partial one.
;   CATEGORY:
;      DMR Image Analysis.
;   CALLING SEQUENCE:
;      FULLMAP,'Partial_Skymap_Filename'
;   INPUTS:
;      Filename = Filename of partial skymap; must be in single-quotes.
;   OUTPUTS:
;      None; but it does create an output file named:  'Filename.Ext_FULL'
;   RESTRICTIONS:
;      Only works with DSK7_SKY.rdl record structures.
;   MODIFICATION HISTORY:
;      Created:  K. Galuk, COBE/DMR, STX, April 1991.
;-
@DMRSKY.inc
DSK_In = P(0)
Bracket = STRPOS(Filename,']')
I = STRPOS(Filename,'.',Bracket+1)
Outfile = STRMID(Filename, Bracket+1, I-1) + '_FULL'
Print,' Input  File: ',Filename
Print,' Output File: ',Outfile
OPENR,Lun1,Filename,/GET_LUN
OPENW,Lun2,Outfile,/GET_LUN
WHILE NOT EOF(Lun1) DO BEGIN
   READU,Lun1,DSK_In
   P(DSK_In.Pixel) = DSK_In
ENDWHILE
CLOSE,Lun1
WRITEU,Lun2,P
CLOSE,Lun2
FREE_LUN,Lun1
FREE_LUN,Lun2
Print,' Full-skymap file ready'
RETURN
END
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


