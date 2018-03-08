;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     CAPTURE copies an image from a workstation window to disk.
;
;DESCRIPTION:
;     Capture COPIES everything on a specified window, first into an
;     IDL array of the total window size, using 'tvrd'. This is a byte 
;     array scaled from 0 to 255. This is then stored on disk in an 
;     unformatted file. 
;
;CALLING SEQUENCE:
;     CAPTURE,outfile,[window]
;
;ARGUMENTS (I = input, O = output, [] = optional):
;
;     OUTFILE     O     str      Name of unformatted file to store 
;                                the output image. 
;     WINDOW	 [I]   int  	 The number identifying the IDL tv 
;                                window to capture. The IDL window '0' 
;                                is captured by default, unless 
;                                otherwise specified.
;                                  
;WARNINGS:
;     1.  The output file is created in the current or default 
;         directory. The user should make sure there is enough space.
;     2.  This procedure can only be called from a Workstation 
;         (X-Windows) session.
;
;EXAMPLE:
;     To capture IDL default display window (0), and store it in a file 
;     named 'output.dat', type:
;
;     CAPTURE, 'output.dat'
;
;     To capture IDL TV display window No. 1, and store it in a file 
;     named 'output.dat', type:
;
;     CAPTURE,1,'output.dat'
;#
; 
;COMMON BLOCKS:
;     None  
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES) :
;     The procedure uses 'tvrd' to read the display on the TV window
;     into a byte array scaled from 0 to 255. The dimensions of this
;     output array and the array itself are stored in an unformatted 
;     form.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     TVRD
; 
; MODIFICATION HISTORY:
;   Written by Mila Mitra, General Sciences Corporation, March 1991
;   SER 9616
;.TITLE
;Routine CAPTURE
;-
pro capture,outfile,window
;
;SETS WINDOW 
if (n_elements(window) eq 0) then window=0
wset,window
;
;USES TVRD TO READ DISPLAY INTO AN ARRAY
print,'Now starting tvrd'
tvwin=tvrd(0,0,!d.x_vsize,!d.y_vsize)
s=size(tvwin)
print,'Dimensions of the TVRD array are:'
print,s(1),s(2)
;
;WRITES TO OUTPUT UNFORMATTED FILE
openw,1,outfile+"/unf"
writeu,1,s(1),s(2)
writeu,1,tvwin
close,1
return
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


