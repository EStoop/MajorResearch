;+ NAME: 
; BTB_MOON_ANGLES calculates and returns the Moon_Re_Bs angles
;
;DESCRIPTION:
; The BTB_MOON_ANGLES routine is an UIDL-based tool function which calculates
; and returns an array of Moon_Re_BS angles with the following specifications:
;
;   1- The first dimension is the angle between moon's position vector
;      and the DIRBE boresight.
;      Units = Degrees, Range = (0, 180) degrees
;
;   2- The second dimension is the Azimuth of moon about the DIRBE boresight.
;      Units = Degrees, Range = (-180, 180) degrees
;      Azimuth = 0 is in the direction normal to the leading edge.
;      Azimuth = 90 points to the spin axis.
;
; This function requires only two input ARRAY parameters of times and pixels.
; The input array of times should be given in VAX ADT format only.
;
; WARNING:
;   This routine is written in FORTRAN and uses COBETRIEVE.
;
; EXAMPLE:
;
; UIDL> FILE = 'CSDR$DIRBE_PASS2:BCO_ACOD'
; UIDL> T1 = '901150000'
; UIDL> T2 = '901160000'
; UIDL> S = READ_TOD(file,'time,pixel_no',t1,t2,TIMES,PIXELS,maxrec=500)
; UIDL> ANGLES = BTB_MOON_ANGLES(TIMES,PIXELS)
; UIDL> help, angles
;
; ANGLES             FLOAT        = Array(2,500)
;
;.TITLE
;Routine BTB_MOON_ANGLES
;-

function btb_moon_angles
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


