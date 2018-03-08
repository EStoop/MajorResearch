function do_meps,data,meps    ;function to apply MEPS to one Mf of 
;SDM data
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    DO_MEPS  applies MUX sequence to one major frame of DIRBE data.
;
;DESCRIPTION:
;    Given an input 16x256 DIRBE science data array arranged in
;    process order, plus the MUX sequence appropriate for that 
;    data array, DO_MEPS produces an output science data array 
;    which has been sorted into detector order (1a,1b,1c,2a...10)
;
;CALLING SEQUENCE:
;    OUT_DATA = DO_MEPS( IN_DATA, MEPS)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    IN_DATA  I  int arr  16-band x 256-sample data array.
;                          (Typically, the user has obtained
;                           this from field DADRBSCI2 of either
;                           NBS_CTOD, BMD_MA, BLI_TOIRSCLD)
;    MEPS     I  byt arr  MEPS sequence (user supplied from
;                          DAMEPS field of either NBS_CTOD,
;                          BMD_MA, BLI_TOIRSCLD).
;    OUT_DATA O  int arr  16x256 output data array, resorted into
;                          in detector order.
;
;WARNINGS:
;    1.  Only works properly if detectors are not sampled redundantly
;        within the MUX.
;
;EXAMPLE:
;    The user has obtained 1 Major Frame's worth of SDM data and has
;    stored in in a 16x256 integer array called SCIDAT.  He/she
;    has also obtained the appropriate MEPS sequence for that major
;    frame (e.g., from BMD_MA.DAMEPS) and stored it in IDL variable 
;    MEPS_SEQ.  The following call to DO_MEPS
;    will re-order the 256 observations for each of the 16 DIRBE 
;    channels into detector order (1a,1b,1c....10), and place them
;    in output array outsci:
;           outsci = do_meps(scidat,meps_seq)
;       
;#
;COMMON BLOCKS:
;    None.
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Matches detector with MUX process and reorders data
;     array.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS,ETC.:
;    None.
;
;MODIFICATION HISTORY:
;    Written: Spiesman (GSC) 1989
;
; SPR 9616
;.TITLE
;Routine DO_MEPS
;-
;

sdm=intarr(16,256)                 ;space for new data

m=bytarr(16)
mm=bytarr(16)

for i=0,15 do mm(i)=meps(2*i+1)

for i=0,15 do begin
  it=where(mm eq i)
  m(i)=it
end

; re- order from 1a 2a 3a... to 1a 1b 1c...
s1=m(1)
s2=m(2)
s3=m(3)
s5=m(5)
s6=m(6)
s7=m(7)

m(1)=s3
m(2)=s6
m(3)=s1
m(5)=s7
m(6)=s2
m(7)=s5

for i=0,15 do for j=0,255 do sdm(i,j)=data(m(i),j)


return,sdm

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


