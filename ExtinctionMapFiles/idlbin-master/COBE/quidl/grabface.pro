        FUNCTION grabface, dataset, facenum
;
;+NAME/ONE-LINE DESCRIPTION:
;     GRABFACE grabs a cube face from an unfolded T or sixpack
;
; DESCRIPTION: 
;    This function extracts one cube face from a dataset which is either
;    in skycube "T" format with a 4:3 aspect ratio or sixpack format
;    with a 3:2 aspect ratio.  The format of the dataset need not be
;    stated in the call to the function.  The input dataset may also be
;    a spectral cube (e.g. FIRAS data).  Any numerical data type can be 
;    accommodated; the output will always be a face with the same data 
;    type as the input.  
;
;
; CALLING SEQUENCE:
;    outface = GRABFACE (dataset, facenum)
;
; ARGUMENTS: (I=input, O=output, []=optional)
;
;    OUTFACE     O   any type   Output one face with same data type as 
;                                 input dataset
;    DATASET     I   any type   Input 2-D or 3-D dataset in either
;                                 sixpack or skycube "T" format
;    FACENUM     I   int        Face number to be extracted
;
; WARNINGS:
;
; EXAMPLE:
;    If you have a skycube called FIRAS_SIGNAL then the call
;
;    UIDL> FACE4 = GRABFACE (FIRAS_SIGNAL, 4)
;
;    will create a square array containing the data values in face 4,
;    the galactic center.
;
;#
;  COMMON BLOCKS:  None
;
;  PROCEDURE AND PROGRAMMING NOTES:
;      The program is written to manipulate 3-D datasets; therefore,
;      2-D datasets are first REFORMed into 3-D arrays with a degenerate 
;      3rd dimension (1) and then the work is done.  The input and 
;      output arrays are REFORMED back to 2-D if that was their original
;      form.  Skycubes and sixpacks are distinguished by their aspect
;      ratios and then the appropriate loop of a CASE calculates the 
;      position coordinates for the faces.
;
;  LIBRARY CALLS:  None
;
;  MODIFICATION HISTORY:
;     Written by Celine Groden, USRA    22 Dec 1992   SPR 10388
;
;.TITLE
;Routine GRABFACE
;-
;
;  Establish error condition to return to calling routine in case of
;  trouble.
;
	ON_ERROR, 2
;  
;  Check if input was in fact supplied and that is it either in 
;  sixpack or skycube format.
;
	dims = SIZE (dataset)
	IF dims(0) EQ 0 THEN MESSAGE, "No input data set provided!"
  	IF (2*dims(1) NE 3*dims(2) AND 3*dims(1) NE 4*dims(2)) THEN $
	     MESSAGE, "Input data must be in skycube or sixpack format"
;
;  Check if a valid face number (0-5) was provided.  Also flag if NO 
;  face number was supplied.
;
	test=size(facenum)
	IF test(1) EQ 0 THEN MESSAGE, "No face number was supplied"
	IF facenum GT 5 OR facenum LT 0 $
	                   THEN MESSAGE, "Face Number must be between 0-5"
;
;  The variable "dim_switch" will identify if the original dataset 
;  is temporarily converted to a 3-D object and will need to be
;  REFORMed to a 2-D array at the end of the routine.  The variable
;  "six_switch" is used to identify whether or not the input dataset
;  is in skycube format (six_switch=0) or sixpack (six_switch=1)
;  format.
;
	six_switch = 0
	dim_switch = 0
;
;  If the input dataset is in sixpack format, turn six_switch on.
;
  	IF (2*dims(1) EQ 3*dims(2)) THEN six_switch = 1
;
;  Define the length of the side of one face.
;
	CASE six_switch OF
	   0: xlen = dims(1) / 4.
	   1: xlen = dims(1) / 3.
	ENDCASE
;
;  Temporarily convert 2-D datasets to 3-D sets.  If data is neither
;  2-D nor 3-D, flag it.
;
	IF (dims(0) NE 2 AND dims(0) NE 3) THEN $
	       MESSAGE, "Input data must be 2-D or 3-D"
;
	IF (dims(0) EQ 2) THEN BEGIN
	       dataset = REFORM (dataset, dims(1), dims(2), 1)
	       dim_switch = 1
	         ;Toggle dim_switch "on" to identify that the dataset 
	         ;needs to be converted back to 2-D at the end of the 
                 ;routine.
	ENDIF
;
;  Set up the coordinates of each of the six face positions within the
;  unfolded right T.  "F0XLO" means "Face 0, low X coord", i.e. the x 
;  coord of the left edge.  Similarly, "F4YLO" is the Y coordinate of the 
;  bottom of face 4, etc.  Each face is XLEN on a side.
;
	CASE six_switch OF
	  1: CASE facenum OF
	      0: BEGIN
		 facexlo = 2 * xlen
		 faceylo = xlen  
		 END
	      1: BEGIN
		 facexlo = 2* xlen 
		 faceylo = 0.
		 END
	      2: BEGIN
		 facexlo = xlen
		 faceylo = 0.
		 END
	      3: BEGIN
		 facexlo = 0.
		 faceylo = 0.
		 END
	      4: BEGIN
		 facexlo = 0.
		 faceylo = xlen 
		 END
	      5: BEGIN
		 facexlo = xlen
		 faceylo = xlen
		 END
	     ENDCASE
	  0: CASE facenum OF
	      0: BEGIN
		 facexlo = 3 * xlen
		 faceylo = 2 * xlen  
		 END
	      1: BEGIN
		 facexlo = 3 * xlen
		 faceylo = xlen  
		 END
	      2: BEGIN	
		 facexlo = 2 * xlen
		 faceylo = xlen
		 END
	      3: BEGIN
		 facexlo = xlen
		 faceylo = xlen
		 END
	      4: BEGIN
		 facexlo = 0.
		 faceylo = xlen
		 END
	      5: BEGIN
		 facexlo = 3 * xlen
		 faceylo = 0.
		 END
	     ENDCASE
	ENDCASE
;
;  Create face array
;
outface = dataset(facexlo:facexlo+xlen-1., faceylo:faceylo+xlen-1., *)
;
;  Switch datasets which were converted to 3-D objects back to 2-D
;  objects.
;
	IF (dim_switch EQ 1) THEN $
	    dataset = REFORM (dataset, dims(1), dims(2))
;
	RETURN, outface
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


