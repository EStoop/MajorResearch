PRO	NGAUSS_FUNCT,X,A,F,PDER
; NAME:
;	NGAUSS_FUNCT
;
; PURPOSE:
;	EVALUATE THE SUM OF A GAUSSIAN AND OPTIONALLY A 2ND ORDER POLYNOMIAL
;	AND OPTIONALLY RETURN THE VALUE OF IT'S PARTIAL DERIVATIVES.
;	NORMALLY, THIS FUNCTION IS USED BY CURVEFIT TO FIT THE
;	SUM OF A LINE AND A VARYING BACKGROUND TO ACTUAL DATA.
;
; CATEGORY:
;	E2 - CURVE AND SURFACE FITTING.
; CALLING SEQUENCE:
;	NGAUSS_FUNCT,X,A,F,PDER
; INPUTS:
;	X = VALUES OF INDEPENDENT VARIABLE.
;	A = PARAMETERS OF EQUATION DESCRIBED BELOW.
; OUTPUTS:
;	F = VALUE OF FUNCTION AT EACH X(I).
;
; OPTIONAL OUTPUT PARAMETERS:
;	PDER = (N_ELEMENTS(X),6) OR (N_ELEMENTS(X),5) or (N_ELEMENTS(X),3) 
;               ARRAY CONTAINING THE PARTIAL DERIVATIVES.  P(I,J) 
;               = DERIVATIVE AT ITH POINT W/RESPECT TO JTH PARAMETER.
; COMMON BLOCKS:
;	NONE.
; SIDE EFFECTS:
;	NONE.
; RESTRICTIONS:
;	NONE.
; PROCEDURE:
;       IF N_ELEMENTS(A) EQ 3 THEN:
;       	F = A(0)*EXP(-Z^2/2) 
;       IF N_ELEMENTS(A) EQ 5 THEN :
;               F = A(0)*EXP(-Z^2/2) + A(3) + A(4)*X 
;       IF N_ELEMENTS(A) EQ 6 THEN :
;               F = A(0)*EXP(-Z^2/2) + A(3) + A(4)*X + A(5)*X^2
;	Z = (X-A(1))/A(2)
; MODIFICATION HISTORY:
;       WRITTEN, 9-aug-1993 J Newmark
;
	ON_ERROR,2                        ;Return to caller if an error occurs
	if a(2) ne 0.0 then Z = (X-A(1))/A(2) $	;GET Z
	else z= 10.

	EZ = EXP(-Z^2/2.)*(ABS(Z) LE 7.) ;GAUSSIAN PART IGNORE SMALL TERMS
        F = A(0) * EZ
        IF (N_ELEMENTS(A) EQ 5) THEN $
         	F = F + A(3) + A(4)*X ;FUNCTIONS.
        IF (N_ELEMENTS(A) EQ 6) THEN $
         	F = F + A(3) + A(4)*X + A(5)*X^2 ;FUNCTIONS.
        IF (N_PARAMS() EQ 3) THEN RETURN  ;NEED PARTIAL?
;
    IF (N_ELEMENTS(A) EQ 6) THEN BEGIN
	PDER = FLTARR(N_ELEMENTS(X),6) ;YES, MAKE ARRAY.
	PDER(0,0) = EZ		;COMPUTE PARTIALS
	if a(2) ne 0. then PDER(0,1) = A(0) * EZ * Z/A(2)
	PDER(0,2) = PDER(*,1) * Z
	PDER(*,3) = 1.
	PDER(0,4) = X
	PDER(0,5) = X^2
    ENDIF 
    IF (N_ELEMENTS(A) EQ 5) THEN BEGIN
	PDER = FLTARR(N_ELEMENTS(X),5) ;YES, MAKE ARRAY.
        PDER(0,0) = EZ		;COMPUTE PARTIALS
	if a(2) ne 0. then PDER(0,1) = A(0) * EZ * Z/A(2)
	PDER(0,2) = PDER(*,1) * Z
        PDER(*,3) = 1.
	PDER(0,4) = X
    ENDIF 
    IF (N_ELEMENTS(A) EQ 3) THEN BEGIN
        PDER = FLTARR(N_ELEMENTS(X),3) ;YES, MAKE ARRAY.
	PDER(0,0) = EZ		;COMPUTE PARTIALS
	if a(2) ne 0. then PDER(0,1) = A(0) * EZ * Z/A(2)
	PDER(0,2) = PDER(*,1) * Z
    ENDIF
    RETURN
END



Function gfit, x, y, weight=weight, a, sigmaa, baserem=baserem
;+
; NAME:
;	GFIT
;
; PURPOSE:
; 	Fit the equation y=f(x) where:
;
; 		F(x) = A0*EXP(-z^2/2) 
; 			and
;		z=(x-A1)/A2
;
;        If a baseline linear fit is requested then :
;  
;               F(x) = A0*EXP(-z^2/2) + A3 + A4*x 
;
;        If a baseline quadratic fit is requested then :
;  
;               F(x) = A0*EXP(-z^2/2) + A3 + A4*x + A5*x^2
;
;	A0 = height of exp, A1 = center of exp, A2 = sigma (the width).
;	A3 = constant term, A4 = linear term, A5 = quadratic term.
; 	The parameters A0, A1, A2, A3 are estimated and then CURVEFIT is 
;	called.
;
; CATEGORY:
;	fitting
;
; CALLING SEQUENCE:
;	Result = GFIT(X, Y [,WEIGHT=WEIGHT] [,A] [,SIGMAA] [,baserem=basetype])
;
; INPUTS:
;	X:	The independent variable.  X must be a vector.
;	Y:	The dependent variable.  Y must have the same number of points
;		as X.
;
; OPTIONAL INPUT PARAMETERS:
;       WEIGHT: A row vector of weights, the same lengthh as x and y.
;               If no weighting then weight(i)=1.0
;       BASEREM: If basetype = quadratic then a quadratic fit is made to a 
;                baseline. If basetype = linear then a linear fit is made to
;                a baseline.
;
; OUTPUTS:
;	The fitted function is returned.
;
; OPTIONAL OUTPUT PARAMETERS:
;	A:	The coefficients of the fit.  A is either a 3, 5 or 6 element 
;               vector as described under PURPOSE.
;       SIGMAA:  A vector of standard deviations for the parameters in A.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	None.
;
; RESTRICTIONS:
;	The peak or minimum of the Gaussian must be the largest
;	or smallest point in the Y vector.
;
; PROCEDURE:
;	If the (MAX-AVG) of Y is larger than (AVG-MIN) then it is assumed
;	that the line is an emission line, otherwise it is assumed there
;	is an absorbtion line.  The estimated center is the MAX or MIN
;	element.  The height is (MAX-AVG) or (AVG-MIN) respectively.
;	The width is found by searching out from the extrema until
;	a point is found less than the 1/e value.
;
; MODIFICATION HISTORY:
;	Written 09-aug-1993 J. Newmark
;-
;
on_error,2                      ;Return to caller if an error occurs
n = n_elements(y)		;# of points.
c = poly_fit(x,y,1,yf)		;fit a straight line.
yd = y-yf			;difference.

ymax=max(yd) & xmax=x(!c) & imax=!c	;x,y and subscript of extrema
ymin=min(yd) & xmin=x(!c) & imin=!c
IF (keyword_set(baserem)) THEN BEGIN
 IF (STRUPCASE(baserem) EQ 'QUADRATIC') THEN a=fltarr(6) ELSE a=FLTARR(5) 
ENDIF ELSE a=fltarr(3)             ;coefficient vector
if abs(ymax) gt abs(ymin) then i0=imax else i0=imin ;emiss or absorp?
i0 = i0 > 1 < (n-2)		;never take edges
dy=yd(i0)			;diff between extreme and mean
del = dy/exp(1.)		;1/e value
i=0
while ((i0+i+1) lt n) and $	;guess at 1/2 width.
	((i0-i) gt 0) and $
	(abs(yd(i0+i)) gt abs(del)) and $
	(abs(yd(i0-i)) gt abs(del)) do i=i+1
IF (keyword_set(baserem)) THEN BEGIN
  IF (STRUPCASE(baserem) EQ 'QUADRATIC') THEN $
  a = [yd(i0), x(i0), abs(x(i0)-x(i0+i)), c(0), c(1), 0.] $
  ELSE a = [yd(i0), x(i0), abs(x(i0)-x(i0+i)), c(0), c(1)]
ENDIF ELSE a = [yd(i0), x(i0), abs(x(i0)-x(i0+i))]              ;estimates
!c=0				;reset cursor for plotting
IF (not keyword_set(weight)) THEN weight=replicate(1.,n)
return,curvefit(x,y,weight,a,sigmaa, $
		function_name = "NGAUSS_FUNCT") ;call curvefit
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


