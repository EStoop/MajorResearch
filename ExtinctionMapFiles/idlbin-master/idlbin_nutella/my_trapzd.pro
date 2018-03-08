;******************************************************************************
;
;MY_TRAPZD Computes the nth stage of refinement of an extended trapezoidal 
;rule.
;
;       func - scalar string giving name of function to be integrated.   This
;               must be a function of one variable.
;       A,B -  scalars giving the limits of the integration
;
;       S -    scalar giving the total sum from the previous interations on 
;               input and the refined sum after the current iteration on output.
;       step - LONG scalar giving the number of points at which to compute the
;               function for the current iteration.   If step is not defined on
;               input, then S is intialized using the average of the endpoints
;               of limits of integration.
;******************************************************************************

PRO MY_TRAPZD,func,a,b,Mstar,alpha,s,step

ON_ERROR,2

IF N_ELEMENTS(step) EQ 0 THEN BEGIN          ;Initialize?

	junk=CHECK_MATH(1) ;If a math error occurs, it is

	s1=CALL_FUNCTION(func,A,Mstar,alpha) ;likely to occur at the endpoints

	IF CHECK_MATH() NE 0 THEN MESSAGE,'ERROR-Illegal lower bound of '+STRTRIM(A,2)+'to function'+STRUPCASE(func)

	s2=CALL_FUNCTION(func,B,Mstar,alpha)

	IF CHECK_MATH() NE 0 THEN MESSAGE,'ERROR-Illegal upper bound of'+STRTRIM(B,2)+'to function'+STRUPCASE(func)

	s=0.5d*(DOUBLE(B)-A)*(s1+s2) ;First approx is average of endpoints

	step=1l

ENDIF ELSE BEGIN

	tnm=FLOAT(step)               

	del=(B-A)/tnm                    ;Spacing of the points to add

	x=A+0.5*del+FINDGEN(step)*del  ;Grid of points @ compute function

	sum=CALL_FUNCTION(func,x,Mstar,alpha)

	S=0.5d*(S+(DOUBLE(B)-A)*TOTAL(sum,/DOUBLE)/tnm)     

	step=2*step

ENDELSE

RETURN

END 
