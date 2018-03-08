;******************************************************************************
;
;MY_QSIMP integrates a function using Simpson's rule to specified accuracy.
;
;       func - scalar string giving name of function of one variable to 
;               be integrated
;
;       A,B  - numeric scalars giving the lower and upper bound of the 
;               integration
;
;       S - Scalar giving the approximation to the integral of the specified
;               function between A and B.
;
;       EPS - scalar specifying the fractional accuracy before ending the 
;               iteration.  Default = 1E-6
;
;       MAX_ITER - Integer specifying the total number iterations at which 
;               QSIMP will terminate even if the specified accuracy has not yet
;               been met.   The maximum number of function evaluations will be
;               2^(MAX_ITER).    Default value is MAX_ITER = 20
;This function was adapted from the built-in procedure QSIMP
;******************************************************************************

pro my_qsimp, func, A, B, mstar, alpha, S, EPS=eps, MAX_ITER = max_iter

 On_error,2                                  ;Return to caller

 if N_params() LT 6 then begin
    print,'Syntax - QSIMP, func, A, B, mstar, alpha, S, [ MAX_ITER = , EPS = ]
    print,' func - scalar string giving function name
    print,' A,B - endpoints of integration, S - output sum'
    return
 endif

 my_zparcheck, 'QSIMP', func, 1, 7, 0, 'Function name'       ;Valid inputs?
 my_zparcheck, 'QSIMP', A, 2, [1,2,3,4,5], 0, 'Lower limit of Integral'
 my_zparcheck, 'QSIMP', B, 3, [1,2,3,4,5], 0, 'Upper limit of Integral'

 if not keyword_set(EPS) then eps = 1.e-6              ;Set defaults
 if not keyword_set(MAX_ITER) then max_iter = 20

 ost = (oS = -1.e30)
 for i = 0,max_iter - 1 do begin
    my_trapzd,func,A,B,mstar,alpha,st,it
    S = (4.*st - ost)/3.
    if ( abs(S-oS) LT eps*abs(oS) ) then return
    os = s
    ost = st
 endfor

 message,/CON, $
        'WARNING - Sum did not converge after '+ strtrim(max_iter,2) + ' steps'
 print,S,oS,eps

 return
 end
