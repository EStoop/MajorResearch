PRO gfunct, X, A, F, pder
	bx = EXP(A[0] * X)
	
;	F = A[0] * bx + A[2]
;	F = A[0] + A[1]*X^2 
;	F = .153*(exp(-0.5*((x - A[1])/A[0])^2))
;	F = A[2] * exp(-0.5*(X- A[0]/A[1])^2)

	F = A[0]*(1 + (X / A[1])^2)^(-A[2])
; I = Ic (1 + (r / alpha)^2)^(-beta)   
	IF N_PARAMS() GE 4 THEN $
    		pder = [[bx], [A[0] * X * bx], [replicate(1.0, N_ELEMENTS(X))]]
END


