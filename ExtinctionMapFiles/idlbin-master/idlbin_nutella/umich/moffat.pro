function moffat, X,P
	return, P(0)*(1 + ((X) / P(1))^2)^(-2.592)
	; I = Ic (1 + (r / alpha)^2)^(-beta)   
end
