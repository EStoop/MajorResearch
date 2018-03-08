function moffat2, X, P
	return,P(2) * (1 + (2^(1/P(0))-1) * (X/P(1))^2)^ (-P(0)) 
end