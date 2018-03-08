function poisson, X, P
	return,P(1)* exp(-P(0))*( P(0)^(X)) / (factorial(X))

end
