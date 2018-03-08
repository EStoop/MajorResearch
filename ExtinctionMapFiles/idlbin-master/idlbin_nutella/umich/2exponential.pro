function exp2, X, P
	return, P(0)*exp(-X/(P(1))) + P(2)*exp(-X/(P(3)))
end
