FUNCTION devauc2, x, P
return, (P(0)) * (exp(-7.67*(((x/P(1))^(1.0/4.0)) - 1.0))) + (P(2)) * (exp(-7.67*(((x/(P(3)))^(1.0/4.0)) - 1.0)))

END 