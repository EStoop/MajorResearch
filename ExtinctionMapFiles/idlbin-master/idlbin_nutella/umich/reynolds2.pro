FUNCTION reynolds2, x, P

return, P(0) / ((1 + (x/P(1)))^2) + P(2) / ((1 + (x/P(3)))^2)

END 
