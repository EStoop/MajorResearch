FUNCTION sersic2, x, P
return, ((P(0)) * (exp(-7.67*(((x/P(1))^(1.0/P(2))) - 1.0))))  + ((P(3)) * (exp(-7.67*(((x/P(4))^(1.0/P(5))) - 1.0))))  

END 
