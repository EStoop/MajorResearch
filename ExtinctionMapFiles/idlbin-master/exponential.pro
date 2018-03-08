function exponential, x, P
  
  ;return, P(0) - P(1)*exp(-P(2)*x)
 return, P(0)*exp(-(x/(P(1))))
end

