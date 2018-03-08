function polynomial_5, x,y, P
  ;fifth order 
  return, P(0) + P(1)*x^2 + P(2)*x*y + P(3)*y^2 $
          + P(4)*x^3 + P(5)*(x^2)*y + P(6) *x*(y^2) + P(7)*y^3 $
          + P(8)*x^4 + P(9)*(x^3)*y + P(10)*(x^2)*(y^2) + P(11)*x*(y^3) + P(12)*y^4 $
          + P(13)*x^5 + P(14)*(x^4)*y + P(15)*(x^3)*(y^2) + P(16)*(x^2)*(y^3) + P(17)*x*(y^4) + P(18)*y^5
  
end

function polynomial_4, x,y, P
  ;fouth order 
  return, P(0) + P(1)*x^2 + P(2)*x*y + P(3)*y^2 $
          + P(4)*x^3 + P(5)*(x^2)*y + P(6) *x*(y^2) + P(7)*y^3 $
          + P(8)*x^4 + P(9)*(x^3)*y + P(10)*(x^2)*(y^2) + P(11)*x*(y^3) + P(12)*y^4
end

function polynomial_3, x,y, P
  ;third order 
  return, P(0) + P(1)*x^2 + P(2)*x*y + P(3)*y^2 $
          + P(4)*x^3 + P(5)*(x^2)*y + P(6) *x*(y^2) + P(7)*y^3 
          
end

function poly_3, x,P
  ;third order 
  return, P(0) + P(1)*x^2 + P(2)*x^3          
end
function poly_4, x,P
  ;third order 
  return, P(0) + P(1)*x^2 + P(2)*x^3  + P(3)*x^4    
end
function polynomial_2, x, P
  ;third order 
  return, P(0) + P(1)*x^2 ;+ P(2)*y^2 ;+ P(3)*x*y 
          
end

function polynomial_hora, x, y, P
  ;what Joe Hora used for the array location dependence correction

  return, P(0) + P(1)*(x-128) + P(2)*(y-128) + P(3)*(x-128)*(y-128) +$
          P(4)*(x-128)^2 + P(5)*(y-128)^2
end
