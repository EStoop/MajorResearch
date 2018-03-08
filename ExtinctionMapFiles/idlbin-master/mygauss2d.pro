FUNCTION mygauss2d, x,y, P
xprime = (x-P(4))*cos(P(6)) - (y - P(5)) *sin(P(6))
yprime = (x-P(4))*sin(P(6)) + (y-P(5))*cos(P(6))
u = (xprime / P(2) )^2 + (yprime / P(3) )^2

;xprime = (x-P(3))*cos(P(5)) - (y - P(4)) *sin(P(5))
;yprime = (x-P(3))*sin(P(5)) + (y-P(5))*cos(P(5))
;u = (xprime / P(2) )^2 + (yprime / P(2) )^2

return, P(0) + P(1)*exp(-u/2.) 


;P[0] = P0 = constant offset
;P[1] = P1 = amplitude
;P[2] = a = s width of Gaussian in the X direction
;;;;P[3] = b = s width of Gaussian in the Y direction
;P[4] = h = X centroid
;P[5] = k = Y centroid
;P[6] = T = Theta, the rotation of the ellipse from the X axis in radians.
 
;http://astro.berkeley.edu/~johnjohn/idl_course/hw5.html
end
