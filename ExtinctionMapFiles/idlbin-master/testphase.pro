pro testphase
tdiff = findgen(100) / 10.
phi_orbit = tdiff; 2*!dPI * tdiff
alpha = ABS(phi_orbit-!dpi)
phase_curve =  (sin(alpha) + (!dpi-alpha)*cos(alpha))/!dpi  ;; Phase curve for a Lambert sphere from Seager book 

;p = plot(phase_curve)
alpha2 = alpha 
y = (sin(alpha) + (!dpi-alpha)*cos(alpha))/!dpi   ;dividing by pi brings the max down to 1 (normalizes)


p = plot(alpha , y)
p = plot(alpha, sin(alpha), color = 'red', overplot = p)
p = plot(alpha, sin(alpha - !dpi/2), color = 'blue', overplot = p)

end
