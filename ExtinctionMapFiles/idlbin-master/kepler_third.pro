pro kepler_third,Pdays, M1 
;P^2 = (4*!dpi^2*(a^3))/(G*(M1+M2))

  ;;convert: Period needs to be in seconds, assuming entered in days
  Psec = Pdays*24.*60.*60.D
  ;;convert: mass needs to be in Kg, assuming entered in solar masses
  M1 = M1* 1.98855D30

  ;;calculate a given P

  G = 6.673D-11                 ; (m^3Kg^-1s^-2)
  M2 = 0.

  top =(G)*(M1+M2)* (Psec^2)
  bottom = (4*!dpi^2)        
  print, top, bottom
  a = top/bottom
  ;a = a / 1.496D11

  print, a, 'not working'
  

end
