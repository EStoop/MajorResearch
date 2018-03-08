pro test_phase 
bmjd = [55109.0,55204.2,55205.2,55206.6,55207.6,55207.7,55207.8,55207.9,55208.,55208.5,55208.8,55209.1,55209.5,55210.2,55213.2,55300.]
period = 3.1
utmjd_center = 55207.8
intended_phase = 0.5

bmjd_dist = bmjd - utmjd_center ; how many UTC away from the transit center

phase =( bmjd_dist / period )- fix(bmjd_dist/period)
print, 'bmjd', bmjd
print, 'bmjd - utmjd', bmjd_dist   
print, 'bmjd_dist / period', bmjd_dist / period
print, ' before phase',  phase;, format = '(A,F0)'    


low = where(phase lt-0.5 and phase ge -1.0)
phase(low) = phase(low) + 1.

high = where(phase gt 0.5 and phase le 1.0)
phase(high) = phase(high) - 1.0

print, 'middle ', phase

if intended_phase gt 0.4 and intended_phase lt 0.6 then begin ;secondary eclipse
   print, 'secondary eclipse intended'
   phase = phase+0.5
endif 
  print, ' after phase',  phase;, format = '(A,F0)'

end
