pro test_trans_dur
;what are some possible transit durations for the hd7924 system, 
;could it be as low as <0.005 phase or ~20 min?

period = 5.3978 ;days

a = .057D;au
a = a *149597871. ; in km
print, 'a', a

Rs = 0.78D; solar radii
Rs = Rs * 695500.
print, 'Rs', Rs

i_deg= 50.
i_rad = i_deg*!dtor


Rp = 1.4 ; earth radii range is 1.5 - 6 earth radii
Rp = Rp*6378.1
print, 'Rp', Rp

insides_top = (1+(Rp/Rs))^2 - ((a/Rs)*cos(i_rad))^2
insides_bottom = 1 - (cos(i_rad))^2
insides = insides_top/ insides_bottom
Tdur = (period/!Pi) * asin( (Rs/a)*sqrt(insides))

print, 'transit duration', Tdur,' days or ', Tdur *24., ' hours' 


end
