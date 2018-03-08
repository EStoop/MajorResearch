pro transit_unc, planetname
 ;; now = 2457558.5               ; June 19, 2016
 ;; now = 2458849.500000          ;Jan 1 2020
  now = 2459945.500000; jan 1 2023
  if planetname eq 'WASP-62b' then begin
     sigma_T0 = .00027          ; days
     sigma_P = 3E-6             ;days
     T0 = 2455855.39
     P =4.411953                ;days
  endif
  

  if planetname eq 'WASP-63b' then begin
     P = 4.37809
     sigma_P = 6E-6
     T0 = 2455921.6527
     sigma_T0 = 0.0005
  endif
  
  if planetname eq 'WASP-100b' then begin
     P = 2.849375
     sigma_P =  8e-06
     T0 = 2456272.3395
     sigma_T0 = 0.0009
  endif
  
  if planetname eq 'HAT-P-11' then begin
     P = 4.887804
     sigma_P = 4e-06
     T0 = 2454605.89132
     sigma_T0 = 0.00032
  endif
  
  if planetname eq 'tres-2b' then begin
     P=2.4706133738
     sigma_P =1.87e-08
     T0 = 2453957.6355
     sigma_T0 = 1.29e-05
  endif
  
  time_elapsed = now - T0
  period_elapsed = time_elapsed / P
  print, 'period_elapsed', period_elapsed
  
  sigma_Tn = sqrt(sigma_T0^2 + (period_elapsed*sigma_P)^2) ;in days
  sigma_Tn = sigma_Tn * 24.*60 ; in minutes
  print, 'sigma timing', sigma_Tn, 'min.'
end
