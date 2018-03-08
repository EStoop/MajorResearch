pro latent_intercept
  redcolor = FSC_COLOR("Red", !D.Table_Size-2)
  bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)

  known_flux = [5.16811E7, 1.2228E8, 3.6333E8, 1.0292E9, 2.5377E9]
  known_flux  = (known_flux / known_flux(0) )* 1E7
  y_int = [570, 830, 1311, 1701, 2027]
 
  plot, known_flux, y_int, xtitle = 'known flux', ytitle = 'y intercept', psym = 2, charsize = 2, /xlog, xrange = [1E4, 1E9]

  
  ;fit with an exponential function.
  start = [2000.,1500.,0.1/1E7]
  noise = fltarr(n_elements(y_int))
  noise[*] = 1.0
  result= MPFITFUN('exponential',known_flux,y_int, noise, start) ; fit a function
  x = findgen(500) * 1E7
  ;print, x
  oplot, x,  result(0) - result(1)*exp(-result(2)*x), color = bluecolor
  ;oplot, known_flux,  result(0) - result(1)*exp(-result(2)*known_flux), color = bluecolor
  
  ;try with a different function
;P(0) + P(1)*alog10(x)
  start = [ -5E3,800.]
  result= MPFITFUN('logarithm',known_flux,y_int, noise, start) ; fit a function
  oplot, known_flux, result(0) + result(1)*alog10(known_flux), color = redcolor
end
