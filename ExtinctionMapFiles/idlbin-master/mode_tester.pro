pro mode_tester

  aordir = '/Volumes/iracdata/flight/IWIC/'
  cd, aordir
  filename = './IRAC026200/cal/0041083136/IRAC.1.0041083136.4516290.1.skydark.fits'
  fits_read, filename, data, header, exten_no = 0
  data = data[*,*,0]
  
;get rid of nans
  a = where(finite(data), n_finite)
  data = data[a]
  mode_calc, data, mode,/verbose
  print, 'mode', mode
  
end
