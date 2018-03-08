pro shiftfinal

;HST ACS data reduction astrometry  
;this code adds togethter shift files which are output from multiple runs of tweakshifts, and makes them useful for multidrizzle.



readcol, '/Users/jkrick/hst/raw/shifts08.combine.txt', name, x, y, rot, mag, format="A"
print, x
openw, outlun, '/Users/jkrick/hst/raw/shifts08.final.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj08refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, format='(A,F10.5,A, F10.5,F10.5,F10.2)',name(i), x(i)+208.336966*0.05, "  ", y(i)-136.278158*0.05, rot(i) + 359.9383014, mag(i)
close, outlun
free_lun, outlun


end
