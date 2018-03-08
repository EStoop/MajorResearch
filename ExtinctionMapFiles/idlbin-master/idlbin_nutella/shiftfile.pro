pro shiftfile

;HST ACS data reduction astrometry  
;this code adds togethter shift files which are output from multiple runs of tweakshifts, and makes them useful for multidrizzle.


readcol, '/Users/jkrick/hst/raw/shifts01.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts01.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts01.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj01refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts02.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts02.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts02.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj01refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts03.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts03.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts03.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj03refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts04.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts04.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts04.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj04refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts05.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts05.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts05.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj05refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts06.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts06.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts06.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj06refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts07.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts07.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts07.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj07refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts08.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts08.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts08.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj08refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts09.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts09.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts09.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj09refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts10.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts10.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts10.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj10refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts11.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts11.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts11.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj11refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts12.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts12.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts12.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj12refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts13.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts13.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts13.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj13refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts14.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts14.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts14.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj14refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts15.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts15.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts15.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj15refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts16.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts16.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts16.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9a16refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts17.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts17.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts17.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj17refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts18.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts18.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts18.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj18refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts19.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts19.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts19.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj19refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts20.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts20.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts20.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj20refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts21.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts21.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts21.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj21refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts22.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts22.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts22.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj22refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts23.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts23.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts23.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj23refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts24.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts24.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts24.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj24refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

readcol, '/Users/jkrick/hst/raw/shifts25.txt', name, x, y, rot, mag, format="A"
readcol, '/Users/jkrick/hst/raw/shifts25.2.txt', name2, x2, y2, rot2, mag2, format="A"
openw, outlun, '/Users/jkrick/hst/raw/shifts25.combine.txt', /get_lun
printf, outlun, "# frame: output"
printf, outlun,"# reference: j9aj25refwcs.2.fits[wcs]"
printf, outlun,"# form: delta "
printf, outlun,"# units: pixels" 
for i = 0, n_elements(x) - 1 do printf, outlun, name(i), x(i)+x2(i), y(i)+y2(i), rot(i) + rot2(i), mag(i)
close, outlun
free_lun, outlun

end
