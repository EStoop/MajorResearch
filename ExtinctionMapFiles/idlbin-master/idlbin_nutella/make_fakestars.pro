pro make_fakestars

;make 100 images with 2 stars in each image.
;each image should have a different separation value and a different set of magnitudes
;separatation values range from 0 to 20 around the center of a 512x512 image
;magnitudes range from 19 to 24  (sextractor thinks it is 2 mags brighter, aka 17 to 22)


sep = randomu(5, 100) * 20.

mag1 = randomu(50, 200)*5. + 19.
mag2 = randomu(500, 200)*5. + 19.


;need to figure out how to plot with mag changing
;will want to measure a delta mag somehow.
;for now
mag1[*] = 21.3
mag2[*] = 20.5

openw, cllun, '/Users/jkrick/spitzer/irac/fakegal/fakestarscript.cl', /get_lun
openw, prolun, '/Users/jkrick/spitzer/irac/fakegal/fakestarlist.txt', /get_lun
for i = 0, n_elements(sep) - 1 do begin
   openw, datlun, strcompress("/Users/jkrick/spitzer/irac/fakegal/stars." + string(i)+ ".dat",/remove_all) , /get_lun
   printf, datlun, '200.5  200.5  ', mag1(i)
   printf, datlun, '200.5 ', 200.5 + sep(i), mag2(i)
   close, datlun
   free_lun, datlun

   printf, cllun, strcompress("mkobjects " + strcompress( " 'stars." + string(i) + ".fits' ", /remove_all) + " obj = " + strcompress("'stars."+string(i)+".dat'",/remove_all))

   printf, prolun, strcompress("/Users/jkrick/spitzer/irac/fakegal/stars."+string(i)+".fits", /remove_all)
endfor

close, cllun
free_lun, cllun
close, prolun
free_lun, prolun


;also need to print out an input file of all the filenames for iracphot.fakegal.pro




end
