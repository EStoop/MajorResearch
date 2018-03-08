pro medcombine, dirname

;quick code to median combine a stack of images
;dirname =  '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/latent/r48724480/ch2/bcd/'
cd, dirname
command  = strcompress( 'ls SPITZER*_bcd.fits > bcdlist.txt')
spawn, command
readcol,'bcdlist.txt',fitsname, format = 'A', /silent

bigim_sub = fltarr(32, 32, 128, n_elements(fitsname))  
bigim_full = fltarr(256,256, n_elements(fitsname))  

for i = 0, n_elements(fitsname) - 1 do begin
;   print, i, fitsname(i)
   fits_read, fitsname(i), data, header
   naxis = sxpar(header, 'NAXIS')
   if naxis eq 3 then bigim_sub(0, 0, 0, i) = data else bigim_full(0, 0, i) = data
endfor
 
if naxis eq 3 then medcombine = median(bigim_sub, dimension = 4) else medcombine = median(bigim_full, dimension = 3)

medname = strcompress(dirname + '/medcombine.fits',/remove_all)
fits_write, medname, medcombine, header

end
