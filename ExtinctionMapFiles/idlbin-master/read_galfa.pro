pro read_galfa
data = mrdfits(  '/Users/jkrick/plumes/coma/galfa/GALFA_HI_RA+DEC_196.00+26.35_N.fits', 0, Header)
naxis1 = sxpar(header, 'NAXIS1')
naxis2 = sxpar(header, 'NAXIS2')

print, data[325,355,490]

;fits_read, '/Users/jkrick/plumes/coma/galfa/GALFA_HI_RA+DEC_196.00+26.35_N.fits', data, header
s = data[325,355,*] ;plume location
;s = data[303,163,*]
x = findgen(n_elements(s))

;make x in velocity space
;ok these are in 0.184 km/s in the z direction from -188 to 188
;0 - 2048 -> -188 to 188
x = (x *0.184) - 188.
;p = plot( x, s, psym = 1)
;p.xtitle = 'Velocity (Km/s)'
;p.ytitle = 'K (Tb)'


;try smoshing the velocity space together by averaging over some of the velocity bins
data2 = fltarr(naxis1, naxis2)
for i = 0, naxis1 -1 do begin
   for j = 0, naxis2 - 1 do begin
      smoosh = mean(data[i,j,800:1240])
      data2[i,j] = smoosh
   endfor
endfor

fits_write, '/Users/jkrick/plumes/coma/galfa/mean_800_1240.fits', data2, Header
;try contour
;first display image
plothist, data2, xrange = [-100, 500]
;im = image(data2[200:400,290:450], min_value = -50, max_value = 500)

;c = contour(sim,min_value = -50, max_value =200,n_levels = 10,/overplot);, c_label_show = 1)

;print, s
end
