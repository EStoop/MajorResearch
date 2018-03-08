pro adjust_flux, photlist, channel


;++++++++++++++++++++++++++++++++++++++
; Give it a file list of phot data and channel and 
;it will compute the delta flux from pixel phase maps, 
;and print them out.
;
INPUT:
;    file: string containing filename of fits data
;    ra, and dec of object
;
;PJL 09/09
;++++++++++++++++++++++++++++++++++++

spawn, 'wc '+photlist, wc
foo=strsplit(wc, /extract)
np=foo(0)


if (channel EQ 1) then begin  
   c='1'
endif  else begin 
   c='2'
endelse

restore, 'pixel_phase_img_ch'+c+'.sav'


readcol, photlist,  ra, dec, xc, yc, flux, fluxerr, sn , NUMLINE=np

xphase = (xc mod 1) - NINT (xc mod 1)  ;;; NINT is nearest integer function
yphase = (yc mod 1) - NINT (yc mod 1)

;; Use interpolation to interpolate over pixel phase map

interp_fractional_flux= interp2d(relative_flux, x_center, y_center, xphase, yphase, /regular)
corr_flux = flux / interp_fractional_flux 


;; Should place this into a file.

print, photlist, ' RA       DEC       xcent  ', '   ycent   ', '   flux   ', '   pixel phase corrected flux ',   ' PHOTERR'
for j = 0, np-1 do begin
     print,  ra(j), dec(j), xc(j), yc(j), flux(j), corr_flux(j) , fluxerr(j)
endfor


print, ' ////// '
result = moment(corr_flux)
print, 'Means: ', 'Variance: ', 'Skewness: ', 'Kurtosis: '
print, result(0), result(1), result(2), result(3)



print, 'you are finished'
end


