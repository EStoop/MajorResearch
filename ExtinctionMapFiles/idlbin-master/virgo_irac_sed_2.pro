pro virgo_irac_sed
!P.multi = [0,1,1]
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)

;plume A
a_length = 668. ;arcsec

ch1signal = .00067 ;Mjy/sr
ch1noise = .00036 ;Mjy/sr

ch2signal = .00065 ;Mjy/sr
ch2noise =.00052  ;Mjy/sr

b_sb = 29.5 ; mag/sq. arcsec
b_noisesb = 29.0 ; mag/sq. arcsec

v_sb = 28.6 ; mag/sq. arcsec
v_noisesb = 28.5 ; mag/sq. arcsec

;halo E1
e_length =310

;ch1signal = .00129
;ch1noise = .00046

;ch2signal = .00111
;ch2noise = .00059

;v_sb = 28.2
;v_noisesb = 28.5

;b_sb = 29.05
;b_noisesb = 29.0 ; mag/sq. arcsec


length = a_length
;convert irac SB to vega mags over the area in question
;ch1-----------------------------------
ch1signal = ch1signal / 3282.8 ; Mjy./sq. degree
ch1signal = ch1signal / ((3.6D3)^2) ; Mjy / sq. arcsec
ch1signal = ch1signal * (length^2) ; Mjy
ch1signal = ch1signal * 1E6 * 1E-23 ; erg/s/cm2/Hz

ch1sigma = ch1noise/ 3282.8 ; Mjy./sq. degree
ch1sigma = ch1sigma / ((3.6D3)^2) ; Mjy / sq. arcsec
ch1sigma = ch1sigma * (length^2) ; Mjy
ch1sigma = ch1sigma * 1E6 * 1E-23 ; erg/s/cm2/Hz

ch1noise = (ch1signal + ch1sigma)  ; erg/s/cm2/Hz

ch1_abmag = flux_to_magab(ch1signal)
ch1_vegamag = ch1_abmag -2.7 
ch1_noise_abmag = flux_to_magab(ch1noise)
ch1_noise_vegamag = ch1_vegamag - (ch1_noise_abmag -2.7 )

print, 'ch1 vega', ch1_vegamag, ' +-', abs(ch1_noise_vegamag)
;print, ch1_noise_abmag - 2.7   ;limiting mag for plume.vega.param

;ch2-----------------------------------
ch2signal = ch2signal / 3282.8 ; Mjy./sq. degree
ch2signal = ch2signal / ((3.6D3)^2) ; Mjy / sq. arcsec
ch2signal = ch2signal * (length^2) ; Mjy
ch2signal = ch2signal * 1E6 * 1E-23 ; erg/s/cm2/Hz

ch2sigma = ch2noise/ 3282.8 ; Mjy./sq. degree
ch2sigma = ch2sigma / ((3.6D3)^2) ; Mjy / sq. arcsec
ch2sigma = ch2sigma * (length^2) ; Mjy
ch2sigma = ch2sigma * 1E6 * 1E-23 ; erg/s/cm2/Hz

ch2noise = (ch2signal + ch2sigma)  ; erg/s/cm2/Hz

ch2_abmag = flux_to_magab(ch2signal)
ch2_vegamag = ch2_abmag -3.25 
ch2_noise_abmag = flux_to_magab(ch2noise)
ch2_noise_vegamag = ch2_vegamag - (ch2_noise_abmag -3.25 )

print, 'ch2 vega', ch2_vegamag, ' +-', abs(ch2_noise_vegamag)
;print, ch2_noise_abmag - 3.25   ;limiting mag for plume.vega.param

;B---------------------------------------
b_vegamag = b_sb -2.5*alog10(length^2)

b_noise = b_noisesb - 2.5*alog10(length^2)

flux = magab_to_flux(b_vegamag)  ; know this is AB and not vega, but I just want the delta
flux_noise = magab_to_flux(b_noise)
flux_total = flux + flux_noise
b_vegamag_noise = b_vegamag -  flux_to_magab(flux_total)

print, 'b vega', b_vegamag, '+-', b_vegamag_noise

;V---------------------------------------
v_vegamag = v_sb -2.5*alog10(length^2)

v_noise = v_noisesb - 2.5*alog10(length^2)

flux = magab_to_flux(v_vegamag)  ; know this is AB and not vega, but I just want the delta
flux_noise = magab_to_flux(v_noise)
flux_total = flux + flux_noise
v_vegamag_noise = v_vegamag -  flux_to_magab(flux_total)

print, 'v_vega', v_vegamag, '+-',v_vegamag_noise

;hyperz-----------------------------------------

openw, outlun, '~/virgo/irac/hyperz_e.txt', /get_lun
printf, outlun,  '1', b_vegamag, v_vegamag, ch1_vegamag, ch2_vegamag, b_vegamag_noise, v_vegamag_noise, abs(ch1_noise_vegamag), abs(ch2_noise_vegamag), format = '(I10, F10.6, F10.6, F10.6, F10.6, F10.6, F10.6, F10.6, F10.6)'
close, outlun
free_lun, outlun

print, '/Users/jkrick/virgo/irac/zphot.plumes.param'
commandline = '~/nutella/bin/hyperz /Users/jkrick/virgo/irac/zphot.plumes.param'
spawn, commandline

;command2 = 'mv ~/idlbin/*.spe /Users/jkrick/virgo/irac/'
;spawn, command2

;look at the output
x = [.4424, .5504, 3.564, 4.5118]
readcol, '/Users/jkrick/virgo/irac/hyperz_e.obs_sed',idz,b, v, ch1, ch2, berr, verr, ch1err, ch2err,format="A",/debug

y = [b,v,ch1, ch2]              ;photometry
yerr = [berr, verr, ch1err, ch2err]

plot,x*1E4, y, psym = 5, thick = 3, charthick = 1, xthick = 1, ythick = 1, $
     xtitle = "Wavelength (Angstrom)", ytitle = "flux(microjansky)", /ylog,/xlog, yrange = [1E2, 1E4]

errplot, x*1E4, y - yerr, y + yerr

readcol,strcompress('/Users/jkrick/virgo/irac/1.spe', /remove_all),x2, y2,format="A",/silent
oplot, x2, y2, thick = 3, color = greencolor


;readcol, '/Users/jkrick/virgo/irac/ell2/hyperz_a.obs_sed',idz,b, v, ch1, ch2, berr, verr, ch1err, ch2err,format="A",/debug


;readcol,strcompress('/Users/jkrick/virgo/irac/ell5/1.spe', /remove_all),x2, y2,format="A",/silent;
;oplot, x2, y2, thick = 3

;readcol,strcompress('/Users/jkrick/virgo/irac/sb/1.spe', /remove_all),x2, y2,format="A",/silent
;oplot, x2, y2, thick = 3, color = redcolor

;readcol,strcompress('/Users/jkrick/virgo/irac/sc/1.spe', /remove_all),x2, y2,format="A",/silent
;oplot, x2, y2, thick = 3, color = bluecolor
  
;readcol,strcompress('/Users/jkrick/virgo/irac/sd/1.spe', /remove_all),x2, y2,format="A",/silent
;oplot, x2, y2, thick = 3, color = greencolor


end
