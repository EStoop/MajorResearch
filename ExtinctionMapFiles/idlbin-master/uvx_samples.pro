pro uvx_samples
ps_open, filename='/Users/jkrick/Virgo/Galex/UVX/all_samples.ps',/portrait,/square,/color

!P.multi = [0,1,1]
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
greencolor = FSC_COLOR("Green", !D.Table_Size-3)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-4)

readcol, '/Users/jkrick/virgo/galex/uvx/uvx_samples.csv', name,	name_other, spiral,	ra,	dec,	FUV_mag_mast,	FUV_mag_err_mast,	NUV_mag_mast,	NUV_mag_err_mast,	FUV_mag_gil,	FUV_mag_err_gil,	NUV_mag_gil,	NUV_mag_err_gil,	V_mag,	V_mag_err,	sdss_g,	sdss_g_err, format = 'A,A,A,F10.5, F10.5, F10.2, F10.5, F10.2, F10.5, F10.2, F10.5, F10.2, F10.5, F10.2, F10.5, F10.2, F10.5'


;spirals
plotsym, 0, /fill
gil = where(FUV_mag_gil gt 0 and spiral eq 'Y' and FUV_mag_gil gt 0 and NUV_mag_gil gt 0 and V_mag gt 0)
print, 'virgo', n_elements(gil)
x = FUV_mag_gil(gil) - V_mag(gil)
xerr = sqrt(FUV_mag_err_gil(gil)^2 + V_mag_err(gil)^2)
y = FUV_mag_gil(gil) - NUV_mag_gil(gil)
yerr =  sqrt(FUV_mag_err_gil(gil)^2 + NUV_mag_err_gil(gil)^2)
ploterror,x , y, 3*xerr, 3*yerr, psym = 8, xrange = [-2, 10], yrange = [4,-8], xtitle = 'FUV - V', ytitle = 'FUV - NUV', thick = 3, xthick = 3, ythick =3, charthick = 3


;uvx ellipticals
plotsym, 4, /fill
ell =  where(FUV_mag_gil gt 0 and spiral eq 'N')
print, 'uvx ellipticals', n_elements(ell)
x = FUV_mag_gil(ell) - V_mag(ell)
xerr = sqrt(FUV_mag_err_gil(ell)^2 + V_mag_err(ell)^2)
y = FUV_mag_gil(ell) - NUV_mag_gil(ell)
yerr =  sqrt(FUV_mag_err_gil(ell)^2 + NUV_mag_err_gil(ell)^2)
oploterror,x , y, 3*xerr, 3*yerr, psym = 8, color = redcolor, errcolor=redcolor, errthick = 3

;ICL plume
plotsym, 8, /fill
x = FUV_mag_mast(30) - V_mag(30)
xerr =  sqrt(FUV_mag_err_mast(30)^2 + V_mag_err(30)^2)
y = FUV_mag_mast(30) - NUV_mag_mast(30)
yerr = sqrt(FUV_mag_err_mast(30)^2 + NUV_mag_err_mast(30)^2)
oploterror,x  , y, xerr, yerr, psym = 8, color = greencolor, errcolor = greencolor, errthick = 3

plotsym, 3, /fill
x = FUV_mag_mast(31) - V_mag(31)
xerr =  sqrt(FUV_mag_err_mast(31)^2 + V_mag_err(31)^2)
y = FUV_mag_mast(31) - NUV_mag_mast(31)
yerr = sqrt(FUV_mag_err_mast(31)^2 + NUV_mag_err_mast(31)^2)
oploterror,x  , y, 3*xerr, 3*yerr, psym = 8, color = greencolor, errcolor = greencolor, errthick = 3


;what about the rest of the gildepaz sample
;it is in a fits file
 ftab_ext, '/Users/jkrick/virgo/galex/uvx/gildepaz.fits', 'Morph,Name,asyFUV,e_asyFUV,asyNUV,e_asyNUV,Vmag,e_Vmag', type, gal_name, FUV, eFUV, NUV, eNUV, V, eV
good = (where(FUV gt 0 and NUV gt 0 and V gt 0))

;spirals
plotsym, 0, /fill
sp = where(type(good) eq 'Sa' or type(good) eq 'Sb' or type(good) eq 'Sc' or type(good) eq 'Sd' or type(good) eq 'Sab' or type(good) eq 'Sbc' or type(good) eq 'Scd')
x =  FUV(good(sp)) - V(good(sp))
xerr =  sqrt(eFUV(good(sp))^2 + eV(good(sp))^2)
y =  FUV(good(sp))  - NUV(good(sp))
yerr =  sqrt(eFUV(good(sp))^2  + eNUV(good(sp))^2)
oploterror,x,y,3*xerr, 3*yerr, psym = 8, errthick = 3

;ellipticals
plotsym, 4, /fill
el = where(type(good) eq 'E' or type(good) eq 'E0' or type(good) eq 'E1' or type(good) eq 'E2' or type(good) eq 'E3' or type(good) eq 'E4' or type(good) eq 'E5' or type(good) eq 'E6' or type(good) eq 'E7', count)
print, 'el', count
x =  FUV(good(el)) - V(good(el))
xerr =  sqrt(eFUV(good(el))^2 + eV(good(el))^2)
y =  FUV(good(el))  - NUV(good(el))
yerr =  sqrt(eFUV(good(el))^2  + eNUV(good(el))^2)
oploterror,x,y,3*xerr, 3*yerr,psym = 8, color = redcolor, errcolor=redcolor, errthick = 3

q = where(FUV(good(el)) - V(good(el)) lt 4)
print, gal_name(good(el(q))), V(good(el(q)))

; add tidal features
;from Neff.csv

FUV_r = [3.12,5.24,2.71,2.11,6.33,3.24,3.93,1.09,2.07,6.36,3.42,3.82,4.68,4.65,1.78,2.46,4.55,5.2,4.59,3.4,4.78,4.64,3.89,5.25]
FUV_V = [2.91,5.03,2.5,1.9,6.12,3.03,3.72,0.88,1.86,6.15,3.21,3.61,4.47,4.44,1.57,2.25,4.34,4.99,4.38,3.19,4.57,4.43,3.68,5.04]
FUV_NUV = [0.26,0.78,0.12,-0.32,0.89,-0.03,0.25,-0.25,-0.24,1.04,0.39,0.34,1,0.12,-0.09,0.1,0.78,0.96,0.71,0.4,1.16,0.63,0.7,1.09]
xerr = [0.212132034,0.122474487,0.157162336,0.651766829,0.15,0.227815715,0.29563491,0.427901858,0.378549865,0.352703842,0.137477271,0.439886349,1.083374358,0.233666429,0.334065862,0.398120585,0.167928556,0.16643317,0.189208879,0.242280829,0.189208879,0.189208879,0.218860686,0.214009346]
yerr = [0.13,0.08,0.1,0.4,0.09,0.14,0.18,0.25,0.23,0.21,0.09,0.27,0.65,0.15,0.21,0.23,0.06,0.06,0.08,0.13,0.08,0.08,0.11,0.1]

;now without the central galaxies
FUV_V=[2.5,1.9,3.03,3.72,0.88,1.86,3.61,4.47,4.44,1.57,2.25,4.57,4.43,3.68,5.04]
FUV_NUV=[0.12,-0.32,-0.03,0.25,-0.25,-0.24,0.34,1,0.12,-0.09,0.1,1.16,0.63,0.7,1.09]
eFUV_NUV=[0.1,0.4,0.14,0.18,0.25,0.23,0.27,0.65,0.15,0.21,0.23,0.08,0.08,0.11,0.1]
eFUV_V=[0.157162336,0.651766829,0.227815715,0.29563491,0.427901858,0.378549865,0.439886349,1.083374358,0.233666429,0.334065862,0.398120585,0.189208879,0.189208879,0.218860686,0.214009346]


oploterror, FUV_V, FUV_NUV, eFUV_V, eFUV_NUV, psym = 2, color = bluecolor, errcolor = bluecolor

ps_close, /noprint,/noid

end


