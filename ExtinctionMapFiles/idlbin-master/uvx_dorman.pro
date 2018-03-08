pro uvx_dorman
!P.multi = [0,0,0]
ps_open, filename='/Users/jkrick/Virgo/Galex/UVX/dorman_arrow.ps',/portrait,/square,/color


;readlargecol, '/Users/jkrick/Virgo/Galex/UVX/dorman_galex_sdss.txt',name,		objid,				specObjID,				distance_arcmin,			IAUName,	                           ra,	                              dec,	                 band,	fuv_mag,	nuv_mag,	fuv_flux,	nuv_flux,	e_bv,	isThereSpect, fuv_fwhm_world,	nuv_fwhm_world, nuv_mag1, nuv_mager, fuv_mag1, fuv_magerr,	nuv_mag_best,	nuv_magerr_best, nuv_mag_auto, nuv_magerr_auto, nuv_kron_radius, fuv_mag_auto, fuv_magerr_auto,	fuv_mag_best,	fuv_magerr_best,	fuv_kron_radius,     name,              objID,	                         ra,	                   dec,	         run,	   rerun, camcol, field,  type,	modelMag_u,	modelMag_g,	modelMag_r,	modelMag_i,	modelMag_z, format = 'A,F,F,F10.5,A,F10.5,F10.5,I,F10.2, F10.2, F10.2, F10.2, F10.2, A, F10.5, F10.5, F10.2, F10.5, F10.2, F10.5, F10.2, F10.5, F10.2, F10.5, F2, F10.2, F10.5, F10.2, F10.2, F10.2, A, L, F10.5, F10.5, I,I,I, A, F10.2, F10.2, F10.2, F10.2

fuv_mag = [17.25,17.77,16.65,18.14,18.37,16.60,16.57,17.68,16.27,15.86,16.50,17.07,16.84,18.53]
fuv_magerr = [.0008,.0255,.0085,.0279,.1049,.0131,.0115,.0702,.0129,.0271,.0116,.0145,.0145,.0914]

nuv_mag =[15.85,15.82,15.10,15.97,15.43,15.47,14.75,14.76,15.06,14.48,15.30,15.46,14.77,15.53]
nuv_magerr =[.0107,.0074,.0041,.0045,.0054,.0054,.0031,.0129,.0128,.0119,.0031,.0051,.0043,.00117]

g_mag = [11.88,11.79,10.51,11.32,11.4,11.51,10.74,14.05,10.98,10.79,11.27,11.21,11.14,11.98]
g_magerr = [1.52,1.58,1.61,1.69,1.57,1.51,1.49,9.71,1.46,1.55,1.45,1.47,1.52,1.68]
g_magerr = g_magerr*1E-3


x =  fuv_mag - g_mag
xerr = sqrt(fuv_magerr^2 + g_magerr^2)

y = fuv_mag - nuv_mag
yerr = sqrt(fuv_magerr^2 + nuv_magerr^2)

;plot, x, y, psym = 2, xtitle = 'FUV - g', ytitle = 'FUV - NUV',$
;      xrange = [0,10], yrange = [4,-2]

ploterror, x, y, 3*xerr, 3*yerr, psym = 2, xtitle = "FUV - g'", ytitle = 'FUV - NUV',$
      xrange = [-1,8], yrange = [4,-7], ystyle = 1, thick = 3, charthick = 3, xthick = 3, ythick = 3

;for the plume
;assume g - V = 0.2 from fukugita 1996 (B-V) = 0.6
hsb_fuv_nuv = 17.7 - 18.9
hsb_fuv_nuv_err = sqrt(.1^2 + .9^2)
hsb_fuv_V = 17.7 - 15.87
hsb_fuv_V_err = sqrt(.1^2 + .1^2)
hsb_fuv_g = hsb_fuv_V +0.2
hsb_fuv_g_err = hsb_fuv_V_err 

lsb_fuv_nuv = 16.7 - 18.23
lsb_fuv_nuv_err = sqrt(.1^2 + 1.6^2)
lsb_fuv_V =16.7 - 16.51
lsb_fuv_V_err = sqrt(.1^2 +.3^2)
lsb_fuv_g = lsb_fuv_V + 0.2
lsb_fuv_g_err=lsb_fuv_V_err

;plotsym, 2, 5
;plots,  lsb_fuv_g, lsb_fuv_nuv, psym = 5
;oploterror,  lsb_fuv_g, lsb_fuv_nuv,3*lsb_fuv_g_err,3* lsb_fuv_nuv_err, psym = 5,  errthick = 3
oploterror,  lsb_fuv_g, lsb_fuv_nuv,3*lsb_fuv_g_err, 0.0, psym = 5, errthick = 3

;plots,  hsb_fuv_g, hsb_fuv_nuv, psym = 6
;oploterror,  hsb_fuv_g, hsb_fuv_nuv,hsb_fuv_g_err, hsb_fuv_nuv_err, psym = 5,/lobar, errthick = 3
;oploterror,  hsb_fuv_g, hsb_fuv_nuv,3*hsb_fuv_g_err, 3*hsb_fuv_nuv_err, psym = 5,errthick = 3
oploterror,  hsb_fuv_g, hsb_fuv_nuv,3*hsb_fuv_g_err, 0.0, psym = 6, errthick = 3


arrow, hsb_fuv_g, hsb_fuv_nuv, hsb_fuv_g  , hsb_fuv_nuv - 3*hsb_fuv_nuv_err,/data, thick = 3, hthick = 3
arrow, lsb_fuv_g, lsb_fuv_nuv, lsb_fuv_g  , lsb_fuv_nuv - 3*lsb_fuv_nuv_err,/data, thick = 3, hthick = 3



;put the directions of the colors
;plotsym, 6,5
;plots, -1,-4,psym=8
;plotsym, 2, 5
;plots,-1,-4, psym = 8
ps_close, /noprint,/noid


readcol, '/Users/jkrick/virgo/galex/uvx/uvx_samples.csv', name,	name_other,	ra,	dec,	FUV_mag_mast,	FUV_mag_err_mast,	NUV_mag_mast,	NUV_mag_err_mast,	FUV_mag_gil,	FUV_mag_err_gil,	NUV_mag_gil,	NUV_mag_err_gil,	V_mag,	V_mag_err,	sdass_g,	sdss_g_err
end


