PRO magcompare

OPENR, lun, "/n/Godiva1/jkrick/A3888/final/SExtractor2.r.cat", /GET_LUN
corflux = fltarr(12000)
flux = fltarr(12000)
i = 0
WHILE (NOT EOF(lun)) DO BEGIN
    ;!!!! change this if you change the daofind.param , now usign godiva7/galfit/dao..
    READF, lun, o, xcenter, ycenter, a, b, e, f, m, isoarea, fwhm, pa, bkgd, apflux,isocorflux,isoflux
    corflux[i] =isocorflux
    flux[i] = isoflux
    i = i + 1
ENDWHILE

corflux = corflux[0:i-1]
flux = flux[0:i-1]

print, flux

close, lun
free_lun, lun


plot, 24.6 -2.5*alog10(flux),  24.6 -2.5*alog10(corflux), psym = 2,yrange = [15,30],xrange=[15,30]
END

