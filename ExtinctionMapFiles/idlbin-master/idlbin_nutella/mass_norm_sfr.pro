function mass_norm_sfr

common sharevariables

rad = 0.25*rvir  ;want 0.5*r200, and here rvir is actually 2*rvir

;for cluster one need to only use half the cluster area around angle 278
acshead = headfits('/Users/jkrick/nutella/hst/raw/wholeacs.fits');
adxy, acshead, objectnew.ra,objectnew.dec , xcenter, ycenter
adxy, acshead, candra[0], canddec[0], x0, y0
adxy, acshead, candra[1], canddec[1], x1, y1
adxy, acshead, candra[2], canddec[2], x2, y2

dx = xcenter - x0
dy = ycenter - y0
theta = lonarr(n_elements(dx))
for i = long(0), n_elements(dx) -1 do begin
   if dx(i) ge 0. and dy(i) ge 0. then theta(i) = 180/!PI*atan(dy(i)/dx(i))
   if dx(i) le 0. and dy(i) ge 0. then theta(i) = 180. + 180/!PI*atan(dy(i)/dx(i))
   if dx(i) le 0. and dy(i) le 0. then theta(i) =  180. + 180/!PI*atan(dy(i)/dx(i))
   if dx(i) ge 0. and dy(i) le 0. then theta(i) = 360. + 180/!PI*atan(dy(i)/dx(i))
endfor


;select objects within 2rvir with cluster redshifts, and mips detections
member1 = where(sphdist(objectnew.ra, objectnew.dec, candra[0], canddec[0], /degrees) lt rad(0) and nflux gt nfluxlimit $
               and objectnew.zphot gt 0.87 and objectnew.zphot lt 1.13  and objectnew.chisq lt 50. $
               and objectnew.mips24flux gt 17.3 and theta ge 98.3 and theta lt 278.3); and objectnew.mips24flux lt 90)
member2 = where(sphdist(objectnew.ra, objectnew.dec, candra[1], canddec[1], /degrees) lt rad(1) and nflux gt nfluxlimit $
               and objectnew.zphot gt 0.87 and objectnew.zphot lt 1.13  and objectnew.chisq lt 50. $
               and objectnew.mips24flux gt 17.3); and objectnew.mips24flux lt 90)
member3 = where(sphdist(objectnew.ra, objectnew.dec, candra[2], canddec[2], /degrees) lt rad(2) and nflux gt nfluxlimit $
               and objectnew.zphot gt 0.87 and objectnew.zphot lt 1.13  and objectnew.chisq lt 50. $
               and objectnew.mips24flux gt 17.3 and ycenter + 0.791*xcenter lt 32807); and objectnew.mips24flux lt 90)



;combine individual clusters into one big cluster array
membersmall = [member1,member2,member3]

;sort and get rid of repeats in overlapping areas.  
membersmall = membersmall[UNIQ(membersmall,sort(membersmall))]
print, '----------------------------------------------------'
print, 'number of membersmalls within 0.5rvir ', n_elements(membersmall)

;make sure not including AGN in there.  based on hyperz fits
notagn = where(objectnew[membersmall].spt eq 1 or objectnew[membersmall].spt eq 2 or objectnew[membersmall].spt eq 3 or objectnew[membersmall].spt eq 4 or objectnew[membersmall].spt eq 5 or objectnew[membersmall].spt eq 8 or objectnew[membersmall].spt eq 9 or objectnew[membersmall].spt eq 10 or objectnew[membersmall].spt eq 11 or objectnew[membersmall].spt eq 12)

;get total lir without AGN
get_lir, objectnew[membersmall[notagn]].zphot, objectnew[membersmall[notagn]].mips24flux, objectnew[membersmall[notagn]].mips24fluxerr, lirnoagn, lirnoagnerr

lirtosfr=1.71E-10
sfr = lirnoagn * lirtosfr
sfrerr = lirnoagnerr * lirtosfr
totalsfrerr = total(sfrerr)
totalsfr  = total(sfr)
print, 'total, mean sfr inside 0.5r200', totalsfr, mean(sfr)

;but I really want mass_normalized sfr.
;how do I do this with three clusters.  I think I need to divide by
;three masses since I have three clusters worth of sfr added together?
;but is the mass also supposed to be within 0.5R200?

mass = [.62,.36,.36]  ;times 1E14
masserr = [.14,.14,.14]
totalmass = total(mass)
totalmasserr = total(masserr)

massnormsfr =  totalsfr / totalmass
massnormsfrsq = massnormsfr^2

;now for the error bars
;mns = massnormalizedsfr

sigmamnssq = massnormsfrsq * ((totalmasserr^2 / totalmass^2) + (totalsfrerr^2 / totalsfr^2))
sigmamns = sqrt(sigmamnssq)

print, 'mass normalized sfr', massnormsfr, sigmamns

;reproduce Bai plot?
; will have it if I want it.

sfr = [558,201,52.5,14,64,307,753,29.1,21.6,80,71,63,90,369]
sfrerr = [279,100.5,26.3,7,17,153.5,376.5,3.1,19.5,28,23,12,19,55]
mass = [11,4.5,14,4.8,19.3,10,5.7,7.1,7.3,13.6,2.3,.55,4.8,9.5]  ;times 1E14
masserr=[1.0,2.7,4.0,.4,2.0,7.0,1.1,1.5,3.15,.7,1.2,.615,1.45,1.8]
z=[.83,.83,.023,.176,.183,.226,.39,.022,.31,.228,.84,.704,.748,.794]

low = [1,3,10, 11, 12] ; masses under 5E14

mnsbai = sfr / mass
sigmamnsbaisq = (mnsbai^2) * ((masserr^2 / mass^2) + (sfrerr^2 /sfr^2))
sigmamnsbai = sqrt(sigmamnsbaisq)

plot, z, mnsbai, /ylog, xrange=[0,1.1], yrange=[1,1000], psym = 5, xtitle = 'redshift', ytitle='Mass-normalized SFR'
errplot, z, mnsbai - sigmamnsbai, mnsbai + sigmamnsbai
oplot, z(low), mnsbai(low), psym = 6

z = 1.0
xyouts, z,  massnormsfr, '*', alignment=0.5
errplot, z, massnormsfr - sigmamns, massnormsfr + sigmamns



return, 0

end



;-------------------------------------------------------------------------------
;how much area is covered by the chosen sample?
;only had to run this once
;-------------------------------------------------------------------------------

;check that upper right quadrant in ACS pixels to see how many are
;covered by the area chosen above.
;x = indgen(21000- 12000) + 12000
;y = indgen(19800 - 10100) + 10100


;npix = long(0)
;for i = 0, n_elements(x) -1 do begin
;   for j = 0, n_elements(y) - 1 do begin
;      true = 0

;      if (x(i) - x0)^2 + (y(j) - y0)^2 le 3491.^2  then true = 1
;      if (x(i) - x1)^2 + (y(j) - y1)^2 le 2879.^2 then true = 1
;      if (x(i) - x2)^2 + (y(j) - y2)^2 le 2879.^2 and  y(j) + 0.791*x(i) lt 32807 then true = 1

;      if true gt 0 then npix = npix + 1
;   endfor
;endfor
;
;print,'total number of pixels', npix

;convert npix into radius of the equivalent circle in degrees.

;npix into square arcseconds
;0.05"/pix means the area of one pixel is .05^2
;area = npix*.05^2  ;in square arcseconds
;area = area * ((1/60)^2)* ((1/60)^2) ;now in square degrees

;print, area
