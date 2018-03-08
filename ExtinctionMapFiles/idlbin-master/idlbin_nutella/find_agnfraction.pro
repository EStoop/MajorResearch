function find_agnfraction
common sharevariables

;what is the AGN fraction for these clusters based on definition of
;eastman2007?


;get 5.8 luminosities
;1jy = 1E-23 erg/s/cm^2/Hz

;first conver flux in microjy to erg/s/cm^2/Hz
jytoerg = 1D-23
fch3 = (objectnew.irac3flux /1D6) * jytoerg

;now convert lumdist into cm.
;luminositydistance = lumdist(1.0)
;luminositydistance = luminositydistance*1D6
;lumdis = lumdis *3.0856D18


;lch3  = 4.*!Pi*lumdis^2.*fch3 / (1+1.)

Lch3 = (((fch3*4*!Pi*luminositydistance^2)/(1 + 1.0)) * ((3.0856D24)^2) ) ; include conversion from cm to Mpc 
Lch3 = Lch3*5.16E13; (Hz)

;print, 'lch3', Lch3[0:200]

lumdis = luminositydistance*1E6
mu = alog10(lumdis)
mu = 5*mu
mu =mu - 5

;print, 'mu', mu
rvir = 0.5*rvir
;print, 'rvir', rvir
Mr = objectnew.rmag - mu

;print, 'mr', objectnew[member].rmag
;print, 'Mr', objectnew[member].rmag  - mu
nfluxlimit =5

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

;select objects within rvir with cluster redshifts with Mr < -20
fracmember1 = where(sphdist(objectnew.ra, objectnew.dec, candra[0], canddec[0], /degrees) lt rvir(0) and nflux gt nfluxlimit $
               and objectnew.zphot gt 0.87 and objectnew.zphot lt 1.13  and objectnew.chisq lt 50. $
               and objectnew.rmag le 24.8 and theta lt 278.3 );and Lch3 gt 2.5E43); and objectnew.mips24flux lt 90)
fracmember2 = where(sphdist(objectnew.ra, objectnew.dec, candra[1], canddec[1], /degrees) lt rvir(1) and nflux gt nfluxlimit $
               and objectnew.zphot gt 0.87 and objectnew.zphot lt 1.13  and objectnew.chisq lt 50. $
               and objectnew.rmag le 24.8); and Lch3 gt 2.5E43); and objectnew.mips24flux lt 90) 
fracmember3 = where(sphdist(objectnew.ra, objectnew.dec, candra[2], canddec[2], /degrees) lt rvir(2) and nflux gt nfluxlimit $
               and objectnew.zphot gt 0.87 and objectnew.zphot lt 1.13  and objectnew.chisq lt 50. $
               and objectnew.rmag le 24.8 and ycenter + 0.791*xcenter lt 32807); and Lch3 gt 2.5E43); and objectnew.mips24flux lt 90) 

fracmember = [fracmember1,fracmember2,fracmember3]
fracmember = fracmember[UNIQ(fracmember,sort(fracmember))]
print, 'number of fracmembers', n_elements(fracmember)

;print, 'lch3 of fracmembers', lch3[fracmember]
a = where(lch3[fracmember] gt 3.5D43)
print, 'a', n_elements(a)
agn = where(objectnew[fracmember].spt eq 6 or objectnew[fracmember].spt eq 7 or objectnew[fracmember].spt eq 13$
or objectnew[fracmember].spt eq 14 or objectnew[fracmember].spt eq 15 and objectnew[fracmember].mips24flux gt 17.3 and lch3[fracmember] gt 3.5D43)

;agn = where(lch3[fracmember] gt 7D43)
fracmemberagn = fracmember[agn]

print, 'number of fracmemberagns', n_elements(fracmemberagn)

openw, outlunred, '/Users/jkrick/nutella/nep/fracmember.reg', /get_lun
printf, outlunred, 'fk5'
for rc=0, n_elements(fracmember) -1 do  printf, outlunred, 'circle( ', objectnew[fracmember[rc]].ra, objectnew[fracmember[rc]].dec, ' 3")'
close, outlunred
free_lun, outlunred

;plothyperz, fracmemberagn, '/Users/jkrick/nutella/nep/clusters/mips24/fracagn.ps'

print, 'mips', objectnew[fracmemberagn].mips24flux
print, 'spt', objectnew[fracmemberagn].spt
; 'irac colors', alog10(objectnew[fracmemberagn].irac3flux / objectnew[fracmemberagn].irac1flux), alog10(objectnew[fracmemberagn].irac4flux / objectnew[fracmemberagn].irac2flux),

;get 5.8 luminosities







;put it back 
rvir = 2*rvir
return, 0
end
