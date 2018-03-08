function cluster_sample

common sharevariables

;for each galaxy, how many flux measurements does it have?
nflux = get_nflux()
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


;select objects within 2rvir with cluster redshifts, and mips detections
member1 = where(sphdist(objectnew.ra, objectnew.dec, candra[0], canddec[0], /degrees) lt rvir(0) and nflux gt nfluxlimit $
               and objectnew.zphot gt 0.87 and objectnew.zphot lt 1.13  and objectnew.chisq lt 50. $
               and objectnew.mips24flux gt 17.3 and theta ge 98.3 and theta lt 278.3); and objectnew.mips24flux lt 90)
member2 = where(sphdist(objectnew.ra, objectnew.dec, candra[1], canddec[1], /degrees) lt rvir(1) and nflux gt nfluxlimit $
               and objectnew.zphot gt 0.87 and objectnew.zphot lt 1.13  and objectnew.chisq lt 50. $
               and objectnew.mips24flux gt 17.3); and objectnew.mips24flux lt 90)
member3 = where(sphdist(objectnew.ra, objectnew.dec, candra[2], canddec[2], /degrees) lt rvir(2) and nflux gt nfluxlimit $
               and objectnew.zphot gt 0.87 and objectnew.zphot lt 1.13  and objectnew.chisq lt 50. $
               and objectnew.mips24flux gt 17.3 and ycenter + 0.791*xcenter lt 32807); and objectnew.mips24flux lt 90)

;select objects within 2rvir with cluster redshifts
clus1_z1 = where(sphdist(objectnew.ra, objectnew.dec, candra[0], canddec[0], /degrees) lt rvir(0) and nflux gt nfluxlimit $
               and objectnew.zphot gt 0.87 and objectnew.zphot lt 1.13  and objectnew.chisq lt 50. and theta ge 98.3 and theta lt 278.3)
clus2_z1 = where(sphdist(objectnew.ra, objectnew.dec, candra[1], canddec[1], /degrees) lt rvir(1) and nflux gt nfluxlimit $
               and objectnew.zphot gt 0.87 and objectnew.zphot lt 1.13  and objectnew.chisq lt 50. )
clus3_z1 = where(sphdist(objectnew.ra, objectnew.dec, candra[2], canddec[2], /degrees) lt rvir(2) and nflux gt nfluxlimit  $
               and objectnew.zphot gt 0.87 and objectnew.zphot lt 1.13  and objectnew.chisq lt 50. and ycenter + 0.791*xcenter lt 32807)

;select objects within 2rvir with mips detections
clus1_24 = where(sphdist(objectnew.ra, objectnew.dec, candra[0], canddec[0], /degrees) lt rvir(0) and nflux gt nfluxlimit $
               and objectnew.mips24flux gt 17.3 and theta ge 98.3 and theta lt 278.3); and objectnew.mips24flux lt 90)
clus2_24 = where(sphdist(objectnew.ra, objectnew.dec, candra[1], canddec[1], /degrees) lt rvir(1) and nflux gt nfluxlimit $
               and objectnew.mips24flux gt 17.3); and objectnew.mips24flux lt 90)
clus3_24 = where(sphdist(objectnew.ra, objectnew.dec, candra[2], canddec[2], /degrees) lt rvir(2) and nflux gt nfluxlimit $
               and objectnew.mips24flux gt 17.3 and ycenter + 0.791*xcenter lt 32807); and objectnew.mips24flux lt 90)

;any object in the field with redshift of clusters and mips detections
;NOT in the cluster area

wholefield = where( nflux gt nfluxlimit and objectnew.zphot gt 0.87 and objectnew.zphot lt 1.13  and objectnew.chisq lt 50. $
                and objectnew.mips24flux gt 17.3 and sphdist(objectnew.ra, objectnew.dec, candra[0], canddec[0], /degrees) gt rvir(0) $
                and sphdist(objectnew.ra, objectnew.dec, candra[0], canddec[0], /degrees) gt rvir(1)  and $
                sphdist(objectnew.ra, objectnew.dec, candra[0], canddec[0], /degrees) gt rvir(2) ) 

;any object in the field at any redshift
;NOT in the cluster area

allzfield = where( nflux gt nfluxlimit $
                and objectnew.mips24flux gt 17.3 and sphdist(objectnew.ra, objectnew.dec, candra[0], canddec[0], /degrees) gt rvir(0) $
                and sphdist(objectnew.ra, objectnew.dec, candra[0], canddec[0], /degrees) gt rvir(1)  and $
                sphdist(objectnew.ra, objectnew.dec, candra[0], canddec[0], /degrees) gt rvir(2) ) 

;combine individual clusters into one big cluster array
member = [member1,member2,member3]
clus_z1 = [clus1_z1,clus2_z1,clus3_z1]
clus_24 = [clus1_24,clus2_24,clus3_24]

;sort and get rid of repeats in overlapping areas.  
member = member[UNIQ(member,sort(member))]
clus_z1 = clus_z1[UNIQ(clus_z1,sort(clus_z1))]
clus_24 = clus_24[UNIQ(clus_24,sort(clus_24))]

print, 'number of members', n_elements(member)
print, 'number of clus_z1', n_elements(clus_z1)
print, 'number of clus_24', n_elements(clus_24)

;print, 'member j mag',mean(objectnew[member].wircjmag), min(objectnew[member].wircjmag)
;plothist, objectnew[member].wircjmag, bin = 0.1, xrange=[15,25], xtitle='J mag (Vega)'

;plothyperz, member, '/Users/jkrick/nep/clusters/mips24/member.ps'

;plot where these galaxies are
openw, outlunred, '/Users/jkrick/nutella/nep/clusters/mips24/all.reg', /get_lun
printf, outlunred, 'fk5'

for rc=0, n_elements(member) -1 do  printf, outlunred, 'circle( ', objectnew[member[rc]].ra, objectnew[member[rc]].dec, ' 2") # color=red'
for rc=0, n_elements(clus_z1) -1 do  printf, outlunred, 'circle( ', objectnew[clus_z1[rc]].ra, objectnew[clus_z1[rc]].dec, ' 3.6") # color=green'

close, outlunred
free_lun, outlunred


;now figure out which ones are likely AGN, and make new samples with
;and without the AGN.

;since I am trusting photz's from hyperz, I should also trust
;SED shape, aka use hyperz to determine if AGN or not.

noagn = where(objectnew[member].spt eq 1 or objectnew[member].spt eq 2 or objectnew[member].spt eq 3 or objectnew[member].spt eq 4 or objectnew[member].spt eq 5 or objectnew[member].spt eq 8 or objectnew[member].spt eq 9 or objectnew[member].spt eq 10 or objectnew[member].spt eq 11 or objectnew[member].spt eq 12)
membernoagn = member[noagn]

agn = where(objectnew[member].spt eq 6 or objectnew[member].spt eq 7 or objectnew[member].spt eq 13$
or objectnew[member].spt eq 14 or objectnew[member].spt eq 15)
memberagn = member[agn]

noagn = where(objectnew[wholefield].spt eq 1 or objectnew[wholefield].spt eq 2 or objectnew[wholefield].spt eq 3 or objectnew[wholefield].spt eq 4 or objectnew[wholefield].spt eq 5 or objectnew[wholefield].spt eq 8 or objectnew[wholefield].spt eq 9 or objectnew[wholefield].spt eq 10 or objectnew[wholefield].spt eq 11 or objectnew[wholefield].spt eq 12)
wholefieldnoagn = wholefield[noagn]

agn = where(objectnew[wholefield].spt eq 6 or objectnew[wholefield].spt eq 7 or objectnew[wholefield].spt eq 13$
or objectnew[wholefield].spt eq 14 or objectnew[wholefield].spt eq 15)
wholefieldagn = wholefield[agn]

noagn = where(objectnew[clus_z1].spt eq 1 or objectnew[clus_z1].spt eq 2 or objectnew[clus_z1].spt eq 3 or objectnew[clus_z1].spt eq 4 or objectnew[clus_z1].spt eq 5 or objectnew[clus_z1].spt eq 8 or objectnew[clus_z1].spt eq 9 or objectnew[clus_z1].spt eq 10 or objectnew[clus_z1].spt eq 11 or objectnew[clus_z1].spt eq 12)
clus_z1noagn = clus_z1[noagn]


print, 'number of membernoagns', n_elements(membernoagn)
print, 'number of memberagns', n_elements(memberagn)
print, 'number of wholefieldnoagns', n_elements(wholefieldnoagn)

print, 'memberagn', memberagn

;print, 'checking fluxes: memberagn', objectnew[memberagn].mips24flux
;print, 'checking fluxes: membernoagn', objectnew[membernoagn].mips24flux

;plothist, objectnew[memberagn].mips24flux, bin = 1, xtitle='memberagn mips24flux', xrange=[0,50]
;plothist, objectnew[membernoagn].mips24flux, bin = 1, xtitle='membernoagn mips24flux', xrange=[0,50]

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
