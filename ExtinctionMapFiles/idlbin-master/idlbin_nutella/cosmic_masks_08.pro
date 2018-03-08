pro cosmic_masks
!P.multi = [0,0,1]
close, /all
openw, outlun, '/Users/jkrick/palomar/cosmic/slitmasks/june08/masklist_junk.txt', /get_lun 
openw, outlun2, '/Users/jkrick/palomar/cosmic/slitmasks/june08/targets_junk.reg',/get_lun
printf, outlun2, 'fk5'
openw, outlun3, '/Users/jkrick/palomar/cosmic/slitmasks/june08/crossref_good.txt',/get_lun
restore, '/Users/jkrick/idlbin/objectnew.sav'

nflux = get_nflux()

hd=headfits('/Users/jkrick/hst/raw/wholeacs.fits')
rhd = headfits('/Users/jkrick/palomar/lfc/coadd_r.fits')
;ihd = headfits('/Users/jkrick/palomar/lfc/coadd_i.fits')



;find the stars
;don't care if they are saturated

;put Morgan's stars in the list
priority = 1
morganstars = [9878,685,11759,10200]
radec, objectnew[morganstars].ra, objectnew[morganstars].dec, ihr, imin, xsec, ideg, imn, xsc

for i = 0, n_elements(morganstars) -1 do begin
   printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, $
           ihr[i], imin[i], xsec[i], ideg[i], imn[i], xsc[i], objectnew[morganstars[i]].rmag
;   printf, outlun2, 'circle( ', objectnew[morganstars[i]].ra, objectnew[morganstars[i]].dec, ' 3")'
   printf, outlun3, priority, morganstars[i], objectnew[morganstars[i]].ra, objectnew[morganstars[i]].dec
   priority = priority + 1
endfor

;all other stars  (will repeat Morgan stars but no matter)
priority = 10
star = where( objectnew.rmag lt 18.5  and objectnew.rmag gt 0 and objectnew.rfwhm lt 8.0 and nflux gt 3)
badstars = [17687,12732,1331,1274,1135,761]
for j = 0, n_elements(badstars) - 1 do begin
   star = star(where(star ne badstars(j)))
endfor

print, 'stars', n_elements(star)
radec, objectnew[star].ra, objectnew[star].dec, ihr, imin, xsec, ideg, imn, xsc

for i = 0, n_elements(star) -1 do begin
   printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, $
           ihr[i], imin[i], xsec[i], ideg[i], imn[i], xsc[i], objectnew[star[i]].rmag
   printf, outlun2, 'circle( ', objectnew[star[i]].ra, objectnew[star[i]].dec, ' 5")'
   printf, outlun3, priority, star[i], objectnew[star[i]].ra, objectnew[star[i]].dec
   priority = priority + 1

endfor
;plothyperz, star, '/Users/jkrick/palomar/cosmic/slitmasks/june08/stars.ps'
;---------------------------------
;add the targeted galaxies to the list
priority = 200
;read in Jason's list.
listjason =[   2358] 
for i = 0, n_elements(listjason) - 1 do begin
   radec, objectnew[listjason(i)].ra, objectnew[listjason(i)].dec, ihr, imin, xsec, ideg, imn, xsc
   printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, ihr, imin, xsec, ideg, imn, xsc, objectnew[listjason(i)].rmag
   printf, outlun2, 'circle( ', objectnew[listjason[i]].ra, objectnew[listjason[i]].dec, ' 3")'
   printf, outlun3, priority, listjason[i], objectnew[listjason[i]].ra, objectnew[listjason[i]].dec

   priority = priority + 1
endfor


;---------------------------------
;find interesting galaxies



;70 micron
priority = 300
seventy = where(objectnew.rmag le 21.0 and objectnew.acsclassstar lt 1. and $
                objectnew.mips70mag gt 0 and objectnew.mips70mag lt 90  and $
                nflux gt 2 and objectnew.rmag gt 0.)
radec, objectnew[seventy].ra, objectnew[seventy].dec, ihr, imin, xsec, ideg, imn, xsc
for i = 0, n_elements(seventy) -1 do begin
   printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, $
           ihr[i], imin[i], xsec[i], ideg[i], imn[i], xsc[i], objectnew[seventy[i]].rmag
   printf, outlun2, 'circle( ', objectnew[seventy[i]].ra, objectnew[seventy[i]].dec, ' 3")'
   printf, outlun3, priority, seventy[i], objectnew[seventy[i]].ra, objectnew[seventy[i]].dec
   priority = priority + 1
    
endfor
;plothyperz, seventy, '/Users/jkrick/palomar/cosmic/slitmasks/june08/mips70.ps'

;24 micron
priority = 400
;use acs fwhm to get rid of some stars, but also want to include those galaxies which do not have acs detections.
mips24 = where(objectnew.rmag le 21.0 and objectnew.acsclassstar lt 1. and $
                objectnew.mips24mag gt 0 and objectnew.mips24mag lt 90  and $
               nflux gt 2 and objectnew.acsfwhm gt 3.8 and objectnew.rmag gt 17.0)

mips24_2 = where(objectnew.rmag le 21.0 and objectnew.acsclassstar lt 1. and $
                objectnew.mips24mag gt 0 and objectnew.mips24mag lt 90  and $
                 nflux gt 2 and objectnew.acsmag lt 0  and objectnew.rmag gt 17.0)
mips24 = [mips24, mips24_2, 2655]

badstars = [390,2200,2247,2285]
for j = 0, n_elements(badstars) - 1 do begin
   mips24 = mips24(where(mips24 ne badstars(j)))
endfor


print,'mips24', n_elements(mips24)

radec, objectnew[mips24].ra, objectnew[mips24].dec, ihr, imin, xsec, ideg, imn, xsc
for i = 0, n_elements(mips24) -1 do begin
   printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, $
           ihr[i], imin[i], xsec[i], ideg[i], imn[i], xsc[i], objectnew[mips24[i]].rmag
   printf, outlun2, 'circle( ', objectnew[mips24[i]].ra, objectnew[mips24[i]].dec, ' 3")'
   printf, outlun3, priority, mips24[i], objectnew[mips24[i]].ra, objectnew[mips24[i]].dec
   priority = priority + 1

endfor
;plothyperz, mips24, '/Users/jkrick/palomar/cosmic/slitmasks/june08/mips24_2.ps'

;Chandra
priority = 600
chandra = where(objectnew.rmag le 21.0 and objectnew.acsclassstar lt 1. and $
                objectnew.xnetrate gt 0 and nflux gt 2 and objectnew.rmag gt 0.)
radec, objectnew[chandra].ra, objectnew[chandra].dec, ihr, imin, xsec, ideg, imn, xsc
for i = 0, n_elements(chandra) -1 do begin
   printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, $
           ihr[i], imin[i], xsec[i], ideg[i], imn[i], xsc[i], objectnew[chandra[i]].rmag
   printf, outlun2, 'circle( ', objectnew[chandra[i]].ra, objectnew[chandra[i]].dec, ' 3")'
   printf, outlun3, priority, chandra[i], objectnew[chandra[i]].ra, objectnew[chandra[i]].dec
   priority = priority + 1

endfor



;---------------------------------
; find the random galaxies
priority = 1000
priority2 = 2000
priority3 = 3000
priority4 = 4000
priority5 = 5000


zarr = objectnew[where(objectnew.rmag le 21.0 and objectnew.rmag gt 17.0 and objectnew.acsclassstar lt 1.)].zphot
plothist,  zarr, xhist, yhist, bin = 0.05, /noprint, xrange=[0,2]

z1 = where(objectnew.rmag le 21.0 and objectnew.rmag gt 17.0 and objectnew.acsclassstar lt 1. $
           and objectnew.zphot gt 0.05 and objectnew.zphot le 0.2 and nflux ge 5 and $
           objectnew.acsfwhm gt 3.8 and objectnew.rfwhm ge 8.0) 
z1_2 = where(objectnew.rmag le 21.0 and objectnew.rmag gt 17.0 and objectnew.acsclassstar lt 1. $
           and objectnew.zphot gt 0.05 and objectnew.zphot le 0.2 and nflux ge 5 and $
             objectnew.acsmag lt 0.and objectnew.rfwhm ge 8.0) 
z1 = [z1, z1_2]
z1bad = [2200,2247,3667,4066,14576]
for j = 0, n_elements(z1bad) - 1 do begin
   z1 = z1(where(z1 ne z1bad(j)))
endfor
print, 'z1', n_elements(z1)

z2 = where(objectnew.rmag le 21.0 and objectnew.rmag gt 17.0 and objectnew.acsclassstar lt 1. $
           and objectnew.zphot gt 0.2 and objectnew.zphot le 0.4 and nflux ge 3 and $
           objectnew.acsfwhm gt 3.8 and objectnew.rfwhm ge 8.0) 
z2_2 = where(objectnew.rmag le 21.0 and objectnew.rmag gt 17.0 and objectnew.acsclassstar lt 1. $
           and objectnew.zphot gt 0.2 and objectnew.zphot le 0.4 and nflux ge 3 and $
             objectnew.acsmag lt 0. and objectnew.rfwhm ge 8.0) 
z2 = [z2,z2_2]
z2bad = [4982,10214]
for j = 0, n_elements(z2bad) - 1 do begin
   z2 = z2(where(z2 ne z2bad(j)))
endfor
print, 'z2', n_elements(z2)

z3 = where(objectnew.rmag le 21.0 and objectnew.rmag gt 17.0 and objectnew.acsclassstar lt 1. $
           and objectnew.zphot gt 0.4 and objectnew.zphot le 0.8 and nflux ge 3 and $
           objectnew.acsfwhm gt 3.8 and objectnew.rfwhm ge 8.0) 
z3_2 = where(objectnew.rmag le 21.0 and objectnew.rmag gt 17.0 and objectnew.acsclassstar lt 1. $
           and objectnew.zphot gt 0.4 and objectnew.zphot le 0.8 and nflux ge 3 and $
             objectnew.acsmag lt 1.0 and objectnew.rfwhm ge 8.0) 
z3 = [z3,z3_2]
z3bad=[4432,4644]
for j = 0, n_elements(z3bad) - 1 do begin
   z3 = z3(where(z3 ne z3bad(j)))
endfor
print, 'z3', n_elements(z3)

z4 = where(objectnew.rmag le 21.0 and objectnew.rmag gt 17.0 and objectnew.acsclassstar lt 1. $
           and objectnew.zphot gt 0.8  and nflux ge 3 and objectnew.acsfwhm gt 3.8 and objectnew.rfwhm ge 8.0) 
z4_2 = where(objectnew.rmag le 21.0 and objectnew.rmag gt 17.0 and objectnew.acsclassstar lt 1. $
           and objectnew.zphot gt 0.8  and nflux ge 3 and objectnew.acsmag lt 0. and objectnew.rfwhm ge 8.0) 
z4 = [z4,z4_2]
z4bad = [31454,35352,35636]
for j = 0, n_elements(z4bad) - 1 do begin
   z4 = z4(where(z4 ne z4bad(j)))
endfor
print, 'z4', n_elements(z4)

;plothyperz, z1, '/Users/jkrick/palomar/cosmic/slitmasks/june08/z1.ps'
;plothyperz, z2, '/Users/jkrick/palomar/cosmic/slitmasks/june08/z2.ps'
;plothyperz, z3, '/Users/jkrick/palomar/cosmic/slitmasks/june08/z3.ps'
;plothyperz, z4, '/Users/jkrick/palomar/cosmic/slitmasks/june08/z4.ps'

radec, objectnew[Z1].ra, objectnew[Z1].dec, ihr, imin, xsec, ideg, imn, xsc
for i = 0, n_elements(z1) - 1 do begin
   printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, $
           ihr[i], imin[i], xsec[i], ideg[i], imn[i], xsc[i], objectnew[z1[i]].rmag
   printf, outlun2, 'circle( ', objectnew[z1[i]].ra, objectnew[z1[i]].dec, ' 3")'
   printf, outlun3, priority, z1[i], objectnew[z1[i]].ra, objectnew[z1[i]].dec
   priority = priority + 1
endfor
 
radec, objectnew[Z2].ra, objectnew[Z2].dec, ihr, imin, xsec, ideg, imn, xsc
for i = 0, n_elements(z2) - 1 do begin
   printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority2, $
           ihr[i], imin[i], xsec[i], ideg[i], imn[i], xsc[i], objectnew[z1[i]].rmag
   printf, outlun2, 'circle( ', objectnew[z2[i]].ra, objectnew[z2[i]].dec, ' 3")'
   printf, outlun3, priority2, z2[i], objectnew[z2[i]].ra, objectnew[z2[i]].dec
   priority2 = priority2 + 1
endfor

radec, objectnew[Z3].ra, objectnew[Z3].dec, ihr, imin, xsec, ideg, imn, xsc
for i = 0, n_elements(z3) - 1 do begin
   printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority3, $
           ihr[i], imin[i], xsec[i], ideg[i], imn[i], xsc[i], objectnew[z3[i]].rmag
   printf, outlun2, 'circle( ', objectnew[z3[i]].ra, objectnew[z3[i]].dec, ' 3")'
   printf, outlun3, priority3, z3[i], objectnew[z3[i]].ra, objectnew[z3[i]].dec
   priority3 = priority3 + 1
endfor

radec, objectnew[Z4].ra, objectnew[Z4].dec, ihr, imin, xsec, ideg, imn, xsc
for i = 0, n_elements(z4) - 1 do begin
   printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority4, $
           ihr[i], imin[i], xsec[i], ideg[i], imn[i], xsc[i], objectnew[z4[i]].rmag
   printf, outlun2, 'circle( ', objectnew[z4[i]].ra, objectnew[z4[i]].dec, ' 3")'
   printf, outlun3, priority4, z4[i], objectnew[z4[i]].ra, objectnew[z4[i]].dec
   priority4 = priority4 + 1
endfor

;add morgan's fainter stars just in case
priority9 = 9000
morganfaint = [752,794,11726]
radec, objectnew[morganfaint].ra, objectnew[morganfaint].dec, ihr, imin, xsec, ideg, imn, xsc
for i = 0, n_elements(morganfaint) - 1 do begin
   printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority9, $
           ihr[i], imin[i], xsec[i], ideg[i], imn[i], xsc[i], objectnew[morganfaint[i]].rmag
   printf, outlun3, priority9, morganfaint[i], objectnew[morganfaint[i]].ra, objectnew[morganfaint[i]].dec
   priority9 = priority9 + 1
endfor
    
 
print, "there are ",priority - 1000, " objectnews with mr < 21 and z < 0.2"
print, "there are ", priority2-2000, " objectnews with mr < 21 and 0.2 < z > 0.4"
print, "there are ", priority3-3000, " objectnews with mr < 21 and 0.4 < z > 0.8"
print, "there are ", priority4-4000, " objectnews with mr < 21 and 0.8 < z "

;--------------------------------------------------



close, outlun
free_lun, outlun
close, outlun2
free_lun, outlun2
close, outlun3
free_lun, outlun3

end



;things with photz gt 1.0
;priority = 300
;distant = where(objectnew.rmag le 21.0 and objectnew.acsclassstar lt 1. and $
;                objectnew.zphot gt 1.0 and objectnew.prob gt 80  and $;
;                nflux gt 2 and objectnew.rmag gt 0.)
;radec, objectnew[distant].ra, objectnew[distant].dec, ihr, imin, xsec, ideg, imn, xsc
;for i = 0, n_elements(distant) -1 do begin
;   printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, $
;           ihr[i], imin[i], xsec[i], ideg[i], imn[i], xsc[i], objectnew[distant[i]].rmag
;   printf, outlun2, 'circle( ', objectnew[distant[i]].ra, objectnew[distant[i]].dec, ' 3")'
;   priority = priority + 1

;endfor

