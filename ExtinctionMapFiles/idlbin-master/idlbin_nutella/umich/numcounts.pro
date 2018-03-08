PRO numcounts

close, /all
vmag = fltarr(10000)
openr, lun2, "/n/Godiva1/jkrick/A3888/cmd/SExtractor.V2.cat", /GET_LUN
i = 0
WHILE (NOT eof(lun2)) DO BEGIN
    READF, lun2,junk, vx, vy, vfluxaper,vmagaper,verraper,vfluxbest,vmagbest,verrbest,vfwhm, vra,vdec,viso    
    
    IF (vfwhm GT 4.0) AND (vfluxbest GT 0.0) AND (viso GT 6.0)THEN BEGIN ;not star or junk
       vmag(i) = vmagbest
        i = i + 1
    ENDIF ELSE BEGIN
        print, junk, vx, vy, vfluxaper,vmagaper,verraper,vfluxbest,vmagbest,verrbest,vfwhm, vra,vdec,viso  
    ENDELSE

ENDWHILE

Vmag = Vmag[0:i - 1]
print, i

sortindex = sort(Vmag)
sortedmag = Vmag[sortindex]

;device, true=24
;device, decomposed=0

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A3888/cmd/numcounts.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plothist, sortedmag, xhist, yhist, bin = 0.5, xtitle = "Vmag", ytitle = "number of galaxies",$
title = "A3888 cluster + background galaxies , z = 0.15, richness class 2", thick = 3, $
charthick = 3, xthick = 3, ythick = 3

print, xhist, yhist
;plot, xhist, yhist

device, /close
set_plot, mydevice

END
