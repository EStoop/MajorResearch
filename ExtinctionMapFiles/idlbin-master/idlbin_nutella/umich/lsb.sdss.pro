PRO lsb


OPENR, lun,'/n/Godiva1/jkrick/A3888/final/lsb.sdss', /GET_LUN  ;galVtab2

;read in the radial profile
color = fltarr(200)
i = 0
WHILE (NOT eof (lun)) DO BEGIN
    readf, lun, c
    color(i) = c
    i = i + 1
ENDWHILE



close, lun
free_lun, lun

color = color[0:i-1]
print, color

plothist, color, xhist, yhist, bin = 0.05
;plot, xhist, yhist


END
