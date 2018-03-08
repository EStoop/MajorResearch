PRO galaxev
close, /all		;close all files = tidiness

device, true=24
device, decomposed=0
colors = GetColor(/load, Start=1)

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/Users/jkrick/umich/icl/color/galaxev.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plot,findgen(10) - findgen(10),findgen(10) - findgen(10), charthick = 3, xtitle = "z",ytitle = "V-r",$
  xthick = 3, ythick = 3, xrange = [0,0.4], xstyle = 1,color = colors.black,$
  yrange = [0,1.3], ystyle = 1, thick = 3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

filelist11 = ['/Users/jkrick/umich/icl/color/m52tf11.5.color',$
              '/Users/jkrick/umich/icl/color/m42tf11.5.color'];,$

FOR n = 0, 1, 1 DO begin
    OPENR, lun,filelist11[n], /GET_LUN

;read in the radial profile
    rows= 108
    zarr= FLTARR(rows)
    agearr= FLTARR(rows)
    cevarr= FLTARR(rows)
    crfarr= FLTARR(rows)
    cnearr = FLTARR(rows)
    karr= FLTARR(rows)
    ekarr= FLTARR(rows)
    
    z = 0.0			
    age = 0.0
    cev = 0.0
    junk = 0.0
    
    j = 0
    FOR j=0,rows - 1 DO BEGIN
        READF, lun, z,junk,age,junk,crf,cne,cev,junk,junk,junk,ek,k
        zarr(j) = z
        agearr(j) = age
        cevarr(j) = cev - 0.3
        cnearr(j) = cne
        crfarr(j) = crf
        karr(j) = k
        ekarr(j) = ek
    ENDFOR
    
    close, lun
    free_lun, lun

    oplot,zarr,cevarr,thick = 3, color = colors.black, linestyle = n

ENDFOR



xyouts, 0.01,0.95,"zf=3.1", charthick = 2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

filelist10 = ['/Users/jkrick/umich/icl/color/m52tf10.color',$
              '/Users/jkrick/umich/icl/color/m42tf10.color'];,$

FOR n = 0, 1, 1 DO begin
    OPENR, lun,filelist10[n], /GET_LUN

;read in the radial profile
    rows= 90
    zarr= FLTARR(rows)
    agearr= FLTARR(rows)
    cevarr= FLTARR(rows)
    crfarr= FLTARR(rows)
    cnearr = FLTARR(rows)
    karr= FLTARR(rows)
    ekarr= FLTARR(rows)
    
    z = 0.0			
    age = 0.0
    cev = 0.0
    junk = 0.0
    
    j = 0
    FOR j=0,rows - 1 DO BEGIN
        READF, lun, z,junk,age,junk,crf,cne,cev,junk,junk,junk,ek,k
        zarr(j) = z
        agearr(j) = age
        cevarr(j) = cev - 0.3
        cnearr(j) = cne
        crfarr(j) = crf
        karr(j) = k
        ekarr(j) = ek
    ENDFOR
    
    close, lun
    free_lun, lun
    
    
    oplot,zarr,cevarr,color = colors.blue, thick = 3, linestyle = n

ENDFOR

xyouts, 0.01, 0.90, "zf=1.9", charthick = 2, color = colors.blue


;------------------------------------------------------------------------
filelist11 = ['/Users/jkrick/umich/icl/color/constm52tf11.5.color',$
              '/Users/jkrick/umich/icl/color/constm42tf11.5.color'];,$

FOR n = 0, 1, 1 DO begin
    OPENR, lun,filelist11[n], /GET_LUN

;read in the radial profile
    rows= 108
    zarr= FLTARR(rows)
    agearr= FLTARR(rows)
    cevarr= FLTARR(rows)
    crfarr= FLTARR(rows)
    cnearr = FLTARR(rows)
    karr= FLTARR(rows)
    ekarr= FLTARR(rows)
    
    z = 0.0			
    age = 0.0
    cev = 0.0
    junk = 0.0
    
    j = 0
    FOR j=0,rows - 1 DO BEGIN
        READF, lun, z,junk,age,junk,crf,cne,cev,junk,junk,junk,ek,k
        zarr(j) = z
        agearr(j) = age
        cevarr(j) = cev - 0.3
        cnearr(j) = cne
        crfarr(j) = crf
        karr(j) = k
        ekarr(j) = ek
    ENDFOR
    
    close, lun
    free_lun, lun

    oplot,zarr,cevarr,thick = 3, color = colors.gray, linestyle = n

ENDFOR




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

filelist10 = ['/Users/jkrick/umich/icl/color/constm52tf10.color',$
              '/Users/jkrick/umich/icl/color/constm42tf10.color'];,$

FOR n = 0, 1, 1 DO begin
    OPENR, lun,filelist10[n], /GET_LUN

;read in the radial profile
    rows= 90
    zarr= FLTARR(rows)
    agearr= FLTARR(rows)
    cevarr= FLTARR(rows)
    crfarr= FLTARR(rows)
    cnearr = FLTARR(rows)
    karr= FLTARR(rows)
    ekarr= FLTARR(rows)
    
    z = 0.0			
    age = 0.0
    cev = 0.0
    junk = 0.0
    
    j = 0
    FOR j=0,rows - 1 DO BEGIN
        READF, lun, z,junk,age,junk,crf,cne,cev,junk,junk,junk,ek,k
        zarr(j) = z
        agearr(j) = age
        cevarr(j) = cev - 0.3
        cnearr(j) = cne
        crfarr(j) = crf
        karr(j) = k
        ekarr(j) = ek
    ENDFOR
    
    close, lun
    free_lun, lun
    
    
    oplot,zarr,cevarr,color = colors.powderblue, thick = 3, linestyle = n

ENDFOR


;------------------------------------------------------------------------
filelist11 = ['/Users/jkrick/umich/icl/color/expm52tf11.5.color',$
              '/Users/jkrick/umich/icl/color/expm42tf11.5.color'];,$

FOR n = 0, 1, 1 DO begin
    OPENR, lun,filelist11[n], /GET_LUN

;read in the radial profile
    rows= 108
    zarr= FLTARR(rows)
    agearr= FLTARR(rows)
    cevarr= FLTARR(rows)
    crfarr= FLTARR(rows)
    cnearr = FLTARR(rows)
    karr= FLTARR(rows)
    ekarr= FLTARR(rows)
    
    z = 0.0			
    age = 0.0
    cev = 0.0
    junk = 0.0
    
    j = 0
    FOR j=0,rows - 1 DO BEGIN
        READF, lun, z,junk,age,junk,crf,cne,cev,junk,junk,junk,ek,k
        zarr(j) = z
        agearr(j) = age
        cevarr(j) = cev - 0.3
        cnearr(j) = cne
        crfarr(j) = crf
        karr(j) = k
        ekarr(j) = ek
    ENDFOR
    
    close, lun
    free_lun, lun

;    oplot,zarr,cevarr,thick = 3, color = colors.forestgreen, linestyle = n

ENDFOR




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

filelist10 = ['/Users/jkrick/umich/icl/color/expm52tf10.color',$
              '/Users/jkrick/umich/icl/color/expm42tf10.color'];,$

FOR n = 0, 1, 1 DO begin
    OPENR, lun,filelist10[n], /GET_LUN

;read in the radial profile
    rows= 90
    zarr= FLTARR(rows)
    agearr= FLTARR(rows)
    cevarr= FLTARR(rows)
    crfarr= FLTARR(rows)
    cnearr = FLTARR(rows)
    karr= FLTARR(rows)
    ekarr= FLTARR(rows)
    
    z = 0.0			
    age = 0.0
    cev = 0.0
    junk = 0.0
    
    j = 0
    FOR j=0,rows - 1 DO BEGIN
        READF, lun, z,junk,age,junk,crf,cne,cev,junk,junk,junk,ek,k
        zarr(j) = z
        agearr(j) = age
        cevarr(j) = cev - 0.3
        cnearr(j) = cne
        crfarr(j) = crf
        karr(j) = k
        ekarr(j) = ek
    ENDFOR
    
    close, lun
    free_lun, lun
    
    
    oplot,zarr,cevarr,color = colors.magenta, thick = 3, linestyle = n

ENDFOR

;-----------------------------------------------------------------------
filelist11 = ['/Users/jkrick/umich/icl/color/dur100m52tf11.5.color',$
              '/Users/jkrick/umich/icl/color/dur100m42tf11.5.color'];,$

FOR n = 0, 1, 1 DO begin
    OPENR, lun,filelist11[n], /GET_LUN

;read in the radial profile
    rows= 108
    zarr= FLTARR(rows)
    agearr= FLTARR(rows)
    cevarr= FLTARR(rows)
    crfarr= FLTARR(rows)
    cnearr = FLTARR(rows)
    karr= FLTARR(rows)
    ekarr= FLTARR(rows)
    
    z = 0.0			
    age = 0.0
    cev = 0.0
    junk = 0.0
    
    j = 0
    FOR j=0,rows - 1 DO BEGIN
        READF, lun, z,junk,age,junk,crf,cne,cev,junk,junk,junk,ek,k
        zarr(j) = z
        agearr(j) = age
        cevarr(j) = cev - 0.3
        cnearr(j) = cne
        crfarr(j) = crf
        karr(j) = k
        ekarr(j) = ek
    ENDFOR
    
    close, lun
    free_lun, lun

;    oplot,zarr,cevarr,thick = 3, color = colors.cyan, linestyle = n

ENDFOR




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

filelist10 = ['/Users/jkrick/umich/icl/color/dur100m52tf10.color',$
              '/Users/jkrick/umich/icl/color/dur100m42tf10.color'];,$

FOR n = 0, 1, 1 DO begin
    OPENR, lun,filelist10[n], /GET_LUN

;read in the radial profile
    rows= 90
    zarr= FLTARR(rows)
    agearr= FLTARR(rows)
    cevarr= FLTARR(rows)
    crfarr= FLTARR(rows)
    cnearr = FLTARR(rows)
    karr= FLTARR(rows)
    ekarr= FLTARR(rows)
    
    z = 0.0			
    age = 0.0
    cev = 0.0
    junk = 0.0
    
    j = 0
    FOR j=0,rows - 1 DO BEGIN
        READF, lun, z,junk,age,junk,crf,cne,cev,junk,junk,junk,ek,k
        zarr(j) = z
        agearr(j) = age
        cevarr(j) = cev - 0.3
        cnearr(j) = cne
        crfarr(j) = crf
        karr(j) = k
        ekarr(j) = ek
    ENDFOR
    
    close, lun
    free_lun, lun
    
    
 ;   oplot,zarr,cevarr,color = colors.red, thick = 3, linestyle = n

ENDFOR
;-----------------------------------------------------------------------
;filelist11 = ['/Users/jkrick/umich/icl/color/dur10m52tf11.5.color',$
;              '/Users/jkrick/umich/icl/color/dur10m42tf11.5.color'];,$
;
;FOR n = 0, 1, 1 DO begin
;    OPENR, lun,filelist11[n], /GET_LUN
;
;;read in the radial profile
;    rows= 108
;    zarr= FLTARR(rows)
;    agearr= FLTARR(rows)
;    cevarr= FLTARR(rows)
;    crfarr= FLTARR(rows)
;    cnearr = FLTARR(rows)
;    karr= FLTARR(rows)
;    ekarr= FLTARR(rows)
;    
;    z = 0.0			
;    age = 0.0
;    cev = 0.0
;    junk = 0.0
;    
;    j = 0
;    FOR j=0,rows - 1 DO BEGIN
;        READF, lun, z,junk,age,junk,crf,cne,cev,junk,junk,junk,ek,k
;        zarr(j) = z
;        agearr(j) = age
;        cevarr(j) = cev - 0.3
;        cnearr(j) = cne
;        crfarr(j) = crf
;        karr(j) = k
;        ekarr(j) = ek
;    ENDFOR
;    
;    close, lun
;    free_lun, lun
;
;    oplot,zarr,cevarr,thick = 3, color = colors.lawngreen, linestyle = n
;
;ENDFOR
;
;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;filelist10 = ['/Users/jkrick/umich/icl/color/dur10m52tf10.color',$
;              '/Users/jkrick/umich/icl/color/dur10m42tf10.color'];,$
;
;FOR n = 0, 1, 1 DO begin
;    OPENR, lun,filelist10[n], /GET_LUN
;
;;read in the radial profile
;    rows= 90
;    zarr= FLTARR(rows)
;    agearr= FLTARR(rows)
;    cevarr= FLTARR(rows)
;    crfarr= FLTARR(rows)
;    cnearr = FLTARR(rows)
;    karr= FLTARR(rows)
;    ekarr= FLTARR(rows)
;    
;    z = 0.0			
;    age = 0.0
;    cev = 0.0
;    junk = 0.0
;    
;    j = 0
;    FOR j=0,rows - 1 DO BEGIN
;        READF, lun, z,junk,age,junk,crf,cne,cev,junk,junk,junk,ek,k
;        zarr(j) = z
;        agearr(j) = age
;        cevarr(j) = cev - 0.3
;        cnearr(j) = cne
;        crfarr(j) = crf
;        karr(j) = k
;        ekarr(j) = ek
;    ENDFOR
;    
;    close, lun
;    free_lun, lun
;    
;    
;    oplot,zarr,cevarr,color = colors.orange, thick = 3, linestyle = n
;
;ENDFOR

;------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

filelist10 = ['/Users/jkrick/umich/icl/color/m52tf13.color',$
              '/Users/jkrick/umich/icl/color/m42tf13.color'];,$

FOR n = 0, 1, 1 DO begin
    OPENR, lun,filelist10[n], /GET_LUN

;read in the radial profile
    rows= 136
    zarr= FLTARR(rows)
    agearr= FLTARR(rows)
    cevarr= FLTARR(rows)
    crfarr= FLTARR(rows)
    cnearr = FLTARR(rows)
    karr= FLTARR(rows)
    ekarr= FLTARR(rows)
    
    z = 0.0			
    age = 0.0
    cev = 0.0
    junk = 0.0
    
    j = 0
    FOR j=0,rows - 1 DO BEGIN
        READF, lun, z,junk,age,junk,crf,cne,cev,junk,junk,junk,ek,k
        zarr(j) = z
        agearr(j) = age
        cevarr(j) = cev - 0.3
        cnearr(j) = cne
        crfarr(j) = crf
        karr(j) = k
        ekarr(j) = ek
    ENDFOR
    
    close, lun
    free_lun, lun
    
    
;   oplot,zarr,cevarr,color = colors.black, thick = 3, linestyle = n

ENDFOR
;-----------------------------------------------------------

;oplot, findgen(100) - findgen(100) + 0.15, findgen(3)
;oplot, findgen(100) - findgen(100) + 0.18, findgen(3)
;oplot, findgen(100) - findgen(100) + 0.23, findgen(3)
;oplot, findgen(100) - findgen(100) + 0.31, findgen(3)
;oplot, findgen(100) - findgen(100) + 0.048, findgen(3), linestyle = 2
;oplot, findgen(100) - findgen(100) + 0.058, findgen(3), linestyle = 2
;oplot, findgen(100) - findgen(100) + 0.062, findgen(3), linestyle = 2
;oplot, findgen(100) - findgen(100) + 0.087, findgen(3), linestyle = 2
;oplot, findgen(100) - findgen(100) + 0.096, findgen(3), linestyle = 2

xcolor = [0.05,0.06,0.06,0.09,0.1,0.15,0.18,0.23,0.31,0.31]
ycolor = [0.58,1.25,1.15,0.86,0.82,0.6,0.0,1.1,0.6 , 0.8]

oplot, xcolor, ycolor, thick = 3, psym = 2
lower = [0.08,0.11,0.03,0.4,0.1,0.2,0.1,0.2,0.2,0.1]
upper = [0.08,0.11,0.03,0.4,0.1,0.2,0.6,0.3,0.2,0.1]
errplot, xcolor,ycolor-lower, ycolor+upper

device, /close
set_plot, mydevice

END

