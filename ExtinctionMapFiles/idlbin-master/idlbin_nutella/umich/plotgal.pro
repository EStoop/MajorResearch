PRO plotgal

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;feb 2003
;makes plots of galaxy profiles given data files with radius and
;counts
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;device, true=24
;device, decomposed=0
colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/pele/jkrick/ind/proftest/idlgal.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

OPENR, lun,'/n/pele/jkrick/ind/proftest/proflist', /GET_LUN
rows = 28
names = strarr(28)
onename = strarr(1)
FOR j=0,rows-1 DO BEGIN
    READF, lun, onename
    names(j) = onename
ENDFOR

;make plot so I can just use oplot in the forloop
xjunk = [0,1]
yjunk = [0,0]
plot, xjunk, yjunk,psym = 3,xrange = [0,20],$
  ytitle = 'normalized counts', title = 'Galaxy profile fun', $
  xtitle = 'radius in pixels'

;first four plots
FOR j=0,3 DO BEGIN
    openr, lun1, names(j), /GET_LUN
    rows= 59
    radius = FLTARR(rows)
    counts = FLTARR(rows)

    r = 0.0                    ;radius
    c = 0.0                    ;mean counts

    FOR k=0,rows-1 DO BEGIN
        READF, lun1, r,c
        radius(k) = r
        counts(k) = c
    ENDFOR

    oplot, radius, counts, thick = 3, color = colors.blue
    xyouts, 10, 0.9, 'B/T ranges from 0.1 - 1,'
    xyouts, 10, 0.85, 'for large re, small rd'
    close, lun1
    free_lun, lun1

ENDFOR

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;make plot so I can just use oplot in the forloop
xjunk = [0,1]
yjunk = [0,0]
plot, xjunk, yjunk,psym = 3,xrange = [0,20],$
  ytitle = 'normalized counts', title = 'Galaxy profile fun', $
  xtitle = 'radius in pixels'

;first four plots
FOR j=4,7 DO BEGIN
    openr, lun1, names(j), /GET_LUN
    rows= 59
    radius = FLTARR(rows)
    counts = FLTARR(rows)

    r = 0.0                    ;radius
    c = 0.0                    ;mean counts

    FOR k=0,rows-1 DO BEGIN
        READF, lun1, r,c
        radius(k) = r
        counts(k) = c
    ENDFOR

    oplot, radius, counts, thick = 3, color = colors.red
    xyouts, 10, 0.9, 'B/T ranges from 0.1 - 1,'
    xyouts, 10, 0.85, 'for small re, large rd'
    close, lun1
    free_lun, lun1

ENDFOR

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;make plot so I can just use oplot in the forloop
xjunk = [0,1]
yjunk = [0,0]
plot, xjunk, yjunk,psym = 3,xrange = [0,20],$
  ytitle = 'normalized counts', title = 'Galaxy profile fun', $
  xtitle = 'radius in pixels'

;first four plots
FOR j=8,11 DO BEGIN
    openr, lun1, names(j), /GET_LUN
    rows= 59
    radius = FLTARR(rows)
    counts = FLTARR(rows)

    r = 0.0                    ;radius
    c = 0.0                    ;mean counts

    FOR k=0,rows-1 DO BEGIN
        READF, lun1, r,c
        radius(k) = r
        counts(k) = c
    ENDFOR

    oplot, radius, counts, thick = 3, color = colors.green
    xyouts, 10, 0.9, 'same galaxy'
    xyouts, 10, 0.85, 'four different gim2d fits'
    close, lun1
    free_lun, lun1

ENDFOR

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;make plot so I can just use oplot in the forloop
xjunk = [0,1]
yjunk = [0,0]
plot, xjunk, yjunk,psym = 3,xrange = [0,20],$
  ytitle = 'normalized counts', title = 'Galaxy profile fun', $
  xtitle = 'radius in pixels'

;first four plots
FOR j=12,15 DO BEGIN
    openr, lun1, names(j), /GET_LUN
    rows= 59
    radius = FLTARR(rows)
    counts = FLTARR(rows)

    r = 0.0                    ;radius
    c = 0.0                    ;mean counts

    FOR k=0,rows-1 DO BEGIN
        READF, lun1, r,c
        radius(k) = r
        counts(k) = c
    ENDFOR

    oplot, radius, counts, thick = 3, color = colors.orange
    xyouts, 10, 0.9, 'rd varies'
    xyouts, 10, 0.85, 'everything else constant'
    close, lun1
    free_lun, lun1

ENDFOR
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;make plot so I can just use oplot in the forloop
xjunk = [0,1]
yjunk = [0,0]
plot, xjunk, yjunk,psym = 3,xrange = [0,20],$
  ytitle = 'normalized counts', title = 'Galaxy profile fun', $
  xtitle = 'radius in pixels'

;first four plots
FOR j=16,19 DO BEGIN
    openr, lun1, names(j), /GET_LUN
    rows= 59
    radius = FLTARR(rows)
    counts = FLTARR(rows)

    r = 0.0                    ;radius
    c = 0.0                    ;mean counts

    FOR k=0,rows-1 DO BEGIN
        READF, lun1, r,c
        radius(k) = r
        counts(k) = c
    ENDFOR

    oplot, radius, counts, thick = 3, color = colors.navy
    xyouts, 10, 0.9, 're varies'
    xyouts, 10, 0.85, 'everything else constant'
    close, lun1
    free_lun, lun1

ENDFOR
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;make plot so I can just use oplot in the forloop
xjunk = [0,1]
yjunk = [0,0]
plot, xjunk, yjunk,psym = 3,xrange = [0,40],$
;  yrange=[0,0.2],$
  ytitle = 'normalized counts', title = 'Galaxy profile fun', $
  xtitle = 'radius in pixels'

;first four plots
FOR j=20,23 DO BEGIN
    openr, lun1, names(j), /GET_LUN
    rows= 59
    radius = FLTARR(rows)
    counts = FLTARR(rows)

    r = 0.0                    ;radius
    c = 0.0                    ;mean counts

    FOR k=0,rows-1 DO BEGIN
        READF, lun1, r,c
        radius(k) = r
        counts(k) = c
    ENDFOR

    oplot, radius, counts, thick = 3, color = colors.pink
    xyouts, 10, 0.9, 'db varies 0.1,0.01,0.001,0.0'
    xyouts, 10, 0.85, 'everything else constant'
    close, lun1
    free_lun, lun1

ENDFOR
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;make plot so I can just use oplot in the forloop
xjunk = [0,1]
yjunk = [0,0]
plot, xjunk, yjunk,psym = 3,xrange = [0,40],$
;  yrange=[0,0.2],$
  ytitle = 'normalized counts', title = 'Galaxy profile fun', $
  xtitle = 'radius in pixels'

;first four plots
FOR j=24,27 DO BEGIN
    openr, lun1, names(j), /GET_LUN
    rows= 59
    radius = FLTARR(rows)
    counts = FLTARR(rows)

    r = 0.0                    ;radius
    c = 0.0                    ;mean counts

    FOR k=0,rows-1 DO BEGIN
        READF, lun1, r,c
        radius(k) = r
        counts(k) = c
    ENDFOR

    oplot, radius, counts, thick = 3, color = colors.green
    xyouts, 10, 0.9, 'lt varies 150,170,190,210'
    xyouts, 10, 0.85, 'everything else constant  *normalixzed*'
    close, lun1
    free_lun, lun1

ENDFOR


device, /close
set_plot, mydevice

END

