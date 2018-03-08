PRO sbhist

close,/all
;device, true=24
;device, decomposed=0
mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A3888/final/sbhist.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

fits_read, "/n/Godiva1/jkrick/A3888/final/fullV2.fits", data, header
;fits_read, "/n/Godiva2/jkrick/A141/original/largeV.fits", data, header
datarr = fltarr(4100000.)
j = float(0.)
FOR x = 800, 2400, 1 DO BEGIN   ; 800 2400
    FOR y = 1200, 2400, 1 DO BEGIN   ;1200 2400
;        IF data[x,y] GT -1 THEN begin 
            datarr[j] = data[x,y] - 2.
            j = j + 1
;        ENDIF

    ENDFOR
ENDFOR

datarr = datarr[0:j-1]
;datarr = 24.3 - 2.5*alog10(datarr/(0.259^2.))

;-------------------------------------------
;make and plot a histogram of data values
;-------------------------------------------
plothist,datarr, xarr, yarr,bin = 0.0001, /noplot

plot, xarr, (yarr)*(0.259^2.), xrange = [-0.01,0.2], thick = 3,$
  charthick = 3, xthick = 3, ythick = 3, xtitle = "couts/pix", $
  ytitle = "area (square arcseconds)", /ylog, yrange=[1E-2,1E4], ystyle = 1


err = dindgen(j) - dindgen(j) + 1
start = [0.001,0.01, 2.E5]
result= MPFITFUN('gauss',xarr,yarr, err, start)   

gaussy = ((result(2)/sqrt(2.*!Pi) * $
  exp(-0.5*((xarr - result(0))/result(1))^2.)))

gaussy2 = result(2)/sqrt(2.*!Pi) * $
  exp(-0.5*((xarr - 0.000299)/0.000666)^2.)

oplot, xarr, (gaussy)*(0.259^2.), thick = 3, linestyle = 2
;oplot, xarr, (yarr - gaussy)*(0.259^2.), thick = 3, linestyle = 3

;-----------------------------------
;want a cumulative distribution
;-----------------------------------
i = 0
sum = fltarr(n_elements(yarr))  ; all flux
sum(0) = 0
sumg = fltarr(n_elements(gaussy))   ;noise (gaussian fit)
sumg(0) = 0

;FOR i = 1, n_elements(yarr) - 1, 1 DO BEGIN
;    sum(i) = sum(i-1) +yarr[i] ;- gaussy2[i]
;    sumg(i) = sumg(i-1) +gaussy2[i]
;ENDFOR

;sum = sum / sum(i-1)    ;normalize
;sumg = sumg / sumg(i-1)    ;normalize

;plot, xarr, sum, thick = 3, yrange = [-.1,1],xrange = [-0.001,0.01]
;oplot, xarr, sumg, thick = 3, linestyle = 2
;oplot, xarr, fltarr(n_elements(xarr))

;oplot, xarr, sum - sumg

;xyouts,0.001, -0.15,'28.9', alignment=0.5
;xyouts,0.003, -0.15,'27.7', alignment=0.5
;xyouts,0.005, -0.15,'27.1', alignment=0.5

;oplot, findgen(10) - findgen(10) + 0.001, findgen(10)
;oplot, findgen(10) - findgen(10) + 0.003, findgen(10)
;oplot, findgen(10) - findgen(10) + 0.005, findgen(10)


device, /close
set_plot, mydevice

END






;plot yarr - gaussian distribution
 ;plot, xarr, yarr -( result(2)/sqrt(2.*!Pi) * $
 ; exp(-0.5*((xarr - result(0))/result(1))^2.)), thick = 3, $
;  linestyle = 3,xrange = [-0.01,0.01], yrange = [1,1E4],/ylog

;oplot, xarr, yarr -( result(2)/sqrt(2.*!Pi) * $
;  exp(-0.5*((xarr - 0.000299)/0.000666)^2.)), thick = 3

;phist, datarr,/cumul, 0.00, 200., 0.0005, xrange=[0.0,0.01],yrange=[0.7,1.0]  ;max(yarr)
;phist, gaussy, /cumul, 0.00, 200., 0.0005,/overplot

