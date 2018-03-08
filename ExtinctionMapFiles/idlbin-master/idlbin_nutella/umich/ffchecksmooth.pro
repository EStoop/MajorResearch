PRO ffchecksmooth
close, /all

;device, true=24
;device, decomposed=0

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'
colors = GetColor(/Load, Start=1)

device, filename = '/Users/jkrick/umich/icl/A3880/noiser.47.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

;sigma = fltarr(13)
;slength = fltarr(13)
;slength = [4,39,49,69,79,109,129,149,169,189,209,249,289];,329,369]

j = 0
;--------------------------------------------------------------------------------------------

;i = long(0)
;npix = 0
;aver = 0.
;averarr = fltarr(200000)
;sigarr = fltarr(13)
;FOR j=0,12,1 DO BEGIN
    i = long(0)
    npix = 0
    aver = 0.
    averarr = fltarr(211000)
;    filename = "/n/Godiva1/jkrick/A3888/background/largeV"+ string(slength(j))+".out"
;    filename = STRCOMPRESS(filename, /REMOVE_ALL)
;    print, filename
    filename = "/Users/jkrick/umich/icl/A3880/med47.out" ;fullr.79
    openr, lun, filename, /GET_LUN

    WHILE (NOT eof(lun)) DO BEGIN
        READF, lun, npix, aver
        IF npix EQ 25 AND aver GT -0.09 THEN BEGIN  ;1.9  AND aver GT -0.1
            averarr(i) = aver ;-2.
            i = i +1 
        ENDIF
    ENDWHILE
    
    close, lun
    free_lun, lun

    averarr = averarr[0:i-1]
    print, "we have ",i," measures of the mean"

    plothist, averarr, xhist, yhist, bin=0.00005,thick = 3,  charthick = 3, xthick = 3, ythick = 3,$
              xtitle = " counts/pix",$;, title = filename, $
              xrange=[-0.01,0.01];,yrange = [0,100];4E4
    yhist = double(yhist)
    err = fltarr(N_ELEMENTS(xhist)) + 1.
    start = [0.001,0.001, 200000]
    result = MPFITFUN('gauss',xhist,yhist, err, start) ;, weights=weight);,PARINFO=pi);, weights=weight)
xarr = (findgen(15000) - 700)/10000 

    oplot, xarr, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - result(0))/result(1))^2.)

;    print, "result",slength(j), "  ", result
;    sigarr(j) = abs(result(1))
;ENDFOR




;plot, slength*0.259, 24.3 - 2.5 * (alog10(sigarr/(0.259^2))), yrange = [28,30],/xlog, thick = 3, $
;  charthick = 3, xthick = 3, ythick = 3, xtitle = "smoothing length(arcseconds)", $
;      ytitle = "stddev(mag/arcsec^2)", xrange = [1,100]
;oplot, slength*0.259, 24.3 - 2.5*alog10(1/(22* (slength)* (0.259^2))), linestyle = 2

device, /close
set_plot, mydevice

END

;plothist, av, xhist, yhist, bin = 0.0001, thick = 3, xrange = [-.01,.01],ythick = 3,peak = 1,$
;charthick = 3, xthick = 3, ytitle = "N", title = "V band ff and background accuracy", xtitle = "counts/s"

;weight = fltarr(N_ELEMENTS(xhist))
;weight[0:42] = 1000.

;start = [0.001,0.001, 1.]
;start = [1.5]
;result = MPFITFUN('gauss',xhist,yhist, err, start);, weights=weight);,PARINFO=pi);, weights=weight)
;print, "gauss result", result


;xaxis = findgen(1000)/100000. 
;xaxis2 = findgen(1000)/(-100000.) 
;test = reverse(xaxis2)
;axis = [test, xaxis]

;oplot, axis, result(2)/sqrt(2.*!Pi) * exp(-0.5*((axis -result(0))/result(1))^2.), thick = 3, color = colors.red
;oplot, axis, result(0)/sqrt(2.*!Pi) * exp(-0.5*((axis )/0.0009)^2.), thick = 3, color = colors.red

;sigma(j) = result(1)

;ENDFOR

;-----------------------------------------------------------------------------------------------

;plot, slength*0.259, 24.3 -2.5 * (alog10(sigma/(0.259^2))), yrange = [28,30],/xlog
;oplot, slength*0.259, 24.3 - 2.5*alog10(1/(22* slength* (0.259^2))), linestyle = 2
;device, /close
;set_plot, mydevice
;openr, lun, "/n/Godiva1/jkrick/A3888/backstatnew115", /GET_LUN

;i = 0
;av115 = fltarr(5000)
;npix = 0
;aver = 0.
;WHILE (NOT eof(lun)) DO BEGIN
;    READF, lun, npix, aver
;    IF npix GT 2000 THEN BEGIN
;        av115(i) = aver - 2.
;;        print, av(i)
;        i = i +1 
;    ENDIF
;ENDWHILE

;av115 = av115[0:i-1]
;close, lun
;free_lun, lun


;plot, xhist, yhist;, xrange = [-0.3,0.3], charthick = 3, xthick = 3, thick = 3,$
;ythick = 3;, xtitle = "Intensity (counts/s)", title = "V-band Pixel Intensity Histograms", ystyle = 5


;plothist, av115, xhist115,yhist115, bin = 0.0001, peak = 1, /noplot

;result115 = MPFITFUN('gauss',xhist115,yhist115, err, start, weights=weight);,PARINFO=pi);, weights=weight)
;print, "gauss115 result", result115
;oplot, axis, result115(2)/sqrt(2.*!Pi) * exp(-0.5*((axis - result115(0))/result115(1))^2.), $
;thick = 3, color = colors.blue
