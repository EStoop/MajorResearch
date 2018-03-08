PRO ffcheckiraf
close, /all

;device, true=24
;device, decomposed=0

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'
colors = GetColor(/Load, Start=1)

device, filename = '/n/Godiva2/jkrick/A141/ffchecksmoothr.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

sigma = fltarr(15)
slength = fltarr(15)
;slength = [39,49,59,69,79,89,109,129,149,169,209,249,289,329,369]

j = 0
;--------------------------------------------------------------------------------------------
;openr, lun, "/n/Godiva1/jkrick/A3888/backlistVlarge", /GET_LUN
i = 0
av = fltarr(5000)
npix = 0
aver = 0.
npixarr = fltarr(20000)
averarr = fltarr(20000)
;FOR j=0,14,1 DO BEGIN
filename = "/n/Godiva2/jkrick/A141/med75.out"
;    filename = "/n/Godiva1/jkrick/A3888/background/backlist"+ string(slength(j))+".large.rout"
;    filename = STRCOMPRESS(filename, /REMOVE_ALL)
;    print, filename
    openr, lun, filename, /GET_LUN

    WHILE (NOT eof(lun)) DO BEGIN
        READF, lun, npix, aver
        IF npix GT 0 AND aver GT -100. THEN BEGIN
;        av(i) = aver - 2.
;        print, av(i)
            npixarr(i) = npix
            averarr(i) = aver -2.
            i = i +1 
        ENDIF
    ENDWHILE
;av = av[0:i-1]
    
    close, lun
    free_lun, lun
    
;ENDFOR
npixarr = npixarr[0:i-1]
averarr = averarr[0:i-1]

print, "i", i

;sort them
sortindex= sort(npixarr)
sortednpix = npixarr[sortindex]
sortedaver = averarr[sortindex]

;bin them, use histogram to reverse engineer the bins
plothist, sortednpix, xhist, yhist, bin = 200, xtitle = "npix", ytitle = "number of boxes"

;make histograms to get at sigma
sigarr = fltarr(N_Elements(yhist))
slength = fltarr(N_Elements(yhist))
;plothist, sortedaver[0:yhist(0)],bin = 0.001, x2,y2,/noplot
;y2 = double(y2)
err = fltarr(N_ELEMENTS(xhist)) + 1.
;result = MPFITFUN('gauss',x2,y2, err, start);, weights=weight);,PARINFO=pi);, weights=weight)
;sigarr(0) = result(1)
;startsum = yhist(0)

;openw, outlun, "/n/Godiva1/jkrick/A3888/background/ffcheckiraf.output", /get_lun
startsum = 0
FOR q = 0, 120, 1 DO BEGIN
    IF yhist(q) GT 20 THEN begin
;        printf, outlun, q, yhist(q),startsum, startsum + yhist(q)
        plothist, sortedaver[startsum:(startsum + yhist(q))],bin = 0.001, x2,y2,/noplot
        y2 = double(y2)
        start = [0.001,0.001, 2*(yhist(q))]
        result = MPFITFUN('gauss',x2,y2, err, start) ;, weights=weight);,PARINFO=pi);, weights=weight)
;        
        sigarr(q) = result(1)
    ENDIF ELSE BEGIN
        sigarr(q) = sigarr(q-1)
    ENDELSE

    startsum = startsum + yhist(q)

    IF q NE N_ElEMENTS(yhist) - 1 THEN BEGIN
        slength(q) = mean(sortednpix[startsum:(startsum + yhist(q))])
        
    ENDIF ELSE BEGIN
        slength(q) = sortednpix[startsum]
;        printf,outlun, "shouldn't be here"
    ENDELSE;

ENDFOR
;close, outlun
;free_lun, outlun

sigarr = abs(sigarr[0:120])
slength = slength[0:120]
print, "sigarr", sigarr
print, "slength", slength

plot, (sqrt(slength))*0.259, 24.3 -2.5 * (alog10(sigarr/(0.259^2))), yrange = [28,30],/xlog, thick = 3, $
  charthick = 3, xthick = 3, ythick = 3, xtitle = "smoothing length(arcseconds)", $
ytitle = "stddev(mag/arcsec^2)"
oplot, (sqrt(slength))*0.259, 24.3 - 2.5*alog10(1/(22* (sqrt(slength))* (0.259^2))), linestyle = 2

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
