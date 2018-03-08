PRO backchecksmooth

close, /all
mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'
colors = GetColor(/Load, Start=1)

device, filename = '/n/Godiva1/jkrick/A3888/background/backchecksmooth.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

openr, rlun, "/n/Godiva1/jkrick/A3888/backlistchecksmoothr", /GET_LUN
openr, Vlun, "/n/Godiva1/jkrick/A3888/backlistchecksmoothV", /GET_LUN

i = 0
npix = 0
aver = 0.
averarr = fltarr(210000)
WHILE (NOT eof(rlun)) DO BEGIN
    READF, rlun, npix, aver
    IF npix EQ 16 AND aver GT -100. THEN BEGIN
        averarr(i) = aver -2.
        i = i +1 
    ENDIF
ENDWHILE

averarr = averarr[0:i-1]
print, "we have ",i," measures of the mean"

plothist, averarr, xhist, yhist, thick = 3,  charthick = 3, xthick = 3, ythick = 3,$
          xtitle = "mean counts/s/pix", ytitle = "number", title = "r histogram" 
yhist = double(yhist)
err = fltarr(N_ELEMENTS(xhist)) + 1.
start = [0.001,0.001, 1000]
result = MPFITFUN('gauss',xhist,yhist, err, start) ;, weights=weight);,PARINFO=pi);, weights=weight)
print, "result", result
;---------------------------------------------------------------------------------
i = 0
npix = 0
aver = 0.
averarr = fltarr(210000)
WHILE (NOT eof(Vlun)) DO BEGIN
    READF, Vlun, npix, aver
    IF npix EQ 16 AND aver GT -100. THEN BEGIN
        averarr(i) = aver -2.
        i = i +1 
    ENDIF
ENDWHILE

averarr = averarr[0:i-1]
print, "we have ",i," measures of the mean"

plothist, averarr, xhist, yhist, thick = 3,  charthick = 3, xthick = 3, ythick = 3,$
          xtitle = "mean counts/s/pix", ytitle = "number", title = "V histogram" 
yhist = double(yhist)
err = fltarr(N_ELEMENTS(xhist)) + 1.
start = [0.001,0.001, 1000]
result = MPFITFUN('gauss',xhist,yhist, err, start) ;, weights=weight);,PARINFO=pi);, weights=weight)
print, "result", result

END
