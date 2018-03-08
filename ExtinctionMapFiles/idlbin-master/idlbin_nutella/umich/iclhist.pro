PRO iclhist

close,/all
mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A3888/final/iclhist.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color


fits_read, "/n/Godiva1/jkrick/A3888/final/centersmoV2.fit", data, header
data = data - 2.0
data = data[800:2300,1000:2400]
plothist, data, xhist, yhist, bin=0.001, xrange=[-0.05,0.05], /ylog, charthick=3, thick = 3,$
 xtitle = "counts/s ", ytitle = "number of pixels",yrange=[1E-0,1E6]


err = dindgen(1000) - dindgen(1000) + 1
start = [0.0,0.05,1E6]
result = MPFITFUN('gauss',xhist[930:1070],yhist[930:1070], err, start)

oplot, xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.)
oplot, xhist, 1E3/sqrt(2.*!Pi) * exp(-0.5*((xhist - 0.0)/0.007)^2.)

phist, data,-0.01,0.01,0.0001,/cumul
oplot, findgen(3) - 1, fltarr(3) + 0.5 


;phist, data,-0.06,0.06,0.0001,/abs, yrange=[0,1E4]
;oplot, xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.)
device, /close
set_plot, mydevice


END
