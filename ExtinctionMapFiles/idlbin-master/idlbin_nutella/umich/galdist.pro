PRO galdist
close, /all
magisocor = fltarr(10000)
i = 0
openr, lun, "/n/Godiva1/jkrick/A3888/final/SExtractor2.V.cat", /get_lun
WHILE (NOT EOF(lun)) DO BEGIN
   
    READF, lun, o, xcenter, ycenter,junk,junk,junk,junk,junk,magiso,mag,junk,junk,junk,flux_isocor,flux_iso
    magisocor[i] = 24.6 -2.5*(alog10(flux_isocor))
    i = i +1
    
ENDWHILE
magisocor = magisocor[0:i-1]

close,lun
free_lun, lun

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A3888/final/galhist.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color


plothist, magisocor, xhist, yhist, /noplot, bin = 0.25
plot, xhist, yhist, psym = 4, thick = 3

device, /close
set_plot, mydevice


END
