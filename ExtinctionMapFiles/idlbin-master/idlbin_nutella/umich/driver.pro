PRO driver
close, /all


i = 0
mag = fltarr(123)
openr, lun,"/n/Godiva4/jkrick/A114/couch01.gal", /get_lun
WHILE (NOT(EOF(lun)) )DO BEGIN
    readf,lun, junk,junk,junk,junk,junk,junk, m
    mag(i) = m
     i = i + 1
ENDWHILE
mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps';

device, filename = '/n/Godiva4/jkrick/A114/redshiftdist.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color;

plothist, mag, bin = 0.01, thick = 3
device, /close
set_plot, mydevice

close, lun
free_lun, lun

END
