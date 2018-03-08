PRO pimbblet_sigma

close,/all
colors = GetColor(/load, Start=1)
mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A3888/pimbblet/dist.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

OPENR, lun, "/n/Godiva1/jkrick/A3888/pimbblet/larcs.cz.txt", /GET_LUN
czarr = dblarr(218)
czerrarr = fltarr(218)
i = 0
cz = 0.D
WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun,rmag,rmagerr,ra,dec,cz,czerr
    czarr(i) = cz
    czerrarr(i) = czerr
    i = i + 1
ENDWHILE
close, lun
free_lun, lun


plothist, czarr,bin = 0.05E4, xrange=[3.8E4,5.2E4]
mean = biweight_mean(czarr, sigma)
print, mean, sigma



OPENR, lun1, "/n/Godiva1/jkrick/A3888/cmd/teaguememcz.txt", /GET_LUN
czarr = dblarr(71)
czerrarr = fltarr(71)
i = 0
cz = 0.D
WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun1,cz,czerr
    czarr(i) = cz
    czerrarr(i) = czerr
    i = i + 1
ENDWHILE
close, lun
free_lun, lun


plothist, czarr,bin = 0.05E4, xrange=[3.8E4,5.2E4]
mean = biweight_mean(czarr, sigma)
print, mean, sigma

device, /close
set_plot, mydevice


END
