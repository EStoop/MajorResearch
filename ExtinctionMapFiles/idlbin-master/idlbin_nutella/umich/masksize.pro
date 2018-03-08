PRO masksize

close, /all

openr, lun, "/n/Godiva1/jkrick/A3888/final/masksize.log", /get_lun
Rarr = fltarr(10000)
smaarr = fltarr(10000)
i = 0
WHILE(NOT eof(lun)) DO BEGIN
    readf, lun, xcenter, ycenter, sma, smb
    Rarr[i] = sqrt((xcenter - 1556.)^2 + (ycenter - 1784.)^2)
    smaarr[i] = sma
    i = i + 1
ENDWHILE

Rarr= Rarr[0:i-1]
smaarr = smaarr[0:i-1]

ps_open, file = "/n/Godiva1/jkrick/A3888/final/masksize.ps", /portrait, /color, ysize = 7

;plot, Rarr, smaarr, thick = 3, charthick=3, xthick=3,ythick=3, psym = 2,$
;  xtitle = "distance from center of cluster", ytitle = "mask radius"
plothist, smaarr, xrange=[0,100], xtitle="semi-major axis of mask",$
  ytitle="number of masks", thick = 3, charthick=3, xthick=3,ythick=3

ps_close, /noprint, /noid

close, lun
free_lun, lun

END
