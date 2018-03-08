PRO revssma
close, /all
openr, lun1, "/n/Godiva1/jkrick/A3888/final/gallist", /get_lun
openw,lun1out, "/n/Godiva1/jkrick/A3888/final/revssma.dat", /get_lun
i = 0
rearr = fltarr(500)
smaarr = fltarr(500)

WHILE (NOT EOF(lun1)) DO BEGIN
    readf, lun1, xcenter, ycenter, sma, smb
    true = 0
    openr, lun2, "/n/Godiva1/jkrick/A3888/final/re.large.log", /get_lun
    WHILE (NOT EOF(lun2))DO BEGIN
        readf, lun2, xcenterre, ycenterre, re

        IF (xcenter GT xcenterre - 1.5) AND (xcenter LT xcenterre + 1.5) AND $  
          (ycenter  GT ycenterre - 1.5) AND (ycenter LT ycenterre + 1.5) THEN BEGIN  
           ;object is the same

            ;skip the galaxies in the very center
            IF (xcenter GT 1236) AND (xcenter LT 1636) AND $
                 (ycenter GT 1324) AND (ycenter LT 1724) THEN true = 1

            ;clipping by hand
;            IF (re GT 11.0) THEN true = 1

            IF (true EQ 0) THEN BEGIN
                printf, lun1out, sma, re
                rearr(i) = re
                smaarr(i) = sma
                i = i +1
            ENDIF

        ENDIF 
    ENDWHILE
    close, lun2
    free_lun, lun2
ENDWHILE
close,lun1
close,lun1out
free_lun, lun1
free_lun, lun1out

rearr = rearr[0:i-1]
smaarr = smaarr[0:i-1]

err = dindgen(i) - dindgen(i) + 1
start = [12.0]

result = MPFITFUN('linear',smaarr,rearr, err, start) 
resultbiweight = ROBUST_LINEFIT(smaarr,rearr, yfit, sig, coeff_sig)

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'
colors = GetColor(/load, Start=1)

device, filename = '/n/Godiva1/jkrick/A3888/final/revssma.ps', /portrait, $
  BITS=8, scale_factor=0.9 

plot, smaarr, rearr, psym = 2, thick = 3,xrange=[0,80], charthick = 3, xthick = 3, ythick = 3,$
   ytitle = "re", xtitle = "sma of masks in pixels"


oplot, (1/result(0))*findgen(80), findgen(80)
oplot, smaarr, result(0)*smaarr + sig, linestyle = 2, color = colors.blue
oplot, smaarr, result(0)*smaarr - sig, linestyle = 2, color = colors.blue

oplot, 5*findgen(80), findgen(80)
oplot, 22*findgen(80), findgen(80)




;oplot, smaarr, yfit, thick = 5
;oplot, smaarr, yfit+sig, thick = 5, linestyle = 1

device, /close
set_plot, mydevice

END
