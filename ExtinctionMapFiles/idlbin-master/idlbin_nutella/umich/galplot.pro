PRO galplot
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;feb 2003
; makes plots of galaxy profiles given simple equations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


device, true=24
device, decomposed=0
colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva2/jkrick/A141/idlgalaxies.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;make radius array
r = fltarr(4001.0)
FOR i = 0.0, 4000.0 DO BEGIN
    r[i] = i/10.0
ENDFOR

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;sersic
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

n = 1.0
k = 1.9992*n - 0.3271
re = 10.0
sige = 23.0

sigr = sige + 1.086*(k*( (r/re)^(1/n)- 1))
    plot, r, sigr,thick = 3, $
	YRANGE = [30,15] ,xrange = [1,100], /xlog,$
	ytitle = 'mag/sqarcsec', charthick = 3, xthick = 3 ,ythick = 3,$
  xtitle = 'radius in arcseconds'

det = 25.45  ; detection threshold FOR A3888 in V
rad = re*( ( (det-sige)/(1.086*k)) +1)^n
print,rad, det, "    A141"
det = 26.0  ; detection threshold FOR A3888 in V
rad = re*( ( (det-sige)/(1.086*k)) +1)^n
print,rad, det, "    A3888"


xyouts, rad, det, "*"
xyouts, 2*rad, det, "*"

;FOR n = 0.5,3.0,0.5 DO BEGIN
;    k = 1.9992*n - 0.3271
;    sigr = sige + 1.086*(k*( (r/re)^(1/n)- 1))
;    
;    oplot, r, sigr, color = colors.blue
;ENDFOR
;n = 0.5
;k = 1.9992*n - 0.3271
;re = 12.0
;sige = 20.0;

;sigr = sige + 1.086*(k*( (r/re)^(1/n)- 1))
;oplot, r, sigr,thick = 3, color = colors.red
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;exponential disk
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
;FOR H = 6.,10.,1 DO BEGIN
;    sigr = 1.086*(r/H) + 20.   ;in magnitude world
;    oplot, r, sigr, thick = 3,color = colors.green
;ENDFOR

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

device, /close
set_plot, mydevice

END

