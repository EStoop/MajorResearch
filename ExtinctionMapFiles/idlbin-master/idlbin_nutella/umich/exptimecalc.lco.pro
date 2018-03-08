;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;April 2004
;extimate exposure time for IMACS on Magellan
;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO exptimecalc

close, /all
;device, true=24
;device, decomposed=0

;colors = GetColor(/load, Start=1)

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = 'exptimeF2.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

;how long of an exptime is required to measure an XXth mag galaxy to
;S/N of 10?
m = findgen(50)/10 + 18.
print, m
;convert mag into flux (erg / (s cm^2 angstrom))
fnaught = 3.75E-9                 ;taken from astrophysical quantities for Vband A0V star
f = fnaught*10^(m/(-2.5))         ; m = -2.5log(f/f0)  => f= f0*10^(m/-2.5)

;convert to photons/(s cm^2 angstrom) by dividing by hnu
hc = 1.989E-16                                 ;erg cm
lamda = 5500.E-8                               ;cm
erg = hc/lamda
f = f  / (erg)       

;put in area of primary ; diameter = 650cm
f = f *3.32E5

;correct for extinction through atmosphere
;of course this is for airmass = 1; tau = 0.2 V band
f = f *0.80

;correct for extinction through telescope system
;including IMACS and grating f/2 with 200l grism
f = f*0.20

;correct for galaxy flux falling outside of slit
;asume re ~1", so slit = 0.5re (in radius)
;0.5re includes ~ 25% of the light of the galaxy
;psf will smear out some of that
;(trujillo)
f = f *1

;correct for final resolution (grating + slit smearing)
res = 2                                                 ;angstroms/pix
f = f *res                                              ;now in units of photons/(s pix)

;adjust for gain   now in units of counts/s pix
;assume it is a factor of a few
;;gain = 2.
gain = 1.                            ;  RAB-RAB-RAB
f = f /gain

;f is source counts


;correct for skynoise
;from laSilla webpage
;days from new moon; V SB
;  0    21.8   
;  3    21.7   
;  7    21.4   
; 10   20.7   
; 14   20.0   

;--------------------------------------------------------------------------
plate = .111                                  ; arcsecs/pixel
swidth = 1                                     ; arcsecond

msky = fltarr(5)
msky = [21.8,21.7,21.4,20.7,20.0]   ;but these are surface brightnesses
fsky = fnaught*10^(msky/(-2.5))  ; now flux /squararcsec or erg / (s cm^2 angstrom squarcsec)
fsky = fsky * (plate^2)*swidth  ;now we are back to flux from point source=  erg / (s cm^2 angstrom)
fsky = fsky  / (erg)                       ; now photons/(s cm^2 angstrom)
fsky = fsky *3.32E5                      ; now electrons/(s angstrom)
fsky = fsky/gain                           ; now counts/(s angstrom)
fsky = fsky *res                             ; now counts/(s pix)
fsky = fsky*0.20                            ;efficiency of the detector

;have not dealt with the extinction due to atmospheric transmission because this is the
;atmospheric transmission


;assume readnoise of 4
;;R = 4  
R = 5                           ;  RAB-RAB-RAB

;see what happens if we change moon
FOR j = 0, 4 DO begin
    bsquared = (100*(f + fsky(j))) ^ 2
    fourac = 4*(f^2)*(R^2)
    bee = 100*(f+fsky(j))
    t = (bee + sqrt(bsquared + fourac) )/(2*f^2)
    IF (j EQ 0) THEN BEGIN
        plot, m, t, xtitle = "Vmagnitude", ytitle= " exposure time (seconds)",/ylog,$
              title='Gain='+string(gain)+'  Res (A/pix)='+string(res)
    ENDIF ELSE BEGIN
        oplot, m, t, linestyle = 2
    ENDELSE
print, j
print,  t

ENDFOR





device, /close
set_plot, mydevice

end
