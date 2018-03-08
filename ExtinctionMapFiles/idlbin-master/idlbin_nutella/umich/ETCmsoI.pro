;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;April 2004      J. Krick
;April 2005      modifications by M. S. Oey
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
m = findgen(50)/10 + 15.
print, m
;convert mag into flux (erg / (s cm^2 angstrom))
;fnaught = 3.75E-9                 ;taken from astrophysical quantities for Vband A0V star
;fnaught = 6.65E-9                 ;**MSO:  taken from Zombeck for B
fnaught = 8.32E-10                 ;**MSO:  taken from Zombeck for I
f = fnaught*10^(m/(-2.5))         ; m = -2.5log(f/f0)  => f= f0*10^(m/-2.5)

;convert to photons/(s cm^2 angstrom) by dividing by hnu
hc = 1.989E-16                                 ;erg cm
; lamda = 5500.E-8                               ;cm
; lamda = 4400.E-8                               ;**MSO:  for B 
lamda = 8600.E-8                               ;**MSO:  for I
erg = hc/lamda
f = f  / (erg)       

;put in area of primary ; diameter = 650cm
f = f *3.32E5

;correct for extinction through atmosphere
;of course this is for airmass = 1; tau = 0.2 V band
;f = f *0.80
;of course this is for airmass = 1; tau = 0.26 B **MSO
; f = f *0.74
;of course this is for airmass = 1; tau = 0.11 I **MSO
f = f *0.89

;correct for extinction through telescope system
;including IMACS and grating
;f = f*0.20
;f = f*0.10               ; at 4000A **MSO
f = f*0.15               ; EST at 8600A FOR 600 l mm-1 **MSO 

;correct for galaxy flux falling outside of slit
;asume re ~1", so slit = 0.5re (in radius)
;0.5re includes ~ 25% of the light of the galaxy
;psf will smear out some of that
;(trujillo)
;f = f *0.25                     ; !!!!!!!  assume all light in slit   TB!!!!!!!!!! 
f = f *0.75                     ; !!!!!!!  most light in slit  **MSO

;correct for final resolution (grating + slit smearing)
;res = 2                       ;angstroms/pix
; res = 0.74                   ;f/4 camera 300 lines/mm
res = 0.8                     ;f/4 camera 600 lines/mm  **MSO
f = f *res                    ;now in units of photons/(s pix)

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
plate = 0.111                                  ; arcsecs/pixel
swidth = 1                                     ; arcsecond

msky = fltarr(5)
; msky = [21.8,21.7,21.4,20.7,20.0]   ;but these are surface brightnesses V
; msky = [22.7,22.4,21.6,20.7,19.5]   ;but these are surface brightnesses B
msky = [19.9,19.9,19.7,19.5,19.2]   ;but these are surface brightnesses I
fsky = fnaught*10^(msky/(-2.5))  ; now flux /squararcsec or erg / (s cm^2 angstrom squarcsec)
fsky = fsky * (plate^2)*swidth  ;now we are back to flux from point source=  erg / (s cm^2 angstrom)
fsky = fsky  / (erg)                       ; now photons/(s cm^2 angstrom)
fsky = fsky *3.32E5                      ; now electrons/(s angstrom)
fsky = fsky/gain                           ; now counts/(s angstrom)
fsky = fsky *res                             ; now counts/(s pix)
; fsky = fsky*0.20                            ;efficiency of the detector V
; fsky = fsky*0.10                            ;efficiency of the detector 4000A **MSO
fsky = fsky*0.15                            ;EST FOR 8600A FOR 600 l mm-1 **MSO

;have not dealt with the extinction due to atmospheric transmission because this is the
;atmospheric transmission


;assume readnoise of 4
;;R = 4  
R = 5                           ;  RAB-RAB-RAB

;see what happens if we change moon
print, "EXPOSURE TIMES FOR 0, 3, 7, 10, 14 DAYS FROM NEW MOON"
print
SNR = 50.
FOR j = 0, 4 DO begin
;    bsquared = (100*(f + fsky(j))) ^ 2
;    fourac = 4*(f^2)*(R^2)
;    bee = 100*(f+fsky(j))
;;    fourac = 4.*(f^2.)*(R^2.)
    fourac = 4.*(f^2.)*(R^2.)*(SNR^2.)
    bee = (SNR^2.)*(f+fsky(j))
    bsquared = bee^2.
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
