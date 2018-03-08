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

device, filename = '/Users/jkrick/palomar/cosmic/exptime.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

;how long of an exptime is required to measure an XXth mag galaxy to
;S/N of 10?
m = findgen(50)/10 + 18.
print, m
;convert mag into flux (erg / (s cm^2 angstrom))
fnaught = 1.75E-9                 ;taken from astrophysical quantities for Rband A0V star
f = fnaught*10^(m/(-2.5))         ; m = -2.5log(f/f0)  => f= f0*10^(m/-2.5)

;convert to photons/(s cm^2 angstrom) by dividing by hnu
hc = 1.989E-16                                 ;erg cm
lamda = 5500.E-8                               ;cm
erg = hc/lamda
f = f  / (erg)       

;put in area of primary ; diameter = 508cm; pi*r^2
f = f *2.03E5

;correct for extinction through atmosphere
;of course this is for airmass = 1; mainly guessing over this lambda range
;http://home.fnal.gov/~annis/astrophys/filters/palomar.html
f = f *0.70

;correct for extinction through telescope system
;from fig 4 Kells 98, 300l/mm grating
f = f*0.24

;correct for galaxy flux falling outside of slit
;asume re ~1", so slit = 0.5re (in radius)
;0.5re includes ~ 25% of the light of the galaxy
;psf will smear out some of that
;(trujillo)
f = f *1

;correct for final resolution (grating + slit smearing)
res = 3.1                                                 ;angstroms/pix for 300l/mm
f = f *res                                              ;now in units of photons/(s pix)

;adjust for gain   now in units of counts/s pix
;assume it is a factor of a few
;;gain = 2.
gain = 3.1                           ;  RAB-RAB-RAB
f = f /gain

;f is source counts


;correct for skynoise
;from various websites, not very exact
;days from new moon; R SB
;fukugita r-R = 0.35 for galaxies, probably not for sky, but closer than 0
;  0    21.0   
;  10   19.5 

;--------------------------------------------------------------------------
plate = .3988                                 ; arcsecs/pixel
swidth = 1.5                                    ; arcsecond

msky = fltarr(2)
msky = [21.0,19.5]   ;but these are surface brightnesses
fsky = fnaught*10^(msky/(-2.5))  ; now flux /squararcsec or erg / (s cm^2 angstrom squarcsec)
fsky = fsky * (plate^2)*swidth  ;now we are back to flux from point source=  erg / (s cm^2 angstrom)
fsky = fsky  / (erg)                       ; now photons/(s cm^2 angstrom)
fsky = fsky *2.03E5                      ; now electrons/(s angstrom)
fsky = fsky/gain                           ; now counts/(s angstrom)
fsky = fsky *res                             ; now counts/(s pix)
fsky = fsky*0.24                            ;efficiency of the detector

;have not dealt with the extinction due to atmospheric transmission because this is the
;atmospheric transmission


;readnoise
;;R = 12 e-
R = 12./3.1                        ;  users manual. ccd characteristics chapter

;see what happens if we change moon
FOR j = 0, 1 DO begin
    bsquared = (100*(f + fsky(j))) ^ 2
    fourac = 4*(f^2)*(R^2)
    bee = 100*(f+fsky(j))
    t = (bee + sqrt(bsquared + fourac) )/(2*f^2)
    IF (j EQ 0) THEN BEGIN
        plot, m, t, xtitle = "Rmagnitude", ytitle= " exposure time (seconds)",/ylog,$
              title="grism 300l"
    ENDIF ELSE BEGIN
        oplot, m, t, linestyle = 2
    ENDELSE
print, j
print,  t

ENDFOR





device, /close
set_plot, mydevice

end
