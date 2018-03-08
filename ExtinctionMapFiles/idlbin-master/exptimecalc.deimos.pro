;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Sep 2011
;extimate exposure time for DEIMOS on Keck
;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO exptimecalc

close, /all
;how long of an exptime is required to measure an XXth mag galaxy to
;S/N of 10?
m = findgen(100)/10 + 20.
print, m
;convert mag into flux (erg / (s cm^2 angstrom))
fnaught = 1.75E-9                 ;taken from astrophysical quantities for Rband A0V star
f = fnaught*10^(m/(-2.5))         ; m = -2.5log(f/f0)  => f= f0*10^(m/-2.5)

;convert to photons/(s cm^2 angstrom) by dividing by hnu
hc = 1.989E-16                                 ;erg cm
lamda = 5100.E-8                               ;cm
erg = hc/lamda
f = f  / (erg)       

;put in area of primary ; diameter = 508cm; pi*r^2
f = f *7.85E5  ;Keck 10mtr diameter

;correct for extinction through atmosphere
;of course this is for airmass = 1; mainly guessing over this lambda range
;http://home.fnal.gov/~annis/astrophys/filters/palomar.html
;f = f *0.70

;correct for extinction through telescope system including atmosphere
;http://www.ucolick.org/~ripisc/results.html
;600l/mm grating/550 nm/ @6000

f = f*0.29

;correct for galaxy flux falling outside of slit
;asume re ~1", so slit = 0.5re (in radius)
;0.5re includes ~ 25% of the light of the galaxy
;psf will smear out some of that
;(trujillo)
f = f *0.9  ; a guess for large galaxies

;correct for final resolution (grating + slit smearing)
res = 0.65                                                 ;angstroms/pix for 600/mm 7500 blaze
f = f *res                                              ;now in units of photons/(s pix)

;adjust for gain   now in units of counts/s pix
gain = 1.0                        
f = f /gain

;f is source counts


;correct for skynoise
;from various websites, not very exact
;days from new moon; R SB
;fukugita r-R = 0.35 for galaxies, probably not for sky, but closer than 0
;  0    21.0   
;  10   19.5 

;--------------------------------------------------------------------------
plate = .119                                ; arcsecs/pixel blue side
swidth = 1.0                                    ; arcsecond

msky = fltarr(2)
;guess for new and full moon from http://www.palomar.caltech.edu:8000/maintenance/darksky/index.tcl
msky = [22.0,21.25]   ;but these are surface brightnesses
fsky = fnaught*10^(msky/(-2.5))  ; now flux /squararcsec or erg / (s cm^2 angstrom squarcsec)
fsky = fsky * (plate^2)*swidth  ;now we are back to flux from point source=  erg / (s cm^2 angstrom)
fsky = fsky  / (erg)                       ; now photons/(s cm^2 angstrom)
fsky = fsky *2.03E5                      ; now electrons/(s angstrom)
fsky = fsky/gain                           ; now counts/(s angstrom)
fsky = fsky *res                             ; now counts/(s pix)
fsky = fsky*0.29                            ;efficiency of the detector

;have not dealt with the extinction due to atmospheric transmission because this is the
;atmospheric transmission


;readnoise
R =3.0 ;e-
R = R/gain                       ;http://www2.keck.hawaii.edu/observing/kecktelgde/ktelinstupdate.pdf

;see what happens if we change moon
FOR j = 0, 1 DO begin
    bsquared = (100*(f + fsky(j))) ^ 2
    fourac = 4*(f^2)*(R^2)
    bee = 100*(f+fsky(j))
    t = (bee + sqrt(bsquared + fourac) )/(2*f^2)
    IF (j EQ 0) THEN BEGIN
        s = plot( m, t, xtitle = "magnitude", ytitle= " exposure time (seconds)",/ylog,yrange = [1E2, 1E7],xrange = [20, 28])
              
    ENDIF ELSE BEGIN
        s = plot( m, t, linestyle = 2,/overplot)
    ENDELSE
;print, j
;print,  t

yline = findgen(31) + 3600.
xline = findgen(31)
s = plot (xline, yline, /overplot, linestyle = 2)

yline2 = [1E2, 1E3, 1E4, 1E5, 1E6, 1E7]
xline2 = fltarr(6) + 25.7
s = plot (xline2, yline2, /overplot, linestyle = 3)

yline2 = [1E2, 1E3, 1E4, 1E5, 1E6, 1E7]
xline2 = fltarr(6) + 25.1
s = plot (xline2, yline2, /overplot, linestyle = 3)


ENDFOR






end
