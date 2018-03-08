pro halpha_exptime

fuv = 3.049D-27   ; erg/s/cm^2/Hz  from virgo_galex.pro
;over an area of 13695 pixels at 1.5" / pix
area = 13695*1.5*1.5  ;in square arcsec


;convert flux to luminosity
mpctocm = 3.08568D24
D = 18.D ;Mpc
D = D * mpctocm
Luv = 4.*(!Pi)*(D^2.)*fuv  ;erg/s/Hz
print, 'luv', luv, 'erg/s/Hz'

;convert UV luminosity to SFR (Kennicutt)
SFR = (1.4E-28)*(Luv)   ; Msun / yr
print, 'sfr', sfr, 'msun/yr'

;convert SFR to Halpha luminosity
Lha = SFR / (7.9E-42)   ; erg/s
print, 'lHa', Lha, 'erg/s'

;then into Halpha flux
fha = Lha/ (4.*!Pi*D^2.)
print, 'fha', fha, 'erg/s/cm2'

;this value is flux over the whole plume = 13695 FUV pixels = 3.08E4 sqarcsec
;so what is the flux per 1 squarcsec = binning ~5x5

fha_sqarcsec = fha / 3.08E4


;need flux density to go to magnitudes
;assume that line is spread over a certain number of Angstroms?
;guess 5 angs = 
lambda = 5.
c =   3.E8 ; m/s
c = c / 1.E10 ;in angstroms/s
hz = c / lambda
fha_sqarcsec = fha_sqarcsec / hz
print, 'fha', fha_sqarcsec, 'erg/s/cm2/Hz'

;convert flux density to magnitudes
;this is then the magnitudes in the Halpha line
AB = -2.5*alog10(fha_sqarcsec) - 48.6
print, 'AB', AB



;what is the rband magnitude of the plume?
;have B-V = 0.69, need to assume something about SED shape to get to r

;fukugita et al 1996, for this B-V,  V - R = 0.55
; and r - R = .35
; and r  = r(AB) +.226
; R = R (AB) -0.055

V = 15.87 ;johnson
R = V -0.55
r = R - .35
rab = r -0.226

;again this is over the whole plume
;want this in say 3" apertures
;easier to go into flux

fr = 10^((rab + 48.6)/(-2.5)) ;erg/s/cm^2/Hz  over whole plume

fr = fr / 3.08E4  ; now erg/s/cm2/Hz over 1 square arcsecond
fr = fr * 9.; 7.068  ; now erg/s/cm2/Hz over 3"diameter aperture

;halpha filter is 100Ang wide
;r filter is 1470 ang wide
;assume flux is uniformly spread out.

fr = fr / 14.7


rab = -2.5*alog10(fr) - 48.6 ;magnitude in 3" diameter aperture
print, 'r mag (AB)', rab

;what is the exposure time it would take to get to this depth?
;Jason's proposals say mR = 26 in 3" diameter aperture in 4 hours
fm = 10^((26 + 48.6)/(-2.5))
fd =  10^((rab + 48.6)/(-2.5))
print, fd, fm
exptime = 4./ ((fd/fm)^2.)
print, 'exptime', exptime
end
