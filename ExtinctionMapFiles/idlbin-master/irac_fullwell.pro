pro irac_fullwell

; what is full well depth in electrons for ch2 0.4s full array?
; http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/11/ gives this in mJy

full_mjy = 980.
gain = 3.7
exptime =  0.2
fluxconv =  0.1469
sbtoe = gain*exptime/fluxconv

;convert to microJy
eflux = full_mjy *1E3
;convert to Mjy/sr
eflux = eflux / 34.98
;convert to e
eflux = eflux *sbtoe

print, 'ch2 0.4 full ', eflux

; what is full well depth in electrons for ch2 0.4s sub array?

full_mjy = 820.
gain = 3.7
exptime =  0.36
fluxconv =  0.1469
sbtoe = gain*exptime/fluxconv

;convert to microJy
eflux = full_mjy *1E3
;convert to Mjy/sr
eflux = eflux / 34.98
;convert to e
eflux = eflux *sbtoe

print, 'ch2 0.4 sub ', eflux


; what is full well depth in electrons for ch2 0.6s full array?

full_mjy = 650.
gain = 3.7
exptime =  0.4
fluxconv =  0.1469
sbtoe = gain*exptime/fluxconv

;convert to microJy
eflux = full_mjy *1E3
;convert to Mjy/sr
eflux = eflux / 34.98
;convert to e
eflux = eflux *sbtoe

print, 'ch2 0.6 full ', eflux

;--------------------------------------------
; what is full well depth in electrons for ch1 0.4s full array?
; http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/11/ gives this in mJy

full_mjy = 950.
gain = 3.7
exptime =  0.2
fluxconv =  0.1253
sbtoe = gain*exptime/fluxconv

;convert to microJy
eflux = full_mjy *1E3
;convert to Mjy/sr
eflux = eflux / 34.98
;convert to e
eflux = eflux *sbtoe

print, 'ch1 0.4 full ', eflux

; what is full well depth in electrons for ch1 0.4s sub array?

full_mjy = 1000.
gain = 3.7
exptime =  0.36
fluxconv =  0.1253
sbtoe = gain*exptime/fluxconv

;convert to microJy
eflux = full_mjy *1E3
;convert to Mjy/sr
eflux = eflux / 34.98
;convert to e
eflux = eflux *sbtoe

print, 'ch1 0.4 sub ', eflux

; what is full well depth in electrons for ch1 0.1s sub array?

full_mjy = 4000.
gain = 3.7
exptime =  0.08
fluxconv =  0.1253
sbtoe = gain*exptime/fluxconv

;convert to microJy
eflux = full_mjy *1E3
;convert to Mjy/sr
eflux = eflux / 34.98
;convert to e
eflux = eflux *sbtoe

print, 'ch1 0.1 sub ', eflux



end
