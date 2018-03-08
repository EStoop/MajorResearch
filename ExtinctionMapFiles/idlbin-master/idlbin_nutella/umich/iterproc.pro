PRO ITERPROC, MYFUNCT, p, iter, fnorm, FUNCTARGS=fcnargs, $
                PARINFO=parinfo, QUIET=quiet, DOF=dof
COMMON sharegal, inten, ierror, length

print, "iteration number ", iter, " paramters", p

;will make my own value for chisq, which is the chisq for the fit
;of the convolution of the psf and the Sersic

;1d
devaucarr = fltarr(length)
moffatarr = fltarr(length)

scale = 0.5  ;need to change this to be a normalizing factor

FOR x = 0,length-1,1 DO BEGIN
    devaucarr[x] = (p(0)) * (exp(-7.67*(((x/p(1))^(1.0/4.0)) - 1.0)))
    moffatarr[x] = scale * (1 + ((x) / 3.218)^2)^(-2.592)
ENDFOR
;print, devaucarr

fd = fft(devaucarr, -1)
fm = fft(moffatarr, -1)

;print, fd
;plot, fd
prod = fd * fm

con = fft(prod, 1)

;print, "did the convolution:"

;ok now I need to recalculate the chisq given this convolution.
blah  = total(((float(con) - inten)/ ierror)^2)

print, "orig chisq:", fnorm, " new chisq = ", float(blah)

fnorm = float(blah)
END



;need to create galaxy data that is an image of the galaxy with this
;iterations parameters

;2d
;FITS_READ, '/n/Truffle1/jkrick/A3888/blank.fits', blankdata, blankheader
;newgal = subdev2(xcen, ycen, maxradius, blankdata, p, PA, ellipticity)
;fits_read, "/n/Truffle1/jkrick/A3888/npsf.fits", psfdata, psfheader
