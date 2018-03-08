PRO fillfact

close,/all
;device, /close
set_plot, 'x'


clusterxcen = 1554.0
clusterycen = 1781.0
rvir = 3121.0   ;r200

n = [0, 150., 300, 450, 600, 750]
i = 0

FOR i = 0, 4, 1 DO BEGIN

    galvolume = 0.
    ngals = 0.

    OPENR, lun, "/n/Godiva1/jkrick/A3888/final/member.V.full.input", /GET_LUN
    WHILE (NOT EOF(lun)) DO BEGIN
        READF, lun, o, xcenter, ycenter, a, b, e, f, m, isoarea, fwhm, pa, bkgd, apflux,isocorflux,isoflux
        dist = sqrt(((clusterxcen - xcenter)^2) + ((clusterycen - ycenter)^2))
        
        IF dist GT n(i) AND dist LT n(i+1) THEN BEGIN
                                ;sum up all the volume which is filled with galaxies
            
            sma = sqrt(isoarea/(!PI*(1-e)))
            smb = isoarea/(!PI*sma)
            rad = (sma+smb)/2.
            rad = 2.*rad ;we mask to roughly twice the isophotal radius
            volume = (4./3.)*!PI*(rad^3.)
            
            galvolume = galvolume + volume
            ngals = ngals + 1
            
        ENDIF
    ENDWHILE
    
    ;find the total volume in a cylinder at the correct radii
    h = sqrt(rvir^2 - (2*n(i+1))^2)
    h1 = sqrt(rvir^2 - (2*n(i))^2)

    totalvolume = !PI*(n(i+1)^2)*h - !PI*(n(i)^2)*h1
;    totalvolume = (4./3.)*!PI*(n(i+1)^3.) - ((4./3.)*!PI*(n(i)^3.))


    fillfactor = galvolume / totalvolume
    
    print, "within ", n(i),"pixels ", fillfactor, ngals, galvolume, totalvolume

ENDFOR

fill = [.328, .095, .013, .006, .009]
xfill = [75, 225, 375, 525, 675]
err = [1,1,1,1,1]
plot, xfill, fill

start = [-0.13,40.0]
result = MPFITFUN('linear',xfill,fill, err, start)   
oplot, xfill, result(0)*xfill + result(1), thick = 3

start = [100., 400.0]
result = MPFITFUN('exponential',xfill,fill, err, start)   
oplot, xfill, result(0)*exp(-xfill/(result(1))), thick = 3
END
