pro test_acs_back

fits_read, '/Users/jkrick/plumes/centaurus/hst/j96f01bbq_flt.fits', data, header, exten_no = 4

bkgd = data[400:1100,600:1200]
bkgd = data[2500:3200,1000:1600]
print, mean(bkgd), median(bkgd)
plothist, bkgd, xhist, yhist, bin = 0.1, xrange = [10,100]

start = [mean(bkgd),stddev(bkgd), 1E3]
noise = fltarr(n_elements(yhist))
noise[*] = 1                    ;equally weight the values
result= MPFITFUN('mygauss',xhist,yhist, noise, start) 
;----------------------------------------------------------



xarr = findgen(100)

oplot, xarr, (result(2))/sqrt(2.*!Pi) * exp(-0.5*((xarr - (result(0)))/(result(1)))^2.), $
       thick=5,  color = green


end
