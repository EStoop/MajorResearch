pro pmap_max

 fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/pmap_ch2_500x500_0043_111129.fits', pmapdata, pmapheader
pmapdata = pmapdata/ max(pmapdata)
fits_write,  '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/pmap_ch2_500x500_0043_maxnorm.fits', pmapdata, pmapheader

end
