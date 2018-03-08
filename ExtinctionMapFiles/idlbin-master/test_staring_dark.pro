pro test_staring_dark

fits_read, '/Users/jkrick/irac_warm/darks/staring/staringdark_s0p4s_ch2.fits', stare0p4sdata, stare0p4sheader


fits_read, '/Users/jkrick/irac_warm/darks/staring/r51667200/ch2/cal/SPITZER_I2_51525888_0000_1_C9972578_sdark.fits', p4sdata, p4sheader


diffim = stare0p4sdata - p4sdata
fits_write, '/Users/jkrick/irac_warm/darks/staring/diffim_ch2_0p4s.fits', diffim, p4sheader



end
