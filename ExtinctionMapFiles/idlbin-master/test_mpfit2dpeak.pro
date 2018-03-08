pro test_mpfit2dpeak

  fits_read, '/Users/jkrick/external/irac_warm/trending/r41591808/ch1/bcd/SPITZER_I1_41591808_0007_0000_2_bcd.fits', data, header
  slice = data[*, *, 10]        ; pick one subarray image

  ;;is there a star on that image
  
  yfit = mpfit2dpeak(slice, A)
  print, A
  if A(1) gt 10*A(0) then print, 'Got one'
end
