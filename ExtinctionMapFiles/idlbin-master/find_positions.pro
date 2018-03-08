pro find_positions

  ra_ref = 133.14756
  dec_ref = 28.330094
  aorname = ['r43981312','r43981568','r43981824','r43981056'] 
 
  for a = 0, n_elements(aorname) -1 do begin 

     ;change these to your actual image locations
   fits_read, '/Users/jkrick/irac_warm/pcrs_planets/55cnc/'+aorname(a)+'/ch2/bcd/SPITZER_I2_'+strmid(aorname(a),1) + '_4500_0000_1_bcd.fits', data, header
   fits_read, '/Users/jkrick/irac_warm/pcrs_planets/55cnc/'+aorname(a)+'/ch2/bcd/SPITZER_I2_'+strmid(aorname(a),1) + '_4500_0000_1_bunc.fits', unc, uncheader

   ;measure centroids, do aperture photometry
   get_centroids_for_calstar_jk,data, header, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                     x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                     x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                     xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                     xfwhm, yfwhm, /WARM
        
   ;convert x & y into ra and dec in degrees
   xyad, header, x3(10), y3(10), ra, dec

   ;concert ra and dec in degrees into hrs/min/sec deg/min/sec.
   radec, ra, dec, ihr, imin, xsec, ideg, imn, xsc

   ;get the date from the headers
   date = sxpar(header, 'DATE_OBS')

   print, aorname(a), '  ', date, x3(10), y3(10), x3s(10), y3s(10), ihr, imin, xsec, ideg, imn, xsc, $
          format = '( A, A, A, F10.3, F10.3, F10.3, F10.3, I10, I10, F10.4, I10, I10, F10.4)'
 
endfor

end
