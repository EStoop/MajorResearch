pro pcrsfaint_v2
;this test was run in July 2011 to test the magnitude range for which pcrs peakup to irac mode can be used.
;test re-run Sep 2011

!P.multi = [0,1,1]
;read in an AOR of data
;measure centroids.

colorarr = [  'dark_red', 'red', 'orange']



faint_aor = ['0044250880', '0044251136', '0044251392']

t0arr = fltarr(n_elements(faint_aor))
deltaxarr =  fltarr(n_elements(faint_aor))
deltayarr =  fltarr(n_elements(faint_aor))
aorname = strarr(n_elements(faint_aor))
totalarr = strarr(n_elements(faint_aor))
framarr = strarr(n_elements(faint_aor))
for a = 0, n_elements(faint_aor) - 1 do begin
   ni = 0
   print, 'working on ', faint_aor(a)
   dir = '/Users/jkrick/iracdata/flight/IWIC/IRAC028700/bcd/' + string(faint_aor(a)) 
   

   CD, dir                      ; change directories to the correct AOR directory
                                ;first = 5
   command  =  " ls *bcd_fp.fits > /Users/jkrick/irac_warm/pcrsfaint/cbcdlist.txt"
   spawn, command
   
   readcol,'/Users/jkrick/irac_warm/pcrsfaint/cbcdlist.txt',fitsname, format = 'A', /silent
   yarr = fltarr(n_elements(fitsname))
   xarr = fltarr(n_elements(fitsname))
   c = 0
   
   for i =1, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
                               ;print, 'working on ', fitsname(i)         
                                ;now were did it really point to?

      if i eq 1 then begin
         header = headfits(fitsname(i)) ;
         aorname[a] = strmid(sxpar(header, 'AORLABEL'),10,5)
         print, 'aorname', aorname[a]
         print, 'naxis1', sxpar(header, 'naxis1')
         ra_ref = sxpar(header, 'RA_REF')
         dec_ref = sxpar(header, 'DEC_REF')
         print, 'ra, dec', ra_ref, dec_ref
         starttime = sxpar(header, 'SCLK_OBS')
         framtime =  sxpar(header, 'FRAMTIME')
         print, 'framtime', framtime
         framarr[a] = framtime

      endif
      
      
      get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5],ra = ra_ref, dec = dec_ref, /silent
      
      xarr[i] = x_center
      yarr[i] = y_center
   endfor                       ; for each fits file in the AOR
   
   print, 'xarr', xarr

   if a eq 0 then begin
      xy = plot(xarr, yarr, '6r1o',xrange = [23.5,24.5], yrange = [231.0,232.0], xtitle = 'X pix',ytitle = 'Y pix',  color = colorarr[a], sym_filled = 1, title = 'pcrsfaint V=12.5 12s ch2full')
  ;    t = plot(xarr, yarr, '6r1s',xrange = [23.5,24.5], yrange = [230.0,233.0], xtitle = 'X pix',ytitle = 'Y pix',  color = colorarr[a], sym_filled = 1)
    endif

   if a gt 0  then begin
      xy.Select
      xy = plot( xarr, yarr, '6r1o',/overplot, color = colorarr[a], sym_filled = 1,/current)
   endif


   print, 'total time ' ,(framtime * n_elements(fitsname)) / 60., ' minutes'
   print, '--------------------------------------'
   
   if a gt 3 then begin         ; fullarray 
      totalarr[a] = (framtime * n_elements(fitsname)) / 60.
   endif else begin             ;subarray 
      totalarr[a] = ((framtime * n_elements(fitsname)) / 60.) * 64.
   endelse
   

   timearr = findgen(n_elements(xarr)) * 12

   if a eq 0 then begin
      st = plot(timearr, yarr,'6r1+',yrange = [231.0,232.0], xrange = [0, 250], xtitle = 'Time(seconds)',ytitle = 'Y pix', color = colorarr[a], title = 'pcrsfaint')
      st2 = plot(timearr, xarr,'6r1+',yrange = [23.5,24.5], xrange = [0, 250], xtitle = 'Time(seconds)',ytitle = 'X pix', color = colorarr[a], title = 'pcrsfaint')
   endif

   if a gt 0 then begin
      st.Select
      st = plot( timearr, yarr, '6r1+',/overplot, color = colorarr[a],/current)
      st2.Select
      st2 = plot( timearr, xarr, '6r1+',/overplot, color = colorarr[a],/current)
  endif



endfor                          ;for each AOR

;for a = 0, n_elements(faint_aor) -1 do begin
                                ;print, 'aorname', aorname[a]
;   ttext = text( 13.7,  (14.40 + (a/15.)), string(aorname[a]) + string(fix(framarr[a])) + string(totalarr[a]),/data,  color = colorarr[a], font_size = 10)

;endfor
xsweet = 15.120 + 9.0
ysweet = 256. - 9. - ( 31. - 15.085)

xy.Select

box_x = [xsweet-0.1, xsweet-0.1, xsweet + 0.1, xsweet + 0.1, xsweet -0.1]
box_y = [ysweet-0.1, ysweet +0.1, ysweet +0.1, ysweet - 0.1,ysweet -0.1]
line4 = polyline(box_x, box_y, thick = 3, color = !color.black,/data,/current)

; save these files
st.Save , '/Users/jkrick/irac_warm/pcrsfaint/yvstime.png'
st2.Save ,  '/Users/jkrick/irac_warm/pcrsfaint/xvstime.png'
xy.Save ,  '/Users/jkrick/irac_warm/pcrsfaint/xy.png'

end
