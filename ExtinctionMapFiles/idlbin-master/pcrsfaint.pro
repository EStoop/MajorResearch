pro pcrsfaint
;this test was run in July 2011 to test the magnitude range for which pcrs peakup to irac mode can be used.
;test re-run Sep 2011

!P.multi = [0,1,1]
;read in an AOR of data
;measure centroids.

colorarr = [  'dark_red', 'dark_red','red', 'red', 'orange',  'orange', 'gold', 'gold', 'spring_green', 'spring_green', 'dark_green','dark_green', 'blue','blue',  'black','black', 'grey','grey']


;faint_aor = ['0042537984','0042538240','0042538496','0042538752','0042577152','0042577408','0042577664','0042577920','0042578176','0042578432','0042578688','0042578944','0042579200','0042579456','0042541824','0042541568','0042542080','0042542336']

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
   ;dir = '/Users/jkrick/irac_warm/pcrsfaint/' + string(faint_aor(a)) 
   dir = '/Users/jkrick/iracdata/flight/IWIC/IRAC028700/bcd/' + string(faint_aor(a)) 
   

   CD, dir                      ; change directories to the correct AOR directory
                                ;first = 5
   command  =  " ls *bcd_fp.fits > /Users/jkrick/irac_warm/pcrsfaint/cbcdlist.txt"
   spawn, command
   
   readcol,'/Users/jkrick/irac_warm/pcrsfaint/cbcdlist.txt',fitsname, format = 'A', /silent
   yarr = fltarr(n_elements(fitsname)*64)
   xarr = fltarr(n_elements(fitsname)*64)
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
      
      if i eq n_elements(fitsname) - 1 then stoptime = sxpar(header, 'SCLK_OBS')
      
      get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5],ra = ra_ref, dec = dec_ref, /silent
      
      deltaxarr[a] = mean(x_center) ;- 15.0
      deltayarr[a] = mean(y_center) ;- 15.0
                                ;print,'mean(x_center, ycenter)', mean(x_center), mean(y_center)
      
      if mean(x_center) gt 20 then begin
         x_center = x_center - 9.0
         y_center = y_center - 216.5
      endif

      if i eq 0 then begin
         xarr = x_center
         yarr = y_center
      endif else begin
         xarr = [xarr, x_center]
         yarr = [yarr, y_center]
      endelse
      

 
  

   endfor                       ; for each fits file in the AOR
   
   if a eq 0 then begin
      t = plot(xarr, yarr, '6r1s',xrange = [13.5,15.5], yrange = [14.0,16.0], xtitle = 'X pix',ytitle = 'Y pix',  color = colorarr[a], sym_filled = 1)
  ;    t = plot(xarr, yarr, '6r1s',xrange = [23.5,24.5], yrange = [230.0,233.0], xtitle = 'X pix',ytitle = 'Y pix',  color = colorarr[a], sym_filled = 1)
    endif

   if a gt 0 and ((a/2.) eq fix(a/2.)) then begin
      t = plot( xarr, yarr, '6r1s',/overplot, color = colorarr[a], sym_filled = 1)
   endif
   if a gt 0 and ((a/2.) ne fix(a/2.)) then begin
      t = plot( xarr, yarr, '6r1o',/overplot, color = colorarr[a], sym_filled = 1)
   endif


print, 'total time ' ,(framtime * n_elements(fitsname)) / 60., ' minutes'
print, '--------------------------------------'

if a gt 3 then begin ; fullarray 
   totalarr[a] = (framtime * n_elements(fitsname)) / 60.
endif else begin ;subarray 
   totalarr[a] = ((framtime * n_elements(fitsname)) / 60.) * 64.
endelse

endfor                          ;for each AOR

for a = 0, n_elements(faint_aor) -1 do begin
   ;print, 'aorname', aorname[a]
   ttext = text( 13.7,  (14.40 + (a/15.)), string(aorname[a]) + string(fix(framarr[a])) + string(totalarr[a]),/data,  color = colorarr[a], font_size = 10)

endfor
mode_x =  14.7456 
mode_y =  15.0683
xsweet = 15.120
ysweet = 15.085

box_x = [mode_x-0.1, mode_x-0.1, mode_x + 0.1, mode_x + 0.1, mode_x -0.1]
box_y = [mode_y-0.1, mode_y +0.1, mode_y +0.1, mode_y - 0.1,mode_y -0.1]
line3 = polyline(box_x, box_y, thick = 3, color = !color.black,/data)

box_x = [xsweet-0.1, xsweet-0.1, xsweet + 0.1, xsweet + 0.1, xsweet -0.1]
box_y = [ysweet-0.1, ysweet +0.1, ysweet +0.1, ysweet - 0.1,ysweet -0.1]
line4 = polyline(box_x, box_y, thick = 3, color = !color.black,/data)

line_x = [13.65, 14.3]
line_y = [14.8, 14.8]
line5 = polyline(line_x, line_y, thick = 2, color = !color.black,/data)

line_x = [13.65, 14.3]
line_y = [15.32, 15.32]
line5 = polyline(line_x, line_y, thick = 2, color = !color.black,/data)

t2 = text(14.5, 14.8, 'Failed', baseline = [0,1.0,0])
end
