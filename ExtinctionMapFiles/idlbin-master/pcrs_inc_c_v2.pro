pro pcrs_inc_c_v2

  colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' , ' red',  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'green', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET' , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU',   ' POWDER_BLUE' ]

ra =     269.72751   
dec=     67.79463

  ;read in all the pixel phase correction data
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/pmap_ch2_500x500_0043_111129.fits', pmapdata, pmapheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/psigma_ch2_500x500_0043_111129.fits', psigmadata, psigmaheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/occu_ch2_500x500_0043_111129.fits', occdata, occheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/xgrid_ch2_500x500_0043_111129.fits', xgriddata, xgridheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/ygrid_ch2_500x500_0043_111129.fits', ygriddata, ygridheader


;second test in Sep

incc_aor = ['0044254976','0044255232','0044255488','0044255744','0044256000','0044256256','0044256512','0044256768','0044257024','0044257280']

;44254976 ,44255232 ,44255488 ,44255744 ,44256000 ,44256256 ,44256512 ,44256768 ,44257024 ,44257280

for a = 0, n_elements(incc_aor) - 1 do begin
   ni = 0
   print, 'working on ', incc_aor(a)
   dir = '/Users/jkrick/iracdata/flight/IWIC/IRAC028700/bcd/' + string(incc_aor(a)) 
   

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
         ra_ref = sxpar(header, 'RA_REF')
         dec_ref = sxpar(header, 'DEC_REF')
         print, 'ra, dec', ra_ref, dec_ref
      endif
      
      
      get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5],ra = ra_ref, dec = dec_ref, /silent     
      if i eq 1 then begin
         xarr = x_center
         yarr = y_center
      endif else begin
         xarr = [xarr, x_center]
         yarr = [yarr, y_center]
      endelse
      
   endfor                       ; for each fits file in the AOR
   
   ;print, 'xarr', xarr

   if a eq 0 then begin
      xy = plot(xarr, yarr, '6r1*',xrange = [14.5,15.5], yrange = [14.5,15.5], xtitle = 'X pix',ytitle = 'Y pix',  color = colorarr[a], title = 'pcrs peakup Inc_c NPM1p67.0536 ')
  ;    t = plot(xarr, yarr, '6r1s',xrange = [23.5,24.5], yrange = [230.0,233.0], xtitle = 'X pix',ytitle = 'Y pix',  color = colorarr[a], sym_filled = 1)
    endif

   if a gt 0  then begin
      xy.Select
      xy = plot( xarr, yarr, '6r1*',/overplot, color = colorarr[a],/current)
   endif
  
   timearr = findgen(n_elements(xarr)) * 0.4

   if a eq 0 then begin
      st = plot(timearr, yarr,'6r1+',yrange = [14.5, 15.5], xrange = [0, 450], xtitle = 'Time(seconds)',ytitle = 'Y pix', color = colorarr[a], title = 'inc_c')
      st2 = plot(timearr, xarr,'6r1+',yrange = [14.5, 15.5], xrange = [0, 450], xtitle = 'Time(seconds)',ytitle = 'X pix', color = colorarr[a], title = 'inc_c')
   endif

   if a gt 0 then begin
      st.Select
      st = plot( timearr, yarr, '6r1+',/overplot, color = colorarr[a],/current)
      st2.Select
      st2 = plot( timearr, xarr, '6r1+',/overplot, color = colorarr[a],/current)
  endif



endfor                          ;for each AOR

mode_x =  14.7456 
mode_y =  15.0683
xsweet = 15.120
ysweet = 15.085
xy.Select
box_x = [mode_x-0.1, mode_x-0.1, mode_x + 0.1, mode_x + 0.1, mode_x -0.1]
box_y = [mode_y-0.1, mode_y +0.1, mode_y +0.1, mode_y - 0.1,mode_y -0.1]
line3 = polyline(box_x, box_y, thick = 3, color = !color.black,/data)

box_x = [xsweet-0.1, xsweet-0.1, xsweet + 0.1, xsweet + 0.1, xsweet -0.1]
box_y = [ysweet-0.1, ysweet +0.1, ysweet +0.1, ysweet - 0.1,ysweet -0.1]
;line4 = polyline(box_x, box_y, thick = 3, color = !color.black,/data)

; save these files
st.Save , '/Users/jkrick/irac_warm/pcrs_inc_c/yvstime.png'
st2.Save ,  '/Users/jkrick/irac_warm/pcrs_inc_c/xvstime.png'
xy.Save ,  '/Users/jkrick/irac_warm/pcrs_inc_c/xy.png'
end
