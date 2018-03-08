pro hd158460_nonlinear

;are we seeing evidence for nonlinearity in the HD158460 data?

restore, '/Users/jkrick/irac_warm/snapshots/snapshots_corr.sav'
 aorname = ['r44497152','r44497664','r44497920','r44498176','r44498432','r44498688','r44498944','r44499200','r44499456','r44500224','r44500480','r44500736', 'r44500992', 'r44501248', 'r44501504' ]
 colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' , ' red',  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'green', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET'  , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU',   ' POWDER_BLUE' ]
mean_peak = fltarr(n_elements(aorname))
mean_corrflux = fltarr(n_elements(aorname))

 for a = 0,n_elements(aorname) - 1 do begin
   mmm, snapshots[a].peakpix, mode, sigma, skew
   mean_peak[a] = mode
   mmm, snapshots[a].corrflux, mode, sigma, skew
   mean_corrflux[a]=mode

    if a eq 0 then begin
      p = plot(snapshots[a].peakpix,snapshots[a].corrflux,  '1o', sym_filled = 1, sym_size = 0.1, color = colorarr[a], ytitle = 'Corrected Aperture Flux (Jy)', xtitle ='Raw Peak Pixel (DN)', title = 'Snapshot observations HD158460 0.4s Ch2 sub', yrange = [1.02, 1.05], xrange = [2.5E4, 2.9E4] )
   endif else begin
       p = plot(snapshots[a].peakpix,snapshots[a].corrflux,  '1o', sym_filled = 1, sym_size = 0.1, color = colorarr[a],/overplot)
endelse
     
endfor

     ;what does this look like for the pmap data?

;restore, '/Users/jkrick/irac_warm/pmap/pmap_corr.sav'
 ;aorname=['r44464128','r44463872','r44463616','r44463360','r44463104','r44462848','r44462592','r44462336','r44462080','r44461824']

;mean_peak = fltarr(n_elements(aorname))
;mean_corrflux = fltarr(n_elements(aorname))
;for a = 0, n_elements(aorname) - 1 do begin
   
;   mmm, snapshots[a].peakpix, mode, sigma, skew
;   mean_peak[a] = mode
;   mmm, snapshots[a].corrflux, mode, sigma, skew
;   mean_corrflux[a]=mode

;   if a eq 0 then begin
;      pm = plot(snapshots[a].peakpix,snapshots[a].corrflux,  '1o', sym_filled = 1, sym_size = 0.1, color = colorarr[a], ytitle = 'Corrected Aperture Flux (Jy)', xtitle ='Raw Peak Pixel (DN)', title = 'Pmap observations BD+67 0.1s Ch2 sub', yrange = [0.41, 0.44], xrange = [2400, 3300])
 ;  endif else begin
;       pm = plot(snapshots[a].peakpix,snapshots[a].corrflux,  '1o', sym_filled = 1, sym_size = 0.1, color = colorarr[a],/overplot)
;    endelse
     
;endfor


mean_corrflux = mean_corrflux(sort(mean_peak))
mean_peak = mean_peak(sort(mean_peak))

   pm = plot(mean_peak, mean_corrflux, '3+',linestyle = 6, symbol = 'star', /sym_filled, sym_size = 2,/overplot)


  end


;aorname = ['r44497152','r44497664','r44497920','r44498176','r44498432','r44498688','r44498944','r44499200','r44499456','r44500224','r44500480','r44500736', 'r44500992', 'r44501248', 'r44501504' ]
;peakarr = dblarr(n_elements(aorname) * 64.*75.)
;count = 0.D
;for a =0, 0 do begin; n_elements(aorname) - 1 do begin
;     print, 'working on ',aorname(a)
;     dir = '/Users/jkrick/irac_warm/snapshots/' + string(aorname(a)) 
;     CD, dir                    ; change directories to the correct AOR directory
;     command  =  " ls ch2/raw/*.fits > /Users/jkrick/irac_warm/snapshots/rawlist.txt"
;     spawn, command
; 
;     readcol,'/Users/jkrick/irac_warm/snapshots/rawlist.txt',fitsname, format = 'A', /silent
;     print, n_elements(fitsname) 
;
;     for j = 0, n_elements(fitsname) - 1 do begin
;        ;find the peak pixel value.
;        print, fitsname(j)
;        fits_read, fitsname(j), data, header
;        mc = reform(data, 32, 32, 64)
;;tvscl, mc(*,*,0)
;        
;        img = mc
;; convert to real
;        img=img*1.
;        
;; fix the problem with unsigned ints
;        fix=where((img LT 0),count)
;        if (count GT 0) then img[fix]=img[fix]+65536.
;        
;; flip the InSb
;        img=65535.-img
;        
;;        fits_write, '/Users/jkrick/irac_warm/snapshots/r44497152/ch2/raw/test_0023.fits', img, header
;;then use ds9 to display
;  ;print, img[15, 16]
;
;        ;for i = 0 , 63 do begin
;        ;   arr = img(*,*,i)
;        ;   print, 'i', i, arr[15,16]
;        ;   peakarr[count] = arr[15,16]
;        ;   count = count + 1
;        ;endfor
;       
;     endfor  ; for each fits file
;  endfor  ;for each AOR
;
