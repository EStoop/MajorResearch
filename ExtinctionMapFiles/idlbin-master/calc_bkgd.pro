function  calc_bkgd, eim, header, ra_ref, dec_ref
  im = eim
  ;;figure out and then convert image to electrons
;  sbtoe = sxpar(header, 'GAIN') * sxpar(header, 'EXPTIME') / sxpar(header,'FLUXCONV')
;  im = im*sbtoe

  ;;make an array to hold the return background values
  bkgdarr = fltarr(3)
;  bkgdarr = fltarr(3,64)

  ;;mask region of the star with nan    
  star_radius = 4
  adxy, header, ra_ref, dec_ref, x_ref, y_ref
  x_ref = round(x_ref)
  y_ref = round(y_ref)
  im[(x_ref-star_radius):(x_ref +star_radius), (y_ref -star_radius):(y_ref+star_radius)] = !Values.F_NAN 
;  im[(x_ref-star_radius):(x_ref +star_radius), (y_ref -star_radius):(y_ref+star_radius), *] = !Values.F_NAN 

  ;;mask the outer row/columns
  im[31,*] = !Values.F_NAN 
  im[0,*]= !Values.F_NAN 
  im[*,0]= !Values.F_NAN 
  im[31,*]= !Values.F_NAN 

  ;;for each of the subarray frames
;  for i = 0, 63 do begin

  ;;make a histogram of the rest of the background values
     plothist, im[*,*], xhist, yhist, bin = .1, /noprint,/noplot,/NAN
;     plothist, im[*,*,i], xhist, yhist, bin = .1, /noprint,/noplot,/NAN
         

  ;;fit a gaussian
     start = [mean(im[*,*],/nan), 0.5, 200]
;     start = [mean(im[*,*,i],/nan), 0.5, 200]
     junkerr = fltarr(n_elements(xhist)) + 1.0
     g_result= MPFITFUN('mygauss',xhist, yhist,junkerr, start, perror = perror, bestnorm = bestnorm,/quiet)  

;     if i lt 5 then begin
        ;;overplot said gaussian and its mean value
;        ph = barplot(xhist, yhist,  xtitle = 'Background', ytitle = 'Number', fill_color = 'sky_blue', xrange = [-2,2] )
;        ph = plot(xhist, g_result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - g_result(0))/g_result(1))^2.), color = 'black', overplot = ph)
;        ph = plot(intarr(150) + g_result(0), indgen(150), linestyle = 2, thick = 2, overplot = ph)
;     endif

     ;;figure how many data points were being used for the fit
     good = where(finite(im[*,*]) gt 0, ngood)
;     good = where(finite(im[*,*,i]) gt 0, ngood)

     bkgdarr[*] = [g_result(0), abs(g_result(1)), ngood]
;     bkgdarr[*,i] = [g_result(0), g_result(1), ngood]
;  endfor
;  print, 'at end of bkgdarr' , bkgdarr
  return, bkgdarr
     
end
