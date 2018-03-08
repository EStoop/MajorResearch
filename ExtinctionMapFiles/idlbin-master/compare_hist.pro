pro compare_hist
;  planetname = 'WASP-14b'
;  savefilename = [ '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/WASP-14b_phot_ch2_2.25000.sav',  '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/WASP-14b_phot_ch2_2.25000sdcorr.sav']
  planetname = 'HD158460'
  savefilename = [ '/Users/jkrick/irac_warm/pcrs_planets/HD158460/HD158460_phot_ch2_2.25000.sav',  '/Users/jkrick/irac_warm/pcrs_planets/HD158460/HD158460_phot_ch2_2.25000sdcorr.sav']

                              ;all the wasp-14b snaps
 planetinfo = create_planetinfo()
  aorname= planetinfo[planetname, 'aorname_ch2'] 
  
  for s = 0, n_elements(savefilename) -1 do begin
     restore, savefilename(s)
     print, savefilename(s)
     for a = 0, n_elements(aorname) - 1 do begin
        if a eq 0 then begin
           corrfluxarr= (planethash[aorname(a),'corrflux']) 
           fluxarr = (planethash[aorname(a),'flux']) 
        endif else begin
           corrfluxarr = [corrfluxarr, (planethash[aorname(a),'corrflux'])]
           fluxarr = [fluxarr, (planethash[aorname(a),'flux'])]
        endelse

     endfor  ; for each AOR
     
; I already throw out the first frame, but want to really remove the first five frames
     index = indgen(n_elements(corrfluxarr))
     goodframe = where(index mod 63 gt 3)
     corrfluxarr = corrfluxarr(goodframe)
     fluxarr = fluxarr(goodframe)


     plothist, corrfluxarr, xhist, yhist, bin = 0.0001, /noplot
     start = [mean(corrfluxarr,/nan), 0.0003, 9000]
     error = fltarr(n_elements(xhist)) + 1.0          ; uniform errors
     result= MPFITFUN('mygauss',xhist,yhist, error, start,/quiet) 
     print, 'sigma', result(1)
     ;normalize
     xhist = xhist / result(0)
     
     if s eq 0 then begin
        h = plot(xhist, yhist, thick = 2, xtitle = 'Aperture flux', ytitle = 'Number', title = planetname, $
           color = 'Navy', name = string('campaign dark corrflux '+string(result(1))), xrange = [0.98, 1.02]);/xlog)
        t= text(1.007, 1010, 'cmpgn corr' + string(result(1)), /data, color = 'Navy')
     endif else begin

        h2 = plot(xhist, yhist, color = 'red',thick = 2,/overplot, name = string('superdark corrflux' + string(result(1))))
        t= text(1.007, 940, 'sperdrk corr' + string(result(1)), /data, color = 'red')

     endelse

     plothist, fluxarr, xhist, yhist, bin = 0.0001, /noplot
     start = [mean(fluxarr,/nan), 0.0003, 6000]
     error = fltarr(n_elements(xhist)) + 1.0             ; uniform errors
     result= MPFITFUN('mygauss',xhist,yhist, error, start,/quiet)
     print, 'sigma', result(1)
   
                                ;normalize
     xhist = xhist / result(0)
     
     if s eq 0 then begin
        h3 = plot(xhist, yhist, color = 'black',thick = 2,/overplot, name = string('campaign flux' + string(result(1))))
        t= text(1.007, 870, 'cmpgn flux' + string(result(1)), /data, color = 'black')

     endif else begin
        h4 = plot(xhist, yhist, color = 'cyan',thick = 2,/overplot, name = string('superdark flux' + string(result(1))))
         t= text(1.007, 800, 'sperdrk flux' + string(result(1)), /data, color = 'cyan')
    endelse

  endfor

;  l = legend(target = [h, h2, h3, h4], position = [0.0595, 5500], /data, /auto_text_color, label = ['1', '2', '3', '4'])

end
