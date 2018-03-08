pro plot_preaor_pdf

  savename =  '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_08.sav'
  restore, savename
  plotnamey = '/Users/jkrick/external/irac_warm/trending/centroiding_ycen_05.pdf'
  plotnamex = '/Users/jkrick/external/irac_warm/trending/centroiding_xcen_05.pdf'
  plotnamef = '/Users/jkrick/external/irac_warm/trending/centroiding_flux_05.pdf'

  aorlist = planethash.keys()
  print, 'n_elements aorlist', n_elements(aorlist)
  ;;loadct, 5
  
  for n =0,30 - 1 do begin
     print, '--------------------------'
     timearr = planethash[aorlist(n)].timearr
     bmjdarr = planethash[aorlist(n)].bmjdarr
     ycen = planethash[aorlist(n)].ycen
     xcen = planethash[aorlist(n)].xcen
     flux = planethash[aorlist(n)].flux

     pid = planethash[aorlist(n)].pid
     prepid=planethash[aorlist(n)].prepid
     preaor = planethash[aorlist(n)].preaor
     np = n_elements(preaor) - 1
     print, 'n, aors ', n, ' ', aorlist(n), ' ',preaor
     print, 'aor min_dur',  planethash[aorlist(n)].min_dur
     print, 'pid', pid, prepid[-1]
     
     if n gt 0 and prepid[-1] eq pid then begin
        if planethash[aorlist(n)].min_dur gt 60. then begin
           
           ;;make sure we are not looking at the pre-aor itself.
           print, 'got a matched pid set', pid, ' ', aorlist(n)
           if planethash.HasKey(preaor(np)) eq 1 then begin
              ppmin_dur = planethash[preaor(np)].min_dur
              
              preaor = preaor[-1]
              alltime = [planethash[preaor].timearr, timearr]
              allycen = [planethash[preaor].ycen, ycen]
              allxcen = [planethash[preaor].xcen, xcen]
              allflux = [planethash[preaor].flux, flux]
              time0 = alltime(0)


              ;;shrink these arrays, don't need every point,
              ;;maybe every fifth point?
              num = findgen(n_elements(allxcen))
              good = where(num mod 5 eq 0)
              alltime = alltime(good)
              allycen = allycen(good)
              allxcen = allxcen(good)
              allflux = allflux(good)
              
              p1y = plot((alltime - time0)/60./60., allycen,'1s', sym_size = 0.3, /sym_filled,color = 'black',$
                         xtitle = 'time(hr)', ytitle = 'ycen', title = aorlist(n), $
                         yrange = [mean(ycen,/nan) -1.1, mean(ycen,/nan) +1.1], buffer = 1)
              p2y = plot((planethash[preaor].timearr - time0)/60./60.,planethash[preaor].ycen,'1s',sym_size = 0.3, $
                        /sym_filled,color = 'red', overplot = p1y)
              p1y.Save, plotnamey,/append,/bitmap 


 ;             p1x = plot((alltime - time0)/60./60., allxcen,'1s', sym_size = 0.3, /sym_filled,color = 'black',$
 ;                        xtitle = 'time(hr)', ytitle = 'xcen', title = aorlist(n), $
 ;                        yrange = [mean(xcen,/nan) -1.1, mean(xcen,/nan) +1.1], buffer = 1)
 ;             p2x = plot((planethash[preaor].timearr - time0)/60./60.,planethash[preaor].xcen,'1s',sym_size = 0.3, $
 ;                       /sym_filled,color = 'red', overplot = p1x)
 ;             p1x.Save, plotnamex,/append;,bitmap = 1
              
              
 ;             p1f = plot((alltime - time0)/60./60., allflux,'1s', sym_size = 0.3, /sym_filled,color = 'black',$
 ;                        xtitle = 'time(hr)', ytitle = 'flux', title = aorlist(n), $
 ;                        yrange = [mean(flux,/nan) -0.02, mean(flux,/nan) +0.02], buffer = 1)
 ;             p2f = plot((planethash[preaor].timearr - time0)/60./60.,planethash[preaor].flux,'1s',sym_size = 0.3, $
 ;                       /sym_filled,color = 'red', overplot = p1f)
 ;             p1f.Save, plotnamef,/append;,bitmap = 1

           endif
        endif
        
     endif else begin

        num = findgen(n_elements(xcen))
        good = where(num mod 5 eq 0)
        timearr = timearr(good)
        ycen = ycen(good)
        xcen = xcen(good)
        flux = flux(good)
        
        
        p1y = plot((timearr - timearr(0))/60./60., ycen,'1s', sym_size = 0.3, /sym_filled,color = 'black',$
                   xtitle = 'time(hr)', ytitle = 'ycen', title = aorlist(n), $
                   yrange = [mean(ycen,/nan) -1.1, mean(ycen,/nan) +1.1], buffer = 1)
        p1y.Save, plotnamey,/append,/bitmap 

;        p1x = plot((timearr - timearr(0))/60./60., xcen,'1s', sym_size = 0.3, /sym_filled,color = 'black',$
;                   xtitle = 'time(hr)', ytitle = 'xcen', title = aorlist(n), $
;                   yrange = [mean(ycen,/nan) -1.1, mean(ycen,/nan) +1.1], buffer = 1)
;        p1x.Save, plotnamex,/append;,bitmap = 1
        
;        p1f = plot((timearr - timearr(0))/60./60., flux,'1s', sym_size = 0.3, /sym_filled,color = 'black',$
;                   xtitle = 'time(hr)', ytitle = 'flux', title = aorlist(n), $
;                  yrange = [mean(flux,/nan) -0.02, mean(flux,/nan) +0.02], buffer = 1)
;        p1f.Save, plotnamef,/append;,bitmap = 1

     endelse
    
  endfor
;;  p1x.Save, plotnamex, /close
  p1y.Save, plotnamey, /close
;;  p1f.Save, plotnamef, /close


  
end


