pro measure_drift
;;XXX need to measure duration and slope of the short term drift with
;;and without pre_aors
  savename = '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval.sav'
  restore, savename
  aorlist = planethash.keys()
  print, 'n_elements aorlist', n_elements(aorlist)
  short_drift = fltarr(n_elements(aorlist))
  loadct, 5
  
  for n =120,125 do begin
     print, '--------------------------'
     timearr = planethash[aorlist(n)].timearr
     bmjdarr = planethash[aorlist(n)].bmjdarr
     ycen = planethash[aorlist(n)].ycen
     pid = planethash[aorlist(n)].pid
     prepid=planethash[aorlist(n)].prepid
     preaor = planethash[aorlist(n)].preaor
     ;;print, 'timearrr 0', n_elements(timearr)
     if prepid[-1] eq pid then begin
        print, 'got a matched pair', n
        preaor = preaor[-1]
        alltime = [planethash[preaor].timearr, timearr]
        allycen = [planethash[preaor].ycen, ycen]
        time0 = alltime(0)
                                ;alltime = (alltime - time0)/60./60. ;
                                ;                     now in hours
                                ;                     instead WINDOW,
                                ;                     0, XSIZE=400,
                                ;                     YSIZE=400of sclk
        ;;print, 'n pre', n_elements(planethash[preaor].timearr), 'n science', n_elements(timearr)
       ;; print, 'alltime 0', n_elements(alltime)
        ;print, 'science y', ycen[0:10]
        ;print, 'pre y', planethash[preaor].ycen
        plot,(alltime - time0)/60./60., allycen,psym = 3,xtitle = 'time(hr)', ytitle = 'ycen', title = aorlist(n), yrange = [mean(ycen,/nan) -1.1, mean(ycen,/nan) +1.1], ystyle = 1
        oplot,(planethash[preaor].timearr - time0)/60./60.,planethash[preaor].ycen,psym=1, color = 80;,color = 'red'
                                ;p1.Save, plotname,/append,bitmap = 1
        CURSOR, x1,y1, /DOWN,/DATA
        print, 'cursor position', x1,y1
        planethash[aorlist(n)].short_drift = x1
     endif else begin
        ;;just plot the single AOR
        time0 = timearr(0)
        plot,(timearr - time0)/60./60., ycen,psym = 3,xtitle = 'time(hr)', ytitle = 'ycen', title = aorlist(n), yrange = [mean(ycen,/nan) -1.1, mean(ycen,/nan) +1.1], ystyle = 1
                CURSOR, x1,y1, /DOWN,/DATA
        print, 'cursor position single AOR', x1,y1
        planethash[aorlist(n)].short_drift = x1

     endelse
     
    
  endfor
  ;p1.Save, plotname, /close
  save, planethash, filename=savename

end


;;       p1 = plot((alltime - time0)/60./60., allycen,'1s', sym_size = 0.3, /sym_filled,color = 'black',xtitle = 'time(hr)', ytitle = 'ycen', title = aorlist(n), yrange = [mean(ycen,/nan) -1.5, mean(ycen,/nan) +1.5])
;;        p2 = plot((planethash[preaor].timearr - time0)/60./60.,planethash[preaor].ycen,'1s',sym_size = 0.3, /sym_filled,color = 'red', overplot = p1)
       ;; plotname = '/Users/jkrick/external/irac_warm/trending/centroiding_ycen.pdf'
