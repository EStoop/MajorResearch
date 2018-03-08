pro plot_preaor

  lastsave = find_newest_file()
 ;; lastsave = '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_07.sav'
  ;;savename = savename + '_pixval.sav'
  savecheck = file_test(lastsave)
  print, 'savecheck', savecheck
  if savecheck gt 0 then begin
     ;;first pull over the short_drift save file and figure out which
     ;;AORs need to be updated
     ;;do some string manipulation to get the correct savename from
     ;;lastsave
     shortsave = lastsave.substring(0, -30) + "short_drift_" + lastsave.substring(-6)
     restore, shortsave
     ;;aorlist variable has the aors which were already processed
     doneaorlist = doneaorlist;list(aorlist,/EXTRACT)
     doneshortdrift = totalshortdrift;short_driftarr
     doneslopedrift = totalslopedrift;slope_driftarr
     
     ;;now pull over the save file
     print, 'restoring big save file', lastsave
     restore, lastsave
     aorlist = planethash.keys()
     allaor = aorlist
     help, doneaorlist, aorlist
     ;;find those AORs in aorlist, but not done
     ;;first make these into strarrays
     stringaorlist = aorlist.ToArray(Type="STRING") 
     stringdone = doneaorlist.ToArray(Type="STRING")
     todo = cmset_op(stringaorlist, 'AND', /NOT2, stringdone)     
     ;;now set aorlist to todo so that we use only those left to be
     ;;done
     aorlist = todo
  endif else begin
     print, 'COULDN"T FIND SAVE FILE'
  endelse

  ;;; starting a new save file
  ;;  print, 'restoring big save file', lastsave
  ;;   restore, lastsave
  ;;   aorlist = planethash.keys()
  ;;   allaor = aorlist
 ;;;;
  
  plotname = '/Users/jkrick/external/irac_warm/trending/centroiding_ycen.pdf'

  ;;aorlist = planethash.keys()
  print, 'n_elements aorlist', n_elements(aorlist)
  ;;loadct, 5
  short_driftarr = fltarr(n_elements(aorlist))
  slope_driftarr = short_driftarr
  
  for n =0,n_elements(aorlist) - 1 do begin
     print, '--------------------------'
     timearr = planethash[aorlist(n)].timearr
     bmjdarr = planethash[aorlist(n)].bmjdarr
     ycen = planethash[aorlist(n)].ycen
     pid = planethash[aorlist(n)].pid
     prepid=planethash[aorlist(n)].prepid
     preaor = planethash[aorlist(n)].preaor
     np = n_elements(preaor) - 1
     print, 'n, aors ', n, ', ', aorlist(n), ', ',preaor
     print, 'aor min_dur',  planethash[aorlist(n)].min_dur
     print, 'pid', pid, prepid[-1]
     help, planethash[aorlist(n)].min_dur
     if n gt 0 and prepid[-1] eq pid then begin
        mindur = planethash[aorlist(n)].min_dur
        if mindur(0) gt 60. then begin
           
           ;;make sure we are not looking at the pre-aor itself.
           print, 'got a matched pid set', pid, ' ', aorlist(n)
           if planethash.HasKey(preaor(np)) eq 1 then begin
              ppmin_dur = planethash[preaor(np)].min_dur
           
              preaor = preaor[-1]
              alltime = [planethash[preaor].timearr, timearr]
              allycen = [planethash[preaor].ycen, ycen]
              time0 = alltime(0)
              
              plot,(alltime - time0)/60./60., allycen,psym = 3,xtitle = 'time(hr)', ytitle = 'ycen', title = aorlist(n), yrange = [mean(ycen,/nan) -1.1, mean(ycen,/nan) +1.1], ystyle = 1, xrange = [0, 5], xstyle = 1
              oplot,(planethash[preaor].timearr - time0)/60./60.,planethash[preaor].ycen,psym=1, color = 80 ;,color = 'red'
              
              ;;p1 = plot((alltime - time0)/60./60., allycen,'1s', sym_size = 0.3, /sym_filled,color = 'black',xtitle = 'time(hr)', ytitle = 'ycen', title = aorlist(n), yrange = [mean(ycen,/nan) -1.1, mean(ycen,/nan) +1.1])
              ;;p2 = plot((planethash[preaor].timearr - time0)/60./60.,planethash[preaor].ycen,'1s',sym_size = 0.3, /sym_filled,color = 'red', overplot = p1)
              ;;p1.Save, plotname,/append,bitmap = 1
              
              CURSOR, X1, Y1, /DOWN,/DATA
              CURSOR, X2, Y2, /DOWN,/DATA
              ;;directions
              ;;negative X1 = pre-aor, or can't measure the
              ;;value due to badish data = value of short and slope = NAN
              ;;negative X2-X1 = no short term drift = 0's
              
              print, 'cursor position pre-aor', X1, Y1, X2, Y2
              if X1 gt 0 then begin  ;;if nagative then is a preaor so keep at nan
                 sd= X2 - X1
                 slope = (Y2-Y1) / (X2 - X1) ; preserve direction - could be negative
                 if X2 - X1 lt 0 then begin
                    sd = 0 
                    slope = 0
                 endif
                 print, 'slope', slope
                 planethash[aorlist(n)].short_drift = sd
                 planethash[aorlist(n)].slope_drift = slope
                 short_driftarr(n) = sd
                 slope_driftarr(n) = slope
                 
              endif else begin
                 
               planethash[aorlist(n)].short_drift = alog10(-1)
                 planethash[aorlist(n)].slope_drift = alog10(-1)
                 short_driftarr(n) = alog10(-1)
                 slope_driftarr(n) = alog10(-1)
              endelse
              
           endif
        endif
        
     endif else begin
        
        ;;just plot single AOR
        plot,(timearr - timearr(0))/60./60., ycen,psym = 3,xtitle = 'time(hr)', ytitle = 'ycen', title = aorlist(n), yrange = [mean(ycen,/nan) -1.1, mean(ycen,/nan) +1.1], ystyle = 1, xrange = [0, 5], xstyle = 1
        CURSOR, X1, Y1, /DOWN,/DATA
        CURSOR, X2, Y2, /DOWN,/DATA
        print, 'cursor position single', X1, Y1, X2, Y2

        ;;p1 = plot((timearr - timearr(0))/60./60., ycen,'1s', sym_size = 0.3, /sym_filled,color = 'black',xtitle = 'time(hr)', ytitle = 'ycen', title = aorlist(n), yrange = [mean(ycen,/nan) -1.1, mean(ycen,/nan) +1.1])
        
        ;;p1.Save, plotname,/append,bitmap = 1
 
        if X1 gt 0 then begin  ;;if nagative then is a preaor so keep at nan
           sd= X2 - X1
           slope = (Y2-Y1) / (X2 - X1) ; preserve direction - could be negative
           if X2 - X1 lt 0 then begin
              sd = 0 
              slope = 0
           endif
           print, 'slope', slope
           planethash[aorlist(n)].short_drift = sd
           planethash[aorlist(n)].slope_drift = slope
        endif
        
     endelse
     
     
  endfor
  ;;p1.Save, plotname, /close
  save, planethash, filename=savename
  ;;oops but this only saves the new ones, not the whole thing.
  ;;concatenate those already done with the new ones
  test = list(aorlist, /extract)
  doneaorlist.add, test, /extract
  totalshortdrift = [doneshortdrift, short_driftarr]
  totalslopedrift = [doneslopedrift, slope_driftarr]

  ;;  starting a new save file
  ;;doneaorlist = list(aorlist, /extract)
  ;;totalshortdrift = short_driftarr
  ;;totalslopedrift =  slope_driftarr

  ;;

  
  save, doneaorlist, totalshortdrift, totalslopedrift, filename = '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/short_drift_08.sav'
end


       ;; plotname = '/Users/jkrick/external/irac_warm/trending/centroiding_ycen.pdf'
