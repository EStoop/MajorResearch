
pro count_dce
  
  savenames = [ '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_01.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_02.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_03.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_04.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_05.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_06.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_07.sav'] ;'/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/track_centroids_pixval.sav',
  
  starts = 0
  stops =  n_elements(savenames) - 1
  totaldcecount = 0L
  for s = starts, stops do begin
     
     restore, savenames(s)
     aorlist = planethash.keys()

     for n = 0,  n_elements(aorlist) - 1 do begin
        print, 'working on, ',n, ', ',  aorlist(n), ', ', n_elements(planethash[aorlist(n)].xcen)
        totaldcecount = totaldcecount + n_elements(planethash[aorlist(n)].xcen)
     endfor

  endfor

  print, 'totaldcecount: ', totaldcecount

end
