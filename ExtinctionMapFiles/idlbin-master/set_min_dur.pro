pro set_min_dur

  savename =  '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/track_centroids_pixval.sav'
  restore, savename
  aorlist = planethash.keys()
  print, 'n_elements aorlist', n_elements(aorlist)

  readcol, '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/out_centroiding_allWarmMissionToDate.txt', aorname,pid,startUTC,campaign,min_dur_cat,RA,Dec,readoutfull,datacollect36,datacollect45, format = '(L10, L10, A, A, F10.4, D10.6, D10.6,A,A,A )', delimiter='|', skipline =7520

  for n =1, n_elements(aorlist) - 1 do begin
     a = where(aorname eq aorlist(n))
     print, aorlist(n), aorname(a)
     premin_dur = min_dur_cat(a-1)
    ;; print, 'a, min_dur, premin_dur', a, min_dur_cat(a), premin_dur
     planethash[aorlist(n)].min_dur = min_dur_cat(a)
  endfor
  
  save, planethash, filename=savename

end
