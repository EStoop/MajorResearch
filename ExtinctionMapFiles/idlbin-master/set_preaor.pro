pro set_preaor
save_names = ['/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_01.sav', '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_02.sav', '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_03.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_04.sav', '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_05.sav', '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_06.sav', '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_07.sav', '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_08.sav']


for s = 0, n_elements(save_names) -1 do begin
   ;;does the file exist?
   ft = file_test(save_names(s))
   if ft gt 0 then begin
      restore, save_names(s)
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
      
      
      ;;save
      ;; savename = '/Users/jkrick/track_centroids.sav'
      print, 'saving', save_names(s)
      save, planethash, filename=save_names(s)
   endif
   
endfor

end
