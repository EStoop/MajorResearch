pro track_starname

  savenames = [ '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_01.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_02.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_03.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_04.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_05.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_06.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_07.sav'] ;'/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/track_centroids_pixval.sav',
  
  starts = 0
  stops = n_elements(savenames) - 1
  totalaorcount = 0
  for s = starts, stops do begin
     
     restore, savenames(s)
     aorlist = planethash.keys()
     
     if s eq starts then begin
        sigmax = fltarr(2* (stops + 1 - starts) *n_elements(aorlist))*alog10(-1) ; two is the fudge factor to make sure arrays are large enough
        print, n_elements(sigmax)
        aorname = lonarr(n_elements(sigmax)) ;strarr(n_elements(sigmax))
        starname = strarr(n_elements(sigmax))
     endif
     for n = 0,  n_elements(aorlist) - 1 do begin
        ;print, '--------------'
        ;print, 'working on ',n, ' ', totalaorcount, ' ', aorlist(n), n_elements(planethash[aorlist(n)].xcen)
        starnamestr = planethash[aorlist(n)].starname
        
        aorname[totalaorcount] = aorlist(n)
        ;print, 'aorname should be long', aorname[totalaorcount]
        
        
        starname[totalaorcount] = starnamestr(0)
        ;print, 'starname', starname[totalaorcount]
        totalaorcount++  ;;keep track of the total number of aors processed
     endfor
  endfor
  starname = starname[0:totalaorcount - 1]

  ;;remove pre-aors
  a = where(starname.Contains('preaor'), complement = b)
  
  print, 'starnames', starname(b)
end


