pro add_earthpoint_centroids

  ;;first read in the csv file with the dates and downlink pitch
  ;;angles from Effertz
  test = read_csv_jim('/Users/jkrick/irac_warm/Spitzer_pitch.csv')
  jd = julday(test.month, test.day, test.year)
  mjd = jd - 2400000.5

  for i = 0, n_elements(mjd) - 1 do print, jd(i), test.month(i), test.day(i), test.year(i)
  
  savenames = [ '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_01.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_02.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_03.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_04.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_05.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_06.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_07.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_08.sav']

  newsavenames = [ '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_01_earth.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_02_earth.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_03_earth.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_04_earth.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_05_earth.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_06_earth.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_07_earth.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_08_earth.sav']
  
  starts = 0
  stops =  0;n_elements(savenames) - 1
  for s = starts, stops do begin
     print, 'restoring ', savenames(s)
     restore, savenames(s)
     aorlist = planethash.keys()
     for n = 0,  n_elements(aorlist) - 1 do begin
        print, '--------------'
        print, 'working on ', aorlist(n), n_elements(planethash[aorlist(n)].xcen)
        bmjdarr = planethash[aorlist(n)].bmjdarr
        
        ;;take first bmjd
        bmjd0 = bmjdarr[0]
        
        ;;which value in csv file has the same mjd?
        bin = Value_Locate(mjd, bmjd0) ;;choose the closest one
        pe = test.pitch_earth(bin)
        
        ;;now I need to get this pitch angle value into the database
        ;;alter the planethash
        planethash[aorlist(n),'pitch_earth'] = pe

        ;;double check
        print, mjd(bin), bmjd0, pe
     endfor
     save, planethash, filename=newsavenames(s)
  endfor
  
  

end
