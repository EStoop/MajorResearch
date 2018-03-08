pro test_read
  ;;figure out where the exoplanet stares are.

 
;;  readcol, '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/weeks456_553_666_749_768.txt', aorname,pid,startUTC,campaign,min_dur,RA,Dec,readoutfull,datacollect36,datacollect45, format = '(L10, I10, A, A, F10.4, D10.6, D10.6,A,A,A )', delimiter='|', skipline = 1

  
  readcol, '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/out_centroiding_allWarmMissionToDate.txt', aorname,pid,startUTC,campaign,min_dur,RA,Dec,readoutfull,datacollect36,datacollect45, format = '(L10, L10, A, A, F10.4, D10.6, D10.6,A,A,A )', delimiter='|', skipline =7520
  start_year = fix(strmid(startUTC, 0,4))
  start_month = fix(strmid(startUTC, 5,2))
  start_day = fix(strmid(startUTC, 8,2))
  start_hr = fix(strmid(startUTC, 11, 2))
  start_min = fix(strmid(startUTC, 14, 2))
  ;;convert times to JD in order to compare later
  start_jd = julday(start_month, start_day, start_year, start_hr, start_min)

   
    ;;convert 3 digit campaign number to string
  campaign_name = strcompress('IRAC0'+ campaign + '00',/remove_all)

  ;;check if full or subarray
  naxis = strarr(n_elements(readoutfull))
  for r = 0, n_elements(readoutfull) - 1 do begin
     if readoutfull[r] eq 't' then naxis[r] ='2' else naxis[r] ='3'
  endfor

  ;;convert aorname from long to a string for compatability with
  ;;phot_exoplanet.pro
  aorname = string(aorname)
  aorname = strtrim(aorname, 1)  ;remove leading whitespaces

  ;;make an array which holds the previous 1 hour of AORnames, ra,
  ;;dec, and start times
  ;;takes too long, just keep the previous 5 AORs.
  preaor = strarr(n_elements(start_jd), 5) ;;list(length=n_elements(start_jd));ptrarr(n_elements(start_jd),/allocate_heap)
  prera = fltarr(n_elements(start_jd), 5)
  predec = fltarr(n_elements(start_jd), 5)
  prejd = fltarr(n_elements(start_jd), 5)
  prepid =  fltarr(n_elements(start_jd), 5)
  pre36= strarr(n_elements(start_jd), 5)
  pre45= strarr(n_elements(start_jd), 5)

  preaor[0,*] = '0'
  prera[0,*] = 0               
  predec[0,*] = 0               
  prejd[0,*] = 0               
  prepid[0,*] = 0               
  pre36[0,*] = 0
  pre45[0,*] = 0
  for j = 5, n_elements(start_jd) - 1 do begin
     preaor[j,*] = aorname[j-5:j-1]
     prera[j,*] = ra[j-5:j-1]
     predec[j,*] = dec[j-5:j-1]
     prejd[j,*] = start_jd[j-5:j-1]
     prepid[j,*] = pid[j-5:j-1]
     pre36[j,*] = datacollect36[j-5:j-1]
     pre45[j,*] = datacollect45[j-5:j-1]

  endfor

 
  ;;pare down to only include those longer than some threshold
  thresh = 130                 ;two hours seems appropriate for now
  longstare = where(min_dur gt thresh, n_longstare, complement = tooshort)
  print, 'number of long AORs', n_longstare
  if n_longstare gt 0 then remove, tooshort, pid, campaign_name, start_jd, aorname, naxis, datacollect36, datacollect45, min_dur
  ;;need to also remove those from the multi-dimensional arrays which
  ;;remove can't handle
  preaor = preaor[longstare, *]
  prera = prera[longstare, *]
  predec = predec[longstare, *]
  prejd = prejd[longstare, *]
  prepid = prepid[longstare, *]
  pre36 = pre36[longstare, *]
  pre45 = pre45[longstare, *]

  ;;and want to make sure they really are stares and not dithers
  ;;so compare to Elena's list of exoplanet and BD pids and
  ;;remove those pid's not on the list
  readcol, '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/staring_pids.txt',  spid, name, format = '(L5, A)'
  starepid = intarr(n_elements(pid))
  for i = 0, n_elements(pid) -1 do begin
     a = where(pid(i) eq spid, good)
     starepid[i] = good
  endfor
 ;; print, starepid
  bad = where(starepid lt 1, nbad, complement = goodstare)
  badpids = pid(bad)
  b = badpids[UNIQ(badpids, SORT(badpids))]
  
  if nbad gt 0 then begin
     remove, bad, pid, campaign_name, start_jd, aorname, naxis, datacollect36, datacollect45, min_dur
     preaor = preaor[goodstare, *]
     prera = prera[goodstare, *]
     predec = predec[goodstare, *]
     prejd = prejd[goodstare, *]
     prepid = prepid[goodstare, *]
     pre36 = pre36[goodstare, *]
     pre45 = pre45[goodstare, *]
     
  endif
  print, 'number of staring, long AORs' ,n_elements(start_jd)
  ;;look at how much total data/time we are talking about here
     ;;do a few calculations
     ;;print, TargetName(longstare)
     ;print, 'nele', n_elements(planets), n_longstare

     total_min = total(min_dur)
     total_sec = total_min *60.
     
     ;;assume 1 2s sub image is 128s and 312K
     total_frame = floor(total_sec / 128.)
     total_size = total_frame *312
     print, 'total time', total_min / 60./24.
     print, 'total size', total_size /1E6, 'G'
     print, 'total number of AORs', n_elements(pid)


;;figure out time distribution
     caldat, start_jd, month, day, year
     for i = 0, n_elements(start_jd) - 1 do print, i, month(i), day(i), year(i)
  
  
end
  ;;readcol, '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/allstares_07_28_16.txt', TargetName, RA,    Dec,   pid,   AOT, min_dur, start_day, start_time,   AOR_key,   AOR_LABEL, format = '(A, A, A, I10, A, F10.2, A, A, L10, A )'


;;     for j = 2, n_elements(start_jd) -1 do begin
;;     if j mod 1000 eq 0 then print, 'j', j
;;     pre = where( min_dur gt 20 and start_jd gt start_jd(j) - 0.1 and start_jd lt start_jd(j), n_pre)

;;    if n_pre gt 0 then begin
;;        preaor[j] = aorname(pre)
;;        prera[j] = ra(pre)
;;        predec[j] = dec(pre)
;;        prejd[j] = start_jd(pre)
;;        prepid[j] = pid(pre)
;;        pre36[j] = datacollect36(pre)
;;        pre45[j] = datacollect45(pre)
;;     endif else begin
        ;;just use the previous single AOR - just means the one before
        ;;                                   it is longer than 1 hour.
;;        preaor[j] = aorname(j-1)
;;        prera[j] = ra(j-1)
;;        predec[j] = dec(j-1)
;;        prejd[j] = start_jd(j-1)
;;        prepid[j] = pid(j-1)
;;        pre36[j] = datacollect36(j-1)
;;        pre45[j] = datacollect45(j-1)
;;        
;;     endelse
     
;;  endfor


 ;;using the webpage that Elena maintains, but this doesn't
  ;;have campaign name which I need to access the data in the archive
;;  planets = webget("http://ssc.spitzer.caltech.edu/warmmission/scheduling/observinglogs/extrasolarplanetsWarm.txt")
;;  planets = planets.text[22:*]
;;  aorname = strmid(planets, 106, 8)  ;; will be screwed if format of this file changes
;;  min_dur = strmid(planets, 71, 6)
;;  start_time = strmid(planets, 79, 21)

