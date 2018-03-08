pro scheduling, wknum , special_cals
;This program reads in csv files from the scheduling products and determines if the normal IRAC cals
;are scheduled approriately


dirloc = '/Users/jkrick/irac_warm/scheduling/'
readcol, strcompress(dirloc + 'wk' + string(wknum) + '_schedule.csv',/remove_all), title ,  start , stop,  duration , ra , dec ,  pid_name,format='(A,A,A, A, A, A, A)', skipline = 3,delim = ',',/silent

;test to make sure the table is getting read in correctly
;for i = 0, n_elements(title)-1  do print, 'duration: ', duration(i)

;now read in the brights:
readcol, strcompress(dirloc + 'wk' + string(wknum) + '_bright.csv',/remove_all),b_num,   b_start , b_key, b_pid, b_title, b_type,b_target, b_hit, b_ra, b_dec, b_dist, b_flux, b_unit, b_flag,format='(I3, A,I10, I6, A,A, A, A, F10.10,F10.12,F10.5, F10.4, A, I)', skipline = 3,delim = ',',count = nbrights, /silent


;test to make sure the table is getting read in correctly
;for i = 0, n_elements(b_num)-1  do print, 'b_ra: ', b_title(i)

;find the cals
calstar = where(strmatch(title, '*calstar*', /fold_case) eq 1, calstarcount)
;for c = 0, calstarcount -1 do print, title(calstar(c))

;seperate the ecliptic from cvz cals
dec_deg = Fix(strmid(dec,1,2))
;for i = 0, n_elements(title)-1  do print, 'dec: ', dec_deg(i)
cvz = where(dec_deg(calstar) gt 57 and dec_deg(calstar) lt 70, cvzcount)
e = where(dec_deg(calstar) le 57 or dec_deg(calstar) ge 70,ecount)
;for i = 0, n_elements(cvz)-1  do print, 'cvz title: ', title(calstar(cvz(i)))
;for i = 0, n_elements(e)-1  do print, 'e title: ', title(calstar(e(i)))
if ecount eq 0 then print, 'no ecliptic cals'
if cvzcount eq 0 then print, 'no cvz cals, must be week 2'
if cvzcount gt 0 then print, 'cvz cals, must be week 1'

print, '-----------------------'
print, 'The ecliptic cal we are using is:', title(calstar(e(0)))
print, 'Check that this is the same as the other week in campaign.'
print, '-----------------------'

;find the darks
skydrk = where(strmatch(title, 'skydrk*', /fold_case) eq 1, skydrkcount)
;for c = 0, skydrkcount -1 do print, title(skydrk(c))
if skydrkcount eq 0 then print, 'no sky darks'

;find the downlinks
downlink = where(strmatch(title, '*downlink*', /fold_case) eq 1, downlinkcount)
;for c = 0, downlinkcount -1 do print, title(downlink(c))

;find the flat
skyflt = where(strmatch(title, 'skyflt*', /fold_case) eq 1, skyfltcount)
;for c = 0, skyfltcount  - 1 do print, title(skyflt(c))
if skyfltcount eq 0 then print, 'no sky flats, must be week 2'
if skyfltcount eq 1 then print, 'skyflt, must be week 1'
if skyfltcount gt 1 then print, 'there is more than one skyflt observation, ', skyfltcount

;find the pmaps
pmap = where(strmatch(title, 'pmap*', /fold_case) eq 1, pmapcount)
;   print, '------------------'
;if pmapcount ne 2 then print, 'wrong number of pmaps: ', pmapcount
;if pmapcount eq 2 then print, 'found 2 pmaps which is good'
;endif
;     print, '------------------'
                         ;re-format duration keyword into floats that can be added together
dur_days = float(strmid(duration,0,1))
dur_hrs = float(strmid(duration, 2, 2))
dur_mins = float(strmid(duration, 5, 2))
dur_secs = float(strmid(duration, 8, 4))

;re-format start/stop keyword into floats 
stop_days = float(strmid(stop, 5, 3))
stop_hrs = float(strmid(stop, 9, 2))
stop_mins = float(strmid(stop, 12, 2))
stop_secs = float(strmid(stop, 15, 4))
start_days = float(strmid(start, 5, 3))
start_hrs = float(strmid(start, 9, 2))
start_mins = float(strmid(start, 12, 2))
start_secs = float(strmid(start, 15, 4))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;work with darks
trim_drk = strmid(title(skydrk), 7)
;for c = 0, n_elements(trim_drk) -1 do print, trim_drk(c)

; the darks really should be right after the downlink
;print, 'downlink, dark', downlink(0), skydrk(0)
if skydrk(0) - downlink(0) gt 1 then begin
   print, 'darks are not right after downlink', skydrk(0) - downlink(0) 
  print, '-----------------------'
endif

;are there the correct number of darks?
if n_elements(skydrk) ne 13 then begin
   print, 'wrong number of skydarks' , n_elements(skydrk)
   print, '-----------------------'
endif else begin
 ;print, 'number of skydrks', n_elements(skydrk)
; are the darks in order?
 if strmid(trim_drk(0),0,6) ne 's0p02s' then print, 'first skydark out of order'
 if strmid(trim_drk(1),0,5) ne 's0p1s' then print, 'second skydark out of order'
 if strmid(trim_drk(2),0,5) ne 's0p4s' then print, 'third skydark out of order'
 if strmid(trim_drk(3),0,3) ne 's2s' then print, 'fourth skydark out of order'
 if strmid(trim_drk(4),0,5) ne 'f0p4s' then print, 'fifth skydark out of order'
 if strmid(trim_drk(5),0,3) ne 'f2s' then print, 'sixth skydark out of order'
 if strmid(trim_drk(6),0,5) ne 'hdr6s' then print, 'seventh skydark out of order'
 if strmid(trim_drk(7),0,3) ne 'f6s' then print, 'eigth skydark out of order'
 if strmid(trim_drk(8),0,6) ne 'hdr12s' then print, 'ninth skydark out of order'
 if strmid(trim_drk(9),0,4) ne 'f12s' then print, 'tenth skydark out of order'
 if strmid(trim_drk(10),0,6) ne 'hdr30s' then print, 'eleventh skydark out of order'
 if strmid(trim_drk(11),0,4) ne 'f30s' then print, 'twelfth skydark out of order'
 ;if strmid(trim_drk(12),0,7) ne 'hdr100s' then print, 'thirteenth skydark out of order'
 if strmid(trim_drk(12),0,5) ne 'f100s' then print, 'thirteenth skydark out of order'

endelse

; there is a subset of darks which absolutely must be there even for dynamic scheduling.
;s0p1, s0p4s, f0p4s, f2s, f12s, f100s

dname = strmid(trim_drk, 0, 5)
d100 = where(dname eq 'f100s', d100count)
if d100count lt 1 then print, '100s dark not included'
d0p4 = where(dname eq 'f0p4s', d0p4count)
if d0p4count lt 1 then print, 'f0p4s dark not included'
s0p4 = where(dname eq 's0p4s', s0p4count)
if s0p4count lt 1 then print, 's0p4s dark not included'
s0p1 = where(dname eq 's0p1s', s0p1count)
if s0p1count lt 1 then print, 's0p1s dark not included'

dname = strmid(trim_drk, 0, 4)
d12 = where(dname eq 'f12s', d12count)
if d12count lt 1 then print, 'f12s dark not included'

dname = strmid(trim_drk, 0, 3)
d2 = where(dname eq 'f2s', d2count)
if d2count lt 1 then print, 'f2s dark not included'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;work with calstars

;what is the time difference between ecliptic calstars?
for i = 0, n_elements(calstar(e)) - 2 do begin
   elapsed_days = total(dur_days(calstar(e(i)):calstar(e(i+1))))
   elapsed_hrs =   total(dur_hrs(calstar(e(i)):calstar(e(i+1))))
   elapsed_mins =  total(dur_mins(calstar(e(i)):calstar(e(i+1))))
;;  print,calstar(i), calstar(i+1)
;  print, elapsed_days+ ( elapsed_hrs+ (elapsed_mins / 60.)) / 24. ,' days between calstar ', i, ' and ', i + 1
   if elapsed_days+ ( elapsed_hrs+ (elapsed_mins / 60.)) / 24.  lt 0.5 or elapsed_days+ ( elapsed_hrs+ (elapsed_mins / 60.)) / 24.  gt 1.5 then begin
      print, 'difference between ecliptic calstars should be between 0.5 and 1.5 days'
      print, elapsed_days+ ( elapsed_hrs+ (elapsed_mins / 60.)) / 24. ,' days between calstar ', i, ' and ', i + 1
      print, '-----------------------'
   endif
endfor

;what is the time difference between the first ecliptic and the skydarks?
;   elapsed_days = total(dur_days(skydrk(n_elements(skydrk) - 1):calstar(e(0))))
;   elapsed_hrs =   total(dur_hrs(skydrk(n_elements(skydrk) - 1):calstar(e(0))))
;   elapsed_mins =  total(dur_mins(skydrk(n_elements(skydrk) - 1):calstar(e(0))))
;   if elapsed_days+ ( elapsed_hrs+ (elapsed_mins / 60.)) / 24.  gt 0.5 then begin
;      print, 'difference between first ecliptic calstar and last dark should be  < 0.5 days'
;      print, elapsed_days+ ( elapsed_hrs+ (elapsed_mins / 60.)) / 24. ,' days between calstar and dark'
;     print, '-----------------------'
;   endif
   
;what is the time difference between the first ecliptic and the first cvz?
if cvzcount gt 0 then begin     ;week 1 (week 2 has no cvzcals)
;   print, 'test', e(0), cvz(0), calstar(e(0)), calstar(cvz(0)), dur_days(calstar(e(0))), dur_days(calstar(cvz(0)))
;   if dur_days(calstar(e(0))) le dur_days(calstar(cvz(0))) then begin
   if calstar(e(0)) le calstar(cvz(0)) then begin
      elapsed_days = total(dur_days(calstar(e(0)):calstar(cvz(0))))
      elapsed_hrs =   total(dur_hrs(calstar(e(0)):calstar(cvz(0))))
      elapsed_mins =  total(dur_mins(calstar(e(0)):calstar(cvz(0))))
   endif else begin
      elapsed_days = total(dur_days(calstar(cvz(0)):calstar(e(0))))
      elapsed_hrs =   total(dur_hrs(calstar(cvz(0)):calstar(e(0))))
      elapsed_mins =  total(dur_mins(calstar(cvz(0)):calstar(e(0))))
   endelse

     if elapsed_days+ ( elapsed_hrs+ (elapsed_mins / 60.)) / 24.  gt 0.5 then begin
         print, 'difference between first ecliptic calstar and first cvz should be  < 0.5 days'
         print, elapsed_days+ ( elapsed_hrs+ (elapsed_mins / 60.)) / 24. ,' days between ecliptic and cvz'
         print, '-----------------------'

      endif
   endif
   
;is there one ecliptic cal per downlink?
if downlinkcount / ecount gt 2 then begin
   print, 'there is less than one ecliptic cal per PAO ', downlinkcount, ecount
   print, '-----------------------'
endif



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;work with overall order

if skyfltcount eq 1 then begin  ;week 1

   print, 'assuming week 1'

;last skydrk should be before first skyflt
   time = diff_time(skyflt,skydrk(n_elements(skydrk) - 1), start_days,start_hrs, start_mins, start_secs)
   if  time gt 0.1 then begin
      print, time, ' days between last skydark and skyflt'
      print, '-----------------------'
   endif
   
;first ecliptic cal should execute within 12 hours of the start of the campaign
   time = diff_time(calstar(e(0)),skydrk(0), start_days,start_hrs, start_mins, start_secs)
   if  time gt 0.5 then begin
      print, time, ' days between start of campaign and first ecliptic cal should be < 0.5'
      print, '-----------------------'
   endif
   
;skyflt should be before ecliptic cal
   time = diff_time(calstar(e(0)),skyflt, start_days,start_hrs, start_mins, start_secs)
   if  time lt 0 then begin
      print, time, ' skyflat should be before ecliptic cal'
;      print, 'stop skyflt',abs(stop_days(skyflt) - start_days(calstar(e(0)))) 
      print, '-----------------------'
   endif

; ecliptic cal should be before cvzcals but within 12 hours
   time = diff_time(calstar(cvz(0)), calstar(e(0)), start_days,start_hrs, start_mins, start_secs)
   if  time gt 0.5 then begin
      print, time, ' days between first ecliptic cal and first cvzcal'
      print, '-----------------------'
   endif

;cvz cals must be in the correct order
   trim_cvz = strmid(title(calstar(cvz)), 13, 12)
;   for i = 0, n_elements(trim_cvz) - 1 do begin
;      print, 'trim_cvz',trim_cvz(i);, title(calstar(cvz(i)))
;   endfor
   
   if trim_cvz(0) ne '1812095_12s' then print, 'first cvzcal out of order'
   if trim_cvz(1) ne 'KF08T3_12s_c' then print, 'second cvzcal out of order'
   if trim_cvz(2) ne 'KF06T2_12s_c' then print, 'third cvzcal out of order'
   if trim_cvz(3) ne 'KF06T1_12s' then print, 'fourth cvzcal out of order'
   if trim_cvz(4) ne 'KF09T1_2s' then print, 'fifth cvzcal out of order'
   if trim_cvz(5) ne 'NPM1p60.0581' then print, 'sixth cvzcal out of order'
   if trim_cvz(6) ne 'NPM1p67.0536' then print, 'seventh cvzcal out of order'
   if trim_cvz(7) ne 'NPM1p68.0422' then print, 'eighth cvzcal out of order'
   if trim_cvz(8) ne 'HD165459_p4_' then print, 'ninth cvzcal out of order'
;   if trim_cvz(9) ne 'HD165459_p4s' then print, 'tenth cvzcal out of order'

endif                           ; end week 1

;--------------------------------------------------------------------
;--------------------------------------------------------------------

if skyfltcount eq 0 then begin   ;week 2
print, 'assuming week 2'
;last skydrk should be before first ecliptic 
;;  time = diff_time( calstar(e(0)),skydrk(13) , start_days,start_hrs, start_mins, start_secs)

;;  if  time gt 0.1 then begin
;;     print, time, ' days between last skydark and first ecliptic'
;;     print, '-----------------------'
;;  endif
   

;first ecliptic cal should execute within 12 hours of the start of the week
 time = diff_time(calstar(e(0)),skydrk(0), start_days,start_hrs, start_mins, start_secs)
   if  time gt 0.5 then begin
      print, time, ' days between start of campaign and first ecliptic cal should be < 0.5'
      print, '-----------------------'
   endif
  
endif ; end week 2



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;work with brights
;print, 'starting brights', b_num
;print, 'n bnum', n_elements(b_num)
;print, 'b_num(n_elements(b_num) - 1)',b_num(n_elements(b_num) - 1)
;if b_num(n_elements(b_num) - 1) ne 0 then begin ; make sure there are brights anyway
if nbrights gt 0 then begin ; make sure there are brights anyway
  print, 'there are brights', n_elements(b_num)
   b_start_days = float(strmid(b_start, 5, 3))
   b_start_hrs = float(strmid(b_start, 9, 2))
   b_start_mins = float(strmid(b_start, 12, 2))
   b_start_secs = float(strmid(b_start, 15, 4))

;flag anything that is brighter than 0.5 mags
   for nb = 0, n_elements(b_num) -1 do begin
      if b_flux(nb) le 0.5 then print, 'very bright source', b_flux(nb), ' at time ', b_start(nb)
   endfor




;check that there are no brights within 12 hours of the end of the week.
;compare last bright start time to last overall start time, difference should be greater than 12hrs
;print, 'working on bright star near end of week'
;print, 'n', n_elements(start) - 3, start_days(n_elements(start) - 1), start_days(358)
   time = diff_b_time(n_elements(start) - 3,n_elements(b_start)-1, start_days,start_hrs, start_mins, start_secs, b_start_days, b_start_hrs, b_start_mins, b_start_secs)

   if time lt 0.5 then begin
      print, 'there is a bright target within 12 hours of the end of the week', time, ' days', b_start(n_elements(b_start) - 1)
      print, '-----------------------'
   endif
   
;check that there are no brights within 4 hours of any cals  = darks, flats, calstars
   allcals = [skydrk, skyflt, calstar,pmap]
;for c = 0, n_elements(allcals) - 1 do print, title(allcals(c))

;somehow compare start time of allcalls with star time of brights
;maybe multiply the whole array of allcals star with brights star and find the minimum time?
   for j = 0, n_elements(b_start) - 1 do begin
      
;      print, 'allcals ', start(allcals)
;      print, 'brights', b_start(j)
;      delta_time =  start_days(allcals) -b_start_days(j)+ ( start_hrs(allcals) -b_start_hrs(j)+ ( start_mins(allcals) -b_start_mins(j) + (start_secs(allcals) -b_start_secs(j)) / 60.) / 60.) / 24.
      delta_time = diff_b_time(allcals,j, start_days,start_hrs, start_mins, start_secs, b_start_days, b_start_hrs, b_start_mins, b_start_secs)
      a = where(delta_time gt 0, dcount)
;      print, 'n_a',n_elements(a), 'a', a, 'dcount', dcount
      if  dcount gt 0 then begin
         if min(delta_time(a)) lt 0.16 then begin
            print, 'bright < 4 hrs before a cal: ', b_title(j), ' ', b_start(j), ' ', min(delta_time(a)), ' ', b_flux(j)
            print, '-----------------------'
         endif
      endif
      
      
                                ; print, min(allcals(j) - brights)
   endfor
   
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;check for special cals

if keyword_set(special_cals) then begin
   print, 'special cals included in this week:'
   readcol, special_cals, calname, format = 'A',/silent
   for i = 0, n_elements(calname) -1 do begin
      ;print, calname(i)
      for j = 0, n_elements(title) - 1 do begin
         if calname(i) eq title(j) then print, title(j), ' ', start(j), ' ', stop(j)
      endfor

   endfor
print, '-------------------'
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end



Function diff_time, first, second, start_day, start_hr, start_min, start_sec

  time =  ((start_day(first) - start_day(second)) * 24*60*60) + ((start_hr(first) -  start_hr(second)) * 60 * 60) + (  (start_min(first) - start_min(second)) * 60) + (start_sec(first) - start_sec(second) )

;now convert from seconds back into fractional days
  time = time / 60./60./24.

  return, time
end

;--------------------------------------------------------

Function diff_b_time, first, second, start_day, start_hr, start_min, start_sec, b_start_day, b_start_hr, b_start_min, b_start_sec

  time =  ((start_day(first) - b_start_day(second)) * 24*60*60) + ((start_hr(first) -  b_start_hr(second)) * 60 * 60) + (  (start_min(first) - b_start_min(second)) * 60) + (start_sec(first) - b_start_sec(second) )
;print, 'details', start_day(first),  b_start_day(second) , ((start_hr(first) -  b_start_hr(second)) * 60 * 60) ,  (  (start_min(first) - b_start_min(second)) * 60), (start_sec(first) - b_start_sec(second) )
;now convert from seconds back into fractional days
  time = time / 60./60./24.

  return, time
end

