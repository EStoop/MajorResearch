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

