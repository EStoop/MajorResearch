pro calc_wasp16_obs

t0 =  2455956.61322D  ;central transit
p = 3.118601D
t0 = t0 + 0.5*p

i = lindgen(20) + 190


obs_jd = t0 + i*p

;now subtract 2.565 hours to get the start of the observation
tstart = 2.565/24.
obs_jd = obs_jd - tstart

;now subtract 15 min to get to the start of the vis window
tw = 0.25/24.
obs_jd = obs_jd - tw

;print, obs_jd, format = '(F0)'

;and convert JD to UT date
caldat, obs_jd, month, day, year, hour, minute, second

for j = 0, n_elements(month) - 1 do begin
   print, month(j), day(j), year(j), hour(j), minute(j), second(j)
endfor



end
