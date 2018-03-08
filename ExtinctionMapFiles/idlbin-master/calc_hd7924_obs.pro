pro calc_hd7924_obs

t0 = 2455866.785194D
p = 5.3978D

i = lindgen(20) + 71

obs_jd = t0 + i*p

;now subtract 4 hours to get the start of the observation
tstart = 4.0/24.
obs_jd = obs_jd - tstart

;now subtract 15 min to get to the start of the vis window
tw = 0.25/24.
obs_jd = obs_jd - tw

;print, obs_jd, format = '(F0)'

;and convert JD to UT date
caldat, obs_jd, month, day, year, hour, minute, second

for j = 0, n_elements(month) - 1 do begin
;   print, month(j), day(j), year(j), hour(j), minute(j), second(j)
endfor

;junk for another planet
 o =  [2456383.65439D, 2456401.85028D, 2456420.04617D,  2456602.00510D,   2456620.20099D,   2456638.39688D,   2456656.59277D]
tw = 12.25/24.
o = o - tw

caldat, o, month, day, year, hour, minute, second

for j = 0, n_elements(month) - 1 do begin
   print, month(j), day(j), year(j), hour(j), minute(j), second(j)
endfor


end
