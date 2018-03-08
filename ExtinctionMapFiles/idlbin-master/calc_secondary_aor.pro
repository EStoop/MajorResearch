pro calc_secondary_aor
;all taken from NSTED transit ephemeris calculator

t0 = 2455855.3919500001D  ;central transit
p = 4.411953D
t0 = t0 + 0.5*p

i = lindgen(40) + 128


obs_jd = t0 + i*p

;now subtract 2.565 hours to get the start of the observation
tstart = 4.07/24.
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
