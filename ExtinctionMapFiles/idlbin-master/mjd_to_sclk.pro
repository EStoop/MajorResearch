function mjd_to_sclk, mjd

;conver MJD to JD
jd = mjd + 2400000.5

;convert to day month year
caldat, jd, month, day, year

;convert to day of year
dayofyear = ymd2dn(year, month, day)

;XXXdon't need to be precise for this particular case
hour = 0 & min = 0 & sec = 0.0

;convert to double precision sclk
sclk = double(0)
sclk = sec + (min*60.) + (hour*60.*60.) + ( dayofyear*24.*60.*60.) + ((year - 1980)*365.*24.*60.*60.)

;print, sclk, format = '(D14.2)'
return, sclk

end
