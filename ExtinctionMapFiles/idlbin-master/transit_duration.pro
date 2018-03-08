function transit_duration, period, r_star, r_planet, b, a
;see below website for figures and derivation
;http://www.paulanthonywilson.com/exoplanets/exoplanet-transits-transmission-spectroscopy/the-transit-light-curve/
t_dur = (period/!Pi) * asin(sqrt((r_star + r_planet)^2 - b^2) / a)

return, t_dur

end
