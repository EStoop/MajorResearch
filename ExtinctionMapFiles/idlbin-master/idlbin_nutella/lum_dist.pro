function lum_dist, z

;COMPUTES THE LUMINOSITY DISTANCE IN PARSECS TO AN OBJECT AT REDSHIFT Z

return, prop_motion_dist(z) * (1.0 + z)

end
