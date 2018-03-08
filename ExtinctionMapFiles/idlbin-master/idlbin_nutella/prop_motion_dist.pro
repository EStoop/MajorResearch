function prop_motion_dist, z

;COMPUTES THE PROPER MOTION DISTANCE IN PARSECS TO AN OBJECT AT REDSHIFT Z

return, (3.0 * 10.0^10.0 * qsimp('prop_motion_dist_func', 1.0/(1.0+z), 1.0)) / (3.0 * 10.0^18.0)

end
