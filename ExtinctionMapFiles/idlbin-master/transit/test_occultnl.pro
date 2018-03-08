pro test_occultnl
plotquery=0
b0=dindgen(401)*0.003d0
rl=0.1d0
;start with these coefficients
;c1=1.d0
;c2=0.d0
;c3=0.d0
;c4=0.d0
;need to find real coefficients for limb darkening; use Claret 2000 model
;what to do about the bands not including Spitzer?
;temperature, gravity, and metallicity taken from nsted and Deeg et al., 2000 for K band
;transitnl, c1, c2, c3, c4, b1=11, teff = 5970, lg1 = 4.3, mh1 = .01

;ok from Claret itself for the same parameters as above (actually the sun) but at 3.6microns
c1=0.491	
c2=-0.441	
c3=0.446	
c4=-0.178

;and at 4.5micron
c1=0.554	
c2=-0.765	
c3=0.812	
c4=-0.314

;and at 8 micron
c1=0.466	
c2=-0.711	
c3=0.670	
c4=-0.230
occultnl,rl,c1,c2,c3,c4,b0,mulimbf,mulimb0,plotquery
plot,b0,mulimb0,ys=1,yr=[0.98d0,1.001d0],ytitle='!4l!3',xtitle='z'
oplot,b0,mulimbf,linestyle=1

end
