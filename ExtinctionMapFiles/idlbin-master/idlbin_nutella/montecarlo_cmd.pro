pro montecarlo_cmd
restore, '/Users/jkrick/idlbin/object.old.sav'
!p.multi = [0, 1,2]


; make a similar cmd for all galaxies with irac and acs detections
all = where( objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.acsmag gt 0 and objectnew.acsmag lt 90 )
;-----------------------------------------------------------
;do some monte carlo calculation of the probability of getting the
;same amount of galaxies on the red sequence by taking 711 random
;galaxies instead of cluster galaxies....is some sort of proof of the
;red sequence and therefore the cluster existance.
;-----------------------------------------------------------
m=-0.12
ball = 5.1
rcswidth = 0.5

totcountarr = fltarr(100)
for niter = 0, n_elements(totcountarr) - 1 do begin
;want 711 random galaxies from this distribution
   randselect = n_elements(all)*randomu(S, 711)


;how many of the randselect are within the red sequence as defined in
;cluster.pro?
   totcount = 0
   for count = 0.D, n_elements(randselect) - 1 do begin
      limithigh = (m*objectnew[all[randselect[count]]].irac1mag) + ball +rcswidth ; 2*sig1;0.5 
      limitlow = (m*objectnew[all[randselect[count]]].irac1mag) + ball - rcswidth ;2*sig1;0.5 
      if (objectnew[all[randselect[count]]].acsmag - objectnew[all[randselect[count]]].irac1mag LT limithigh ) AND ( objectnew[all[randselect[count]]].acsmag - objectnew[all[randselect[count]]].irac1mag GT limitlow) then begin
                                ; add this member galaxy to a larger array of all member galaxies
         totcount = totcount + 1
      endif
   endfor
   totcountarr[niter]= totcount
endfor

   print, "mean, stddev", mean(totcountarr), stddev(totcountarr)

   nsigma = (275. - mean(totcountarr)) / stddev(totcountarr)
   print, "nsigma", nsigma
   PRINT, 'probability', (GAUSS_PDF( (275. - mean(totcountarr))/(stddev(totcountarr)) )  ) 

end
