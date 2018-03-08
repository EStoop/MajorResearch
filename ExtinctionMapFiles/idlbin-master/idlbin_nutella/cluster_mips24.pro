pro cluster_mips24

common sharevariables, candra, canddec, sep, objectnew, cluster24, clusterall, dist24, distall, allmips, allgals, allcluster, distallmips, distallgals, distallcluster


restore, '/Users/jkrick/idlbin/objectnew.sav'

!P.multi=[0,1,1]
ps_open, filename='/Users/jkrick/nep/clusters/mips24/distances.ps',/portrait,/square,/color, xsize=4, ysize=4
!P.charthick=3
!P.thick=3
!X.thick=3
!Y.thick=31
;-------------------------------------------------------------------------------
;where are the 24 micron sources in the cluster areas
;-------------------------------------------------------------------------------

canddec = [69.045052,  69.06851,  69.087766]
candra = [264.68365, 264.89228, 264.83053]
sep = 0.019


;pick the sample
junk = cluster_mips24_sample()

  
;plot cumulative distribution of distance from center for mips detections
;need to sort on distance
dist24 = dist24*60.
sortdist24 = dist24(sort(dist24))
N1 = n_elements(sortdist24)
f1 = (findgen(N1) + 1.)/ N1
plot, sortdist24, f1, xtitle='Distance from Cluster Center in Arcminutes', ytitle='Cumulative Fraction';, psym = 10

;look at cluster members
distall = distall*60.
sortdistall = distall(sort(distall))
N1all = n_elements(sortdistall)
f1all = (findgen(N1all) + 1.)/ N1all
;oplot, sortdistall, f1all, psym = 10, linestyle = 5

;-------------------------------------------------------------------------------
;get together some random background regions
;-------------------------------------------------------------------------------
nrand = 30
seed = 423
randx = randomu(seed, nrand)
randy = randomu(seed, nrand)

randx = randx * 9500. + 6800.
randy = randy * 13300. + 3100.

acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits') ;
xyad, acshead, randx, randy, randra, randdec
allranddist = fltarr(nrand*100)
d = 0
for cand=0, n_elements(randdec) - 1 do begin ;for each random area
   
   goodrand = where(sphdist(objectnew.ra, objectnew.dec, randra(cand), randdec(cand), /degrees) lt sep and objectnew.mips24mag gt 0 and objectnew.mips24mag lt 90 )
   
   for c = 0, n_elements(goodrand) - 1 do begin
      allranddist[d] = sphdist(objectnew[goodrand[c]].ra, objectnew[goodrand[c]].dec, randra(cand), randdec(cand), /degrees) 
      d = d + 1
   endfor


endfor
allranddist = allranddist[0:d-1]*60.

;sort allrand on dist
allranddist = allranddist(sort(allranddist))
N2 = n_elements(allranddist)
f2 = (findgen(N2) + 1.)/ N2
oplot, allranddist, f2, psym = 10, linestyle=1

;determine the ks-test distance and probability
kstwo, allranddist, sortdist24, ksdist, prob
print, 'ksdist, prob', ksdist, prob
;xyouts, .003*60.,0.85, string(ksdist) + string(prob)

;now bin distances
;5 bins, < .0038, .0038 < .0076, .0076<.0114, .0114 < .0152,
;.0152<.019

x = [.0019,.0057,.0095,.0133,.0171]
y = [n_elements(where(dist24 lt .0038)),n_elements(where(dist24 ge .0038 and dist24 lt .0076)),n_elements(where(dist24 ge .0076 and dist24 lt .0114)),n_elements(where(dist24 ge .0114 and dist24 lt .0152)),n_elements(where(dist24 ge .152 and dist24 lt .019))]

y2 = [n_elements(where(allranddist lt .0038)),n_elements(where(allranddist ge .0038 and allranddist lt .0076)),n_elements(where(allranddist ge .0076 and allranddist lt .0114)),n_elements(where(allranddist ge .0114 and allranddist lt .0152)),n_elements(where(allranddist ge .152 and allranddist lt .019))]

;plot, x, y, psym = 10
;oplot, x, y2,psym=10, linestyle=1

;plothist, dist, bin=0.019*60./10., xtitle='distance from cluster center in arcminutes', ytitle='N', xrange=[0,0.02*60.], xstyle=1
;plothist, allranddist, bin=0.019*60./10., /overplot, linestyle=1
print, n_elements(dist24), n_elements(allranddist) ,'nelements dist, rand'


;-------------------------------------------------------------------------------
;what is the luminosity of those sources
;-------------------------------------------------------------------------------

ps_close, /noprint,/noid


end

;plothist, objectnew[good].rmag, xrange=[10,30]

;plothyperz, good, strcompress( '/Users/jkrick/nep/clusters/mips24/mips24_'+string(i)+'.ps')


;openw, outlun, strcompress('/Users/jkrick/nep/clusters/mips24/cluster'+string(i)+'.reg'), /get_lun
;printf, outlun, 'fk5'
;for rc=0, n_elements(good) -1 do  printf, outlun, 'circle( ', objectnew[good[rc]].ra, objectnew[good[rc]].dec, ' 3")'
;close, outlun
;free_lun, outlun


;openw, outlun2, strcompress('/Users/jkrick/nep/clusters/mips24/cluster'+string(i)+'.txt'), /get_lun
;printf, outlun, 'number, ra, dec, zphot, rmag, acsmag, irac1mag, mips24flux, acsmatchdist, comment'
;for rc=0, n_elements(good) -1 do  printf, outlun2,  format='(I10,F10.5,F10.5,F10.2,F10.2,F10.2,F10.2,F10.2)', good[rc], objectnew[good[rc]].ra, objectnew[good[rc]].dec, objectnew[good[rc]].zphot, objectnew[good[rc]].rmag,objectnew[good[rc]].acsmag, objectnew[good[rc]].irac1mag,objectnew[good[rc]].mips24flux
;close, outlun2
;free_lun, outlun2
