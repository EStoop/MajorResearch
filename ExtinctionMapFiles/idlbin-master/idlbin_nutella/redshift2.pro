pro redshift2

restore, '/Users/jkrick/idlbin/objectnew.sav'

!p.multi = [0,0, 1]

;from Matt Ashby

;    Have you thought about the following exercise.  If you look at sources which are *RED* in [3.6]-[4.5] and blue in [5.8]-[8.0], i.e., [3.6]-[4.5] *>* 0 and [5.8]-[8.0] < 0 (where [] denotes Vega-relative magnitude) that are also detected by MIPS 24 (a band which you're clearly very familiar with), you should be detecting primarily z=2 objects by virtue of their combined redshifted 8uqm PAH emission + redshifted OH opacity minimum.  I really wonder what you'll find if you try it -- it works well in the EGS dataset, with low-z interloper contamination pretty durn low.  It may be that a selection of this kind would ID high-z clusters.  Just grasping at straws, but the talk I saw last week on this topic was very convincing.  Maybe that's a good topic for a subsequent paper?

;vega = ab + abtovega
;ch1	    -2.7
;ch2	    -3.25
;ch3	    -3.73
;ch4 	    -4.37

candidate = where((objectnew.irac1mag -2.7) - (objectnew.irac2mag - 3.25) gt 0.5 and (objectnew.irac3mag -3.73) - (objectnew.irac4mag - 4.37) lt 0 and objectnew.mips24flux gt 0)

print, n_elements(candidate), "candidate"

plothist, objectnew[candidate].zphot, bin =0.1

;make *.reg file for viewing objcts in ds9
openw, outlunred, '/Users/jkrick/nep/z=2_candidate.reg', /get_lun
printf, outlunred, 'fk5'
for rc=0, n_elements(candidate) -1 do  printf, outlunred, 'circle( ', objectnew[candidate[rc]].ra, objectnew[candidate[rc]].dec, ' 3")'
close, outlunred
free_lun, outlunred

;plothyperz, candidate, '/Users/jkrick/nep/z=2.ps


;--------------------------------
;what about just selecting things with photz =2, and at least 5
;detections or something.

nflux = get_nflux()
nfluxlimit =5

cand2 = where(objectnew.zphot gt 1.8 and objectnew.zphot lt 2.2 and nflux ge nfluxlimit )
;make *.reg file for viewing objcts in ds9
openw, outlunred, '/Users/jkrick/nep/z=2_cand2.reg', /get_lun
printf, outlunred, 'fk5'
for rc=0, n_elements(cand2) -1 do  printf, outlunred, 'circle( ', objectnew[cand2[rc]].ra, objectnew[cand2[rc]].dec, ' 4") # color = red'
close, outlunred
free_lun, outlunred


end
