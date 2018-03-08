pro cluster_flamex

;z=0.17
dm = 39.6
dm = dm -0.17

;

;----------------------------------------------------------------
;read in the 'matched' cluster data
;----------------------------------------------------------------

ftab_ext, '/Users/jkrick/nep/clusters/flamex/match_k_33_34_35.fits',[7,8,11,13,1,2,15,16] ,ra, dec, kmag, kmagerr, bmag1, bmagerr1, bmag2, bmagerr2;'col8, col9, col16,col35,mag_auto_1, magerr_auto_1,mag_auto_2, magerr_auto_2'

;combine bmag arrays into one array
bmag = fltarr(n_elements(bmag1) )
for i = 0l, n_elements(bmag1) -1 do begin
   if bmag1(i) gt 0 and bmag1(i) lt 99 then bmag(i) = bmag1(i)
   if bmag2(i) gt 0 and bmag2(i) lt 99  then bmag(i) = bmag2(i)
endfor


;vega to AB
kmag = kmag + 1.87
bmag = bmag - .16


;define cluster area
racluster=217.4978562 
deccluster =34.0273418

sep = .072;0.085
start = where( sphdist(ra, dec, racluster, deccluster, /degrees)lt sep and kmag lt 40 and bmag gt 0 and bmag lt 40)   ;   

;----------------------------------------------------------------
;plot cmd and  k band distribution
;----------------------------------------------------------------

!p.multi=[0,1,2]
ps_open, filename= '/Users/jkrick/nep/clusters/flamex/cmd.ps', /portrait,  /color
vsym, /polygon, /fill

plot, kmag(start), bmag(start) - kmag(start), psym = 8, xrange=[12,24], xstyle = 9, ythick=3,xthick=3, $
      charthick=3, xtitle='K (AB)', ytitle ='B-K (AB)'

    ;fit cmd by eye.
m = -0.12
b = 5.0
sig1 = 1.0
oplot, findgen(40), (m*findgen(40)) + b, thick = 2 ;, color = orangecolor
oplot, findgen(40), (m*findgen(40)) + b + sig1, thick = 2 ;, color = orangecolor
oplot, findgen(40), (m*findgen(40)) + b - sig1, thick = 2 ;, color = orangecolor

 axis, 0, 8, xaxis=1, xrange=[12-dm, 24-dm], xstyle=1, xthick=3, charthick=3

;histogram of all galaxies in the area
;plothist, kmag(start), bin = 1, xstyle=9, xthick=3, ythick=3, thick=3, charthick=3
; axis, 0, 200, xaxis=1, xrange=[12-dm, 24-dm], xstyle=1, xthick=3, charthick=3

;----------------------------------------------------------------
;use cmd color requirement for membership
;----------------------------------------------------------------
 singlecount = 0
 kmemberarr = fltarr(n_elements(start))
 bmemberarr = fltarr(n_elements(start))
 for count = 0, n_elements(start) - 1 do begin
    limithigh = (m*kmag(start(count))) + b + sig1
    limitlow = (m*kmag(start(count)))+ b - sig1
    
    if (bmag(start(count)) - kmag(start(count))) LT limithigh  AND  (bmag(start(count)) - kmag(start(count))) GT limitlow then begin
      kmemberarr(singlecount) = kmag(start(count))
      bmemberarr(singlecount) = bmag(start(count))
      singlecount = singlecount  + 1
    endif
 endfor
 
 kmemberarr=kmemberarr[0:singlecount-1]
 bmemberarr=bmemberarr[0:singlecount-1]
 
 print, "start, singlecount ", n_elements(start), singlecount

;historgram of only members
plothist, kmemberarr, bin = 1, xstyle=9, xthick=3, ythick=3, thick=3, charthick=3, yrange=[0,100], xtitle='K (AB)'
 axis, 0, 100, xaxis=1, xrange=[12-dm, 24-dm], xstyle=1, xthick=3, charthick=3

oplot, [17.4,17.4], [0,150], linestyle=2
oplot, [19.4,19.4], [0,150], linestyle=2
oplot, [21.4,21.4], [0,150], linestyle=2

;-----------------------------------------------------------------------
;what is the average number of galaxies in 100 random regions of the same size as above?
;-----------------------------------------------------------------------
seed = 230
nrand = 100
randra = randomu(seed, nrand)
randdec = randomu(seed, nrand)

print, max(ra), min(ra)
print, max(dec), min(dec)

deltara = max(ra) - min(ra) -0.2
deltadec = max(dec) - min(dec) - 0.2
randra = randra*(deltara) + min(ra) + 0.1
randdec = randdec*(deltadec) + min(dec) + 0.1

gc = 0
goodbrightarr = fltarr(nrand)
gooddimarr = fltarr(nrand)

 
for cand=0, n_elements(randdec) - 1 do begin ;for each random area

   bc = 0
   dc = 0
   
   print, 'cand', cand
   for s = 0l, n_elements(ra) - 1 do begin
      limithigh = (m*kmag(s)) + b + sig1
      limitlow = (m*kmag(s))+ b- sig1
      
      if sphdist(ra(s), dec(s), randra(cand), randdec(cand), /degrees)  lt sep  and kmag(s) le 19.4 and kmag(s) gt 17.4 $
         and bmag(s) - kmag(s) LT limithigh  AND  (bmag(s) - kmag(s)) GT limitlow then begin
         bc = bc + 1
      endif
      if sphdist(ra(s), dec(s), randra(cand), randdec(cand), /degrees) lt sep  and kmag(s) lt 21.4 and kmag(s) gt 19.4 $
         and bmag(s) - kmag(s) LT limithigh  AND  (bmag(s) - kmag(s)) GT limitlow then begin
         dc = dc + 1
      endif
      
   endfor
   goodbrightarr(cand) = bc
   gooddimarr(cand) = dc
endfor 
goodbrightarr = goodbrightarr(where(goodbrightarr gt 0))
gooddimarr = gooddimarr(where(gooddimarr gt 0))

print, goodbrightarr
goodcolorbright = mean(goodbrightarr)
goodcolordim = mean(gooddimarr)
print, "goodcolorbright, goodcolordim " , goodcolorbright, goodcolordim
print, "stddev", stddev(goodbrightarr), stddev(gooddimarr)

  ;----------------------------------------------------------------
;make a measurement
;----------------------------------------------------------------
;break between bright and dim is Mk = -20, or mk = 18.2  /pm 2mags on each side
 bright = where(kmemberarr lt 19.4 and kmemberarr gt 17.4)
 dim = where(kmemberarr ge 19.4 and kmemberarr lt 21.4)

print, 'results, bright, dim (uncorrected): ', n_elements(bright), n_elements(dim)
print, 'results, bright, dim : ', n_elements(bright) - goodcolorbright, n_elements(dim)-goodcolordim


  ps_close, /noprint,/noid
  
end


;----------------------------------------------------------------
;read in the cluster data
;----------------------------------------------------------------
;  readlargecol, '/Users/jkrick/nep/clusters/flamex/smallcat.txt', GALAXY_ID    ,NUMBER          ,XPEAK_IMAGE    ,YPEAK_IMAGE    ,X_IMAGE        ,Y_IMAGE         ,ALPHAPEAK_J2000 ,DELTAPEAK_J2000 ,FLAGS           , FWHM_IMAGE   $
;                ,FLUX_RADIUS     ,MAG_ISO         ,MAG_ISOCOR      ,MAG_AUTO        ,MAG_BEST        ,MAG_APER1,MAG_APER2,MAG_APER3,MAG_APER4,MAG_APER5,$
;                MAG_APER6,MAG_APER7,MAG_APER8,MAG_APER9,MAG_APER10,MAG_APER11,MAG_APER12,MAG_APER13,MAGERR_ISO      ,MAGERR_ISOCOR  $
;                ,MAGERR_AUTO     ,MAGERR_BEST   ,format="A"


  ;assume catalog is in Vega, convert into AB.
;  mag_best = mag_best + 1.87

;get those objects near to the cluster

;  good = where(sqrt((alphapeak_j2000 - 217.4978562)^2 + (deltapeak_J2000-34.0273418)^2) lt sep  )
