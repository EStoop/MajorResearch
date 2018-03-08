pro cluster

restore, '/Users/jkrick/idlbin/object.old.sav'
!p.multi = [0, 1,2]
;!p.multi = [0, 1, 2]
;ps_open, file = "/Users/jkrick/nep/clusters/cluster.cmd.ps", /portrait, /color,/square, xsize=4, ysize=4
!P.charthick = 1
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

openw, outlunred, '/Users/jkrick/nep/clusters/redder.new.reg', /get_lun
printf, outlunred, 'fk5'

openw, outlunback, '/Users/jkrick/nep/clusters/backfaint.reg', /get_lun
printf, outlunback, 'fk5'
openw, outlunclus, '/Users/jkrick/nep/clusters/clusfaint.reg', /get_lun
printf, outlunclus, 'fk5'

memberarr = fltarr(50000)
totcount = 0.D
ck=0
a24 = fltarr(500)
;setup variable for the colors of all clusters
iall = fltarr(210000)
iracall = fltarr(210000)
mips24all = fltarr(210000)
  j = 0l
iallz = fltarr(210000)
iracallz = fltarr(210000)
m2 = 0l 
;I know that the photometry is bad for these objects and they end up with wrong,very red colors, so don't even let them into the sample
;this happens because they are not detections in irac, but are near a bright source, so they will get flux in the aperture as
;determined by the acs photometry.
badnum=[21735,21820,21854,22133,22169,22707,22761,23225]
objectnew[badnum].irac1mag = -99.0

m=-0.12
ball = 5.1
rcswidth = 0.5
;-----------------------------------------------------------------------
;what is the average redshift distribution in 50 random regions of the same size as below?
;-----------------------------------------------------------------------
sep =0.017 ;(55")  0.0062; (720kpc)
;sep = 0.01672
;first find 50 random regions which are not near the edge.
;look at central 13000x13000 region of hst image by pixel numbers                                   
;need 100 pairs of numbers between 1 and 11000
seed = 420
nrand = 10
randx = randomu(seed, nrand)
randy = randomu(seed, nrand)

randx = randx * 9500. + 6800.
randy = randy * 13300. + 3100.
;randx = randx * 11000. + 5500.
;randy = randy * 15000. + 2500.

acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits') ;
xyad, acshead, randx, randy, randra, randdec

randz = fltarr(1000000)
r = 0.D
gc = 0
goodbrightarr = fltarr(nrand)
backmemarr = fltarr(nrand)
gooddimarr = fltarr(nrand)
backcount = fltarr(nrand)
singarr=fltarr(nrand)
randmemberarr = fltarr(nrand*1000.)
totrand = 0
for cand=0, n_elements(randdec) - 1 do begin ;for each random area
   good = where(sphdist(objectnew.ra, objectnew.dec, randra(cand), randdec(cand), /degrees) lt sep and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.acsmag gt 0 and objectnew.acsmag lt 90 )

   dimcount  = 0
   brightcount = 0
   singlecountfield = 0

   for c =0, n_elements(good) -1 do begin
      randz(r) = objectnew[good[c]].zphot
      r = r + 1

      limithigh = (m*objectnew[good[c]].irac1mag) + ball +rcswidth; 2*sig1;rcswidth 
      limitlow = (m*objectnew[good[c]].irac1mag) + ball - rcswidth ;2*sig1;rcswidth 

      if (objectnew[good[c]].acsmag - objectnew[good[c]].irac1mag LT limithigh ) AND ( objectnew[good[c]].acsmag - objectnew[good[c]].irac1mag GT limitlow) and objectnew[good[c]].irac1mag lt 24.6 and objectnew[good[c]].irac1mag gt 22.6 then dimcount = dimcount + 1
      if (objectnew[good[c]].acsmag - objectnew[good[c]].irac1mag LT limithigh )AND ( objectnew[good[c]].acsmag - objectnew[good[c]].irac1mag GT limitlow) and objectnew[good[c]].irac1mag lt 22.6 and objectnew[good[c]].irac1mag gt 20.6 then brightcount = brightcount + 1
      
      if (objectnew[good[c]].acsmag - objectnew[good[c]].irac1mag LT limithigh )AND ( objectnew[good[c]].acsmag - objectnew[good[c]].irac1mag GT limitlow) then begin
         singlecountfield = singlecountfield + 1
         randmemberarr(totrand) = good[c]
         totrand = totrand + 1
      endif

   endfor
   backcount[cand] = n_elements(good)
;   goodcolorbright = where(sphdist(objectnew.ra, objectnew.dec, randra(cand), randdec(cand), /degrees) lt sep and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.acsmag - objectnew.irac1mag gt 2.0 and objectnew.acsmag - objectnew.irac1mag lt 3.5 and objectnew.irac1mag lt 22.6 and objectnew.irac1mag gt 21.1, brightcount)

;   goodcolordim = where(sphdist(objectnew.ra, objectnew.dec, randra(cand), randdec(cand), /degrees) lt sep and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.acsmag - objectnew.irac1mag gt 2.0 and objectnew.acsmag - objectnew.irac1mag lt 3.5 and objectnew.irac1mag lt 24.1 and objectnew.irac1mag gt 22.6,dimcount)

;   backmem = where(sphdist(objectnew.ra, objectnew.dec, randra(cand), randdec(cand), /degrees) lt sep and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.acsmag gt 0 and objectnew.acsmag lt 90 and objectnew.acsmag - objectnew.irac1mag gt 2.0 and objectnew.acsmag - objectnew.irac1mag lt 3.5 and objectnew.irac1mag lt 21.4 , backmemcount)

   goodbrightarr(gc) = brightcount; n_elements(goodcolorbright)
   gooddimarr(gc) = dimcount; n_elements(goodcolordim)
;   backmemarr(gc) =  backmemcount
   gc = gc + 1

;   print, 'field ',n_elements(good),  singlecountfield, brightcount, dimcount
   singarr[cand] = singlecountfield
endfor 

free_lun, outlunback
print, 'singlearrcount', mean(singarr)
print, 'backcount', mean(backcount)

goodbrightarr = goodbrightarr(where(goodbrightarr gt 0))
gooddimarr = gooddimarr(where(gooddimarr gt 0))
backmemarr = backmemarr[0:gc-1]

randmemberarr=randmemberarr[0:totrand - 1]
save, randmemberarr, filename='/Users/jkrick/nep/clusters/randmemberarr.txt'

goodcolorbright = mean(goodbrightarr)
goodcolordim = mean(gooddimarr)
print, 'backmem', mean(backmemarr)
print, "goodcolorbright, goodcolordim " , goodcolorbright, goodcolordim
print, "stddev", stddev(goodbrightarr), stddev(gooddimarr)
randz = randz[0:r-1]

plothist, randz, xhist, yhist, bin = 0.05, /noplot
yhist = yhist / nrand

;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;look at all objects within sep
;-----------------------------------------------------------------------
;candidate clusters
;;canddec = [69.04579,  69.06851,  69.087766];,  68.98017,  69.019103 ]
;;candra = [264.69266, 264.89228, 264.83053];, , 264.66337, 265.27136 ]
canddec = [69.04481,  69.06851,  69.087766];,  68.98017,  69.019103 ]
candra = [264.68160, 264.89228, 264.83053];, , 264.66337, 265.27136 ]

;delete immediately
canddec = [69.04481,    69.087766];,  68.98017,  69.019103 ]
candra = [264.68160, 264.83053];, , 264.66337, 265.27136 ]

openw, outlun23, '/Users/jkrick/nep/clusters/rcslt24.0.reg', /get_lun
printf, outlun23, 'fk5'
printf, outlun23, 'global color=red font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source'

for cand=0,1 do begin
;for cand=0, n_elements(canddec) - 1 do begin
   print, 'working on cluster', candra(cand), canddec(cand)
   
   ;objects within sep of  the center and detected in irac1
   good = where(sphdist(objectnew.ra, objectnew.dec, candra(cand), canddec(cand), /degrees) lt sep and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.acsmag gt 0 and objectnew.acsmag lt 90); and objectnew.nearestdist gt 0.00055) ;and objectnew.nearestdist*60*60/0.6 gt 4.0 
   print, 'n_elements(good), cluster', n_elements(good)
   ;put the starting members in an array of all clusters
     for k = 0l, n_elements(good) - 1 do begin
        iall[j] = objectnew[good[k]].acsmag
        iracall[j] =objectnew[good[k]].irac1mag
        mips24all[j] = objectnew[good[k]].mips24flux
        j = j + 1
     endfor

     goodz = where(sphdist(objectnew.ra, objectnew.dec, candra(cand), canddec(cand), /degrees) lt sep and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90  and objectnew.zphot gt 0.8 and objectnew.zphot lt 1.1 and objectnew.prob gt 20) 
     for l = 0l, n_elements(goodz) - 1 do begin
        iallz[m2] = objectnew[goodz[l]].acsmag
        iracallz[m2] =objectnew[good[l]].irac1mag
        m2 =m2 + 1
     endfor

;-------------------------------------------------------------------
; individual redshift distributions
;---------------------------------------------------------------------

   plothist, objectnew[good].zphot, xhist, yhist, bin = 0.05, xrange=[0,1.5], /xstyle, thick = 3, xthick = 3, ythick = 3, $
             charthick = 3, xtitle='Redshift', ystyle =1, ytitle='Number'
;title=strcompress(string(candra(cand)) + string(canddec(cand)))  yrange=[0,10],

;what is the actual phot-z predicted redshift of the cluster?
;   print, 'xhist', xhist[15:27]
;   print, 'yhist', yhist[15:27]
   start = [1.0, 0.3, 50.0]
   err = fltarr(n_elements(xhist)) +1.0 
   result1= MPFITFUN('gauss',xhist[15:27],yhist[15:27], err, start,/quiet)   
   print, 'gaussian fit', result1

   oplot, xhist, result1(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result1(0))/result1(1))^2.), linestyle = 2

   ;add the average distribution
;;   oplot, xhist, yhist, linestyle = 1,thick=3
;   oplot, xhist, yhist + 3*sqrt(yhist), linestyle = 2,thick=3
;   plothist, randz/nrand, xhist2, yhist2, bin = 0.05, /overplot


;-------------------------------------------------------------------
;cmd fitting section
;---------------------------------------------------------------------
;;   plot, objectnew[good].irac1mag, objectnew[good].acsmag - objectnew[good].irac1mag,psym = 2, yrange=[0,8], xtitle = 'acs', ytitle = 'acs - ch1', xrange=[15,25], thick = 3, xthick = 3, ythick = 3, charthick = 3, xstyle=1

;sort them in rmag space to only fit to the better data
   sortindex= sort(objectnew[good].irac1mag)
   irac1sort= objectnew[good].irac1mag[sortindex]
   rsort = objectnew[good].rmaga[sortindex]
   acssort = objectnew[good].acsmag[sortindex]
   zsort = objectnew[good].zmagbest[sortindex]
   
   a = where(acssort - irac1sort gt 2.0 and acssort - irac1sort lt 4.0 )
;a = where(irac1sort lt 22)
;biweight fit (robust for non-gaussian distributions)
;coeff = ROBUST_LINEFIT( objectnew.rmag, objectnew.vr, yfit, sig, coeff_sig)
   if n_elements(a) gt 1 then begin
      coeff1 = ROBUST_LINEFIT(irac1sort(a),acssort(a) - irac1sort(a), yfit1, sig1, coeff_sig1)
      coeff2 = ROBUST_LINEFIT(irac1sort,rsort - irac1sort, yfit2, sig2, coeff_sig2)
      coeff3 = ROBUST_LINEFIT(irac1sort,zsort - irac1sort, yfit3, sig3, coeff_sig3)
      
;      print, "sig1, sig2", sig1,sig2,sig3
      err = dindgen(n_elements(good)) - dindgen(n_elements(good)) + 1
      start = [-0.02,2.0]
      startr = [-0.02, 4.0]
      sortindex= sort(objectnew[good].irac1mag)
      irac1sort= objectnew[good].irac1mag[sortindex]
      rsort = objectnew[good].rmaga[sortindex]
      acssort = objectnew[good].acsmag[sortindex]
      zsort = objectnew[good].zmagbest[sortindex]
      yfit1sort = yfit1[sortindex]
      yfit2sort = yfit2[sortindex]
      yfit3sort = yfit3[sortindex]
      
   endif
   
   
;instead just pick a rcs by eye and with bc03 idea of where the color should be given the redshift
;want to fix the slope at -0.09 de lucia etal.
;can come back and fix this later
   xvals = findgen(22) +13
   m = -0.12
   b = fltarr(5)
   b = [4.9,4.5,4.5,4.5,5.2]
   
;;   oplot, xvals,  m*xvals + b(cand), thick = 3 , color = orangecolor
;;   oplot, xvals, m*xvals + b(cand) + rcswidth, thick = 3 , color = orangecolor
;;   oplot, xvals,  m*xvals + b(cand) - rcswidth, thick = 3 , color = orangecolor
;;   oplot, xvals,  m*xvals + ball, thick = 3 , color = orangecolor
;;   oplot, xvals, m*xvals + ball + rcswidth, thick = 3 , color = orangecolor
;;   oplot, xvals,  m*xvals + ball - rcswidth, thick = 3 , color = orangecolor


;---------------------------------------------------------------------
;determine which objects are 'members' of the RCS
;---------------------------------------------------------------------

 

   singlecount = 0
   singlememberarr = fltarr(n_elements(good))
   brightcountclus = 0
   dimcountclus = 0
   for count = 0.D, n_elements(good) - 1 do begin
      limithigh = (m*objectnew[good[count]].irac1mag) + ball +rcswidth; 2*sig1;0.5 
      limitlow = (m*objectnew[good[count]].irac1mag) + ball - rcswidth;2*sig1;0.5 
;      limithigh = (m*objectnew[good[count]].irac1mag) + b(cand) +rcswidth; 2*sig1;0.5 
;      limitlow = (m*objectnew[good[count]].irac1mag) + b(cand) - rcswidth;2*sig1;0.5 
      if (objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag LT limithigh )AND ( objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag GT limitlow) then begin
       ; add this member galaxy to a larger array of all member galaxies
         memberarr(totcount) = good[count]
         totcount = totcount + 1
         ;also want to know members of this cluster only
         singlememberarr(singlecount) =good[count]
         singlecount = singlecount  + 1
        printf, outlunclus, 'circle( ', objectnew[good[count]].ra, objectnew[good[count]].dec, ' 3")'

     endif

;      if (objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag LT limithigh ) $
;         AND ( objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag GT limitlow) $
;         and objectnew[good[count]].irac1mag lt 24.6 and objectnew[good[count]].irac1mag gt 22.6 then begin
;         dimcountclus = dimcountclus + 1
         
;      endif

;      if (objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag LT limithigh )AND ( objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag GT limitlow) and objectnew[good[count]].irac1mag lt 22.6 and objectnew[good[count]].irac1mag gt 20.6 then brightcountclus = brightcountclus + 1


;remove the set of statements, reinstate those above

      if (objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag LT limithigh )AND ( objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag GT limitlow) and objectnew[good[count]].irac1mag lt 22.6 and objectnew[good[count]].rmaga lt 24.0 then begin
         brightcountclus = brightcountclus + 1

         printf, outlun23, 'circle( ', objectnew[good[count]].ra, objectnew[good[count]].dec, ' 3")'

      endif

      if (objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag LT limithigh )AND ( objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag GT limitlow) and objectnew[good[count]].irac1mag lt 22.6 and objectnew[good[count]].rmaga lt 24.0 then dimcountclus = dimcountclus + 1

   endfor
   print, 'cluster', singlecount, brightcountclus, dimcountclus, totcount

;what are the objects redder than the RCS?
   redder = where(objectnew[good].acsmag - objectnew[good].irac1mag gt 3.5)
   for rc=0, n_elements(redder) -1 do  printf, outlunred, 'circle( ', objectnew[good[redder[rc]]].ra, objectnew[good[redder[rc]]].dec, ' 3")'
;   plothyperz, good[redder], '/Users/jkrick/nep/clusters/redder.ps'
;print, 'redder', good

;---------------------------------------------------------------------
;histograms of member galaxies
;---------------------------------------------------------------------
   singlememberarr = singlememberarr[0:singlecount - 1]
   print, 'number in rcs in cluster;  brighter than 22.1', n_elements(singlememberarr),  n_elements(singlememberarr(where (objectnew(singlememberarr).irac1mag lt 21.4)))

;;   plothist, objectnew[singlememberarr].irac1mag,  thick = 3, xthick = 3, ythick = 3, charthick = 3, $
;;             xtitle='irac1 mag', xrange=[15,25], xstyle = 9;, yrange=[0,12]
;;   axis, 0,12, xaxis=1,xrange=[-24.1,-16.1], xstyle=1, xthick=3, charthick=3

;   plothyperz, singlememberarr, '/Users/jkrick/nep/clusters/member.ps'

endfor 
close, outlunred
free_lun, outlunred

close, outlun23
free_lun, outlun23

free_lun, outlunclus
iall = iall[0:j-1]
iracall = iracall[0:j-1]
mips24all = mips24all[0:j-1]
memberarr = memberarr[0:totcount - 1]
iallz = iallz[0:m2-1]
iracallz = iracallz[0:m2-1]

save, memberarr, filename = '/Users/jkrick/nep/clusters/memberarr.txt'
print, 'n_red',n_elements(where( iall - iracall gt 3.5))
;print,'memberarr',  n_elements(memberarr)


;----------------------------------------------------------------
;plot combined cmd and histogram of all cand  member galaxies
;----------------------------------------------------------------
dm = 44.1
dm = dm - 0.75  ;for filter bandpass shrinking
dm = dm - 0.7  ;evolution correction
kcorr = 1.2
kcorr = 0.

kcorr = -1.6
kcorr = 0.

ball = 5.2

;;plot, iracall, iall - iracall,psym = 2, $
;;      yrange=[0,6], xtitle = '[3.6] (AB)', ytitle = 'ACS F814W - [3.6]', xrange=[18,25], thick = 3, xthick = 3, $
;;      ythick = 3, charthick = 3, xstyle=9;, title='composite'
;;axis, 0,6, xaxis=1,xrange=[18-dm+kcorr,25-dm+kcorr], xstyle=1, xthick=3, charthick=3

;;oplot, findgen(100),  m*findgen(100) + ball, thick = 3 ;, color = orangecolor
;;oplot, findgen(100), m*findgen(100) + ball + rcswidth, thick = 3 ;, color = orangecolor
;;oplot, findgen(100),  m*findgen(100) +ball - rcswidth, thick = 3 ;, color = orangecolor

;oplot, iracallz, iallz - iracallz, psym = 2, thick=3, color = orangecolor
;plot magnitude detection limit
;;oplot, [25.1, 22], [2.9,6.0], linestyle = 2, thick=3

ps_open, file = "/Users/jkrick/nep/clusters/cmd.allgals.ps", /portrait, /color,/square, xsize=4, ysize=4

im_hessplot, iracall, iall -iracall, yrange=[0,6], xtitle = '[3.6] (AB)', ytitle = 'ACS F814W - [3.6]', xrange=[18,25],/notsquare, nbin2d=15, position = [0.1 ,0.1 , 0.9,0.9], charsize=1, charthick=3, xthick=3, ythick=3, thick=3,xstyle=9, ystyle = 1
axis, 0,6, xaxis=1,xrange=[18-dm+kcorr,25-dm+kcorr], xstyle=1, xthick=3, charthick=3


oplot, findgen(100),  m*findgen(100) + ball, thick = 3 ;, color = orangecolor
oplot, findgen(100), m*findgen(100) + ball + rcswidth, thick = 3 ;, color = orangecolor
oplot, findgen(100),  m*findgen(100) +ball - rcswidth, thick = 3 ;, color = orangecolor
oplot, [25.1, 22], [2.9,6.0], linestyle = 2, thick=3

; make a similar cmd for all galaxies with irac and acs detections
all = where( objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.acsmag gt 0 and objectnew.acsmag lt 90 and objectnew.acsclassstar lt 1.)
;-----------------------------------------------------------
;do some monte carlo calculation of the probability of getting the
;same amount of galaxies on the red sequence by taking 711 random
;galaxies instead of cluster galaxies....is some sort of proof of the
;red sequence and therefore the cluster existance.

;this is done in montecarlo_cmd.pro
;-----------------------------------------------------------

;want 711 random galaxies from this distribution
randselect = 50900*randomu(S, 711)
im_hessplot, objectnew[randselect].irac1mag, objectnew[randselect].acsmag - objectnew[randselect].irac1mag, yrange=[0,6], xtitle = '[3.6] (AB)', ytitle = 'ACS F814W - [3.6]', xrange=[18,25],/notsquare, nbin2d=15, position = [0.1 ,0.1 , 0.9,0.9], charsize=1, charthick=3, xthick=3, ythick=3, thick=3,xstyle=1, ystyle = 1
;axis, 0,6, xaxis=1,xrange=[18-dm+kcorr,25-dm+kcorr], xstyle=1, xthick=3, charthick=3


oplot, findgen(100),  m*findgen(100) + ball, thick = 3 ;, color = orangecolor
oplot, findgen(100), m*findgen(100) + ball + rcswidth, thick = 3 ;, color = orangecolor
oplot, findgen(100),  m*findgen(100) +ball - rcswidth, thick = 3 ;, color = orangecolor
oplot, [25.1, 22], [2.9,6.0], linestyle = 2, thick=3

ps_close, /noprint, /noid


;plothist, iracall,  thick = 3, xthick = 3, ythick = 3, charthick = 3, xtitle='ch1 mag', $
 ;         xrange=[18,25], xstyle = 9, yrange=[0,100]
;axis, 0,100, xaxis=1,xrange=[18-dm+kcorr,25-dm+kcorr], xstyle=1, xthick=3, charthick=3
;oplot, [24.6, 24.6], [0,60], linestyle=2
;oplot, [22.6, 22.6], [0,60], linestyle=2
;oplot, [20.6, 20.6], [0,60], linestyle=2

;plot, objectnew[memberarr].irac1mag, objectnew[memberarr].acsmag - objectnew[memberarr].irac1mag,psym = 2, $
;      yrange=[0,6], xtitle = '[3.6] (AB)', ytitle = 'ACS F814W - [3.6]', xrange=[18,25], thick = 3, xthick = 3, $
;      ythick = 3, charthick = 3, xstyle=9
;axis, 0,6, xaxis=1,xrange=[18-dm+kcorr,25-dm+kcorr], xstyle=1, xthick=3, charthick=3

;ps_close, /noprint, /noid

;ps_open, file = "/Users/jkrick/nep/clusters/cluster.dist.ps", /portrait, /color,/square, xsize=4, ysize=4

;;plothist, objectnew[memberarr].irac1mag,  thick = 3, xthick = 3, ythick = 3, charthick = 3, $
;;          xtitle='[3.6] (AB)', ytitle='Number', xrange=[18,25], xstyle = 9, yrange=[0,80]
;;axis, 0,80, xaxis=1,xrange=[18-dm+kcorr,25-dm+kcorr], xstyle=1, xthick=3, charthick=3
;;oplot, [24.6, 24.6], [0,80], linestyle=2
;;oplot, [22.6, 22.6], [0,80], linestyle=2
;;oplot, [20.6, 20.6], [0,80], linestyle=2

;----------------------------------------------------------------
;want a member array chosen from the composite CMD
;----------------------------------------------------------------
compmemi = fltarr(n_elements(iall))
compmemirac = fltarr(n_elements(iall))
totcompcount = 0
for compcount = 0.D, n_elements(iall) - 1 do begin
   limithigh = (m*iracall(compcount)) + ball +rcswidth    ; 2*sig1;rcswidth 
   limitlow = (m*iracall(compcount)) + ball - rcswidth    ;2*sig1;rcswidth 
   if (iall(compcount) - iracall(compcount) LT limithigh )AND ( iall(compcount) - iracall(compcount)  GT limitlow) then begin
      compmemi(totcompcount) = iall(compcount)
      compmemirac(totcompcount) = iracall(compcount)
      totcompcount = totcompcount + 1
   endif
endfor

compmemi = compmemi[0:totcompcount - 1]
compmemirac = compmemirac[0:totcompcount - 1]


;----------------------------------------------------------------
;make a statistical subtraction of field galaxies
;----------------------------------------------------------------
nclusters = n_elements(candra)
print, 'nclusters', nclusters
;count all galaxies outside of the cluster areas
rad = 0.05
field = where(sphdist(objectnew.ra, objectnew.dec, candra(0), canddec(0), /degrees)  gt rad and sphdist(objectnew.ra, objectnew.dec, candra(1), canddec(1), /degrees)  gt rad  and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90  and objectnew.acsmag gt 0 and objectnew.acsmag lt 90 and objectnew.acsmag - objectnew.irac1mag gt 2.0 and objectnew.acsmag - objectnew.irac1mag lt 3.5) ;and sphdist(objectnew.ra, objectnew.dec, candra(2), canddec(2), /degrees)  gt rad

;plot, objectnew[field].irac1mag, objectnew[field].acsmag - objectnew[field].irac1mag, psym = 8, xtitle= 'ch1', ytitle='acs-ch1', title = 'all field', symsize = 0.2

;bin into the two bins I want to eventually measure.
fieldbright = where(objectnew[field].irac1mag lt 22.6 and objectnew[field].irac1mag gt 20.6)
fielddim = where(objectnew[field].irac1mag ge 22.6 and objectnew[field].irac1mag lt 24.6)

;how much area is that?
;irac detection and acs detection limits it to acs area - area in clusters.

;area in clusters is a circle with areadegree radius
;assume clusters are not overlapping
clusarea = !PI*rad^2
clusarea = nclusters * clusarea    ;using nclusters clusters

;area covered by acs
fits_read, '/Users/jkrick/hst/raw/mosaic_lowres.fits', acsdata, acshead
naxis1=sxpar(acshead,'NAXIS1')-2
naxis2=sxpar(acshead,'NAXIS2')-2
no = where(acsdata eq 0)         ;no valid data
acsarea = naxis1*naxis2 - n_elements(no)   ;number of pixels with valid data.
acsarea = acsarea * .16   ;area in square arcseconds
acsarea = acsarea / 3600. / 3600. - clusarea;area in square degrees.
print, 'acsarea', acsarea    ;total area outside of clusters

;how many field galaxies are in the same area as the cluster area
nfieldbright = n_elements(fieldbright)*!PI*sep^2/ acsarea
nfielddim = n_elements(fielddim)*!PI*sep^2/ acsarea
print, 'nfieldbright, nfielddim', nfieldbright, nfielddim

;----------------------------------------------------------------
;make a measurement
;----------------------------------------------------------------
;break between bright and dim is Mk = -20, or m3.6 = 22.6  /pm 2mags on each side
 bright = where(objectnew[memberarr].irac1mag lt 22.6 and objectnew[memberarr].irac1mag gt 20.6)
 dim = where(objectnew[memberarr].irac1mag ge 22.6 and objectnew[memberarr].irac1mag lt 24.6)

print, 'results, bright, dim (uncorrected): ', n_elements(bright) , n_elements(dim) 
print, 'results, bright, dim : ', n_elements(bright) -nclusters* goodcolorbright, n_elements(dim) - nclusters*goodcolordim

brightcomp = where(compmemirac lt 22.6 and compmemirac gt 20.6)
faintcomp = where(compmemirac ge 22.6 and compmemirac lt 24.6)

print, 'composite uncorrected ', n_elements(brightcomp), n_elements(faintcomp)

;ps_close, /noprint, /noid


end




;first need to make a hyperz input with just those objectnews.
;openw, outlunh, "/Users/jkrick/ZPHOT/hyperz_cat_cluster.txt",/get_lun
;for num = 0, n_elements(good) - 1 do begin
;
;   if objectnew[good(num)].flamjmag gt 0 and objectnew[good(num)].flamjmag ne 99 then begin
;      fab = 1594.*10^(objectnew[good(num)].flamjmag/(-2.5))
;      jmagab = -2.5*alog10(fab) +8.926
;   endif else begin
;      jmagab = objectnew[good(num)].flamjmag
;   endelse
;
;   if objectnew[good(num)].wircjmag gt 0 and objectnew[good(num)].wircjmag ne 99 then begin
;      wircjab = 1594.*10^(objectnew[good(num)].wircjmag/(-2.5))
;      wircjmagab = -2.5*alog10(wircjab) +8.926
;   endif else begin
;      wircjmagab = objectnew[good(num)].wircjmag
;   endelse
;
;   if objectnew[good(num)].wirchmag gt 0 and objectnew[good(num)].wirchmag ne 99 then begin
;      wirchab = 1024.*10^(objectnew[good(num)].wirchmag/(-2.5))
;      wirchmagab = -2.5*alog10(wirchab) +8.926
;   endif else begin
;      wirchmagab = objectnew[good(num)].wirchmag
;   endelse
;
;   if objectnew[good(num)].wirckmag gt 0 and objectnew[good(num)].wirckmag ne 99 then begin
;      wirckab = 666.8*10^(objectnew[good(num)].wirckmag/(-2.5))
;      wirckmagab = -2.5*alog10(wirckab) +8.926
;   endif else begin
;      wirckmagab = objectnew[good(num)].wirckmag
;   endelse
;
;   if objectnew[good(num)].irac1flux lt 0 then err1 = -1. else err1 = 0.05
;   if objectnew[good(num)].irac2flux lt 0 then err2 = -1. else err2 =  0.05;8.926 - 2.5*alog10(0.05E-6*objectnew[good(num)].irac2)
;   if objectnew[good(num)].irac3flux lt 0 then err3 = -1. else err3 =  0.05;8.926 - 2.5*alog10(0.05E-6*objectnew[good(num)].irac3)
;   if objectnew[good(num)].irac4flux lt 0 then err4 = -1. else err4 =  0.05;8.926 - 2.5*alog10(0.05E-6*objectnew[good(num)].irac4)
;   
;   if objectnew[good(num)].imagerra gt 1000. then objectnew[good(num)].imagerra = 1000.
;   if objectnew[good(num)].gmagerra gt 1000. then objectnew[good(num)].gmagerra = 1000.
;   if objectnew[good(num)].rmagerra gt 1000. then objectnew[good(num)].rmagerra = 1000.
;   if objectnew[good(num)].umagerra gt 1000. then objectnew[good(num)].umagerra = 1000.
;   
;
;    printf, outlunh, format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.;;;;;;;;;;;;;2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',$
;                 good(num), objectnew[good(num)].umaga, objectnew[good(num)].gmaga, objectnew[good(num)].rmaga, $
;                 objectnew[good(num)].imaga,  objectnew[good(num)].acsmag,  jmagab, wircjmagab, $
;                 objectnew[good(num)].tmassjmag,  wirchmagab, objectnew[good(num)].tmasshmag, wirckmagab, $
;                 objectnew[good(num)].tmasskmag, objectnew[good(num)].irac1mag,objectnew[good(num)].irac2mag,$
;                 objectnew[good(num)].irac3mag,objectnew[good(num)].irac4mag,   $
;                 objectnew[good(num)].umagerra, objectnew[good(num)].gmagerra, $
;                 objectnew[good(num)].rmagerra, objectnew[good(num)].imagerra, objectnew[good(num)].acsmagerr, objectnew[good(num)].flamjmagerr, $
;                 objectnew[good(num)].wircjmagerr,0.02, objectnew[good(num)].wirchmagerr, 0.02, objectnew[good(num)].wirckmagerr,$
;                 0.02,err1,err2,err3,err4
; endfor
;close, outlunh
;free_lun, outlunh
;
;;run hyperz
;;~/bin/hyperz  zphot.cluster.param
;
;
;
;;plot SED's of keeper objectnews
;
;readcol,'/Users/jkrick/ZPHOT/hyperz_swire_cluster.z_phot',idz, zphota, chia, proba, specta, nagea, agea, $
;        ava, ba, zinf99a,zsup99a,zinf90a,zsup90a,zinf68a,zsup68a,zwma,probwma,Mabsa,$
;        zphot2a,prob2a,format="A"
;
;
;readlargecol, '/Users/jkrick/ZPHOT/hyperz_swire_cluster.obs_sed',idz, u,g,r,i,acs,j,wj,tj,wh,th,wk,tk,one,two,three, four,$
;         uerr, gerr, rerr, ierr,acserr, jerr, wjerr,tjerr,wherr,therr,wkerr,tkerr,oneerr, twoerr, threeerr, fourerr,format="A",/debug;
;
;print, n_elements(zphota), n_elements(u)
;
;;spectral type key
;sptkey=[" Arp220 " ," Ell13  " ," Ell2  " ," I19254 " ," M82  " ," Mrk231" ," QSO1 " ," QSO2 " ," S0   " ," Sa  " ," Sb  " ," Sc  " ," Sd  " ," Sey2  " ," Torus" ]     
;;---------------------------------------------------------------------
;
;
;acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits');
;adxy,  acshead,objectnew[good].ra, objectnew[good].dec, xcenter, ycenter
;
;
;openw, outlun, '/Users/jkrick/nep/cluster.reg', /get_lun
;for count = 0, n_elements(good) - 1 do begin
;   if objectnew[good(count)].zphot gt 0.9 and objectnew[good(count)].zphot lt 1.3 then begin
;      print, "got one"
;      printf, outlun, 'circle( ', xcenter(count), ycenter(count), ' 40)'
;   endif
;
;endfor
;close, outlun
;free_lun, outlun
;
;
;size = 10
;x = [.3540,.4660,.6255,.7680,.8330,1.196, 1.25,1.24,1.635,1.65,2.15,2.16,3.6000,4.5000,5.8,8.]    ;wavelengths of
;
;y = x
;yerr = x
;
;interest =indgen(n_elements(idz)) 
;for n = 0,n_elements(idz) - 1 do begin
;
;;output the acs image of the objectnew
;   print, "working on n ", n
;
;   if objectnew[good(n)].acsmag gt 0 and objectnew[good(n)].acsmag lt 90 then begin
;;     print,    xcenter(n) - size/0.05, xcenter(n)+ size/0.05  ,ycenter(n) -size/0.05, ycenter(n)+size/0.05
;      acsdata = mrdfits('/Users/jkrick/hst/raw/wholeacs.fits', range=[ycenter(n) -size/0.05, ycenter(n)+size/0.05])
;      plotimage, xrange=[xcenter(n) - size/0.05, xcenter(n)+ size/0.05],$
;;                 yrange=[ycenter(n) -size/0.05, ycenter(n)+size/0.05], $
;                 bytscl(acsdata, min = -0.01, max = 0.1),$
;                 /preserve_aspect, /noaxes, ncolors=60
;      xyouts, xcenter(n)- 0.6*size/0.05, -10., strcompress(string(good(n))+ '    ACS' + string(objectnew[good(n)].acsmag)),charthick = 3
;;      xyouts, xcenter(n)- 0.6*size/0.05, -10., string(good(n)),charthick = 3
;   endif
;
;
;
;;plot the fitted SED along with photometry of an individual galaxy
;
;   y = [u(interest(n)),g(interest(n)),r(interest(n)),i(interest(n)),acs(interest(n)), j(interest(n)),wj(interest(n)),tj(interest(n)),wh(interest(n)),th(interest(n)),wk(interest(n)),t;;;;;k(interest(n)),one(interest(n)),two(interest(n)),three(interest(n)), four(interest(n))] ;photometry
;   yerr = [uerr(interest(n)),gerr(interest(n)),rerr(interest(n)),ierr(interest(n)),acserr(interest(n)),jerr(interest(n)),wjerr(interest(n)),tjerr(interest(n)), wherr(interest(n)),therr(interest(n)),wkerr(interest(n)), tkerr(interest(n)),oneerr(interest(n)),twoerr(interest(n)),threeerr(interest(n)), fourerr(interest(n))]

;   plot,alog10(x), alog10(y), psym = 5, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
;        xtitle = "log(microns)", ytitle = "log(flux(microjansky))", title=strcompress("objectnew " + string(idz(n)) +string(zphota(n)) + string(proba(n)) + string(sptkey(specta(n)-1)))

;   errplot, alog10(x), alog10(y - yerr), alog10(y + yerr)
;;
;   readcol,strcompress('/Users/jkrick/ZPHOT/' + string(idz(n)) + '.spe', /remove_all),x2, y2,format="A",/silent
;   oplot, alog10((x2/1E4)), alog10(y2), thick = 3
;   
;endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;canddec = [69.04481,  69.06851,  69.069154,  68.871674,  69.108844,  68.98017,  69.090313,  69.021478,  69.075394, 69.127795, 	69.036086, 68.997594, 	 68.986163,  68.969327,68.99948]

;candra = [264.68160, 264.89228, 265.35535, 265.39991, 264.95902, 264.66337, 264.83102, 265.27616, 265.46662, 265.23219 ,264.91446, 264.7641, 265.34753, 265.17806,265.04271]

 ;  plot, objectnew[good].irac1mag, objectnew[good].rmaga - objectnew[good].irac1mag,psym = 2, yrange=[0,8], xtitle = 'ch 1', ytitle = 'r - ch1', xrange=[16,24], thick = 3, xthick = 3, ythick = 3, charthick = 3

 ;  oplot, xvals, ((resultr(0))*xvals) + resultr(1), thick = 3    ;, color = colors.orange
 ;  oplot, xvals, ((resultr(0))*xvals) + resultr(1)+ sig2 , thick = 3 ;, color = colors.orange
 ;  oplot, xvals, ((resultr(0))*xvals) + resultr(1)- sig2 , thick = 3 ;, color = colors.orange

 ;  oplot, xvals, 27.4 - xvals, thick = 3 , linestyle = 2

;   plot, objectnew[good].irac1mag, objectnew[good].zmagbest - objectnew[good].irac1mag,psym = 2, yrange=[0,8], xtitle = 'ch 1', ytitle = 'z - ch1', xrange=[16,24], thick = 3, xthick = 3, ythick = 3, charthick = 3

;   oplot, xvals, ((resultz(0))*xvals) + resultz(1), thick = 3        ;, color = colors.orange
;   oplot, xvals, ((resultz(0))*xvals) + resultz(1)+ sig3 , thick = 3 ;, color = colors.orange
;   oplot, xvals, ((resultz(0))*xvals) + resultz(1)- sig3 , thick = 3 ;, color = colors.orange
   
;   endif
;   if cand eq 1 then begin
;      print,'working on ',  candra(cand), canddec(cand)
;      oplot, xvals,  (m)*xvals + 4.5, thick = 3 , color = orangecolor
;      oplot, xvals,  (m)*xvals + 5.0, thick = 3 , color = orangecolor
;      oplot, xvals,  (m)*xvals + 4.0, thick = 3 , color = orangecolor
;   endif
;   if cand eq 2 then begin
;;      print,'working on ',  candra(cand), canddec(cand)
;      oplot, xvals,  (m)*xvals + 4.5, thick = 3 , color = orangecolor
;      oplot, xvals,  (m)*xvals + 5.0, thick = 3 , color = orangecolor
;      oplot, xvals,  (m)*xvals + 4.0, thick = 3 , color = orangecolor
;   endif
;   if cand eq 3 then begin
;;      print,'working on ',  candra(cand), canddec(cand)
;      oplot, xvals,  (m)*xvals + 4.5, thick = 3 , color = orangecolor
;      oplot, xvals,  (m)*xvals + 5.0, thick = 3 , color = orangecolor
;      oplot, xvals,  (m)*xvals + 4.0, thick = 3 , color = orangecolor
;   endif
;  if cand eq 4 then begin
;;      print,'working on ',  candra(cand), canddec(cand)
;      oplot, xvals,  (m)*xvals + 5.2, thick = 3 , color = orangecolor
;      oplot, xvals,  (m)*xvals + 5.7, thick = 3 , color = orangecolor
;      oplot, xvals,  (m)*xvals + 4.7, thick = 3 , color = orangecolor
;   endif




;   for i = 0 , n_elements(good) -1 do begin
;      if objectnew[good[i]].acsmag - objectnew[good[i]].irac1mag gt 5 then print, "red", good[i]
;;   endfor
;   b = where( objectnew[good].acsmag - objectnew[good].irac1mag ge 5)
;plothyperz, good[b]
;--------------------------------------------------------
;openw, outlun, '/Users/jkrick/nep/clusters/rand.reg', /get_lun
;for i = 0,n_elements(randx) - 1 do printf, outlun, 'circle( ',randx(i), randy(i),' 1200)'
;close, outlun
;free_lun, outlun


 ;     resultacs = MPFITFUN('linear',irac1sort, yfit1sort,err, start)
;      print, "resultacs", resultacs
;      resultr= MPFITFUN('linear',irac1sort, yfit2sort,err, startr)
;      print, "resultr", resultr
;      resultz= MPFITFUN('linear',irac1sort, yfit3sort,err, start)
;      print, "resultz", resultz

;     oplot, xvals, ((resultacs(0))*xvals) + resultacs(1), thick = 3    ;, color = colors.orange
 ;     oplot, xvals, ((resultacs(0))*xvals) + resultacs(1)+ 0.5 , thick = 3 ;, color = colors.orange
 ;     oplot, xvals, ((resultacs(0))*xvals) + resultacs(1)- 0.5 , thick = 3 ;, color = colors.orange


;plothist, objectnew[memberarr].irac1mag,  thick = 3, xthick = 3, ythick = 3, charthick = 3, xtitle='irac 3.6 flux'
;print, (objectnew[memberarr].irac1mag)
;plothist, objectnew[memberarr].mips24flux, thick = 3, xthick = 3, ythick = 3, charthick = 3, xtitle='mips 24 flux'
;plothist, objectnew[memberarr].rmaga, thick = 3, xthick = 3, ythick = 3, charthick = 3, xtitle='r mag', xrange=[20,30], bin=0.5
;plot, objectnew[memberarr].irac1mag , objectnew[memberarr].wircHmag -  objectnew[memberarr].irac1mag, thick =3, psym = 2, yrange=[-5,5]

;wirckab = 666.8*10^(objectnew[memberarr].wirckmag/(-2.5))
;wirckmagab = -2.5*alog10(wirckab) +8.926
;wircjab = 1594.*10^(objectnew[memberarr].wircjmag/(-2.5))
;wircjmagab = -2.5*alog10(wircjab) +8.926

;plothist, objectnew[memberarr].rmaga - wirckmagab, xrange=[0,10], bin=0.1, charthick = 1, thick = 1
;plothist, objectnew[memberarr].rmaga - wircjmagab, xrange=[0,10], bin=0.5, charthick = 1, thick = 1
;print, n_elements(where(objectnew[memberarr].wircjmag gt 0 and objectnew[memberarr].wircjmag lt 90))

;help, /memory


;----------------------------------------------------------------
;plot combined cmd of cand cluster 1 and 2
;----------------------------------------------------------------
;cluster1
;good1 = where(sqrt((objectnew.ra - candra(0))^2 + (objectnew.dec-canddec(0))^2) lt .015 and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90  ) ;and objectnew.nearestdist*60*60/0.6 gt 4.0 and objectnew.extended lt 1
;good1z = where(sqrt((objectnew.ra - candra(0))^2 + (objectnew.dec-canddec(0))^2) lt .015 and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90   and objectnew.zphot gt 0.7 and objectnew.zphot lt 1.3) ;and objectnew.nearestdist*60*60/0.6 gt 4.0 and objectnew.extended lt 1

;cluster2
;good2 = where(sqrt((objectnew.ra - candra(1))^2 + (objectnew.dec-canddec(1))^2) lt .015 and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90  ) ;and objectnew.nearestdist*60*60/0.6 gt 4.0 and objectnew.extended lt 1
;good2z = where(sqrt((objectnew.ra - candra(1))^2 + (objectnew.dec-canddec(1))^2) lt .015 and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90   and objectnew.zphot gt 0.7 and objectnew.zphot lt 1.3) ;and objectnew.nearestdist*60*60/0.6 gt 4.0 and objectnew.extended lt 1


;good1 = [good1, good2]
;goodz =[good1z, good2z]

;plot, objectnew[good].acsmag, objectnew[good].acsmag - objectnew[good].irac1mag,psym = 2, $
;      yrange=[0,6], xtitle = 'acs', ytitle = 'acs - ch1', xrange=[20,28], thick = 3, xthick = 3, $
;      ythick = 3, charthick = 3, xstyle=9

;oplot, objectnew[goodz].acsmag, objectnew[goodz].acsmag - objectnew[goodz].irac1mag,psym = 2, thick=3, color = yellowcolor

;plot magnitude detection limit
;oplot, [25.1, 28], [0,2.9], linestyle = 2

;add absolute magnitude to top xaxis of plot
;axis, 0,6, xaxis=1,xrange=[-24.1,-16.1], xstyle=1, xthick=3, charthick=3

;-----

;change xaxis
;plot, objectnew[good].irac1mag, objectnew[good].acsmag - objectnew[good].irac1mag,psym = 2, $
;      yrange=[0,6], xtitle = 'ch 1', ytitle = 'acs - ch1', xrange=[18,25.1], thick = 3, xthick = 3, $
;      ythick = 3, charthick = 3, xstyle=9

;oplot, objectnew[goodz].irac1mag, objectnew[goodz].acsmag - objectnew[goodz].irac1mag,psym = 2, thick=3, color = yellowcolor

;plot magnitude detection limit
;oplot, [20.0, 28], [8.0,0], linestyle = 2

;add absolute magnitude to top xaxis of plot
;axis, 0,6, xaxis=1,xrange=[-26.1,-19.3], xstyle=1, xthick=3, charthick=3

;not so good candidates
;notsogoodra = [265.35535, 265.40978,265.34753,265.17806,  264.95557, 265.46662, 265.23219 ,264.91446, 264.7641]
;notsogooddec = [ 69.069154, 68.871061, 68.986163,68.969327,   69.11266,  69.075394, 69.127795,69.036086, 68.997594]

;;plot, iall, iall - iracall,psym = 2, $
;;      yrange=[0,6], xtitle = 'acs', ytitle = 'acs - ch1', xrange=[20,28], thick = 3, xthick = 3, $
;;      ythick = 3, charthick = 3, xstyle=9, title='composite'
;;axis, 0,6, xaxis=1,xrange=[20-dm-kcorr,28-dm-kcorr], xstyle=1, xthick=3, charthick=3

;plot magnitude detection limit
;;oplot, [25.1, 28], [0,2.9], linestyle = 2

;;  plothist, iall,  thick = 3, xthick = 3, ythick = 3, charthick = 3, xtitle='acs mag', xrange=[20,28], xstyle = 9, yrange=[0,100]
;;   axis, 0,100, xaxis=1,xrange=[20-dm-kcorr,28-dm-kcorr], xstyle=1, xthick=3, charthick=3

;----------------------------------------------------------------
;;plot, objectnew[memberarr].acsmag, objectnew[memberarr].acsmag - objectnew[memberarr].irac1mag,psym = 2, $
;;      yrange=[0,6], xtitle = 'acs', ytitle = 'acs - ch1', xrange=[20,28], thick = 3, xthick = 3, $
;;      ythick = 3, charthick = 3, xstyle=9, title='composite'
;;axis, 0,6, xaxis=1,xrange=[20-dm-kcorr,28-dm-kcorr], xstyle=1, xthick=3, charthick=3
;plot magnitude detection limit
;;oplot, [25.1, 28], [0,2.9], linestyle = 2

;;  plothist, objectnew[memberarr].acsmag,  thick = 3, xthick = 3, ythick = 3, charthick = 3, $
;;            xtitle='acs mag', xrange=[20,28], xstyle = 9, yrange=[0,40]
;;   axis, 0,40, xaxis=1,xrange=[20-dm-kcorr,28-dm-kcorr], xstyle=1, xthick=3, charthick=3

;;plothist, objectnew[memberarr].irac1flux, thick = 3, xthick = 3, ythick = 3, charthick = 3, $
;;          xtitle='[3.6] (AB)', bin=0.5, xrange=[0,50]
;ncon = where(objectnew[memberarr].irac1flux lt 1.5) 
;print, 'fraction of galaxies confused, ch1 flux lt 1.5micro', n_elements(ncon) , n_elements(memberarr)

;----------------------------------------------------------------
;change xaxis
   ;point source objects within sep of the center, detected in irac1 and with a redshift ~1
;;   goodz = where(sphdist(objectnew.ra, objectnew.dec, candra(cand), canddec(cand), /degrees) lt sep and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90  and objectnew.extended lt 1 and objectnew.zphot gt 0.7 and objectnew.zphot lt 1.3) ;and objectnew.nearestdist*60*60/0.6 gt 4.0 

   ;overplot those galaxies with redshifts in range of 1
 ;  oplot, objectnew[goodz].acsmag, objectnew[goodz].acsmag - objectnew[goodz].irac1mag,psym = 2, color = yellowcolor

;----------------------------------------------------------------
;----------------------------------------------------------------
;----------------------------------------------------------------
;----------------------------------------------------------------
;----------------------------------------------------------------
;----------------------------------------------------------------
;----------------------------------------------------------------
;----------------------------------------------------------------
;----------------------------------------------------------------
;calculate star formation rates from 24micron sources
;----------------------------------------------------------------

;start over with finding the 24 micron sources
;but memberarr is only for one cluster, the last one specifically.

;source24 = where(objectnew[memberarr].mips24flux gt 0)
;print, 'source24', n_elements(source24)
;print, memberarr[source24]
;convert flux to luminosity
;;wave = findgen(n_elements(a24)) - findgen(n_elements(a24)) + 24.
;;redshift = findgen(n_elements(a24)) - findgen(n_elements(a24)) + 1.0
;;lum24 = nulnu(wave, objectnew[a24].mips24flux*1E-6, redshift)
;convert 24 micron luminosity to L(IR)  ;chary elbaz
;;lumIR = 0.89*lum24^(1.094)

;;plothist, alog10(lumir), bin=0.1, xtitle = 'log(Lir) of all mips sources in cluster areas', thick=3, xthick=3,ythick=3, charthick=3

;error propagation
;;siglum24 = (objectnew[a24].mips24fluxerr/objectnew[a24].mips24flux) * lum24

;;siglumIR = sqrt( lumIR^2*( (0.33^2/0.89^2) + ( siglum24^2/lum24^2)  ) )

;convert L(IR) to SFR   ;chary elbaz
;;sfr = 1.71E-10*lumIR

;;sigsfr = 1.71E-10*siglumIR

;print, 'mean lum24', mean(lum24)
;print, 'mean siglum24', mean(siglum24)

;print, 'mean lumir', mean(lumIR)
;print, 'mean siglumir', mean(siglumir)

;;totalsfr = total(sfr)
;;sigtotalsfr = total(sigsfr)
;print,'total, toterr, mean, meanerr', totalsfr, sigtotalsfr, mean(sfr), mean(sigsfr)



;what about calculating for only those galaxies on the red sequence.


;----------------------------------------------------------------

;24 detected objects
;a24 = a24[0:ck -1]
; plothyperz, a24, '/Users/jkrick/nep/clusters/gal24.ps'
;print, 'a24', a24
;print, 'nelements(a24)', n_elements(a24)
;print, 'a24', a24;, 'ra', objectnew[[a24]].ra, 'dec', objectnew[[a24]].dec


;-------------------------------------------------------------------
; what about 24micron?
;---------------------------------------------------------------------
;   plothist, objectnew[good].mips24flux, bin=0.1, title='mips24 flux in cluster area', xrange=[0,1000]
;;   for cj = 0, n_elements(good) - 1 do begin
;;      if objectnew[good[cj]].mips24flux gt 0 then begin
;;         a24[ck] = good[cj]
;;         ck = ck + 1
 ;        print, 'ck', good[cj], objectnew[good[cj]].ra, objectnew[good[cj]].dec
;;      endif

;;   endfor

;add in the locations of 24 micron detections
;;b24 = where(mips24all gt 0)
;print, n_elements(b24), '  b24'

;oplot, iracall(b24), iall(b24) - iracall(b24), psym=2, color=redcolor, thick=3

;print,'b24', b24

;b24 = where(objectnew[memberarr].mips24flux gt 0)
;print, n_elements(b24), '  b24'
;oplot, objectnew[memberarr(b24)].irac1mag, objectnew[memberarr(b24)].acsmag - objectnew[memberarr(b24)].irac1mag, psym=2, color=redcolor, thick=3

