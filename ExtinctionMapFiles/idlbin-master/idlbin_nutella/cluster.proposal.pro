pro cluster

restore, '/Users/jkrick/idlbin/object.old.sav'
; where is that big nice far away cluster?
!p.multi = [0, 0, 1]

ps_open, file = "/Users/jkrick/noao/gemini/cluster.z.ps", /portrait, /square

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

memberarr = fltarr(500)
totcount = 0
;-----------------------------------------------------------------------
;what is the average redshift distribution in 100 random regions of the same size as below?
;-----------------------------------------------------------------------

;first find 100 random regions which are not near the edge.
;look at central 13000x13000 region of hst image by pixel numbers
;need 100 pairs of numbers between 1 and 11000
seed = 230
randx = randomu(seed, 50)
randy = randomu(seed, 50)

randx = randx * 11000. + 5500.
randy = randy * 15000. + 2500.


acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits') ;
xyad, acshead, randx, randy, randra, randdec


randz = fltarr(10000)
r = 0
for cand=0, n_elements(randdec) - 1 do begin

   good = where(sqrt((object.ra - randra(cand))^2 + (object.dec-randdec(cand))^2) lt .015 and object.irac1mag gt 0 and object.irac1mag lt 90  )

   for c =0, n_elements(good) -1 do begin
      randz(r) = object[good[c]].zphot
      r = r + 1
   endfor

endfor 

randz = randz[0:r-1]

plothist, randz, xhist, yhist, bin = 0.1, /noplot
yhist = yhist / 50.

;-----------------------------------------------------------------------
;look at all objects within 40" or 800 pixels
;-----------------------------------------------------------------------

;candidate clusters
canddec = [  69.06851 ]

candra = [ 264.89228]
;canddec = [69.04481,  69.06851,  68.98017,  69.087766,  69.019103 ]

;candra = [264.68160, 264.89228, 264.66337, 264.83053, 265.27136 ]

notsogoodra = [265.35535, 265.40978,265.34753,265.17806,  264.95557, 265.46662, 265.23219 ,264.91446, 264.7641]
notsogooddec = [ 69.069154, 68.871061, 68.986163,68.969327,   69.11266,  69.075394, 69.127795, 	69.036086, 68.997594]

m = -0.05
b = fltarr(5)
b = [4.9,4.5,4.5,4.5,5.2]
b = 3.7

for cand=0, n_elements(canddec) - 1 do begin
   print, 'working on cluster', candra(cand), canddec(cand)

   good = where(sqrt((object.ra - candra(cand))^2 + (object.dec-canddec(cand))^2) lt .015 and object.irac1mag gt 0 and object.irac1mag lt 90  )


   print,'n_elements(good)',  n_elements(good)

   plothist, object[good].zphot, bin = 0.1, xrange=[0,2], /xstyle, thick = 3, xthick = 3, ythick = 3, charthick = 3, xtitle='z', yrange=[0,15], ystyle =1

;   plothist, object[good2].zphot, bin = 0.1
   oplot, xhist, yhist, linestyle = 1,thick=3
   oplot, xhist, yhist + 2*sqrt(yhist), linestyle = 2,thick=3

   a = where(object[good].mips24mag ge 0 and object[good].mips24mag lt 90 )
;   print,'gals within area with mips24'
;   print, object[good[a]].rmaga
;   plothyperz, good, strcompress('/Users/jkrick/nep/clusters/cluster'+string(candra(cand))+'.ps')

ps_close, /noprint, /noid

;-------------------------------------------------------------------
;cmd fitting section
;---------------------------------------------------------------------
ps_open, filename='/Users/jkrick/noao/gemini/cluster.cmd.ps',/portrait,/square
   plot, object[good].irac1mag, object[good].acsmag - object[good].irac1mag,psym = 2, yrange=[1,4.5], xtitle = 'IRAC ch1 3.6 micron', ytitle = 'ACS I - ch1 3.6', xrange=[16,24], thick = 3, xthick = 3, ythick = 3, charthick = 3, ystyle=1

;sort them in rmag space to only fit to the better data
   sortindex= sort(object[good].irac1mag)
   irac1sort= object[good].irac1mag[sortindex]
   rsort = object[good].rmaga[sortindex]
   acssort = object[good].acsmag[sortindex]
    
   a = where(acssort - irac1sort gt 2.0 and acssort - irac1sort lt 4.0 )
;a = where(irac1sort lt 22)
;biweight fit (robust for non-gaussian distributions)
;coeff = ROBUST_LINEFIT( object.rmag, object.vr, yfit, sig, coeff_sig)
   if n_elements(a) gt 1 then begin
      coeff1 = ROBUST_LINEFIT(irac1sort(a),acssort(a) - irac1sort(a), yfit1, sig1, coeff_sig1)
      coeff2 = ROBUST_LINEFIT(irac1sort,rsort - irac1sort, yfit2, sig2, coeff_sig2)
      

      print, "sig1, sig2", sig1,sig2
      err = dindgen(n_elements(good)) - dindgen(n_elements(good)) + 1
      start = [-0.02,2.0]
      startr = [-0.02, 4.0]
      sortindex= sort(object[good].irac1mag)
      irac1sort= object[good].irac1mag[sortindex]
      rsort = object[good].rmaga[sortindex]
      acssort = object[good].acsmag[sortindex]
      yfit1sort = yfit1[sortindex]
      yfit2sort = yfit2[sortindex]
      
 ;     resultacs = MPFITFUN('linear',irac1sort, yfit1sort,err, start)
;      print, "resultacs", resultacs
;      resultr= MPFITFUN('linear',irac1sort, yfit2sort,err, startr)
;      print, "resultr", resultr
;      resultz= MPFITFUN('linear',irac1sort, yfit3sort,err, start)
;      print, "resultz", resultz
      
      xvals = findgen(12) +13
 ;     oplot, xvals, ((resultacs(0))*xvals) + resultacs(1), thick = 3    ;, color = colors.orange
 ;     oplot, xvals, ((resultacs(0))*xvals) + resultacs(1)+ 0.5 , thick = 3 ;, color = colors.orange
 ;     oplot, xvals, ((resultacs(0))*xvals) + resultacs(1)- 0.5 , thick = 3 ;, color = colors.orange
      
   endif
   
;  if cand eq 0 then begin
;      print,'working on ',  candra(cand), canddec(cand)

   
;instead just pick a rcs by eye and with bc03 idea of where the color should be given the redshift
;want to fix the slope at -0.09 de lucia etal.
;can come back and fix this later

      oplot, xvals,  m*xvals + b(cand), thick = 3 ;, color = orangecolor
      oplot, xvals, m*xvals + b(cand) + 0.3, thick = 3 ;, color = orangecolor
      oplot, xvals,  m*xvals + b(cand) - 0.3, thick = 3 ;, color = orangecolor

      for count = 0, n_elements(good) - 1 do begin
         limithigh = (m*object[good[count]].irac1mag) + b(cand) + 0.5 
         limitlow = (m*object[good[count]].irac1mag) + b(cand) - 0.5 
         if (object[good[count]].acsmag - object[good[count]].irac1mag LT limithigh )AND ( object[good[count]].acsmag - object[good[count]].irac1mag GT limitlow) then begin
           ; add this member galaxy to a larger array of all member galaxies
           memberarr(totcount) = good[count]
           totcount = totcount + 1
;           print, object[good[count]].rmaga
         endif
      endfor


endfor 

memberarr = memberarr[0:totcount - 1]
print,'memberarr',  n_elements(memberarr)

;plothist, object[memberarr].irac1mag,  thick = 3, xthick = 3, ythick = 3, charthick = 3, xtitle='irac 3.6 flux'
;print, (object[memberarr].irac1mag)
;plothist, object[memberarr].mips24flux, thick = 3, xthick = 3, ythick = 3, charthick = 3, xtitle='mips 24 flux'
;plothist, object[memberarr].rmaga, thick = 3, xthick = 3, ythick = 3, charthick = 3, xtitle='r mag', xrange=[20,30], bin=0.5
;plot, object[memberarr].irac1mag , object[memberarr].wircHmag -  object[memberarr].irac1mag, thick =3, psym = 2, yrange=[-5,5]
ps_close, /noprint, /noid

wirckab = 666.8*10^(object[memberarr].wirckmag/(-2.5))
wirckmagab = -2.5*alog10(wirckab) +8.926
wircjab = 1594.*10^(object[memberarr].wircjmag/(-2.5))
wircjmagab = -2.5*alog10(wircjab) +8.926
      
;plothist, object[memberarr].rmaga - wirckmagab, xrange=[0,10], bin=0.1, charthick = 1, thick = 1
;plothist, object[memberarr].rmaga - wircjmagab, xrange=[0,10], bin=0.5, charthick = 1, thick = 1
;print, n_elements(where(object[memberarr].wircjmag gt 0 and object[memberarr].wircjmag lt 90))

;help, /memory

end




;first need to make a hyperz input with just those objects.
;openw, outlunh, "/Users/jkrick/ZPHOT/hyperz_cat_cluster.txt",/get_lun
;for num = 0, n_elements(good) - 1 do begin
;
;   if object[good(num)].flamjmag gt 0 and object[good(num)].flamjmag ne 99 then begin
;      fab = 1594.*10^(object[good(num)].flamjmag/(-2.5))
;      jmagab = -2.5*alog10(fab) +8.926
;   endif else begin
;      jmagab = object[good(num)].flamjmag
;   endelse
;
;   if object[good(num)].wircjmag gt 0 and object[good(num)].wircjmag ne 99 then begin
;      wircjab = 1594.*10^(object[good(num)].wircjmag/(-2.5))
;      wircjmagab = -2.5*alog10(wircjab) +8.926
;   endif else begin
;      wircjmagab = object[good(num)].wircjmag
;   endelse
;
;   if object[good(num)].wirchmag gt 0 and object[good(num)].wirchmag ne 99 then begin
;      wirchab = 1024.*10^(object[good(num)].wirchmag/(-2.5))
;      wirchmagab = -2.5*alog10(wirchab) +8.926
;   endif else begin
;      wirchmagab = object[good(num)].wirchmag
;   endelse
;
;   if object[good(num)].wirckmag gt 0 and object[good(num)].wirckmag ne 99 then begin
;      wirckab = 666.8*10^(object[good(num)].wirckmag/(-2.5))
;      wirckmagab = -2.5*alog10(wirckab) +8.926
;   endif else begin
;      wirckmagab = object[good(num)].wirckmag
;   endelse
;
;   if object[good(num)].irac1flux lt 0 then err1 = -1. else err1 = 0.05
;   if object[good(num)].irac2flux lt 0 then err2 = -1. else err2 =  0.05;8.926 - 2.5*alog10(0.05E-6*object[good(num)].irac2)
;   if object[good(num)].irac3flux lt 0 then err3 = -1. else err3 =  0.05;8.926 - 2.5*alog10(0.05E-6*object[good(num)].irac3)
;   if object[good(num)].irac4flux lt 0 then err4 = -1. else err4 =  0.05;8.926 - 2.5*alog10(0.05E-6*object[good(num)].irac4)
;   
;   if object[good(num)].imagerra gt 1000. then object[good(num)].imagerra = 1000.
;   if object[good(num)].gmagerra gt 1000. then object[good(num)].gmagerra = 1000.
;   if object[good(num)].rmagerra gt 1000. then object[good(num)].rmagerra = 1000.
;   if object[good(num)].umagerra gt 1000. then object[good(num)].umagerra = 1000.
;   
;
;    printf, outlunh, format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.;;;;;;;;;;;;;2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',$
;                 good(num), object[good(num)].umaga, object[good(num)].gmaga, object[good(num)].rmaga, $
;                 object[good(num)].imaga,  object[good(num)].acsmag,  jmagab, wircjmagab, $
;                 object[good(num)].tmassjmag,  wirchmagab, object[good(num)].tmasshmag, wirckmagab, $
;                 object[good(num)].tmasskmag, object[good(num)].irac1mag,object[good(num)].irac2mag,$
;                 object[good(num)].irac3mag,object[good(num)].irac4mag,   $
;                 object[good(num)].umagerra, object[good(num)].gmagerra, $
;                 object[good(num)].rmagerra, object[good(num)].imagerra, object[good(num)].acsmagerr, object[good(num)].flamjmagerr, $
;                 object[good(num)].wircjmagerr,0.02, object[good(num)].wirchmagerr, 0.02, object[good(num)].wirckmagerr,$
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
;;plot SED's of keeper objects
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
;adxy,  acshead,object[good].ra, object[good].dec, xcenter, ycenter
;
;
;openw, outlun, '/Users/jkrick/nep/cluster.reg', /get_lun
;for count = 0, n_elements(good) - 1 do begin
;   if object[good(count)].zphot gt 0.9 and object[good(count)].zphot lt 1.3 then begin
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
;;output the acs image of the object
;   print, "working on n ", n
;
;   if object[good(n)].acsmag gt 0 and object[good(n)].acsmag lt 90 then begin
;;     print,    xcenter(n) - size/0.05, xcenter(n)+ size/0.05  ,ycenter(n) -size/0.05, ycenter(n)+size/0.05
;      acsdata = mrdfits('/Users/jkrick/hst/raw/wholeacs.fits', range=[ycenter(n) -size/0.05, ycenter(n)+size/0.05])
;      plotimage, xrange=[xcenter(n) - size/0.05, xcenter(n)+ size/0.05],$
;;                 yrange=[ycenter(n) -size/0.05, ycenter(n)+size/0.05], $
;                 bytscl(acsdata, min = -0.01, max = 0.1),$
;                 /preserve_aspect, /noaxes, ncolors=60
;      xyouts, xcenter(n)- 0.6*size/0.05, -10., strcompress(string(good(n))+ '    ACS' + string(object[good(n)].acsmag)),charthick = 3
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
;        xtitle = "log(microns)", ytitle = "log(flux(microjansky))", title=strcompress("object " + string(idz(n)) +string(zphota(n)) + string(proba(n)) + string(sptkey(specta(n)-1)))

;   errplot, alog10(x), alog10(y - yerr), alog10(y + yerr)
;;
;   readcol,strcompress('/Users/jkrick/ZPHOT/' + string(idz(n)) + '.spe', /remove_all),x2, y2,format="A",/silent
;   oplot, alog10((x2/1E4)), alog10(y2), thick = 3
;   
;endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;canddec = [69.04481,  69.06851,  69.069154,  68.871674,  69.108844,  68.98017,  69.090313,  69.021478,  69.075394, 69.127795, 	69.036086, 68.997594, 	 68.986163,  68.969327,68.99948]

;candra = [264.68160, 264.89228, 265.35535, 265.39991, 264.95902, 264.66337, 264.83102, 265.27616, 265.46662, 265.23219 ,264.91446, 264.7641, 265.34753, 265.17806,265.04271]

 ;  plot, object[good].irac1mag, object[good].rmaga - object[good].irac1mag,psym = 2, yrange=[0,8], xtitle = 'ch 1', ytitle = 'r - ch1', xrange=[16,24], thick = 3, xthick = 3, ythick = 3, charthick = 3

 ;  oplot, xvals, ((resultr(0))*xvals) + resultr(1), thick = 3    ;, color = colors.orange
 ;  oplot, xvals, ((resultr(0))*xvals) + resultr(1)+ sig2 , thick = 3 ;, color = colors.orange
 ;  oplot, xvals, ((resultr(0))*xvals) + resultr(1)- sig2 , thick = 3 ;, color = colors.orange

 ;  oplot, xvals, 27.4 - xvals, thick = 3 , linestyle = 2

;   plot, object[good].irac1mag, object[good].zmagbest - object[good].irac1mag,psym = 2, yrange=[0,8], xtitle = 'ch 1', ytitle = 'z - ch1', xrange=[16,24], thick = 3, xthick = 3, ythick = 3, charthick = 3

;   oplot, xvals, ((resultz(0))*xvals) + resultz(1), thick = 3        ;, color = colors.orange
;   oplot, xvals, ((resultz(0))*xvals) + resultz(1)+ sig3 , thick = 3 ;, color = colors.orange
;   oplot, xvals, ((resultz(0))*xvals) + resultz(1)- sig3 , thick = 3 ;, color = colors.orange
   
;   endif
;   if cand eq 1 then begin
;      print,'working on ',  candra(cand), canddec(cand)
;      oplot, xvals,  (-0.09)*xvals + 4.5, thick = 3 , color = orangecolor
;      oplot, xvals,  (-0.09)*xvals + 5.0, thick = 3 , color = orangecolor
;      oplot, xvals,  (-0.09)*xvals + 4.0, thick = 3 , color = orangecolor
;   endif
;   if cand eq 2 then begin
;;      print,'working on ',  candra(cand), canddec(cand)
;      oplot, xvals,  (-0.09)*xvals + 4.5, thick = 3 , color = orangecolor
;      oplot, xvals,  (-0.09)*xvals + 5.0, thick = 3 , color = orangecolor
;      oplot, xvals,  (-0.09)*xvals + 4.0, thick = 3 , color = orangecolor
;   endif
;   if cand eq 3 then begin
;;      print,'working on ',  candra(cand), canddec(cand)
;      oplot, xvals,  (-0.09)*xvals + 4.5, thick = 3 , color = orangecolor
;      oplot, xvals,  (-0.09)*xvals + 5.0, thick = 3 , color = orangecolor
;      oplot, xvals,  (-0.09)*xvals + 4.0, thick = 3 , color = orangecolor
;   endif
;  if cand eq 4 then begin
;;      print,'working on ',  candra(cand), canddec(cand)
;      oplot, xvals,  (-0.09)*xvals + 5.2, thick = 3 , color = orangecolor
;      oplot, xvals,  (-0.09)*xvals + 5.7, thick = 3 , color = orangecolor
;      oplot, xvals,  (-0.09)*xvals + 4.7, thick = 3 , color = orangecolor
;   endif




;   for i = 0 , n_elements(good) -1 do begin
;      if object[good[i]].acsmag - object[good[i]].irac1mag gt 5 then print, "red", good[i]
;;   endfor
;   b = where( object[good].acsmag - object[good].irac1mag ge 5)
;plothyperz, good[b]
;--------------------------------------------------------
;openw, outlun, '/Users/jkrick/nep/clusters/rand.reg', /get_lun
;for i = 0,n_elements(randx) - 1 do printf, outlun, 'circle( ',randx(i), randy(i),' 1200)'
;close, outlun
;free_lun, outlun
