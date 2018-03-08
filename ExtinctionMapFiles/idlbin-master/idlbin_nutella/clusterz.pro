pro cluster

restore, '/Users/jkrick/idlbin/object.sav'
; where is that big nice far away cluster?
!p.multi = [0, 0, 1]

ps_open, file = "/Users/jkrick/nep/clusters/cluster.zl.ps", /portrait, xsize = 6, ysize = 6,/color

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

cmdacs = [3.2,2.4,2.8,2.8,1.5,3.7,2.8,3.5,0.9,2.4]
cmdr = [5.2,4.3,4.7,4.7,3.0,6.2,4.7,5.7,2.0,4.3]


;-----------------------------------------------------------------------
;what is the average redshift distribution in 100 random regions of the same size as above.
;-----------------------------------------------------------------------

;first find 100 random regions which are not near the edge.
;look at central 13000x13000 region of hst image by pixel numbers
;need 100 pairs of numbers between 1 and 11000
seed = 230
randx = randomu(seed, 50)
randy = randomu(seed, 50)

randx = randx * 11000. + 5500.
randy = randy * 15000. + 2500.

openw, outlun, '/Users/jkrick/nep/clusters/rand.reg', /get_lun
for i = 0,n_elements(randx) - 1 do printf, outlun, 'circle( ',randx(i), randy(i),' 1200)'
close, outlun
free_lun, outlun

acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits');
xyad, acshead, randx, randy, randra, randdec

;randra=[265.45171,265.19677,264.90424,265.27899,264.76746,265.00615,265.3371]
;randdec=[69.103101,69.054803,69.086835,68.971012,68.994038,68.918891,68.853382]

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

;plot, xhist, yhist, xtitle = 'z', title = "random"


;plot, object[good].irac1mag, object[good].acsmag - object[good].irac1mag,psym = 2, yrange=[0,8], xtitle = 'ch 1', ytitle = 'acs - ch1', xrange=[16,24]

;-----------------------------------------------------------------------
;look at all objects within 40" or 800 pixels
;candidate clusters
;-----------------------------------------------------------------------

canddec = [69.04481,  69.06851,   69.11266,  68.98017,  69.087766,  69.021478,  69.075394, 69.127795, 	69.036086, 68.997594	   ]

candra = [264.68160, 264.89228,  264.95557, 264.66337, 264.83053, 265.27616, 265.46662, 265.23219 ,264.91446, 264.7641 ]

notsogoodra = [265.35535, 265.40978,265.34753,265.17806]
notsogooddec = [ 69.069154, 68.871061, 68.986163,68.969327]


acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits') ;
adxy, acshead, candra,canddec , xcenter, ycenter

fits_read, '/Users/jkrick/spitzer/IRAC/ch1/mosaic.fits',iracdata, iracheader
adxy, iracheader, candra, canddec, iracxcenter, iracycenter

size = 0.015
for cand=0, n_elements(canddec) - 1 do begin

   good = where(sqrt((object.ra - candra(cand))^2 + (object.dec-canddec(cand))^2) lt size and object.irac1mag gt 0 and object.irac1mag lt 90  )

   print, n_elements(good)

   plothist, object[good].zphot, bin = 0.1, xrange=[0,2], xtitle = 'z', title = strcompress(string(candra(cand)) + string(canddec(cand))), yrange=[0,15], ystyle =1, charthick = 3, xthick=3,ythick=3,thick=3
   oplot, xhist, yhist, linestyle = 1,thick=3
   oplot, xhist, yhist + 2*sqrt(yhist), linestyle = 2,thick=3

;-------------------------------------------------------------------
;output the acs image of the object
;-------------------------------------------------------------------

   if  ycenter(cand) - size*60*60/0.05 gt 0 and ycenter(cand) + size*60*60/0.05 lt 20300 then begin;and object[keeper(cand)].acsmag lt 90 
      acsdata = mrdfits('/Users/jkrick/hst/raw/wholeacs.fits', range=[ycenter(cand) -size*60*60/0.05, ycenter(cand)+size*60*60/0.05])
      plotimage, xrange=[xcenter(cand) - size*3600/0.05, xcenter(cand)+ size*3600/0.05],$
;                 yrange=[ycenter(cand) -size/0.05, ycenter(cand)+size/0.05], $
                 bytscl(acsdata, min = -0.01, max = 0.1),$
                 /preserve_aspect, /noaxes, ncolors=60
;      xyouts, xcenter(cand)- 0.6*size/0.05, -10., string(keeper(cand)),charthick = 3
   endif

;   if object[keeper(cand)].irac1mag gt 0 then begin;and object[keeper(cand)].acsmag lt 90 
      plotimage, xrange=[iracxcenter(cand) - size*3600/0.6, iracxcenter(cand)+ size*3600/0.6],$
                 yrange=[iracycenter(cand) -size*3600/0.6, iracycenter(cand)+size*3600/0.6], $
                 bytscl(iracdata, min =0.04, max = 0.1),$
                 /preserve_aspect, /noaxes, ncolors=60
;      xyouts, xcenter(cand)- 0.6*size/0.05, -10., string(keeper(cand)),charthick = 3
;   endif

;-------------------------------------------------------------------
;CMD and fitting section
;---------------------------------------------------------------------

   plot, object[good].irac1mag, object[good].acsmag - object[good].irac1mag,psym = 2, yrange=[0,8], xtitle = 'ch 1', ytitle = 'acs - ch1', xrange=[16,25], charthick = 3, xthick=3,ythick=3,thick=3

;sort them in rmag space to only fit to the better data
   sortindex= sort(object[good].irac1mag)
   irac1sort= object[good].irac1mag[sortindex]
   rsort = object[good].rmaga[sortindex]
   acssort = object[good].acsmag[sortindex]

   a = where(acssort - irac1sort gt 2.0 and acssort - irac1sort lt 4.0 )

;biweight fit (robust for non-gaussian distributions)
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
      
      
      resultacs = MPFITFUN('linear',irac1sort, yfit1sort,err, start)
      print, "resultacs", resultacs
      resultr= MPFITFUN('linear',irac1sort, yfit2sort,err, startr)
      print, "resultr", resultr
      
      xvals = findgen(27)
      oplot, xvals, ((resultacs(0))*xvals) + resultacs(1), thick = 3     ;, color = colors.orange
      oplot, xvals, ((resultacs(0))*xvals) + resultacs(1)+ sig1 , thick = 3 ;, color = colors.orange
      oplot, xvals, ((resultacs(0))*xvals) + resultacs(1)- sig1 , thick = 3 ;, color = colors.orange
      oplot, xvals, 29. - xvals, thick = 3 , linestyle = 2

   endif
   
   oplot, findgen(27), -0.09*findgen(27) + cmdacs(cand) + 1.89, color = redcolor
   if cand eq 4 then oplot, findgen(27), -0.09*findgen(27) + 2.8 + 1.89, color = redcolor
   if cand eq 5 then oplot, findgen(27), -0.09*findgen(27) + 4.7+ 1.89, color = redcolor

   for i = 0 , n_elements(good) -1 do begin
      if object[good[i]].acsmag - object[good[i]].irac1mag gt 5 then print, "red", good[i]
   endfor
   b = where( object[good].acsmag - object[good].irac1mag ge 5)
;plothyperz, good[b]
;--------------------------------------------------------

   plot, object[good].irac1mag, object[good].rmaga - object[good].irac1mag,psym = 2, yrange=[0,8], xtitle = 'ch 1', ytitle = 'r - ch1', xrange=[16,24], charthick = 3, xthick=3,ythick=3,thick=3
   oplot, xvals, ((resultr(0))*xvals) + resultr(1), thick = 3        ;, color = colors.orange
   oplot, xvals, ((resultr(0))*xvals) + resultr(1)+ sig2 , thick = 3 ;, color = colors.orange
   oplot, xvals, ((resultr(0))*xvals) + resultr(1)- sig2 , thick = 3 ;, color = colors.orange
   
   oplot, xvals, 27.4 - xvals, thick = 3 , linestyle = 2
   
   oplot, findgen(27), -0.09*findgen(27) + cmdr(cand) + 1.89, color = redcolor
   if cand eq 4 then oplot, findgen(27), -0.09*findgen(27) + 4.7 + 1.89, color = redcolor
   if cand eq 5 then oplot, findgen(27), -0.09*findgen(27) + 7.5+ 1.89, color = redcolor

;--------------------------------------------------------


;   good2 = where(sqrt((object.ra - candra(cand))^2 + (object.dec-canddec(cand))^2) lt size  )

;   plot, object[good2].rmaga, object[good2].gmaga - object[good2].rmaga,psym = 2, yrange=[-2,4], xtitle = 'r', ytitle = 'g - r', xrange=[18,28], charthick = 3, xthick=3,ythick=3,thick=3
endfor 



ps_close, /noprint, /noid

end
