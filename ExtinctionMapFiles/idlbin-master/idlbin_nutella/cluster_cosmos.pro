pro cluster_cosmos

  !p.multi=[0,1,1]
  ps_open, filename= '/Users/jkrick/nep/clusters/cosmos/dist.ps', /portrait,  /color,/square, xsize=4, ysize=4
  vsym, /polygon, /fill
  
  redcolor = FSC_COLOR("Red", !D.Table_Size-2)
  bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
  greencolor = FSC_COLOR("Green", !D.Table_Size-4)
  yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
  cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
  orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
  purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

;color k corrections
;ignoring cluster 0
  colorkcor=[-0.1,-0.1,-0.1,0.06]
  kcorg=[-.06,-.06,-.06,.04]
  kcork=[.04,.04,.04,-.02]
  kcorb = [-.12,-.12,-.12,0.07]
  colorkcorb=[0,0,0,0]

;  colorkcor=[0,-0.1,-0.1,-0.21,-0.1,0.12,0.06]
;  kcorg=[0,-.06,-.06,-.15,-.06,.08,.04]
;  kcork=[0,.04,.04,.06,.04,-.04,-.02]
;distance modulus
  dm=[38.73,38.73,38.73,38.07]
;  dm=[39.70,38.73,38.73,38.92,38.73,37.80,38.07]

;setup variable for the colors of all clusters
  gall = fltarr(210000)
  kall = fltarr(210000)
  ball = fltarr(210000)
  j = 0l
  allcount =0l
  allmemberg = fltarr(210000)
  allmemberk = fltarr(210000)
  allmemberb = fltarr(210000)
;----------------------------------------------------------------
;read in the cluster data
;----------------------------------------------------------------
  readcol, '/Users/jkrick/nep/clusters/cosmos/cluster.wcs.txt', targetra, targetdec, cosclusternum, z, filename, format="F,F,I,F,A", /silent
;  print, 'n elements z', n_elements(z)
  nclusters = n_elements(z)
  for i = 0, n_elements(z) -1 do begin
                                ;skip the first cluster at z=0.18
     readlargecol, filename(i),      id42, image_tile42,       ra42,      dec42,    x_pix42,    y_pix42,   i_fwhm42,    i_max42,   i_star42,         i_mag_auto42,        auto_offset42,auto_flag42,    u_mag42,u_mag_err42,    b_mag42,b_mag_err42,    v_mag42,v_mag_err42,    g_mag42,g_mag_err42,    r_mag42,r_mag_err42,    i_mag42,i_mag_err42,    z_mag42,z_mag_err42,    k_mag42,k_mag_err42,         i_cfht_mag42,     i_cfht_mag_err42,         u_sdss_mag42,     u_sdss_mag_err42,         g_sdss_mag42,     g_sdss_mag_err42,         r_sdss_mag42,     r_sdss_mag_err42,         i_sdss_mag42,     i_sdss_mag_err42,         z_sdss_mag42,     z_sdss_mag_err42,f814w_mag42,      f814w_mag_err42,nb816_mag42,      nb816_mag_err42,      ebv42,    zphot42,        zerr_68_min42,        zerr_68_max42,        zerr_95_min42,        zerr_95_max42,    tphot42,    rphot42,    chisq42,   nbands42,       mv42,      d9542,  logmass42,     star42,   b_mask42,   v_mask42,   i_mask42,   z_mask42, blend_mask42, format="A", /silent

;     print, 'gmag42', n_elements(g_mag42)


;how deep is the K band, no H band photometry?  

;!p.multi=[0,2,1]
;plothist, k_mag42, bin=0.1, xrange=[10,30]
;plothist, k_mag104, bin=0.1, xrange=[10,30]

;what about B instead of g?
;will get the balmer break

 ;  plothist, b_mag42, bin=0.1, xrange=[10,30]
 ;----------------------------------------------------------------
;select likely cluster members
;----------------------------------------------------------------

;must have k and g band detections 
;go for the first, simplest, with photz in the range of projected cluster redshift near to the center
     sep = 0.0754;0.19
;     start = where( b_mag42 gt 0 and b_mag42 lt 90 and k_mag42 gt 0 and k_mag42 lt 90 $;and zphot42 gt 0 and zphot42 lt 3*z(i) $
;                    and sqrt((ra42 - targetra(i))^2 + (dec42-targetdec(i))^2) lt sep)   ;   
     start= where( b_mag42 gt 0 and b_mag42 lt 90 and k_mag42 gt 0 and k_mag42 lt 90 $;and zphot42 gt 0 and zphot42 lt 3*z(i) $
                    and sphdist(ra42, dec42, targetra(i), targetdec(i), /degrees)  lt sep)   ;   


;put the starting members in an array of all clusters
     for k = 0l, n_elements(start) - 1 do begin
        gall[j] = g_mag42(start(k)) + kcorg(i)
        kall[j] = k_mag42(start(k)) + kcork(i)
        ball[j] = b_mag42(start(k)) + kcorb(i)
        
        j = j + 1
     endfor
     
;what do we expect the color to be (BC03)?  2.3-2.45

;----------------------------------------------------------------
;make CMD for each cluster individually 
;----------------------------------------------------------------
;use color k correction to put them all at the same redshift z=0.1
;;     plot, k_mag42(start) + kcork(i) , b_mag42(start) - k_mag42(start) + colorkcor(i), psym = 8, symsize=0.3,xrange=[14,25], $
;;           xtitle='K mag (AB)', ytitle='B - K @z=0.1', xstyle=9, yrange=[-1,5], ystyle = 1, thick=3, xthick=3, ythick=3, charthick=3
;;     axis, 0, 5, xaxis=1, xrange=[14-dm(i), 25-dm(i)], xstyle=1, xthick=3, charthick=3
     
     
;----------------------------------------------------------------
;fit the cmd
;----------------------------------------------------------------
;sort them in rmag space to only fit to the better data
     sortindex= sort(k_mag42(start))
     ksort= k_mag42(start[sortindex])
     gsort = g_mag42(start[sortindex])
     bsort = b_mag42(start[sortindex])
     
     a = where(bsort - ksort gt 1.0 and bsort - ksort lt 4.0 )
;a = where(ksort lt 22)
;biweight fit (robust for non-gaussian distributions)
;coeff = ROBUST_LINEFIT( objectnew.rmag, objectnew.vr, yfit, sig, coeff_sig)
     if n_elements(a) gt 1 then begin
        coeff1 = ROBUST_LINEFIT(ksort(a),bsort(a) - ksort(a), yfit1, sig1, coeff_sig1)
        coeff2 = ROBUST_LINEFIT(ksort,bsort - ksort, yfit2, sig2, coeff_sig2)
        
        print, "sig1, sig2", sig1,sig2
        err = dindgen(n_elements(ksort)) - dindgen(n_elements(ksort)) + 1
        startcoeff = [-0.02,2.0]
        startr = [-0.02, 4.0]
;   sortindex= sort(objectnew[good].irac1mag)
;   ksort= objectnew[good].irac1mag[sortindex]
;   gsort = objectnew[good].acsmag[sortindex]
   
        yfit1sort = yfit1[sortindex]
        yfit2sort = yfit2[sortindex]
        result = MPFITFUN('linear',ksort, yfit1sort,err, startcoeff,/quiet)
        result2 = MPFITFUN('linear',ksort, yfit2sort,err, startcoeff,/quiet)
        
     endif
   
;     oplot, ksort, ((result(0))*ksort) + result(1), thick = 2  ;, color = colors.orange
;     oplot, ksort, ((result(0))*ksort) + result(1)+ sig1 , thick = 2 ;, color = colors.orange
;     oplot, ksort, ((result(0))*ksort) + result(1)- sig1 , thick = 2 ;, color = colors.orange


     ;fit by eye.
     m = -0.12
     b = [4.8,4.8,4.8,4.6]
     interceptall = 4.7
;;     oplot, ksort, (m*ksort) + b(i), thick = 2 , color = orangecolor
;;     oplot, ksort, (m*ksort) + b(i) + sig1, thick = 2 , color = orangecolor
;;     oplot, ksort, (m*ksort) + b(i) - sig1, thick = 2 , color = orangecolor
     

;----------------------------------------------------------------
;use cmd color requirement for membership
;----------------------------------------------------------------
     singlecount = 0
     singlememberarr = fltarr(n_elements(start))
     for count = 0, n_elements(gsort) - 1 do begin
        limithigh = (m*ksort(count)) + b(i) + 0.5
        limitlow = (m*ksort(count))+ b(i) - 0.5
;        limithigh = (result(0)*ksort(count)) + result(1) + sig1
;        limitlow = (result(0)*ksort(count))+ result(1) - sig1

        if (bsort(count) - ksort(count)) LT limithigh  AND  (bsort(count) - ksort(count)) GT limitlow then begin
           
   ; add this member galaxy to a larger array of all member galaxies
   ;don't bother with indexes, save the actual values, is easier to keep track of
           allmemberg(allcount) = gsort(count) + kcorg(i)
           allmemberk(allcount) = ksort(count) + kcork(i)
           allmemberb(allcount) = bsort(count) + kcorb(i)
           allcount = allcount + 1
         ;also want to know members of this cluster only
           singlememberarr(singlecount) =[count]
           singlecount = singlecount  + 1

        endif
     endfor

     singlememberarr=singlememberarr[0:singlecount-1]
     
     print, "start, singlecount ", n_elements(start), singlecount
;----------------------------------------------------------------
;make color distribution for each cluster individually of cmd member galaxies
;----------------------------------------------------------------

;;     plothist, k_mag42(start), xrange=[14,25], xtitle='K mag (AB)', ytitle='N', yrange=[0,200], ystyle=1,$
;;               xstyle=9,bin=0.5, charthick=3, xthick=3, ythick=3, thick=3
;;     axis, 0, 200, xaxis=1, xrange=[14-dm(i), 25-dm(i)], xstyle=1, xthick=3, charthick=3



  endfor


;  print, ' total number ', j
  gall = gall[0:j-1]
  kall = kall[0:j-1]
  ball = ball[0:j-1]
  allmemberg=allmemberg[0:allcount-1]
  allmemberk=allmemberk[0:allcount-1]
  allmemberb = allmemberb[0:allcount-1]

save, allmemberk, filename= '/Users/jkrick/nep/clusters/cosmos_allmemberk.txt'
;----------------------------------------------------------------
;make this work for all low-z clusters combined
;make combined cmd
;make combined rcs magnitude distribution
;----------------------------------------------------------------
dm = 38.3
dm = dm - 0.1   ;filter size correction

;;plot, kall, ball - kall, psym = 8, symsize=0.3,xrange=[14,25], $ ;title = 'composite',$
;;      xtitle='K (AB)', ytitle='B - K ', xstyle=9, yrange=[-1,5], ystyle = 1, thick=3, xthick=3, ythick=3, charthick=3
;;axis, 0, 5, xaxis=1, xrange=[14-dm, 25-dm], xstyle=1, xthick=3, charthick=3

;;oplot, findgen(100), (m*findgen(100)) + interceptall, thick = 2          ;, color = orangecolor
;;oplot, findgen(100), (m*findgen(100)) + interceptall + 0.5, thick = 2   ;, color = orangecolor
;;oplot, findgen(100), (m*findgen(100)) + interceptall - 0.5, thick = 2   ;, color = orangecolor

;;im_hessplot, kall, ball-kall,/notsquare, nbin2d=15 ,xrange=[14,25],position = [0.1 ,0.1 , 0.9,0.9], charsize=1, $
;;      xtitle='K (AB)', ytitle='B - K ', xstyle=9, yrange=[-1,5], ystyle = 1, thick=3, xthick=3, ythick=3, charthick=3;
;;axis, 0, 5, xaxis=1, xrange=[14-dm, 25-dm], xstyle=1, xthick=3, charthick=3

;;oplot, findgen(100), (m*findgen(100)) + interceptall, thick = 2          ;, color = orangecolor
;;oplot, findgen(100), (m*findgen(100)) + interceptall + 0.5, thick = 2   ;, color = orangecolor
;;oplot, findgen(100), (m*findgen(100)) + interceptall - 0.5, thick = 2   ;, color = orangecolor

;  plothist, kall, xrange=[14,25], xtitle='K mag (AB)', ytitle='N', yrange=[0,900], ystyle=1,$
;               xstyle=9,bin=0.5, charthick=3, xthick=3, ythick=3, thick=3
;  axis, 0, 900, xaxis=1, xrange=[14-dm, 25-dm], xstyle=1, xthick=3, charthick=3
;dm is the dm of z=0.1

;------
;  plot, allmemberk, allmemberb - allmemberk, psym = 8, symsize=0.3,xrange=[14,25], title = 'composite',$
;        xtitle='K mag (AB)', ytitle='B - K @z=0.1', xstyle=9, yrange=[-1,5], ystyle = 1, thick=3, xthick=3, ythick=3, charthick=3
 ; axis, 0, 5, xaxis=1, xrange=[14-dm, 25-dm], xstyle=1, xthick=3, charthick=3


  plothist, allmemberk, xrange=[14,25], xtitle='K  (AB)', ytitle='Number', yrange=[0,200], ystyle=1,$
              xstyle=9,bin=0.5, charthick=3, xthick=3, ythick=3, thick=3
  axis, 0, 200, xaxis=1, xrange=[14-dm, 25-dm], xstyle=1, xthick=3, charthick=3

oplot, [16.7,16.7], [0,2000], linestyle=2
oplot, [18.2,18.2], [0,2000], linestyle=2
oplot, [19.7,19.7], [0,2000], linestyle=2

;----------------------------------------------------------------
;want a member array chosen from the composite CMD
;----------------------------------------------------------------
compmemb = fltarr(n_elements(kall))
compmemk = fltarr(n_elements(kall))
totcompcount = 0
for compcount = 0.D, n_elements(kall) - 1 do begin
   limithigh = (m*kall(compcount)) + interceptall +0.5    ; 2*sig1;0.5 
   limitlow = (m*kall(compcount)) + interceptall - 0.5    ;2*sig1;0.5 
;   print, "lh, lw", limithigh,limitlow, kall(compcount), ball(compcount) - kall(compcount) 
   if (ball(compcount) - kall(compcount) LT limithigh ) AND ( ball(compcount) - kall(compcount)  GT limitlow) then begin
      compmemb(totcompcount) = ball(compcount)
      compmemk(totcompcount) = kall(compcount)
      totcompcount = totcompcount + 1
;      print, 'totcount', totcompcount
   endif
endfor

compmemb = compmemb[0:totcompcount - 1]
compmemk = compmemk[0:totcompcount - 1]



;-----------------------------------------------------------------------
;what is the average number of galaxies in 100 random regions of the same size as above?
;-----------------------------------------------------------------------

;look at all cosmos galaxies with the correct k band fluxes
 readlargecol, '/Users/jkrick/nep/clusters/cosmos/all.txt',      id,     ra,      dec,   b_mag,b_mag_err,     k_mag,k_mag_err,  zphot,        zerr_68_min,        zerr_68_max,        zerr_95_min,        zerr_95_max,    tphot,     chisq,   nbands,       mv,     logmass,     star,   format="A", /silent


seed = 230
nrand = 10
randra = randomu(seed, nrand)
randdec = randomu(seed, nrand)

;print, max(ra), min(ra)
;print, max(dec), min(dec)

randra = randra*1.3 + 149.5
randdec = randdec*1.3 + 1.5
gc = 0
goodbrightarr = fltarr(nrand)
gooddimarr = fltarr(nrand)
totback = 0
comosbackmemberarr = fltarr(nrand*1000.)
for cand=0, n_elements(randdec) - 1 do begin   ;for each random area
   bc = 0
   dc = 0
  for s = 0, n_elements(ra) - 1 do begin
     limithigh = (m*k_mag(s)) + interceptall+ 0.5
     limitlow = (m*k_mag(s))+ interceptall - 0.5

     if sphdist(ra(s), dec(s), randra(cand), randdec(cand), /degrees)  lt sep  and k_mag(s) le 18.2 and k_mag(s) gt 16.2 $
        and b_mag(s) - k_mag(s) LT limithigh  AND  (b_mag(s) - k_mag(s)) GT limitlow then begin
        bc = bc + 1
     endif
     if sphdist(ra(s), dec(s), randra(cand), randdec(cand), /degrees)  lt sep  and k_mag(s) lt 20.2 and k_mag(s) gt 18.2 $
        and b_mag(s) - k_mag(s) LT limithigh  AND  (b_mag(s) - k_mag(s)) GT limitlow then begin
        dc = dc + 1
     endif

     if sphdist(ra(s), dec(s), randra(cand), randdec(cand), /degrees)  lt sep  $
        and b_mag(s) - k_mag(s) LT limithigh  AND  (b_mag(s) - k_mag(s)) GT limitlow then begin

        comosbackmemberarr[totback] = k_mag(s)
        totback = totback + 1
     endif

  endfor

   goodbrightarr(cand) = bc
   gooddimarr(cand) = dc

endfor 

goodcolorbright = mean(goodbrightarr)
goodcolordim = mean(gooddimarr)
print, "goodcolorbright, goodcolordim " , goodcolorbright, goodcolordim
print, "stddev", stddev(goodbrightarr), stddev(gooddimarr)

comosbackmemberarr = comosbackmemberarr[0:totback - 1]
save, comosbackmemberarr, filename = '/Users/jkrick/nep/clusters/cosmos/backmemberarr.txt'
;----------------------------------------------------------------
;make a measurement
;----------------------------------------------------------------
;break between bright and dim is Mk = -20, or mk = 18.2  /pm 2mags on each side
 bright = where(allmemberk lt 18.2 and allmemberk gt 16.2)
 dim = where(allmemberk ge 18.2 and allmemberk lt 20.2)

print, 'results, bright, dim (uncorrected): ', n_elements(bright), n_elements(dim)
print, 'results, bright, dim : ', n_elements(bright) - nclusters*goodcolorbright, n_elements(dim)-nclusters*goodcolordim

brightcomp = where(compmemk lt 18.2 and compmemk gt 16.2)
faintcomp = where(compmemk ge 18.2 and compmemk lt 20.2)

print, 'composite uncorrected ', n_elements(brightcomp), n_elements(faintcomp)

ps_close, /noprint,/noid



end


;----------------------------------------------------------------
;change x-axis
;  plot, ball, ball - kall, psym = 8, symsize=0.3,xrange=[16,28], title = 'composite',$
;        xtitle='B mag (AB)', ytitle='B - K @z=0.1', xstyle=9, yrange=[-1,5], ystyle = 1, thick=3, xthick=3, ythick=3, charthick=3
;  axis, 0, 5, xaxis=1, xrange=[16-dm, 28-dm], xstyle=1, xthick=3, charthick=3


;  plothist, ball, xrange=[16,28], xtitle='B mag (AB)', ytitle='N', yrange=[0,600], ystyle=1,$
;               xstyle=9,bin=0.5, charthick=3, xthick=3, ythick=3, thick=3
;  axis, 0, 600, xaxis=1, xrange=[16-dm, 28-dm], xstyle=1, xthick=3, charthick=3
;dm is the dm of z=0.1

;------
;  plot, allmemberb, allmemberb - allmemberk, psym = 8, symsize=0.3,xrange=[16,28], title = 'composite',$
;        xtitle='B mag (AB)', ytitle='B - K @z=0.1', xstyle=9, yrange=[-1,5], ystyle = 1, thick=3, xthick=3, ythick=3, charthick=3
;  axis, 0, 5, xaxis=1, xrange=[16-dm, 28-dm], xstyle=1, xthick=3, charthick=3


;  plothist, allmemberb, xrange=[16,28], xtitle='B mag (AB)', ytitle='N', yrange=[0,400], ystyle=1,$
;               xstyle=9,bin=0.5, charthick=3, xthick=3, ythick=3, thick=3
;  axis, 0, 400, xaxis=1, xrange=[16-dm, 28-dm], xstyle=1, xthick=3, charthick=3
