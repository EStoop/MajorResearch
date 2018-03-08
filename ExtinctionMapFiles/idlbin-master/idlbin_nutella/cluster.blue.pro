pro clusterblue

restore, '/Users/jkrick/idlbin/object.old.sav'

!p.multi = [0, 1, 1]
;!p.multi = [0, 1, 2]

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)


memberarr = fltarr(50000)
totcount = 0.D
ck=0
a24 = fltarr(500)
;setup variable for the colors of all clusters
iall = fltarr(210000)
iracall = fltarr(210000)
mips24all = fltarr(210000)
  j = 0l
  
;I know that the photometry is bad for these objects and they end up with wrong,very red colors, so don't even let them into the sample
;this happens because they are not detections in irac, but are near a bright source, so they will get flux in the aperture as
;determined by the acs photometry.
badnum=[21735,21820,21854,22133,22169,22707,22761,23225]
objectnew[badnum].irac1mag = -99.0


;-----------------------------------------------------------------------
;what is the average redshift distribution in 50 random regions of the same size as below?
;-----------------------------------------------------------------------
sep =0.017 ;(55")  0.0062; (720kpc)

;first find 50 random regions which are not near the edge.
;look at central 13000x13000 region of hst image by pixel numbers                                   
;need 100 pairs of numbers between 1 and 11000
seed = 420
nrand = 50
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
for cand=0, n_elements(randdec) - 1 do begin   ;for each random area
   good = where(sphdist(objectnew.ra, objectnew.dec, randra(cand), randdec(cand), /degrees) lt sep and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90  )
   for c =0, n_elements(good) -1 do begin
      randz(r) = objectnew[good[c]].zphot
      r = r + 1
   endfor
   backcount[cand] = n_elements(good)
   goodcolorbright = where(sphdist(objectnew.ra, objectnew.dec, randra(cand), randdec(cand), /degrees) lt sep and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.irac1mag lt 22.6 and objectnew.irac1mag gt 20.6 and objectnew.acsmag gt 0 and objectnew.acsmag lt 90, brightcount)

   goodcolordim = where(sphdist(objectnew.ra, objectnew.dec, randra(cand), randdec(cand), /degrees) lt sep and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.irac1mag lt 24.6 and objectnew.irac1mag gt 22.6 and objectnew.acsmag gt 0 and objectnew.acsmag lt 90,dimcount)

   backmem = where(sphdist(objectnew.ra, objectnew.dec, randra(cand), randdec(cand), /degrees) lt sep and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.acsmag gt 0 and objectnew.acsmag lt 90 and objectnew.acsmag - objectnew.irac1mag gt 2.0 and objectnew.acsmag - objectnew.irac1mag lt 3.5 and objectnew.irac1mag lt 21.4 , backmemcount)

   goodbrightarr(gc) = brightcount; n_elements(goodcolorbright)
   gooddimarr(gc) = dimcount; n_elements(goodcolordim)
   backmemarr(gc) =  backmemcount
   gc = gc + 1

 

endfor 


print, 'backcount', mean(backcount)

goodbrightarr = goodbrightarr(where(goodbrightarr gt 0))
gooddimarr = gooddimarr(where(gooddimarr gt 0))
backmemarr = backmemarr[0:gc-1]

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
canddec = [69.04481,  69.06851,  69.087766];,  68.98017,  69.019103 ]
candra = [264.68160, 264.89228, 264.83053];, , 264.66337, 265.27136 ]


;for cand=0, 2 do begin
for cand=0, n_elements(canddec) - 1 do begin
   print, 'working on cluster', candra(cand), canddec(cand)
   
   ;objects within sep of  the center and detected in irac1
   clusterbright = where(sphdist(objectnew.ra, objectnew.dec, candra(cand), canddec(cand), /degrees) lt sep and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.acsmag gt 0 and objectnew.acsmag lt 90 and objectnew.irac1mag lt 22.6 and objectnew.irac1mag gt 20.6, nclusterbright) 


   clusterfaint = where(sphdist(objectnew.ra, objectnew.dec, candra(cand), canddec(cand), /degrees) lt sep and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.acsmag gt 0 and objectnew.acsmag lt 90 and objectnew.irac1mag lt 24.6 and objectnew.irac1mag gt 22.6, nclusterfaint) 

   print,    'nclusterbright nclusterfaint', nclusterbright, nclusterfaint
;
endfor


;;----------------------------------------------------------------
;read in the cosmos  cluster data
;----------------------------------------------------------------

kcork=[.04,.04,.04,-.02]

sep = 0.0754

readcol, '/Users/jkrick/nep/clusters/cosmos/cluster.wcs.txt', targetra, targetdec, cosclusternum, z, filename, format="F,F,I,F,A"
print, 'n elements z', n_elements(z)
nclusters = n_elements(z)
for i = 0, n_elements(z) -1 do begin
   readlargecol, filename(i),      id42, image_tile42,       ra42,      dec42,    x_pix42,    y_pix42,   i_fwhm42,    i_max42,   i_star42,         i_mag_auto42,        auto_offset42,auto_flag42,    u_mag42,u_mag_err42,    b_mag42,b_mag_err42,    v_mag42,v_mag_err42,    g_mag42,g_mag_err42,    r_mag42,r_mag_err42,    i_mag42,i_mag_err42,    z_mag42,z_mag_err42,    k_mag42,k_mag_err42,         i_cfht_mag42,     i_cfht_mag_err42,         u_sdss_mag42,     u_sdss_mag_err42,         g_sdss_mag42,     g_sdss_mag_err42,         r_sdss_mag42,     r_sdss_mag_err42,         i_sdss_mag42,     i_sdss_mag_err42,         z_sdss_mag42,     z_sdss_mag_err42,f814w_mag42,      f814w_mag_err42,nb816_mag42,      nb816_mag_err42,      ebv42,    zphot42,        zerr_68_min42,        zerr_68_max42,        zerr_95_min42,        zerr_95_max42,    tphot42,    rphot42,    chisq42,   nbands42,       mv42,      d9542,  logmass42,     star42,   b_mask42,   v_mask42,   i_mask42,   z_mask42, blend_mask42, format="A"
   
   k = k_mag42 + kcork(i)
   cosmosbright= where( b_mag42 gt 0 and b_mag42 lt 90 and k_mag42 gt 0 and k_mag42 lt 90 $ ;and zphot42 gt 0 and zphot42 lt 3*z(i) $
                        and sphdist(ra42, dec42, targetra(i), targetdec(i), /degrees)  lt sep and k lt 18.2 and k gt 16.2, ncosmosbright) ;   
   cosmosfaint= where( b_mag42 gt 0 and b_mag42 lt 90 and k_mag42 gt 0 and k_mag42 lt 90 $ ;and zphot42 gt 0 and zphot42 lt 3*z(i) $
                        and sphdist(ra42, dec42, targetra(i), targetdec(i), /degrees)  lt sep and k lt 20.2 and k gt 18.2, ncosmosfaint) ;   
   
   print, 'ncosmosbright, ncosmosfaint', ncosmosbright, ncosmosfaint
endfor

end

