pro cluster_mips24


common sharevariables, candra, canddec, sep, objectnew, rvir, member, membernoagn, memberagn, clus_z1, clus_24,nfluxlimit,nflux, beta, luminositydistance, wholefield, wholefieldnoagn, wholefieldagn, clus_z1noagn


restore, '/Users/jkrick/idlbin/objectnew.sav'

!P.multi=[0,1,1]
ps_open, filename='/Users/jkrick/nep/clusters/mips24/sfr.ps',/portrait,/square,/color, xsize=5, ysize=5
!P.charthick=3
!P.thick=3
!X.thick=3
!Y.thick=3
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
blackcolor = FSC_COLOR("Black", !D.Table_Size-4)


;calculate distances in Mpc for reference
z = 1.0    ; maybe one day calculate this more exactly
luminositydistance = lumdist(z)
angulardistance = luminositydistance / ((1+z)^2)    

;what is the conversion between arcminutes and Mpc?
mpcperarcmin = (tan((1./60.)*(!Pi/180.))) * angulardistance 
;-------------------------------------------------------------------------------
;where are the 24 micron sources in the cluster areas
;-------------------------------------------------------------------------------

canddec = [69.045052,  69.06851,  69.087766]
candra = [264.68365, 264.89228, 264.83053]
sep = 0.019

rvir = [.02425,.020,.020]  ;~0.6Mpc
;want rvirial or 2*rvirial
rvir = 2* rvir

;first find a sample of those objects with phot-z's in the
;right range which are within rvir or 2*rvir.
;junk = cluster_mips24_sample()
junk = cluster_sample()

radofsample =0.053 ;radius in degrees of circle with equivalent size of the sample area
 
;-------------------------------------------------------------------------------
;plot cumulative distribution of distance from center for all members
;with mips detections (member), all galaxies in general, and then all cluster
;members regardless of mips detections (clus_z1)
;-------------------------------------------------------------------------------
;measure distance of each member object with mips detection from it's cluster's center.

memberdist0 = sphdist(objectnew[member].ra, objectnew[member].dec, candra(0), canddec(0), /degrees) 
memberdist1 = sphdist(objectnew[member].ra, objectnew[member].dec, candra(1), canddec(1), /degrees) 
memberdist2 = sphdist(objectnew[member].ra, objectnew[member].dec, candra(2), canddec(2), /degrees) 
memberdist = fltarr(n_elements(member))

;in the case of overlaps, keep the smallest distance.
for memnum = 0, n_elements(member) - 1 do begin
   if memberdist0(memnum) lt memberdist1(memnum) and memberdist0(memnum) lt memberdist2(memnum) $
   then memberdist(memnum) = memberdist0(memnum)
   if memberdist1(memnum) lt memberdist0(memnum) and memberdist1(memnum) lt memberdist2(memnum) $
   then memberdist(memnum) = memberdist1(memnum)
   if memberdist2(memnum) lt memberdist1(memnum) and memberdist2(memnum) lt memberdist0(memnum) $
   then memberdist(memnum) = memberdist2(memnum)
endfor

;second test, in the case of overlaps, ignore cluster 3, but keep
;smallest distance, either cluster 1 or 2
memberdist_noclus3 = fltarr(n_elements(member))
for memnum = 0, n_elements(member) - 1 do begin
   if memberdist0(memnum) lt memberdist1(memnum) then memberdist_noclus3(memnum) = memberdist0(memnum)
   if memberdist1(memnum) lt memberdist0(memnum) then memberdist_noclus3(memnum) = memberdist1(memnum)
endfor

;need to sort on distance
memberdist = memberdist*60.  ;now in arcmin
sortmemberdist= memberdist(sort(memberdist))
N1 = n_elements(sortmemberdist)
f1 = (findgen(N1) + 1.)/ N1
plot, sortmemberdist, f1, xtitle='Distance from Cluster Center (arcmin)', ytitle='Cumulative Fraction', xrange=[0,3], xstyle=9;, psym = 10
axis, 0, 1.0, xaxis=1, xrange=[0,3*mpcperarcmin], xstyle = 1;, xthick = 3,charthick = 3

;-------------------------------------------------------------------------------
;measure distance of each z=1 object from it's cluster's center.

clus_z1dist0 = sphdist(objectnew[clus_z1].ra, objectnew[clus_z1].dec, candra(0), canddec(0), /degrees) 
clus_z1dist1 = sphdist(objectnew[clus_z1].ra, objectnew[clus_z1].dec, candra(1), canddec(1), /degrees) 
clus_z1dist2 = sphdist(objectnew[clus_z1].ra, objectnew[clus_z1].dec, candra(2), canddec(2), /degrees) 
clus_z1dist = fltarr(n_elements(clus_z1))

;in the case of overlaps, keep the smallest distance.
for memnum = 0, n_elements(clus_z1) - 1 do begin
   if clus_z1dist0(memnum) lt clus_z1dist1(memnum) and clus_z1dist0(memnum) lt clus_z1dist2(memnum) $
   then clus_z1dist(memnum) = clus_z1dist0(memnum)
   if clus_z1dist1(memnum) lt clus_z1dist0(memnum) and clus_z1dist1(memnum) lt clus_z1dist2(memnum) $
   then clus_z1dist(memnum) = clus_z1dist1(memnum)
   if clus_z1dist2(memnum) lt clus_z1dist1(memnum) and clus_z1dist2(memnum) lt clus_z1dist0(memnum) $
   then clus_z1dist(memnum) = clus_z1dist2(memnum)
endfor

;second test, in the case of overlaps, ignore cluster 3, but keep
;smallest distance, either cluster 1 or 2
clus_z1dist_noclus3 = fltarr(n_elements(clus_z1))
for memnum = 0, n_elements(clus_z1) - 1 do begin
   if clus_z1dist0(memnum) lt clus_z1dist1(memnum) then clus_z1dist_noclus3(memnum) = clus_z1dist0(memnum)
   if clus_z1dist1(memnum) lt clus_z1dist0(memnum) then clus_z1dist_noclus3(memnum) = clus_z1dist1(memnum)
endfor

;need to sort on distance
clus_z1dist = clus_z1dist*60.   ;now in Mpc
sortclus_z1dist= clus_z1dist(sort(clus_z1dist))
N1 = n_elements(sortclus_z1dist)
f1 = (findgen(N1) + 1.)/ N1
oplot, sortclus_z1dist, f1, psym = 10, linestyle = 1


legend, ['member w/ mips24', 'all members', 'random z=1 mips24'], linestyle = [0, 1, 3]
;-------------------------------------------------------------------------------
;get together some random background regions to look at their spatial distribution
;-------------------------------------------------------------------------------
nrand = 30
seed = 45
randx = randomu(seed, nrand)
randy = randomu(seed, nrand)

randx = randx * 9500. + 6800.
randy = randy * 13300. + 3100.

acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits') ;
xyad, acshead, randx, randy, randra, randdec
allranddist = fltarr(nrand*1000)
d = 0
for cand=0, n_elements(randdec) - 1 do begin ;for each random area
   ;don't include stars
   goodrand = where(sphdist(objectnew.ra, objectnew.dec, randra(cand), randdec(cand), /degrees) lt rvir(0) $
                    and objectnew.zphot gt 0.87 and objectnew.zphot lt 1.13  and objectnew.chisq lt 50. $
                    and objectnew.mips24flux gt 17.3 and objectnew.mips24mag lt 90  and nflux gt nfluxlimit)

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
oplot, allranddist, f2, psym = 10, linestyle=3

;determine the ks-test distance and probability for each population
kstwo, allranddist, sortmemberdist, ksdist, prob
print, '-------------------------------------------------------------------------------'
print, 'ksdist for member with 24, prob', ksdist, prob

kstwo, allranddist, sortclus_z1dist, ksdistgals, prob
print, 'ksdist for members, prob', ksdistgals, prob

kstwo, sortmemberdist, sortclus_z1dist, ksdistgals, prob
print, 'ksdist between members w/ and w/o 24 , prob', ksdistgals, prob


;-------------------------------------------------------------------------------
;Is there an overdensity of Mips detections in the cluster compared to
;the background?
;-------------------------------------------------------------------------------
;first, is the spatial density of mips sources over all redshifts higher in the
;clusters than mips sources in the background?
randmips = fltarr(n_elements(randdec))
randmipsz1 = fltarr(n_elements(randdec))

; how many mips detections are within radofsample of randra, randdec?
for rand=0, n_elements(randdec) - 1 do begin ;for each random area
   randmips(rand) = n_elements(where(sphdist(objectnew.ra, objectnew.dec, randra(rand), randdec(rand), /degrees) lt radofsample $
                    and objectnew.mips24flux gt 17.3 and objectnew.mips24mag lt 90  and nflux gt nfluxlimit))


;second, at redshift one, is the spatial density of mips sources in
;the clusters higher than mips sources at redshift one in the background?

   randmipsz1(rand) = n_elements(where(sphdist(objectnew.ra, objectnew.dec, randra(rand), randdec(rand), /degrees) lt radofsample $
                      and objectnew.zphot gt 0.87 and objectnew.zphot lt 1.13  and objectnew.chisq lt 50. $
                      and objectnew.mips24flux gt 17.3 and objectnew.mips24mag lt 90  and nflux gt nfluxlimit))
endfor

print, '-------------------------------------------------------------------------------'
print, 'mean mips sources in background area', mean(randmips), stddev(randmips)
print, 'mean mips sources at z=1 in background area', mean(randmipsz1), stddev(randmipsz1)
print, '-------------------------------------------------------------------------------'



;-------------------------------------------------------------------------------
;what is the luminosity of those sources
;-------------------------------------------------------------------------------


;use ranga's code to match SED's to the 24micron flux given the photz
;of the galaxy.
get_lir, objectnew[member].zphot, objectnew[member].mips24flux, objectnew[member].mips24fluxerr, lir, lirerr
;print, 'total, mean lir, mean lirerr', total(lir), mean(lir) , mean(lirerr)
plothist, lir, bin = 5E9, /xlog, xtitle = 'Lir in solar luminosities', ytitle = 'Number', xrange=[1E9, 5E12], xstyle=9
oplot, findgen(25) - findgen(25) + 1E11, findgen(25), linestyle = 3
oplot, findgen(25) - findgen(25) + 1E12, findgen(25), linestyle = 3

;how many LIRGs, and ULIRGs?
lirg = where(lir ge 1E11 and lir lt 1E12)
print, 'LIRGS #, %', n_elements(lirg), float(n_elements(where(lir ge 1E11 and lir lt 1E12))) / float(n_elements(member))
print, 'ULIRGS #, %', n_elements(where(lir ge 1E12 )), float(n_elements(where(lir ge 1E12 ))) / float(n_elements(member))

;now convert lir into sfr.
lirtosfr=1.71E-10
sfr = lir * 1.71E-10
print, 'total, mean sfr', total(sfr), mean(sfr)
axis, 1E9,8,  xaxis=1, xrange=[1E9*lirtosfr,5E12*lirtosfr], xstyle = 1

;plothist, sfr, xtitle = 'SFR in solar masses/ year', ytitle = 'Number'

;----------------
;compare with all mips sources at z=1 across the field
get_lir, objectnew[wholefield].zphot, objectnew[wholefield].mips24flux, objectnew[wholefield].mips24fluxerr, lir_field, lirerr_field
;normalize to the total number of the member sample
a = float(n_elements(member))
b =  float(n_elements(wholefield))
normfactor = a / b

plothist,  lir_field , xhist_field, yhist_field, bin = 5E9, /xlog, /noprint, /noplot;, xtitle = 'Lir in solar luminosities', ytitle = 'Number', xrange=[1E9, 5E12], xstyle=9
;oplot, xhist_field, yhist_field, psym = 10, linestyle = 1
plot, xhist_field , yhist_field*normfactor,  xtitle = 'Lir in solar luminosities', ytitle = 'Number', xrange=[1E9, 5E12], xstyle=9, psym = 10,/xlog
oplot, findgen(25) - findgen(25) + 1E11, findgen(25), linestyle = 3
oplot, findgen(25) - findgen(25) + 1E12, findgen(25), linestyle = 3
axis, 1E9,8,  xaxis=1, xrange=[1E9*lirtosfr,5E12*lirtosfr], xstyle = 1


;-------------------------
;which are the ulirgs?
ulirg = where(lir ge 1E12)
;print, 'ulirg', member[ulirg], memberdist[ulirg]*mpcperarcmin, objectnew[member[ulirg]].ra, objectnew[member[ulirg]].dec



;-------------------------------------------------------------------------------
;what is the morphology of those sources
;-------------------------------------------------------------------------------

;start with just what type of galaxy was the best fit in hyperz 

;spectral type key
sptkey = strarr(15)
;sptkey=["  "," Arp220 " ," Ell13  " ," Ell2  " ," M82  " ," Mrk231" ," QSO1 " ," QSO2 " ," S0   " ," Sa  " ," Sb  " ," Sc  " ," Sd  " ," Sey1.8", " Sey2  " ," Torus" ]     
sptkey=[" Arp220 " ," Ell13  " ," Ell2  " ," M82  " ," Mrk231" ," QSO1 " ," QSO2 " ," S0   " ," Sa  " ," Sb  " ," Sc  " ," Sd  " ," Sey1.8", " Sey2  "  ]     

nothing = replicate(' ', 17)
plothist, objectnew[member].spt, xhist, yhist, /noprint, /noplot; xrange=[0,16], xstyle = 1, xtickname = nothing
plothist, objectnew[wholefield].spt, xhist_field, yhist_field, /noprint, /noplot; xrange=[0,16], xstyle = 1, xtickname = nothing
yhist_field = yhist_field*normfactor
;want to sort on the most popular kind of object
;sptkey = sptkey[sort(yhist)]
;yhist = yhist[sort(yhist)]

;try sorting on type of objects.
index = [1,2,7,8,9,10,11,3,0,4,5,6,12,13]
sptkey = sptkey(index)
yhist = yhist(index)
yhist_field = yhist_field(index)

plot, xhist, yhist, xrange=[0,17], xstyle = 1, xtickname = nothing, psym = 10, ytitle="Number"
oplot, xhist_field, yhist_field, psym = 10, linestyle = 2

xlabel = findgen(16) 
ylabel = findgen(16) - findgen(16) -0.1
;sptkey=[" Torus", sptkey]
sptkey =[" ", sptkey, " Torus"]

xyouts, xlabel, ylabel, sptkey, orientation=310

;make a vertical line at the end of the histogram
;psym = 10 does not fully get the histogram thing correct
x = [1,1]
y  = [11,0]
oplot, x, y
x = [1,2]
y  = [11,11]
oplot, x, y
x=[15,15]
y = [13,0]
oplot, x, y
x=[14,15]
y = [13,13]
oplot, x, y

x = [15,15]
y = [18,0]
oplot, x, y, linestyle = 2
x = [14,15]
y = [18,18]
oplot, x, y, linestyle = 2
x = [1,2]
y = [7,7]
oplot, x, y, linestyle = 2
x = [1,1]
y = [7,0]
oplot, x, y, linestyle = 2

;-------------------------------------------------------------------------------
;make thumbnail image of the morphology by eye stuff
;-------------------------------------------------------------------------------

;interest = [1196, 1250, 1036, 1026, 2536]
;acsthumbnails, interest, '/Users/jkrick/nep/clusters/mips24/morphthumb.ps'

;-------------------------------------------------------------------------------
;make irac color color plot to see where all of these objects are.
;-------------------------------------------------------------------------------
vsym, /polygon, /fill

plot, alog10(objectnew[clus_z1].irac3flux / objectnew[clus_z1].irac1flux), alog10(objectnew[clus_z1].irac4flux / objectnew[clus_z1].irac2flux), $
      psym = 8, xrange = [-1.5, 1.], yrange=[-1.0, 1.0],/xst,/yst, xtitle='log(S5.8/S3.6)', ytitle = 'log(S8.0/S4.5)', symsize=0.5 

oplot, [-0.2,1.3], [-0.2,-0.2]
oplot, [-0.2,-0.2], [-0.2,0.3]
oplot, [-0.2,1.3], [0.3,1.5]

oplot,alog10(objectnew[member].irac3flux / objectnew[member].irac1flux), alog10(objectnew[member].irac4flux / objectnew[member].irac2flux), $
      psym = 8,symsize=0.5, color=redcolor

; now pick out the ones that hyperz thinks are AGN
agn = where(objectnew[member].spt eq 6 or objectnew[member].spt eq 7 or objectnew[member].spt eq 13$
or objectnew[member].spt eq 14 or objectnew[member].spt eq 15)


oplot,  alog10(objectnew[member[agn]].irac3flux / objectnew[member[agn]].irac1flux), $
        alog10(objectnew[member[agn]].irac4flux / objectnew[member[agn]].irac2flux),  psym = 8,symsize=0.5, color=bluecolor
print, '-------------'
print, 'number of agn within member based on hyperz', n_elements(AGN)

legend, ['all redshift 1', 'all redshift 1 with 24', 'all redshift 1 with 24 with AGN hyperz fit'], psym = [8,8,8], colors=[blackcolor, redcolor, bluecolor]

;-------------------------------------------------------------------------------
;recalculate LIR without the ones measured as AGN

;assume AGNness from hyperz
notagn = where(objectnew[member].spt eq 1 or objectnew[member].spt eq 2 or objectnew[member].spt eq 3 or objectnew[member].spt eq 4 or objectnew[member].spt eq 5 or objectnew[member].spt eq 8 or objectnew[member].spt eq 9 or objectnew[member].spt eq 10 or objectnew[member].spt eq 11 or objectnew[member].spt eq 12)

get_lir, objectnew[member[notagn]].zphot, objectnew[member[notagn]].mips24flux, objectnew[member[notagn]].mips24fluxerr, lirnoagn, lirnoagnerr
;plothist, lirnoagn, bin = 5E9, /xlog, xtitle = 'Lirnoagn in solar luminosities', ytitle = 'Number', xrange=[1E9, 5E12], xstyle=9, yrange=[0,8]
;oplot, findgen(25) - findgen(25) + 1E11, findgen(25), linestyle = 3
;oplot, findgen(25) - findgen(25) + 1E12, findgen(25), linestyle = 3

;how many LIRNOAGNGs, and ULIRNOAGNGs?
print, 'LIRNOAGNGS #, %', n_elements(where(lirnoagn ge 1E11 and lirnoagn lt 1E12)), float(n_elements(where(lirnoagn ge 1E11 and lirnoagn lt 1E12))) / float(n_elements(member[notagn]))
print, 'ULIRNOAGNGS #, %', n_elements(where(lirnoagn ge 1E12 )), float(n_elements(where(lirnoagn ge 1E12 ))) / float(n_elements(member[notagn]))

;now convert lirnoagn into sfr.
lirnoagntosfr=1.71E-10
sfr = lirnoagn * 1.71E-10
print, 'total, mean sfr', total(sfr), mean(sfr)
axis, 1E9,8,  xaxis=1, xrange=[1E9*lirnoagntosfr,5E12*lirnoagntosfr], xstyle = 1


;assume AGNness from Lacy wedge
;something not working right here..******
notagn1 = where(alog10(objectnew[member].irac3flux / objectnew[member].irac1flux) lt -0.2 )


notagn2 = where(alog10(objectnew[member].irac3flux / objectnew[member].irac1flux) ge -0.2 and alog10(objectnew[member].irac4flux / objectnew[member].irac2flux) lt -0.2)

notagn=[notagn1,notagn2]
print, '-------------'
print, 'number of AGN within member based on hyperz', n_elements(member) - n_elements(notagn)

get_lir, objectnew[member[notagn]].zphot, objectnew[member[notagn]].mips24flux, objectnew[member[notagn]].mips24fluxerr, lirnoagn, lirnoagnerr
;plothist, lirnoagn, bin = 5E9, /xlog, xtitle = 'Lirnoagn in solar luminosities', ytitle = 'Number', xrange=[1E9, 5E12], xstyle=9, yrange=[0,8]
;oplot, findgen(25) - findgen(25) + 1E11, findgen(25), linestyle = 3
;oplot, findgen(25) - findgen(25) + 1E12, findgen(25), linestyle = 3

;how many LIRNOAGNGs, and ULIRNOAGNGs?
print, 'LIRNOAGNGS #, %', n_elements(where(lirnoagn ge 1E11 and lirnoagn lt 1E12)), float(n_elements(where(lirnoagn ge 1E11 and lirnoagn lt 1E12))) / float(n_elements(member[notagn]))
print, 'ULIRNOAGNGS #, %', n_elements(where(lirnoagn ge 1E12 )), float(n_elements(where(lirnoagn ge 1E12 ))) / float(n_elements(member[notagn]))

;now convert lirnoagn into sfr.
lirnoagntosfr=1.71E-10
sfr = lirnoagn * 1.71E-10
print, 'total, mean sfr', total(sfr), mean(sfr)
axis, 1E9,8,  xaxis=1, xrange=[1E9*lirnoagntosfr,5E12*lirnoagntosfr], xstyle = 1

;------------------------------------------------------------------------------
;compare our measurements to Bai et al. 2007
junk = mass_norm_sfr()

;------------------------------------------------------------------------------
;what are the UV derived SFR's for the member galaxies with
;mips detections?

;junk = uv_sfr()
;------------------------------------------------------------------------------
;put the member galaxies with 24 detections on the CMD's from paper I

;same plot as paper 1
;arg ; unfortunately need to use the old catalog to make the plots match...
restore, '/Users/jkrick/idlbin/object.old.sav'

iall = fltarr(210000)
iracall = fltarr(210000)
j = 0
sep_paper1 =0.017
for cand=0, n_elements(canddec) - 1 do begin
   good = where(sphdist(objectnew.ra, objectnew.dec, candra(cand), canddec(cand), /degrees) lt sep_paper1 and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.acsmag gt 0 and objectnew.acsmag lt 90 );objectnew.nearestdist gt 0.00055) ;and objectnew.nearestdist*60*60/0.6 gt 4.0 ) ;

   for k = 0l, n_elements(good) - 1 do begin
      iall[j] = objectnew[good[k]].acsmag
      iracall[j] =objectnew[good[k]].irac1mag
      j = j + 1
   endfor
endfor
iall = iall[0:j-1]
iracall = iracall[0:j-1]

;print,'n_elem,ents(iall)', n_elements(iall)


dm = 44.1
dm = dm - 0.75  ;for filter bandpass shrinking
dm = dm - 0.7  ;evolution correction
kcorr = 0.
m = -0.12
ball = 5.2
rcswidth = 0.5

im_hessplot, iracall, iall -iracall, yrange=[0,6], xtitle = '[3.6] (AB)', ytitle = 'ACS F814W - [3.6]', xrange=[18,25],/notsquare, nbin2d=15, position = [0.1 ,0.1 , 0.9,0.9], charsize=1, charthick=3, xthick=3, ythick=3, thick=3,xstyle=9, ystyle = 1
axis, 0,6, xaxis=1,xrange=[18-dm+kcorr,25-dm+kcorr], xstyle=1, xthick=3, charthick=3

oplot, findgen(100),  m*findgen(100) + ball, thick = 3 ;, color = orangecolor
oplot, findgen(100), m*findgen(100) + ball + rcswidth, thick = 3 ;, color = orangecolor
oplot, findgen(100),  m*findgen(100) +ball - rcswidth, thick = 3 ;, color = orangecolor
oplot, [25.1, 22], [2.9,6.0], linestyle = 2, thick=3

;now overplot member galaxies
;switch back to new catalog
restore, '/Users/jkrick/idlbin/objectnew.sav'

oplot, objectnew[member].irac1mag, objectnew[member].acsmag - objectnew[member].irac1mag, psym = 2

;now indicate which ones are AGN based on their SED fits
oplot, objectnew[member[agn]].irac1mag, objectnew[member[agn]].acsmag - objectnew[member[agn]].irac1mag, psym = 6


redmember = where(objectnew[member[notagn]].acsmag - objectnew[member[notagn]].irac1mag gt 2.5)
;plothyperz, member[notagn[redmember]], '/Users/jkrick/nep/clusters/mips24/redmips24.ps'

;maybe instead or inaddition should plot color distribution.
;otherwise it is noticeable that there are galaxies in the member
;sample in places where there are none in the CMD.  damn.  also a few
;junky things up by the detection limit.


plothist, objectnew[member].acsmag - objectnew[member].irac1mag, xtitle =  'ACS F814W - [3.6]', xrange=[1,5], xstyle = 1 $
          , ytitle = 'Number', bin = 0.2


;how many of the mips sources are on the red sequence, roughly
red = where(objectnew[member].acsmag - objectnew[member].irac1mag gt 2.5 )
print, 'number on red sequence', n_elements(red)
print, 'number of nonagn red', n_elements(redmember) 

;what does the cmd of all mips sources at photz=1 look like?

plot, objectnew[wholefield].irac1mag, objectnew[wholefield].acsmag - objectnew[wholefield].irac1mag, psym = 2, yrange=[0,6], xtitle = '[3.6] (AB)', ytitle = 'ACS F814W - [3.6]', xrange=[18,25], xstyle = 1

oplot, findgen(100),  m*findgen(100) + ball, thick = 3 ;, color = orangecolor
oplot, findgen(100), m*findgen(100) + ball + rcswidth, thick = 3 ;, color = orangecolor
oplot, findgen(100),  m*findgen(100) +ball - rcswidth, thick = 3 ;, color = orangecolor
oplot, [25.1, 22], [2.9,6.0], linestyle = 2, thick=3


;------------------------------------------------------------------------------


ps_close, /noprint,/noid


end


;-------------------------------------------------------------------------------------------------


;plothist, objectnewnew[good].rmag, xrange=[10,30]

;plothyperz, good, strcompress( '/Users/jkrick/nep/clusters/mips24/mips24_'+string(i)+'.ps')


;openw, outlun, strcompress('/Users/jkrick/nep/clusters/mips24/cluster'+string(i)+'.reg'), /get_lun
;printf, outlun, 'fk5'
;for rc=0, n_elements(good) -1 do  printf, outlun, 'circle( ', objectnewnew[good[rc]].ra, objectnewnew[good[rc]].dec, ' 3")'
;close, outlun
;free_lun, outlun


;openw, outlun2, strcompress('/Users/jkrick/nep/clusters/mips24/cluster'+string(i)+'.txt'), /get_lun
;printf, outlun, 'number, ra, dec, zphot, rmag, acsmag, irac1mag, mips24flux, acsmatchdist, comment'
;for rc=0, n_elements(good) -1 do  printf, outlun2,  format='(I10,F10.5,F10.5,F10.2,F10.2,F10.2,F10.2,F10.2)', good[rc], objectnewnew[good[rc]].ra, objectnewnew[good[rc]].dec, objectnewnew[good[rc]].zphot, objectnewnew[good[rc]].rmag,objectnewnew[good[rc]].acsmag, objectnewnew[good[rc]].irac1mag,objectnewnew[good[rc]].mips24flux
;close, outlun2
;free_lun, outlun2


;-------------------------------------------------------------------------------------------------
;now bin distances
;5 bins, < .0038, .0038 < .0076, .0076<.0114, .0114 < .0152,
;.0152<.019

;x = [.0019,.0057,.0095,.0133,.0171]
;x = x * 60.
;y = [n_elements(where(sortdistallcluster lt .0038* 60.)),n_elements(where(sortdistallcluster ge .0038* 60. and sortdistallcluster lt .0076* 60.)),n_elements(where(sortdistallcluster ge .0076* 60. and sortdistallcluster lt .0114* 60.)),n_elements(where(sortdistallcluster ge .0114* 60. and sortdistallcluster lt .0152* 60.)),n_elements(where(sortdistallcluster ge .152* 60. and sortdistallcluster lt .019* 60.))]

;y2 = [n_elements(where(allranddist lt .0038* 60.)),n_elements(where(allranddist ge .0038* 60. and allranddist lt .0076* 60.)),n_elements(where(allranddist ge .0076* 60. and allranddist lt .0114* 60.)),n_elements(where(allranddist ge .0114* 60. and allranddist lt .0152* 60.)),n_elements(where(allranddist ge .152* 60. and allranddist lt .019* 60.))]

;plot, x, y, psym = 2;, psym = 10
;oplot, x, y2, psym = 1; linestyle = 1;,psym=10

;plothist, sortdistallcluster, bin=0.019*60./10., xtitle='distance from cluster center in arcminutes', ytitle='N', xrange=[0,0.02*60.], xstyle=1
;plothist, allranddist, bin=0.019*60./10., /overplot, linestyle=1
;print, n_elements(dist24), n_elements(allranddist) ,'nelements dist, rand'

