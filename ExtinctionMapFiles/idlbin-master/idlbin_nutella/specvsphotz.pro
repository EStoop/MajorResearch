pro specvsphot

restore, '/Users/jkrick/idlbin/objectnew.sav'
!P.multi=[0,1,1]
!P.charthick=3
!P.thick=3
!X.thick=3
!Y.thick=3

ps_open, filename='/Users/jkrick/palomar/cosmic/specvsphotz_smaller.ps',/portrait,/square,/color, xsize=5, ysize=5

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)

;read in the excell spreadsheet of spec z's and photz's from the
;catalog.  which now should be the same as any new run of
;photz's since I have just rerun the photz's in the catalog.  08/22.

readcol, '/Users/jkrick/palomar/cosmic/01june2008/spectra/specz.txt', masknum, catnum, photz, xccz, xcczerr,	zxc,	xcr,	xcgrade,	emcz,	emczerr,	zem,	emgrade,	whichmethod, format="(I1,L5,F10.3, F10.2,F10.2, F10.3, F10.2, I1, F10.2,F10.2,F10.3,I1,A)"

catarr = fltarr(n_elements(catnum))
photzarr = fltarr(n_elements(catnum))
speczarr = fltarr(n_elements(catnum))
sptarr = fltarr(n_elements(catnum))
chisqarr = fltarr(n_elements(catnum))
count = 0

for i = 0, n_elements(catnum) - 1 do begin
   if whichmethod(i) eq 'em' and emgrade(i) gt 1 then begin
      catarr(count) = catnum(i)
      photzarr(count) = photz(i)
      speczarr(count) = zem(i)
      count = count + 1
   endif
   if whichmethod(i) eq 'xc' and xcgrade(i) gt 1 then begin
      catarr(count) = catnum(i)
      photzarr(count) = photz(i)
      speczarr(count) = zxc(i)
      count = count + 1
   endif

   
endfor

catarr = catarr[0:count - 1]
photzarr = photzarr[0:count - 1]
speczarr = speczarr[0:count - 1]
sptarr = sptarr[0:count - 1]
chisqarr = chisqarr[0:count - 1]


chisqarr = objectnew[catarr].chisq
sptarr = objectnew[catarr].spt
;but wait if the catalog does get updated in terms of photz then I
;want to have the most recent photz's.
photzarr = objectnew[catarr].zphot


;print, 'test', catarr(12), photzarr(12), chisqarr(12), sptarr(12), objectnew[catarr(12)].zphot, objectnew[catarr(12)].chisq, objectnew[catarr(12)].spt

star = where(xcgrade lt 1)
print, 'n stars', n_elements(star)



;plot, findgen(10), findgen(10),  xrange=[0,1], yrange=[-3,3], xtitle='spec z', ytitle='(photz - specz) / (1 + specz)',/nodata, $
;      charthick = 3, xthick = 3, ythick = 3
;oplot, findgen(10), findgen(10) - findgen(10)
;oplot, speczarr, (photzarr - speczarr) / (1 + speczarr),psym=2,  thick = 3 ;color=bluecolor,

dz1plusz =(photzarr - speczarr) / (1 + speczarr)

sigmaorig = stddev(dz1plusz)
print, 'sigma original', sigmaorig
print, 'n_elements(speczarr)', n_elements(speczarr)
;ps_close, /noprint,/noid

;-----------------------------------------------------------------
;rerun hyperz through running plothyperz.
;plothyperz, catarr, '/Users/jkrick/palomar/cosmic/01june2008/spectra/testhyperz.ps'

;then read in the resulting photz's
;readcol,'/Users/jkrick/ZPHOT/hyperz_swire_target.z_phot',idz, zphota, chia, proba, specta, nagea, agea, $
;        ava, ba, zinf99a,zsup99a,zinf90a,zsup90a,zinf68a,zsup68a,zwma,probwma,Mabsa,$
;        zphot2a,prob2a,format="A"


;for i = 0, n_elements(speczarr) - 1 do begin
;   print, 'old, new', catarr(i), photzarr(i), zphota(i)
;endfor
;-----------------------------------------------------------------
;
;get rid of objects which are agn based on photz template
;agngood = where(specta ne 4 and specta ne 5 and specta ne 6 and specta ne 7 and specta ne 13 and specta ne 14 and specta ne 15)
;agngood = where(sptarr ne 4 and sptarr ne 5 and sptarr ne 6 and sptarr ne 7 and sptarr ne 13 and sptarr ne 14 and sptarr ne 15)
;print, 'n_elements emgood after agn removed', n_elements(agngood)


;instead of agn, try getting rid of things with very high chisq.
;agngood= where(chisqarr lt 50)

;need to disinclude those with close close pairs, < few arcseconds.
;563,1739,2477,2017,2503,1917,1921
good = where(chisqarr lt 50 and catarr ne 563 and catarr ne 1739 and catarr ne 2477 and catarr ne 2017 and catarr ne 2503 and catarr ne 1917 and catarr ne 1921)
;get rid of objects which are agn based on irac colors care of Mark
;Lacy
; have already shown with iraccolorcolor.pro that actually the AGN
;determination within hyperz is very similar to that determined by
;Mark, and what we really want to exclude is those AGN which dominate
;the spectra, if they are hidden, it doesn't matter too much for photoz's
;agngood = where(alog10(objectnew[catarr].ch4flux / objectnew[catarr].ch2flux) lt -0.2 and alog10(objectnew[catarr].ch3flux / objectnew[catarr].ch1flux) lt -0.2))


;replot the sigma plot
plot, findgen(10), findgen(10),  xrange=[0,1], yrange=[-0.5,0.5], xtitle='spec z', $
      ytitle='(photz - specz) / (1 + specz)',/nodata, $
      charthick = 3, xthick = 3, ythick = 3, ystyle = 1
oplot, findgen(10), findgen(10) - findgen(10)
;oplot, speczarr, (zphota(agngood)- speczarr(agngood)) / (1 + speczarr(agngood)),psym=2, color=bluecolor, thick = 3

;dz1plusz = (zphota(agngood)- speczarr(agngood)) / (1 + speczarr(agngood))
dz1plusz =(photzarr(good) - speczarr(good)) / (1 + speczarr(good))

oplot, speczarr(good), dz1plusz,psym=2, thick = 3 ;, color=bluecolor



for count = 0, n_elements(dz1plusz) - 1 do begin
  if dz1plusz(count) ge 0.2 then print, dz1plusz(count), chisqarr(good(count)), "what",  catarr(good(count))
endfor


sigmaafter = stddev(dz1plusz)
print, 'sigmaafter', sigmaafter
print, 'n_elements(good)', n_elements(good)
ps_close, /noprint,/noid

;problem = where((zphota(agngood) -speczarr(agngood)) / (1+speczarr(agngood)) gt 1.0 )
;print, 'problem', idz(agngood(problem))

end

;plot, findgen(10), findgen(10),  xrange=[0,1], yrange=[0,1], xtitle='spec z', ytitle='phot z',/nodata, $
;      charthick = 3, xthick = 3, ythick = 3
;oplot, findgen(10), findgen(10)

;oplot,  zem(emgood), zphota[0:n_elements(emgood) - 1],psym=2, color=bluecolor, thick = 3
;oplot,  zxc(xcgood), zphota[n_elements(emgood):*],psym=5, color=bluecolor, thick = 3

;separate out the galaxies which have specz's that are at least
;ok, but are likely not to be AGN based on the photz fitting.
;embest = where(whichmethod eq 'em' and emgrade gt 2); and objectnew[catnum].spt ne 3 and objectnew[catnum].spt ne 4 and objectnew[catnum].spt ne 5 and objectnew[catnum].spt ne 6 and objectnew[catnum].spt ne 7 and objectnew[catnum].spt ne 13 and objectnew[catnum].spt ne 14 );and objectnew[catnum].prob gt 0)
;print, 'n_elements(embest)', n_elements(embest)
;xcbest = where(whichmethod eq 'xc' and xcgrade gt 2 );and objectnew[catnum].spt ne 3 and objectnew[catnum].spt ne 4 and objectnew[catnum].spt ne 5 and objectnew[catnum].spt ne 6 and objectnew[catnum].spt ne 7 and objectnew[catnum].spt ne 13 and objectnew[catnum].spt ne 14 );and objectnew[catnum].prob gt 0)
;print, 'n_elements(xcbest)', n_elements(xcbest)
;emgood = where(whichmethod eq 'em' and emgrade gt 1 );and objectnew[catnum].spt ne 3 and objectnew[catnum].spt ne 4 and objectnew[catnum].spt ne 5 and objectnew[catnum].spt ne 6 and objectnew[catnum].spt ne 7 and objectnew[catnum].spt ne 13 and objectnew[catnum].spt ne 14 );and objectnew[catnum].prob gt 0)
;print, 'n_elements(emgood)', n_elements(emgood)
;xcgood = where(whichmethod eq 'xc' and xcgrade gt 1 );and objectnew[catnum].spt ne 3 and objectnew[catnum].spt ne 4 and objectnew[catnum].spt ne 5 and objectnew[catnum].spt ne 6 and objectnew[catnum].spt ne 7 and objectnew[catnum].spt ne 13 and objectnew[catnum].spt ne 14 );and objectnew[catnum].prob gt 0)
;print, 'n_elements(xcgood)', n_elements(xcgood)
;good = [emgood, xcgood]
;speczgood = [zem(emgood), zxc(xcgood)]

;print, 'speczgood', speczgood

;oplot,  zem(emgood), (photz(emgood) - zem(emgood)) / (1+zem(emgood)),psym=2, color=bluecolor, thick = 3
;oplot,  zxc(xcgood), (photz(xcgood) - zxc(xcgood)) / (1+zxc(xcgood)),psym=5, color=bluecolor, thick = 3

;oplot,  zem(embest), (photz(embest) - zem(embest)) / (1+zem(embest)),psym = 2, color = redcolor, thick = 3
;oplot,  zxc(xcbest),(photz(xcbest) - zxc(xcbest)) / (1+zxc(xcbest)), psym=5, color = redcolor, thick=3

;oplot,  zem(emgood), (zphota[0:n_elements(emgood) - 1] - zem(emgood)) / (1 +  zem(emgood)),psym=2, color=bluecolor, thick = 3
;oplot,  zxc(xcgood), (zphota[n_elements(emgood):*] - zxc(xcgood)) / (1 + zxc(xcgood)),psym=5, color=bluecolor, thick = 3

