pro check_mips24


close, /all

ps_open, filename='/Users/jkrick/spitzer/mips/mips24/prfvsap.ps',/portrait,/square
numoobjects = 34000

;FITS_READ, '/Users/jkrick/palomar/lfc/coadd_g.fits',gdata, gheader
;FITS_READ, '/Users/jkrick/spitzer/irac/ch1/mosaic.fits',irdata, irheader

iracobject = replicate({iracob, irxcenter:0D,irycenter:0D,tmassJ:0D,tmassh:0D,tmassk:0D,irac1:0D,irac2:0D,irac3:0D, irac4:0D,mips24:0D},numoobjects)
mipsobject = replicate({mipsob, ra24:0D, dec24:0D, xcenter24:0D,ycenter24:0D,flux24:0D, fluxerr:0D, snr24:0D, chi24:0D, status24:0D, ap124:0D, ap224:0D,ap324:0D, ap424:0D,ap524:0D, ap624:0D, bad24:0D},numoobjects)
mipsobject44 = replicate({mipsob44, ra24:0D, dec24:0D, xcenter24:0D,ycenter24:0D,flux24:0D, snr24:0D, chi24:0D, status24:0D, ap124:0D, ap224:0D, bad24:0D},numoobjects)
daoobject = replicate({daoob, id:0D, xcenter:0D, ycenter:0D, ra:0D, dec:0D, flux:0D, fluxerr:0D, bckgrnd:0D, iter:0D},numoobjects)
p = 0
openr, luni, "/Users/jkrick/palomar/LFC/catalog/combined_catalog.prt", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, wcsra, wcsdec, xcenter, ycenter, tmassJa,tmassha,tmasska,irac1a,irac2a,irac3a, irac4a,mips24a

   iracobject[p] ={iracob, wcsra, wcsdec, tmassJa,tmassha,tmasska,irac1a,irac2a,irac3a, irac4a,mips24a}
   p = p + 1
   
endwhile
iracobject = iracobject[0:p - 1]
print, "there are ",p," ir objects"
close, luni
free_lun, luni

;----------------------------------------------------------------

m = 0

openr, luni, "/Users/jkrick/Spitzer/mips/mips24/apex_step2_13/mosaic_extract.tbl", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, num, ra, dec, xcenter, ycenter, junk, junk, flux, ferr, snr, chi2, status, ap1, ap2, ap3,ap4,ap5,ap6,bad

   if status gt 0 then begin
      mipsobject[m] ={mipsob, ra, dec, xcenter, ycenter, flux, ferr, snr, chi2, status, ap1, ap2, ap3,ap4,ap5,ap6,bad}
      m = m + 1
   endif
endwhile
mipsobject = mipsobject[0:m - 1]
print, "there are ",m," ir objects", n_elements(mipsobject.ra24)
close, luni
free_lun, luni

;----------------------------------------------------------------

;----------------------------------------------------------------

m = 0

openr, luni, "/Users/jkrick/Spitzer/mips/mips24/apex_step2/mosaic_extract.tbl", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, num, ra, dec, xcenter, ycenter, junk, junk, flux, snr, chi2, status, ap1, ap2, bad

   if status gt 0 then begin
      mipsobject44[m] ={mipsob44, ra, dec, xcenter, ycenter, flux, snr, chi2, status, ap1, ap2, bad}
      m = m + 1
   endif
endwhile
mipsobject44 = mipsobject44[0:m - 1]
print, "there are ",m," ir objects", n_elements(mipsobject44.ra24)
close, luni
free_lun, luni

;----------------------------------------------------------------
;check apex prf vs. Jason's apertures

; create initial arrays
ir=n_elements(iracobject.irxcenter)
m = n_elements(mipsobject.xcenter24)

irmatch=fltarr(ir)
mmatch=fltarr(m)
irmatch[*]=-999
mmatch[*]=-999


print,'Matching  g to ir'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,m-1 do begin

   dist=sphdist(mipsobject[q].ra24,mipsobject[q].dec24,iracobject.irxcenter,iracobject.irycenter,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0003) then begin
      mmatch[q]=ind
;      mobject[ind].match = 1
  endif
endfor

print,"Finished at "+systime()

matched=where(mmatch GT 0)
nonmatched = where(mmatch lt 0)

print, n_elements(matched), " objects in matched"
print, n_elements(mmatch), "objects in mmatch"
print, n_elements(nonmatched), " objects not matched"
;ratio=iracobject(matched).irac1/robject(gmatch[matched]).rfluxa
;plot,iracobject(matched).irac1,ratio,/xlog,psym=3

;print, "n_elements(irmatch[nonmatched])", n_elements(irmatch[nonmatched])
;print, "gobject(irmatch[matched[547]]).gmaga", gobject(irmatch[matched[547]]).gmaga
;print, "gobject(irmatch[nonmatched[547]]).gmaga", gobject(irmatch[nonmatched[547]]).gmaga


openw, outlun, "/Users/jkrick/spitzer/mips/mips24/matched.txt",/get_lun
;print the matched objects 
for num=0, n_elements(matched) - 1 do begin
   
 printf, outlun, mipsobject(matched[num]).ra24,mipsobject(matched[num]).dec24,iracobject(mmatch[matched[num]]).irxcenter, iracobject(mmatch[matched[num]]).irycenter, mipsobject(matched[num]).flux24, iracobject(mmatch[matched[num]]).mips24

endfor

close, outlun
free_lun, outlun

plot, mipsobject(matched).flux24, iracobject(mmatch[matched]).mips24, xtitle = "prf fitting", ytitle = "JS apertures"$
, psym = 2, thick = 3, xthick=3,ythick=3, xrange=[0,1200],yrange=[0,1200]

oplot, findgen(2000), findgen(2000)
;oplot, findgen(2000), 1.1 * findgen(2000)


;------------------------
; check mips 13pix vs 44 pixel normalization radius
m13=n_elements(mipsobject.ra24)
m44 = n_elements(mipsobject44.ra24)

m13match=fltarr(m13)
m44match=fltarr(m44)
m13match[*]=-999
m44match[*]=-999


print,'Matching  13 to 44'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,m13-1 do begin

   dist=sphdist(mipsobject[q].ra24,mipsobject[q].dec24,mipsobject44.ra24,mipsobject44.dec24,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0003) then begin
      m13match[q]=ind
;      mobject[ind].match = 1
  endif
endfor

print,"Finished at "+systime()

matched=where(m13match GT 0)
nonmatched = where(m13match lt 0)

print, n_elements(matched), " objects in matched"
print, n_elements(mmatch), "objects in mmatch"
print, n_elements(nonmatched), " objects not matched"
;ratio=iracobject(matched).irac1/robject(gmatch[matched]).rfluxa
;plot,iracobject(matched).irac1,ratio,/xlog,psym=3

;print, "n_elements(irmatch[nonmatched])", n_elements(irmatch[nonmatched])
;print, "gobject(irmatch[matched[547]]).gmaga", gobject(irmatch[matched[547]]).gmaga
;print, "gobject(irmatch[nonmatched[547]]).gmaga", gobject(irmatch[nonmatched[547]]).gmaga


plothist, mipsobject44(m13match[matched]).flux24 / mipsobject(matched).flux24, xhist, yhist, bin = 0.0001, /noplot
plot, xhist, yhist, thick = 3, xtitle = "ratio of 44 normalized to 13 normalized", xrange= [1.0,1.02]
plot, mipsobject44(m13match[matched]).flux24, mipsobject(matched).flux24,psym = 2,$
xtitle="prf normalized at 44", ytitle = "prf normalized at 13"
oplot, findgen(2000), findgen(2000)


;compare prf to ap flux in same apex run

plot, mipsobject.flux24, mipsobject.ap124, psym = 2, xrange=[0,1000],yrange=[0,1000],$
xtitle= "prf fitted flux", ytitle = "aperture flux r = 6.5"
oplot, findgen(2000), findgen(2000)
plot, mipsobject.flux24, mipsobject.ap224, psym = 2, xrange=[0,1000],yrange=[0,1000],$
xtitle= "prf fitted flux", ytitle = "aperture flux r = 22"
oplot, findgen(2000), findgen(2000)
plot, mipsobject.ap124, mipsobject.ap224, psym = 2, xrange=[0,1000],yrange=[0,1000],$
xtitle= "aperture flux r = 6.5", ytitle = "aperture flux r = 22"

x = [6.5,14.25,22,29,35,40]
y1 = [mipsobject(8).ap124, mipsobject(8).ap224,mipsobject(8).ap324, mipsobject(8).ap424,mipsobject(8).ap524, mipsobject(8).ap624]
y2 = [mipsobject(11).ap124, mipsobject(11).ap224,mipsobject(11).ap324, mipsobject(11).ap424,mipsobject(11).ap524, mipsobject(11).ap624]
y3 = [mipsobject(13).ap124, mipsobject(13).ap224,mipsobject(13).ap324, mipsobject(13).ap424,mipsobject(13).ap524, mipsobject(13).ap624]
plot,x,y1
oplot, x, y2
oplot, x, y3



;------------------------------------------------------------------
;look at daophot psf fitted fluxes

m = 0

openr, luni, "/Users/jkrick/Spitzer/mips/mips24/dao/mips24.phot", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, num, xcen,ycen, mag, magerr, back, niter, junk,junk,junk
   ;need to convert from magnitudes back to fluxes.
   daoobject[m] = {daoob, num, xcen, ycen, 0.0,0.0, 10^((25-mag)/2.5),  10^((25-magerr)/2.5), back, niter}
   m = m + 1

endwhile
daoobject = daoobject[0:m - 1]
print, "there are ",m," dao objects"
close, luni
free_lun, luni

;need to go from x,y in pixels to ra, dec
FITS_READ, '/Users/jkrick/spitzer/mips/mips24/daophot/mosaic.fits',data, header
xyad, header, daoobject.xcenter, daoobject.ycenter, ra, dec
daoobject.ra = ra
daoobject.dec = dec


;match daophot to apex

; create initial arrays
dao=n_elements(daoobject.ra)
apex = n_elements(mipsobject.ra24)

daomatch=fltarr(dao)
apexmatch=fltarr(apex)
daomatch[*]=-999
apexmatch[*]=-999


print,'Matching  dao to apex'
print,"Starting at "+systime()
dist=daomatch
dist[*]=0

for q=0,apex-1 do begin

   dist=sphdist(mipsobject[q].ra24,mipsobject[q].dec24,daoobject.ra,daoobject.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0005) then begin
      apexmatch[q]=ind
;      mobject[ind].match = 1
  endif
endfor

print,"Finished at "+systime()

matched=where(apexmatch GT 0)
nonmatched = where(apexmatch lt 0)

print, n_elements(matched), " objects in matched"
print, n_elements(apexmatch), "objects in mmatch"
print, n_elements(nonmatched), " objects not matched"

plot, mipsobject(matched).flux24, daoobject(apexmatch[matched]).flux, xtitle = "apex prf fitting", ytitle = "daophot psf fitting"$
, psym = 2, thick = 3, xthick=3,ythick=3;, xrange=[0,1200],yrange=[0,1200]

oplot, findgen(2000), findgen(2000)


;--------------------------------------------------------------------
;check daophot vs. Jason's apertures

; create initial arrays
dao=n_elements(daoobject.ra)
js = n_elements(iracobject.irxcenter)

daomatch=fltarr(dao)
jsmatch=fltarr(js)
daomatch[*]=-999
jsmatch[*]=-999


print,'Matching  dao to js'
print,"Starting at "+systime()
dist=daomatch
dist[*]=0

for q=0, dao-1 do begin

   dist=sphdist( daoobject[q].ra, daoobject[q].dec,iracobject.irxcenter,iracobject.irycenter,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0007) then begin
      daomatch[q]=ind
;      mobject[ind].match = 1
  endif
endfor

print,"Finished at "+systime()

matched=where(daomatch GT 0)
nonmatched = where(daomatch lt 0)

print, n_elements(matched), " objects in matched"
print, n_elements(daomatch), "objects in mmatch"
print, n_elements(nonmatched), " objects not matched"

plot,1.43*daoobject(matched).flux,iracobject(daomatch[matched]).mips24, ytitle = "js  apertures", xtitle = "daophot psf fitting"$
, psym = 2, thick = 3, xthick=3,ythick=3, xrange=[0,1200],yrange=[0,1200]

oplot, findgen(2000), findgen(2000)







undefine, iracobject
undefine, mipsobject
undefine, mipsobject44
undefine, daoobject
ps_close, /noprint,/noid

end


