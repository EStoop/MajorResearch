pro make_catalog


close, /all

!P.multi = [0,2,2]

numoobjects = 74000
colorarr = fltarr(numoobjects)
magarr = fltarr(numoobjects)
FITS_READ, '/Users/jkrick/palomar/lfc/coadd_g.fits',gdata, gheader
FITS_READ, '/Users/jkrick/spitzer/irac/ch1/mosaic.fits',irdata, irheader

uobject = replicate({uob, uxcenter:0D,uycenter:0D,ufluxa:0D,umaga:0D,umagerra:0D,ufwhm:0D,uisoarea:0D,uellip:0D,ufluxmax:0D, uclassstar:0D, uflags:0D},numoobjects)
gobject = replicate({gob, gxcenter:0D,gycenter:0D,gfluxa:0D,gmaga:0D,gmagerra:0D,gfwhm:0D,gisoarea:0D,gellip:0D,gfluxmax:0D,match:0D, gclassstar:0D, gflags:0D},numoobjects)
robject = replicate({rob, rxcenter:0D,rycenter:0D,rfluxa:0D,rmaga:0D,rmagerra:0D,rfwhm:0D,risoarea:0D, rellip:0D,rfluxmax:0D, rclassstar:0D,rflags:0D},numoobjects)
iobject = replicate({iob, ixcenter:0D,iycenter:0D,ifluxa:0D,imaga:0D,imagerra:0D,ifwhm:0D,iisoarea:0D, iellip:0D,ifluxmax:0D, iclassstar:0D, iflags:0D},numoobjects)

acsobject =replicate({acsob,acsconcentration:0D,acsgini:0D,acsmag:0D, acsmu:0D,acsellip:0D,acsflag:0D,acsra:0D,acsdec:0D,acsisoarea:0D,acstheta:0D,acssnr:0D,acsflux:0D, acssky:0D,acsskyrms:0D,acscentralmu:0D,acsassym:0D,acshlightrad:0D,acsm20:0D,acsclassstar:0D,acsmagerr:0D,acsfwhm:0D, acsfluxmax:0D},80000)
;-------
flamjobject = replicate({flamjob,flamjxcenter:0D, flamjycenter:0D, flamjra:0D, flamjdec:0D, flamjmag:0D, flamjmagerr:0D, flamjfwhm:0D, flamjisoarea:0D, flamjellip:0D },numoobjects)
wircjobject = replicate({wircjob, wircjxcenter:0D, wircjycenter:0D, wircjra:0D, wircjdec:0D, wircjmag:0D, wircjmagerr:0D, wircjfwhm:0D, wircjisoarea:0D, wircjellip:0D },numoobjects)
wirchobject = replicate({wirchob, wirchxcenter:0D, wirchycenter:0D, wirchra:0D, wirchdec:0D, wirchmag:0D, wirchmagerr:0D, wirchfwhm:0D, wirchisoarea:0D, wirchellip:0D },numoobjects)
wirckobject = replicate({wirckob, wirckxcenter:0D, wirckycenter:0D, wirckra:0D, wirckdec:0D, wirckmag:0D, wirckmagerr:0D, wirckfwhm:0D, wirckisoarea:0D, wirckellip:0D },numoobjects)

mips24object = replicate({mips24ob, mips24xcenter:0D, mips24ycenter:0D, mips24ra:0D, mips24dec:0D, mips24flux:0D, mips24fluxerr:0D, mips24mag:0D, mips24magerr:0D, mips24bckgrnd:0D, iter:0D,sharp:0D},numoobjects)

mips70object = replicate({mips70ob, mips70xcenter:0D, mips70ycenter:0D, mips70ra:0D, mips70dec:0D, mips70flux:0D, mips70fluxerr:0D, mips70mag:0D, mips70magerr:0D, mips70bckgrnd:0D, mips70isoarea:0D,mips70fwhm:0D, mips70flags:0D},numoobjects)

iracobject = replicate({iracob, irxcenter:0D,irycenter:0D,tmassJflux:0D,tmasshflux:0D,tmasskflux:0D,irac1flux:0D,irac2flux:0D,irac3flux:0D, irac4flux:0D,tmassJmag:0D,tmasshmag:0D,tmasskmag:0D,irac1mag:0D,irac2mag:0D,irac3mag:0D, irac4mag:0D,mips24:0D},numoobjects)

infraredobject = replicate({spitzob, INHERITS iracob,INHERITS mips24ob,INHERITS mips70ob},numoobjects)
;=======
object= replicate({ob, ra:0D, dec:0D, INHERITS uob, INHERITS gob, INHERITS rob, INHERITS iob, INHERITS flamjob, INHERITS wircjob, INHERITS wirchob, INHERITS wirckob, INHERITS spitzob, zphot:0D, chisq:0D, prob:0D, spt:0D, age:0D, av:0D, Mabs:0D, mass:0D, plusmasserr:0D, minusmasserr:0D, masschi:0D, massprob:0D, massage:0D, model:' ',INHERITS acsob}, numoobjects)

;----------------------------------------------------------------
;u catalog has the smallest number of junk obs at the end of the catalog, so use it for size.
nu = 0.
openr, lung, "/Users/jkrick/palomar/LFC/catalog/SExtractor.u.cat", /get_lun
WHILE (NOT EOF(lung)) DO BEGIN
    READF, lung, o, xwcsg, ywcsg, xcenterg, ycenterg, fluxautog, magautog, magerrg, magaper, magapererr, fwhmg, isoareag, fluxmaxg, ellipu,gnum, class_star, flags

   if (magaper gt 21.0)  then begin
;   if (magautog gt 22.5)  then begin
       magautog = magaper
       magerrg = magapererr
    endif

   if gnum ne 0 then begin

;    if (fwhmg gt 0 ) and fluxautog gt 0   then begin
       uobject[gnum-1] ={uob, xwcsg, ywcsg,fluxautog, magautog, magerrg,fwhmg, isoareag, ellipu,fluxmaxg, class_star, flags}
       nu = nu + 1
    ENDIF

 endwhile
uobject = uobject[0:nu - 1]
print, "there are ",nu," u objects"
close, lung
free_lun, lung
;----------------------------------------------------------------

lg = 0.
openr, lung, "/Users/jkrick/palomar/LFC/catalog/SExtractor.g.cat", /get_lun
WHILE (NOT EOF(lung)) DO BEGIN
    READF, lung, o, xwcsg, ywcsg, xcenterg, ycenterg, fluxautog, magautog, magerrg, magaper, magapererr,fwhmg, isoareag, fluxmaxg, ellipg, class_star, flags
;    if o eq 20389 then print, xwcsg, ywcsg, fluxmaxg
  ;  if (fluxautog gt 0 )   then begin
   if (magaper gt 21.0)  then begin
;   if (magautog gt 22.5)  then begin
       magautog = magaper
       magerrg = magapererr
    endif

;   if fwhmg gt 5 then begin
       gobject[lg] ={gob, xwcsg, ywcsg,fluxautog, magautog, magerrg,fwhmg, isoareag, ellipg, fluxmaxg,-1, class_star, flags}
;       if o eq 20389 then print, "found g object", l, gobject[l].gxcenter, gobject[l].gycenter
       lg = lg + 1
 ;   ENDIF

 endwhile
gobject = gobject[0:nu - 1]
print, "there are ",lg," g objects"
close, lung
free_lun, lung

;set them all as unmatched objects
;gobject.match = -1
   
;----------------------------------------------------------------
i = 0.
openr, lunr, "/Users/jkrick/palomar/LFC/catalog/SExtractor.r.cat", /get_lun
WHILE (NOT EOF(lunr)) DO BEGIN
    READF, lunr, o, xwcsr, ywcsr, xcenterr, ycenterr, fluxautor, magautor, magautoerrr, magaper, magapererr,fwhmr, isoarear, fluxmaxr, ellipr,gnum, class_star, flags

    if (magaper gt 21.0)  then begin
;    if (magautor gt 22.5)  then begin
       magautor = magaper
       magautoerrr = magapererr
    endif

   if gnum ne 0 then begin
    
 ;   if fwhmr gt 5 then begin
       robject[gnum-1] ={rob, xwcsr, ywcsr,fluxautor, magautor, magautoerrr,fwhmr, isoarear, ellipr,fluxmaxr, class_star, flags}
      i = i + 1
    endif


 endwhile
robject = robject[0:nu-1]
print, "there are ",i," r objects"
close, lunr
free_lun, lunr
;--------------------------------------------------------
;am changing the zeropoints here, after the catalogs are made.

m = 0.
openr, lung, "/Users/jkrick/palomar/LFC/catalog/SExtractor.i.cat", /get_lun
WHILE (NOT EOF(lung)) DO BEGIN
    READF, lung, o, xwcsg, ywcsg, xcenterg, ycenterg, fluxautog, magautog, magerrg, magaper, magapererr, fwhmg, isoareag, fluxmaxg, ellipi,gnum, class_star, flags

   if (magaper gt 21.0)  then begin
;   if (magautog gt 22.5)  then begin
       magautog = magaper
       magerrg = magapererr
    endif

   if gnum ne 0 then begin

;   if (fwhmg gt 5 )   then begin
       iobject[gnum-1] ={iob, xwcsg, ywcsg,fluxautog, magautog, magerrg,fwhmg, isoareag, ellipi,fluxmaxg, class_star, flags}
       m = m + 1
    ENDIF

 endwhile
iobject = iobject[0:nu - 1]
print, "there are ",m," i objects"
close, lung
free_lun, lung

;----------------------------------------------------------------

p = 0
openr, luni, "/Users/jkrick/nep/flamingos/jband_noheader.txt", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, junk,x_image, y_image,junk,junk,junk,junk,junk,junk,mag_best, magerr_best, junk,junk, junk,isoareaimage, ra,dec,junk,junk,junk, junk,fwhm_image, junk,ell, junk 

   flamjobject[p] ={flamjob, x_image, y_image, ra, dec, mag_best, magerr_best, fwhm_image, isoareaimage, ell}
   p = p + 1
   
 endwhile
flamjobject = flamjobject[0:p - 1]
print, "there are ",p,"flam  j objects"
close, luni
free_lun, luni
;----------------------------------------------------------------


p = 0
openr, luni, "/Users/jkrick/palomar/wirc/jband_wirc_nohead.txt", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, junk,x_image, y_image,junk,junk,junk,junk,junk,junk,mag_best, magerr_best, junk,junk, junk,isoareaimage, ra,dec,junk,junk,junk, junk,fwhm_image, junk,ell, junk 

   wircjobject[p] ={wircjob, x_image, y_image, ra, dec, mag_best, magerr_best, fwhm_image, isoareaimage, ell}
   p = p + 1
   
 endwhile
wircjobject = wircjobject[0:p - 1]
print, "there are ",p," wirc j objects"
close, luni
free_lun, luni
;----------------------------------------------------------------

p = 0
openr, luni, "/Users/jkrick/palomar/wirc/hband_wirc_nohead.txt", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, junk,x_image, y_image,junk,junk,junk,junk,junk,junk,mag_best, magerr_best, junk,junk, junk,isoareaimage, ra,dec,junk,junk,junk, junk,fwhm_image, junk,ell, junk 

   wirchobject[p] ={wirchob, x_image, y_image, ra, dec, mag_best, magerr_best, fwhm_image, isoareaimage, ell}
   p = p + 1
   
 endwhile
wirchobject = wirchobject[0:p - 1]
print, "there are ",p," wirc h objects"
close, luni
free_lun, luni
;----------------------------------------------------------------

p = 0
openr, luni, "/Users/jkrick/palomar/wirc/kband_wirc_nohead.txt", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, junk,x_image, y_image,junk,junk,junk,junk,junk,junk,mag_best, magerr_best, junk,junk, junk,isoareaimage, ra,dec,junk,junk,junk, junk,fwhm_image, junk,ell, junk 

   wirckobject[p] ={wirckob, x_image, y_image, ra, dec, mag_best, magerr_best, fwhm_image, isoareaimage, ell}
   p = p + 1
   
 endwhile
wirckobject = wirckobject[0:p - 1]
print, "there are ",p," wirc k objects"
close, luni
free_lun, luni
;----------------------------------------------------------------
;ACS
mag = 0.D
p = 0l
;openr, luni, "/Users/jkrick/hst/acscat_merged.txt", /get_lun
openr, luni, "/Users/jkrick/hst/raw/wholeacs.cat", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
;   READF, luni, c,g,mag,mu,x,y,s,a,b,flag,ra,dec,isoarea,theta,snr,flux, sky,skyrms,cmu,asym,hlr,m20,class,magerr
   readf, luni, NUMBER,X_WORLD,Y_WORLD,X_IMAGE,Y_IMAGE,MAG_AUTO,FLUX_BEST,FLUXERR_BEST,MAG_BEST,MAGERR_BEST,BACKGROUND,FLUX_MAX,ISOAREA_IMAGE,ALPHA_J2000,DELTA_J2000,THETA_IMAGE,MU_THRESHOLD,MU_MAX,FLAGS,FWHM_IMAGE,CLASS_STAR,ELLIPTICITY

   if magerr_best eq 0. then magerr_best = 0.03
                                ;  acsobject[p] ={acsob,c,g,mag, mu,1-b/a,flag,ra,dec,isoarea,theta,snr,flux, sky,skyrms,cmu,asym,hlr,m20,class,magerr}
   acsobject[p]={acsob,0.,0., mag_best, 0., ellipticity, flags, x_world, y_world, isoarea_image, theta_image, 0., flux_best, background, 0., 0., 0., 0., 0., class_star, magerr_best,fwhm_image,flux_max}
   p = p + 1
   
 endwhile
acsobject = acsobject[0:p-1]
print, "there are ",p," acs objects"

close, luni
free_lun, luni


;----------------------------------------------------------------
p = 0
openr, luni, "/Users/jkrick/palomar/LFC/catalog/combined_catalog.prt", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, wcsra, wcsdec, xcenter, ycenter, tmassJa,tmassha,tmasska,irac1a,irac2a,irac3a, irac4a,mips24a
  
  if tmassja le 0 then begin
      tmassja = -99
      tmassjab = -99
   endif else begin
      tmassjab = -2.5*alog10(tmassja) + 8.926 
   endelse

  if tmassha le 0 then begin
      tmassha = -99
      tmasshab = -99
   endif else begin
      tmasshab = -2.5*alog10(tmassha) + 8.926 
   endelse

  if tmasska le 0 then begin
      tmasska = -99
      tmasskab = -99
   endif else begin
      tmasskab = -2.5*alog10(tmasska) + 8.926 
   endelse

   if irac1a lt 0 then begin
      mag1 = -99. 
      irac1a=-99.
   endif else begin
      mag1 = 8.926 - 2.5*alog10(1E-6*irac1a) 
   endelse

   if irac2a lt 0 then begin
      mag2 = -99. 
      irac2a=-99.
   endif else begin
      mag2 = 8.926 - 2.5*alog10(1E-6*irac2a) 
   endelse

   if irac3a lt 0 then begin
      mag3 = -99. 
      irac3a=-99.
   endif else begin
      mag3 = 8.926- 2.5*alog10(1E-6*irac3a)
   endelse

   if irac4a lt 0 then begin
      mag4 = -99. 
      irac4a=-99.
   endif else begin
      mag4 = 8.926 - 2.5*alog10(1E-6*irac4a)
   endelse


   iracobject[p] ={iracob, wcsra, wcsdec, tmassJa,tmassha,tmasska,irac1a,irac2a,irac3a, irac4a,tmassjab,tmasshab,tmasskab,mag1,mag2,mag3,mag4,mips24a}
   p = p + 1
   
 endwhile
iracobject = iracobject[0:p - 1]
print, "there are ",p," ir objects, plothist"

plothist, iracobject.tmasskflux, title = 'after reading file'

close, luni
free_lun, luni

;----------------------------------------------------------------
;mips24
w = 0
num = long(0)
xcen = 0.
ycen = 0.
mag = 0.
magerr = 0.
back = 0.
niter = 0.
sh = 0.
junk = 0.

openr, lunw, "/Users/jkrick/Spitzer/mips/mips24/dao/mips24.phot", /get_lun
WHILE (NOT EOF(lunw)) DO BEGIN
;   print, w
   READF, lunw, num, xcen,ycen, mag, magerr, back, niter, sh, junk, junk
;   print, xcen
   ;need to convert from magnitudes back to fluxes.
   mips24object[w] = {mips24ob, xcen, ycen, 0.0,0.0, 1.432*(10^((25-mag)/2.5)),  10^((25-(mag -magerr))/2.5) -(10^((25-mag)/2.5)) , mag,magerr,back, niter,sh}
   w = w +1
endwhile
mips24object = mips24object[0:w - 1]
print, "there are ",w," mips24 objects"
close, lunw
free_lun, lunw

;need to go from x,y in pixels to ra, dec
FITS_READ, '/Users/jkrick/spitzer/mips/mips24/daophot/mosaic.fits',data, header
xyad, header, mips24object.mips24xcenter, mips24object.mips24ycenter, ra, dec
mips24object.mips24ra = ra
mips24object.mips24dec = dec
;----------------------------------------------------------------
;mips70
mag = 0.D
num = 0
ra = 0.
dec = 0.
x=0.
y=x
fluxiso=x
fluxerriso=x
fluxauto=x
fluxerrauto=x
kron=x
back=x
isoarea=x
a=x
b=x
theta=x
flag=x
fwhm=x
elon=x
ellip=x
p = 0l
openr, luni, "/Users/jkrick/spitzer/mips/mips70/mosaic/Combine/mosaic.asc", /get_lun, error=err
if (err ne 0.) then print, "file didn't open", ﻿!ERROR_STATE.MSG  

WHILE (NOT EOF(luni)) DO BEGIN

   READF, luni, num,ra, dec, x, y, fluxiso, fluxerriso,fluxauto, fluxerrauto,kron,back,isoarea,a, b, theta, flag, fwhm,elon, ellip
; I do not know what the magnitudes or the correct fluxes are.  there is an area issue.

   mips70object[p] ={mips70ob, x,y, ra, dec, fluxauto*23.5045, fluxerrauto*23.5045, fluxauto, fluxerrauto,back, isoarea, fwhm, flag}

   p = p + 1
   
 endwhile
mips70object = mips70object[0:p-1]
print, "there are ",p," mips70 objects"

close, luni
free_lun, luni


;-----------------------------------------------
;now match mips24 to irac before I even deal with the optical.
;there shouldn't be any mips detections which are not detected in irac.
; create initial arrays
m=n_elements(mips24object.mips24ra)
ir=n_elements(iracobject.irxcenter)

irmatch=fltarr(ir)
mmatch=fltarr(m)
irmatch[*]=-999
mmatch[*]=-999


print,'Matching  mips24 to ir'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,ir-1 do begin

   dist=sphdist(iracobject[q].irxcenter,iracobject[q].irycenter, mips24object.mips24ra,mips24object.mips24dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008) then begin
      irmatch[q]=ind
   endif 
endfor

print,"Finished at "+systime()

matched=where(irmatch GT 0)
nonmatched = where(irmatch lt 0)

print, n_elements(matched),"matched"
print, n_elements(nonmatched),"nonmatched"

;catalog  the matched and nonmatched objects 
for num=0, n_elements(matched) - 1 do begin
   
   infraredobject[num] = {spitzob, iracobject(matched[num]).irxcenter, iracobject(matched[num]).irycenter,iracobject(matched[num]).tmassJflux, iracobject(matched[num]).tmasshflux, iracobject(matched[num]).tmasskflux, iracobject(matched[num]).irac1flux, iracobject(matched[num]).irac2flux,iracobject(matched[num]).irac3flux, iracobject(matched[num]).irac4flux, iracobject(matched[num]).tmassJmag, iracobject(matched[num]).tmasshmag, iracobject(matched[num]).tmasskmag, iracobject(matched[num]).irac1mag, iracobject(matched[num]).irac2mag,iracobject(matched[num]).irac3mag, iracobject(matched[num]).irac4mag,iracobject(matched[num]).mips24, mips24object(irmatch[matched[num]]).mips24xcenter , mips24object(irmatch[matched[num]]).mips24ycenter , mips24object(irmatch[matched[num]]).mips24ra , mips24object(irmatch[matched[num]]).mips24dec , mips24object(irmatch[matched[num]]).mips24flux , mips24object(irmatch[matched[num]]).mips24fluxerr ,  mips24object(irmatch[matched[num]]).mips24mag , mips24object(irmatch[matched[num]]).mips24magerr ,mips24object(irmatch[matched[num]]).mips24bckgrnd , mips24object(irmatch[matched[num]]).iter ,mips24object(irmatch[matched[num]]).sharp,-99.,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,99.0,99.0,-99.0,-99.0,-99.0}

;print, "matched", num, infraredobject[num].irxcenter, infraredobject[num].irycenter, infraredobject[num].ra, infraredobject[num].dec

endfor
count = n_elements(matched)
for num= 0 , n_elements(nonmatched) -1. do begin
   
   infraredobject[count] = {spitzob, iracobject(nonmatched[num]).irxcenter, iracobject(nonmatched[num]).irycenter,iracobject(nonmatched[num]).tmassJflux,iracobject(nonmatched[num]).tmasshflux,iracobject(nonmatched[num]).tmasskflux,iracobject(nonmatched[num]).irac1flux,iracobject(nonmatched[num]).irac2flux,iracobject(nonmatched[num]).irac3flux, iracobject(nonmatched[num]).irac4flux,iracobject(nonmatched[num]).tmassJmag,iracobject(nonmatched[num]).tmasshmag,iracobject(nonmatched[num]).tmasskmag,iracobject(nonmatched[num]).irac1mag,iracobject(nonmatched[num]).irac2mag,iracobject(nonmatched[num]).irac3mag, iracobject(nonmatched[num]).irac4mag,iracobject(nonmatched[num]).mips24, -99.0,-99.0,-99.0,-99.0,-99.0,-99.0,99.0,99.0,-99.0,-99.0,-99.0,-99.,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,99.0,99.0,-99.0,-99.0,-99.0}
count = count + 1

endfor

infraredobject = infraredobject[0:n_elements(matched) + n_elements(nonmatched)]
plothist, infraredobject.tmasskflux, title = 'after infraredobject'
;----------------------------------------------------------------


;print, "filling object structure"
;print, "object[0].uxcenter", object[0].uxcenter


; create initial arrays
r=n_elements(robject.rxcenter)
g=n_elements(gobject.gxcenter)
u=n_elements(uobject.uxcenter)
i=n_elements(iobject.ixcenter)
ir=n_elements(infraredobject.irxcenter)
print, ir, "ir"
irmatch=fltarr(ir)
gmatch=fltarr(g)
irmatch[*]=-999
gmatch[*]=-999


print,'Matching  g to ir'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,ir-1 do begin

   dist=sphdist(infraredobject[q].irxcenter,infraredobject[q].irycenter,gobject.gxcenter,gobject.gycenter,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008) then begin
      irmatch[q]=ind
      gobject[ind].match = 1
  endif
endfor

print,"Finished at "+systime()

matched=where(irmatch GT 0)
nonmatched = where(irmatch lt 0)

print, n_elements(matched), " objects in matched"
print, n_elements(irmatch), "objects in irmatch"
print, n_elements(nonmatched), " objects not matched"
;ratio=infraredobject(matched).irac1/robject(gmatch[matched]).rfluxa
;plot,infraredobject(matched).irac1,ratio,/xlog,psym=3

;print, "n_elements(irmatch[nonmatched])", n_elements(irmatch[nonmatched])
;print, "gobject(irmatch[matched[547]]).gmaga", gobject(irmatch[matched[547]]).gmaga
;print, "gobject(irmatch[nonmatched[547]]).gmaga", gobject(irmatch[nonmatched[547]]).gmaga


openw, outlun, "/Users/jkrick/palomar/lfc/catalog/opticalir.txt", /get_lun
openw, outlun2, "/Users/jkrick/palomar/lfc/catalog/hyperz_cat.txt", /get_lun
openw, outlun3, "/Users/jkrick/palomar/lfc/catalog/junk", /get_lun
openw, outlun4, "/Users/jkrick/palomar/lfc/catalog/junkir", /get_lun

;print the matched objects 
for num=0, n_elements(matched) - 1 do begin
   
   object[num] = {ob, infraredobject(matched[num]).irxcenter,infraredobject(matched[num]).irycenter, uobject(irmatch[matched[num]]).uxcenter,uobject(irmatch[matched[num]]).uycenter,uobject(irmatch[matched[num]]).ufluxa,uobject(irmatch[matched[num]]).umaga,uobject(irmatch[matched[num]]).umagerra,uobject(irmatch[matched[num]]).ufwhm,uobject(irmatch[matched[num]]).uisoarea,uobject(irmatch[matched[num]]).uellip,uobject(irmatch[matched[num]]).ufluxmax, uobject(irmatch[matched[num]]).uclassstar,uobject(irmatch[matched[num]]).uflags, gobject(irmatch[matched[num]]).gxcenter,gobject(irmatch[matched[num]]).gycenter,gobject(irmatch[matched[num]]).gfluxa,gobject(irmatch[matched[num]]).gmaga,gobject(irmatch[matched[num]]).gmagerra,gobject(irmatch[matched[num]]).gfwhm,gobject(irmatch[matched[num]]).gisoarea,gobject(irmatch[matched[num]]).gellip,gobject(irmatch[matched[num]]).gfluxmax , gobject(irmatch[matched[num]]).match ,gobject(irmatch[matched[num]]).gclassstar,gobject(irmatch[matched[num]]).gflags, robject(irmatch[matched[num]]).rxcenter,robject(irmatch[matched[num]]).rycenter,robject(irmatch[matched[num]]).rfluxa,robject(irmatch[matched[num]]).rmaga,robject(irmatch[matched[num]]).rmagerra,robject(irmatch[matched[num]]).rfwhm,robject(irmatch[matched[num]]).risoarea, robject(irmatch[matched[num]]).rellip, robject(irmatch[matched[num]]).rfluxmax,robject(irmatch[matched[num]]).rclassstar,robject(irmatch[matched[num]]).rflags,iobject(irmatch[matched[num]]).ixcenter,iobject(irmatch[matched[num]]).iycenter,iobject(irmatch[matched[num]]).ifluxa,iobject(irmatch[matched[num]]).imaga,iobject(irmatch[matched[num]]).imagerra,iobject(irmatch[matched[num]]).ifwhm,iobject(irmatch[matched[num]]).iisoarea, iobject(irmatch[matched[num]]).iellip,iobject(irmatch[matched[num]]).ifluxmax,iobject(irmatch[matched[num]]).iclassstar,iobject(irmatch[matched[num]]).iflags, -99.0,-99.0,-99.0,-99.0,99.0,99.0,-99.0,-99.0,-99.0, -99.0,-99.0,-99.0,-99.0,99.0,99.0,-99.0,-99.0,-99.0, -99.0,-99.0,-99.0,-99.0,99.0,99.0,-99.0,-99.0,-99.0, -99.0,-99.0,-99.0,-99.0,99.0,99.0,-99.0,-99.0,-99.0,infraredobject(matched[num]).irxcenter,infraredobject(matched[num]).irycenter,infraredobject(matched[num]).tmassJflux,infraredobject(matched[num]).tmasshflux,infraredobject(matched[num]).tmasskflux,infraredobject(matched[num]).irac1flux,infraredobject(matched[num]).irac2flux,infraredobject(matched[num]).irac3flux, infraredobject(matched[num]).irac4flux,infraredobject(matched[num]).tmassJmag,infraredobject(matched[num]).tmasshmag,infraredobject(matched[num]).tmasskmag,infraredobject(matched[num]).irac1mag,infraredobject(matched[num]).irac2mag,infraredobject(matched[num]).irac3mag, infraredobject(matched[num]).irac4mag,infraredobject(matched[num]).mips24, infraredobject(matched[num]).mips24xcenter , infraredobject(matched[num]).mips24ycenter , infraredobject(matched[num]).mips24ra , infraredobject(matched[num]).mips24dec , infraredobject(matched[num]).mips24flux , infraredobject(matched[num]).mips24fluxerr , infraredobject(matched[num]).mips24mag, infraredobject(matched[num]).mips24magerr , infraredobject(matched[num]).mips24bckgrnd , infraredobject(matched[num]).iter ,infraredobject(matched[num]).sharp,-99.,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,99.0,99.0,-99.0,-99.0,-99.0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., "none", 0.,0.,99.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-99.0}

 

endfor

;print the ir objects which don't have an optical match 
for num2=0, n_elements(nonmatched) - 1 do begin
   
   adxy, gheader, infraredobject(nonmatched[num2]).irxcenter,infraredobject(nonmatched[num2]).irycenter, xcen, ycen
   printf,outlun3, xcen,ycen
   adxy, irheader, infraredobject(nonmatched[num2]).irxcenter,infraredobject(nonmatched[num2]).irycenter, xcen2, ycen2
   printf,outlun4, xcen2,ycen2


   object[num]={ob, infraredobject(nonmatched[num2]).irxcenter ,infraredobject(nonmatched[num2]).irycenter ,-99.0,-99.0,-99.0,99.0,-99.0,-99.0,-99.0,-99.0, -99.0,-99.0,-99.0,-99.0,99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,99.0, -99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,99.0, -99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,99.0,-99.0,-99.0,-99.0,-99.0, -99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,infraredobject(nonmatched[num2]).irxcenter ,infraredobject(nonmatched[num2]).irycenter ,infraredobject(nonmatched[num2]).tmassJflux ,infraredobject(nonmatched[num2]).tmasshflux ,infraredobject(nonmatched[num2]).tmasskflux ,infraredobject(nonmatched[num2]).irac1flux ,infraredobject(nonmatched[num2]).irac2flux ,infraredobject(nonmatched[num2]).irac3flux , infraredobject(nonmatched[num2]).irac4flux,infraredobject(nonmatched[num2]).tmassJmag ,infraredobject(nonmatched[num2]).tmasshmag ,infraredobject(nonmatched[num2]).tmasskmag ,infraredobject(nonmatched[num2]).irac1mag ,infraredobject(nonmatched[num2]).irac2mag ,infraredobject(nonmatched[num2]).irac3mag , infraredobject(nonmatched[num2]).irac4mag,infraredobject(nonmatched[num2]).mips24, infraredobject(nonmatched[num2]).mips24xcenter , infraredobject(nonmatched[num2]).mips24ycenter , infraredobject(nonmatched[num2]).mips24ra , infraredobject(nonmatched[num2]).mips24dec , infraredobject(nonmatched[num2]).mips24flux , infraredobject(nonmatched[num2]).mips24fluxerr ,  infraredobject(nonmatched[num2]).mips24mag , infraredobject(nonmatched[num2]).mips24magerr , infraredobject(nonmatched[num2]).mips24bckgrnd , infraredobject(nonmatched[num2]).iter ,infraredobject(nonmatched[num2]).sharp,-99.,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,99.0,99.0,-99.0,-99.0,-99.0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,"none", 0.,0.,99.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-99.}
 
   num = num+1

endfor

print, "num between", num
num = num - 1
;;print the optical objects which are detected in 2 optical bands and not in the IR catalog

for num3 = 0l, n_elements(gobject.match) -1 do begin
   if gobject[num3].match lt 0 then begin
      ;this object did not match to an irac object
      ;is it a real object?

      if (gobject[num3].gfluxa GT 0) and (gobject[num3].gmaga lt 30)$  
         and (gobject[num3].gmagerra lt 2) $
         and (robject[num3].rfluxa GT 0) and (robject[num3].rmaga lt 30)$  
         and (robject[num3].rmagerra lt 2) $ ;and (iobject[num3].ifluxa GT 0) and (iobject[num3].imaga lt 30)$  
 ;        and (iobject[num3].imagerra lt 2) $
      then begin
         ;it is both in r ,g, i catalogs and apparently not noise
         ;so add it to the catalog
 ;         adxy, gheader, gobject[num3].gxcenter,gobject[num3].gycenter, xcen, ycen
 ;       printf,outlun3, xcen,ycen
 ;         adxy, irheader, gobject[num3].gxcenter,gobject[num3].gycenter, xcen2, ycen2
 ;       printf,outlun4, xcen2,ycen2

         object[num]={ob, gobject[num3].gxcenter,gobject[num3].gycenter, uobject[num3].uxcenter,uobject[num3].uycenter,uobject[num3].ufluxa ,uobject[num3].umaga ,uobject[num3].umagerra ,uobject[num3].ufwhm ,uobject[num3].uisoarea ,uobject[num3].uellip ,uobject[num3].ufluxmax ,uobject[num3].uclassstar ,uobject[num3].uflags,gobject[num3].gxcenter,gobject[num3].gycenter,gobject[num3].gfluxa ,gobject[num3].gmaga ,gobject[num3].gmagerra ,gobject[num3].gfwhm ,gobject[num3].gisoarea ,gobject[num3].gellip, gobject[num3].gfluxmax,gobject[num3].match,gobject[num3].gclassstar ,gobject[num3].gflags,robject[num3].rxcenter,robject[num3].rycenter,robject[num3].rfluxa ,robject[num3].rmaga ,robject[num3].rmagerra ,robject[num3].rfwhm ,robject[num3].risoarea ,robject[num3].rellip, robject[num3].rfluxmax,robject[num3].rclassstar ,robject[num3].rflags,iobject[num3].ixcenter,iobject[num3].iycenter,iobject[num3].ifluxa ,iobject[num3].imaga ,iobject[num3].imagerra ,iobject[num3].ifwhm ,iobject[num3].iisoarea ,iobject[num3].iellip, iobject[num3].ifluxmax,iobject[num3].iclassstar ,iobject[num3].iflags,-99.0,-99.0,-99.0,-99.0,99.0,99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,99.0,99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,99.0,99.0,-99.0,-99.0,-99.0,-99.0, -99.0, -99.0, -99.0, 99.0,99.0, -99.0, -99.0, -99.0,-99.0,-99.0,-99.0,-99.0, -99.0,-99.0,-99.0,-99.0,-99.0,99.0,99.0,99.0,99.0, 99.0,99.0,99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,99.0,-99.0,-99.0,-99.0,-99.0,-99.,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,99.0,99.0,-99.0,-99.0,-99.0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,"none",0.,0.,99.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-99.}


         num = num + 1
      endif
   endif
endfor

print, "object ends with ", num
plothist, object.tmasskflux, title = 'after g matching'
;----------------------------------------------------------
;match J to object


; create initial arrays
m=n_elements(flamjobject.flamjra)
ir=n_elements(object.ra)

irmatch=fltarr(ir)
mmatch=fltarr(m)
irmatch[*]=-999
mmatch[*]=-999


print,'Matching flam J to object'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,m-1 do begin

   dist=sphdist( flamjobject[q].flamjra, flamjobject[q].flamjdec,object.ra,object.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin
      mmatch[q]=ind
      object[ind].flamjra = flamjobject[q].flamjra
      object[ind].flamjdec = flamjobject[q].flamjdec
      object[ind].flamjxcenter = flamjobject[q].flamjxcenter
      object[ind].flamjycenter = flamjobject[q].flamjycenter
      object[ind].flamjmag = flamjobject[q].flamjmag
      object[ind].flamjmagerr = flamjobject[q].flamjmagerr
      object[ind].flamjfwhm = flamjobject[q].flamjfwhm
      object[ind].flamjisoarea = flamjobject[q].flamjisoarea
      object[ind].flamjellip = flamjobject[q].flamjellip

;      print, "matched a J", jobject[q].jra, object[ind].ra, q
   endif 
endfor

print,"Finished at "+systime()



;---------------------------------------------------------
;match wirc J and H and K to object

; create initial arrays
m=n_elements(wircjobject.wircjra)
ir=n_elements(object.ra)

irmatch=fltarr(ir)
mmatch=fltarr(m)
irmatch[*]=-999
mmatch[*]=-999


print,'Matching wirc J to object'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,m-1 do begin

   dist=sphdist( wircjobject[q].wircjra, wircjobject[q].wircjdec,object.ra,object.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin
      mmatch[q]=ind
      object[ind].wircjra = wircjobject[q].wircjra
      object[ind].wircjdec = wircjobject[q].wircjdec
      object[ind].wircjxcenter = wircjobject[q].wircjxcenter
      object[ind].wircjycenter = wircjobject[q].wircjycenter
      object[ind].wircjmag = wircjobject[q].wircjmag
      object[ind].wircjmagerr = wircjobject[q].wircjmagerr
      object[ind].wircjfwhm = wircjobject[q].wircjfwhm
      object[ind].wircjisoarea = wircjobject[q].wircjisoarea
      object[ind].wircjellip = wircjobject[q].wircjellip

 ;    print, "matched a J", wircjobject[q].wircjra, object[ind].ra, q
   endif else begin
      print, "did not match a J", wircjobject[q].wircjra,wircjobject[q].wircjdec
   endelse

endfor
print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))

print,"Finished at "+systime()

;----------
; create initial arrays
m=n_elements(wirchobject.wirchra)
ir=n_elements(object.ra)

irmatch=fltarr(ir)
mmatch=fltarr(m)
irmatch[*]=-999
mmatch[*]=-999


print,'Matching wirch to object'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,m-1 do begin

   dist=sphdist( wirchobject[q].wirchra, wirchobject[q].wirchdec,object.ra,object.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin
      mmatch[q]=ind
      object[ind].wirchra = wirchobject[q].wirchra
      object[ind].wirchdec = wirchobject[q].wirchdec
      object[ind].wirchxcenter = wirchobject[q].wirchxcenter
      object[ind].wirchycenter = wirchobject[q].wirchycenter
      object[ind].wirchmag = wirchobject[q].wirchmag
      object[ind].wirchmagerr = wirchobject[q].wirchmagerr
      object[ind].wirchfwhm = wirchobject[q].wirchfwhm
      object[ind].wirchisoarea = wirchobject[q].wirchisoarea
      object[ind].wirchellip = wirchobject[q].wirchellip

;     print, "matched a H", wirchobject[q].wirchra, object[ind].ra, q
   endif 
endfor
print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))

print,"Finished at "+systime()
;----------
; create initial arrays
m=n_elements(wirckobject.wirckra)
ir=n_elements(object.ra)

irmatch=fltarr(ir)
mmatch=fltarr(m)
irmatch[*]=-999
mmatch[*]=-999


print,'Matching wirck to object'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,m-1 do begin

   dist=sphdist( wirckobject[q].wirckra, wirckobject[q].wirckdec,object.ra,object.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin
      mmatch[q]=ind
      object[ind].wirckra = wirckobject[q].wirckra
      object[ind].wirckdec = wirckobject[q].wirckdec
      object[ind].wirckxcenter = wirckobject[q].wirckxcenter
      object[ind].wirckycenter = wirckobject[q].wirckycenter
      object[ind].wirckmag = wirckobject[q].wirckmag
      object[ind].wirckmagerr = wirckobject[q].wirckmagerr
      object[ind].wirckfwhm = wirckobject[q].wirckfwhm
      object[ind].wirckisoarea = wirckobject[q].wirckisoarea
      object[ind].wirckellip = wirckobject[q].wirckellip

;     print, "matched a K", wirckobject[q].wirckra, object[ind].ra, q
   endif 
endfor
print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))

print,"Finished at "+systime()

;---------------------------------------------------------

;match mips70
; create initial arrays
m=n_elements(mips70object.mips70ra)
ir=n_elements(object.ra)

irmatch=fltarr(ir)
mmatch=fltarr(m)
irmatch[*]=-999
mmatch[*]=-999


print,'Matching mips70 to object'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,m-1 do begin

   dist=sphdist( mips70object[q].mips70ra, mips70object[q].mips70dec,object.ra,object.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin
      mmatch[q]=ind
      object[ind].mips70ra = mips70object[q].mips70ra
      object[ind].mips70dec = mips70object[q].mips70dec
      object[ind].mips70xcenter = mips70object[q].mips70xcenter
      object[ind].mips70ycenter = mips70object[q].mips70ycenter
      object[ind].mips70mag = mips70object[q].mips70mag
      object[ind].mips70magerr = mips70object[q].mips70magerr
      object[ind].mips70fwhm = mips70object[q].mips70fwhm
      object[ind].mips70isoarea = mips70object[q].mips70isoarea
      object[ind].mips70flags = mips70object[q].mips70flags
      object[ind].mips70bckgrnd = mips70object[q].mips70bckgrnd

   endif 
endfor
print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))

print,"Finished at "+systime()

;---------------------------------------------------------

; create initial arrays
m=n_elements(acsobject.acsmag)
ir=n_elements(object.ra)

irmatch=fltarr(ir)
mmatch=fltarr(m)
irmatch[*]=-999
mmatch[*]=-999


print,'Matching morpheus to object'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

print, "m", m
for q=long(0),m-1 do begin
   dist=sphdist( acsobject[q].acsra, acsobject[q].acsdec,object.ra,object.dec,/degrees)
;   dist=sphdist( ra(q), dec(q),object.ra,object.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0006)  then begin
      mmatch[q]=ind

      object[ind].acsra=acsobject(q).acsra
      object[ind].acsdec=acsobject(q).acsdec
      object[ind].acsisoarea=acsobject(q).acsisoarea
      object[ind].acsellip= acsobject(q).acsellip
      object[ind].acstheta=acsobject(q).acstheta
      object[ind].acssnr=acsobject(q).acssnr
      object[ind].acsmag=acsobject(q).acsmag
      object[ind].acsmagerr=acsobject(q).acsmagerr
      object[ind].acsflux=acsobject(q).acsflux
      object[ind].acssky=acsobject(q).acssky
      object[ind].acsskyrms=acsobject(q).acsskyrms
      object[ind].acsmu=acsobject(q).acsmu
      object[ind].acscentralmu=acsobject(q).acscentralmu
      object[ind].acsgini=acsobject(q).acsgini
      object[ind].acsassym=acsobject(q).acsassym
      object[ind].acshlightrad=acsobject(q).acshlightrad
      object[ind].acsconcentration=acsobject(q).acsconcentration
      object[ind].acsm20=acsobject(q).acsm20
      object[ind].acsclassstar =acsobject(q).acsclassstar
      object[ind].acsfwhm = acsobject(q).acsfwhm
      object[ind].acsfluxmax = acsobject(q).acsfluxmax
   endif
endfor

matched=where(mmatch GT 0)
nonmatched = where(mmatch lt 0)

print, n_elements(matched),"matched"
print, n_elements(nonmatched),"nonmatched"



;catalog  the nonmatched objects 
num = num - 1.
print, "object stars with", num
for num2 = 0.D, n_elements(nonmatched) -1 do begin
;   print,"inside num,num2", num,num2
   object[num] = {ob, acsobject(nonmatched[num2]).acsra, acsobject(nonmatched[num2]).acsdec,-99.,-99.,-99.,99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,99.,99.,99.,99.,99.,99.,99.,-99.,-99.,-99.,-99.,-99.,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,99.0,99.0,-99.0,-99.0,-99.0,-99.,-99.,-99.,99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,-99.,"none",  acsobject(nonmatched[num2]).acsconcentration, acsobject(nonmatched[num2]).acsgini,acsobject(nonmatched[num2]).acsmag,acsobject(nonmatched[num2]).acsmu,acsobject(nonmatched[num2]).acsellip,acsobject(nonmatched[num2]).acsflag,acsobject(nonmatched[num2]).acsra,acsobject(nonmatched[num2]).acsdec,acsobject(nonmatched[num2]).acsisoarea,acsobject(nonmatched[num2]).acstheta,acsobject(nonmatched[num2]).acssnr,acsobject(nonmatched[num2]).acsflux,acsobject(nonmatched[num2]).acssky,acsobject(nonmatched[num2]).acsskyrms,acsobject(nonmatched[num2]).acscentralmu,acsobject(nonmatched[num2]).acsassym,acsobject(nonmatched[num2]).acshlightrad,acsobject(nonmatched[num2]).acsm20,acsobject(nonmatched[num2]).acsclassstar, acsobject(nonmatched[num2]).acsmagerr, acsobject(nonmatched[num2]).acsfwhm,acsobject(nonmatched[num2]).acsfluxmax}

   num = num + 1
endfor

print,"Finished at "+systime()

;--------------------------------------------------------


object = object[0:num - 1]


for newcount = long(0), n_elements(object.mips24flux) -1 do begin
   if object[newcount].mips24flux eq 0 then begin
      object[newcount].mips24flux = -99.0
      object[newcount].mips24fluxerr = -99.0
   endif

endfor

plothist, iracobject.tmasskflux, title = 'at the end'
save, object, filename='/Users/jkrick/idlbin/object.test.sav'


close, outlun
free_lun, outlun
close, outlun2
free_lun, outlun2
close, outlun3
free_lun, outlun3
close, outlun4
free_lun, outlun4


undefine, robject
undefine, gobject
undefine, object
undefine, iobject
undefine, uobject
undefine, infraredobject
undefine, mips24object
undefine, flamjobject
undefine, wircjobject
undefine, wirchobject
undefine, wirckobject
end


;ftab_ext, '/Users/jkrick/hst/total_catalog.fits', 'Concentration,Gini,mag_auto,SurfaceBrightness,x_image,y_image,class_star,a_image,b_image', c,g,mag,mu,x,y,s,a,b
;ftab_ext, '/Users/jkrick/hst/total_catalog.fits', 'flags,x_world,y_world,isoarea_image,theta_world,signaltonoiseratio,fluxcounts,sky,skyrms',flag,ra,dec,isoarea,theta,snr,flux, sky,skyrms
;ftab_ext, '/Users/jkrick/hst/total_catalog.fits', 'centralsurfacebrightness,rotationalasymmetry,halflightradiusinpixels,m20,class_star,magerr_auto',cmu,asym,hlr,m20,class,magerr

;for p = 0l, n_elements(c) -1 do begin
;acsobject[p] ={acsob,c(p),g(p),mag(p), mu(p),1-b(p)/a(p),flag(p),ra(p),dec(p),isoarea(p),theta(p),snr(p),flux(p), sky(p),skyrms(p),cmu(p),asym(p),hlr(p),m20(p),class(p),magerr(p)}
;endfor
;acsobject = acsobject[0:p-1]
;print, "there are ",p," acs objects"
