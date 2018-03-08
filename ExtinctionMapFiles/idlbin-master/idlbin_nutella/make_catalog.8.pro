pro make_catalog


close, /all

!P.multi = [0,2,2]

numoobjects = 100000
colorarr = fltarr(numoobjects)
magarr = fltarr(numoobjects)
FITS_READ, '/Users/jkrick/palomar/lfc/coadd_g.fits',gdata, gheader
FITS_READ, '/Users/jkrick/spitzer/irac/ch1/mosaic.fits',irdata, irheader

uobject = replicate({uob, ura:0D,udec:0D,uflux:0D,umag:0D,umagerr:0D,ufwhm:0D,uisoarea:0D,uellip:0D,ufluxmax:0D, uclassstar:0D, uflags:0D, unndist:0D, umatchdist:0D},numoobjects)
gobject = replicate({gob, gra:0D,gdec:0D,gflux:0D,gmag:0D,gmagerr:0D,gfwhm:0D,gisoarea:0D,gellip:0D,gfluxmax:0D,gmatch:0D, gclassstar:0D, gflags:0D, gnndist:0D, gmatchdist:0D},numoobjects)
robject = replicate({rob, rra:0D,rdec:0D,rflux:0D,rmag:0D,rmagerr:0D,rfwhm:0D,risoarea:0D, rellip:0D,rfluxmax:0D, rclassstar:0D,rflags:0D, rnndist:0D, rmatchdist:0D},numoobjects)
iobject = replicate({iob, ira:0D,idec:0D,iflux:0D,imag:0D,imagerr:0D,ifwhm:0D,iisoarea:0D, iellip:0D,ifluxmax:0D, iclassstar:0D, iflags:0D, inndist:0D, imatchdist:0D},numoobjects)

zobject = replicate({zob, zra:0D,zdec:0D,zfluxauto:0D,zmagauto:0D,zmagerrauto:0D,zfluxaper:0D,zmagaper:0D,zmagerraper:0D,zfluxbest:0D,zmagbest:0D,zmagerrbest:0D,zfwhm:0D,zisoarea:0D, zellip:0D,zfluxmax:0D, zclassstar:0D, zflags:0D, znndist:0D, zmatchdist:0D,zmatch:0D},numoobjects)

acsobject =replicate({acsob,acsconcentration:0D,acsgini:0D,acsmag:0D, acsmu:0D,acsellip:0D,acsflags:0D,acsra:0D,acsdec:0D,acsisoarea:0D,acstheta:0D,acssnr:0D,acsflux:0D, acssky:0D,acsskyrms:0D,acscentralmu:0D,acsassym:0D,acshlightrad:0D,acsm20:0D,acsclassstar:0D,acsmagerr:0D,acsfwhm:0D, acsfluxmax:0D, acsnndist:0D, acsmatchdist:0D, acsmatch:0D},80000)
;-------
flamjobject = replicate({flamjob,flamjxcenter:0D, flamjycenter:0D, flamjra:0D, flamjdec:0D, flamjmag:0D, flamjmagerr:0D, flamjfwhm:0D, flamjisoarea:0D, flamjellip:0D, flamjnndist:0D, flamjmatchdist:0D, flamjmatch:0D },numoobjects)
wircjobject = replicate({wircjob, wircjxcenter:0D, wircjycenter:0D, wircjra:0D, wircjdec:0D, wircjmag:0D, wircjmagerr:0D, wircjfwhm:0D, wircjisoarea:0D, wircjellip:0D, wircjnndist:0D, wircjmatchdist:0D, wircjmatch:0D },numoobjects)
wirchobject = replicate({wirchob, wirchxcenter:0D, wirchycenter:0D, wirchra:0D, wirchdec:0D, wirchmag:0D, wirchmagerr:0D, wirchfwhm:0D, wirchisoarea:0D, wirchellip:0D , wirchnndist:0D, wirchmatchdist:0D, wirchmatch:0D},numoobjects)
wirckobject = replicate({wirckob, wirckxcenter:0D, wirckycenter:0D, wirckra:0D, wirckdec:0D, wirckmag:0D, wirckmagerr:0D, wirckfwhm:0D, wirckisoarea:0D, wirckellip:0D, wircknndist:0D, wirckmatchdist:0D , wirckmatch:0D},numoobjects)

;nearirobject = replicate({nearirob, nearirra:0D, nearirdec:0D, INHERITS flamjob, INHERITS wircjob, INHERITS wirchob, INHERITS wirckob}, numoobjects)

mips24object = replicate({mips24ob, mips24xcenter:0D, mips24ycenter:0D, mips24ra:0D, mips24dec:0D, mips24flux:0D, mips24fluxerr:0D, mips24mag:0D, mips24magerr:0D, mips24bckgrnd:0D, iter:0D,sharp:0D, mips24match:0D, mips24nndist:0D, mips24matchdist:0D},numoobjects)

mips70object = replicate({mips70ob, mips70xcenter:0D, mips70ycenter:0D, mips70ra:0D, mips70dec:0D, mips70flux:0D, mips70fluxerr:0D, mips70mag:0D, mips70magerr:0D, mips70bckgrnd:0D, mips70isoarea:0D,mips70fwhm:0D, mips70flags:0D, mips70nndist:0D, mips70matchdist:0D, mips70match:0D},numoobjects)

iracobject = replicate({iracob, iracra:0D,iracdec:0D,tmassJflux:0D,tmasshflux:0D,tmasskflux:0D,irac1flux:0D,irac2flux:0D,irac3flux:0D, irac4flux:0D,tmassJmag:0D,tmasshmag:0D,tmasskmag:0D,irac1mag:0D,irac2mag:0D,irac3mag:0D, irac4mag:0D,mips24:0D, iracnndist:0D, iracmatchdist:0D, iracmatch:0D},numoobjects)

infraredobject = replicate({spitzob, INHERITS iracob,INHERITS mips24ob,INHERITS mips70ob},numoobjects)
;=======
object= replicate({ob, ra:0D, dec:0D, INHERITS uob, INHERITS gob, INHERITS rob, INHERITS iob, INHERITS zob, INHERITS flamjob, INHERITS wircjob, INHERITS wirchob, INHERITS wirckob, INHERITS spitzob, zphot:0D, chisq:0D, prob:0D, spt:0D, age:0D, av:0D, Mabs:0D, mass:0D, plusmasserr:0D, minusmasserr:0D, masschi:0D, massprob:0D, massage:0D, model:' ',INHERITS acsob}, numoobjects)

;----------------------------------------------------------------
;template objects
;object= {ob, ra, dec, ura,udec,uflux,umag,umagerr,ufwhm,uisoarea,uellip,ufluxmax, uclassstar, uflags , unndist, umatchdist, gra,gdec,gflux,gmag,gmagerr,gfwhm,gisoarea,gellip,gfluxmax,gmatch, gclassstar, gflags, gnndist, gmatchdist,  rra,rdec,rflux,rmag,rmagerr,rfwhm,risoarea, rellip,rfluxmax, rclassstar,rflags, rnndist, rmatchdist, ira,idec,iflux,imag,imagerr,ifwhm,iisoarea, iellip,ifluxmax, iclassstar, iflags, inndist, imatchdist, zra,zdec,zfluxauto,zmagauto,zmagerrauto,zfluxaper,zmagaper,zmagerraper,zfluxbest,zmagbest,zmagerrbest,zfwhm,zisoarea, zellip,zfluxmax, zclassstar, zflags, znndist, zmatchdist,zmatch,flamjxcenter, flamjycenter, flamjra, flamjdec, flamjmag, flamjmagerr, flamjfwhm, flamjisoarea, flamjellip, flamjnndist, flamjmatchdist,flamjmatch , wircjxcenter, wircjycenter, wircjra, wircjdec, wircjmag, wircjmagerr, wircjfwhm, wircjisoarea, wircjellip, wircjnndist, wircjmatchdist,wircjmatch,wirchxcenter, wirchycenter, wirchra, wirchdec, wirchmag, wirchmagerr, wirchfwhm, wirchisoarea, wirchellip, wirchnndist, wirchmatchdist, wirchmatch,wirckxcenter, wirckycenter, wirckra, wirckdec, wirckmag, wirckmagerr, wirckfwhm, wirckisoarea, wirckellip, wircknndist, wirckmatchdist, wirckmatch, iracra,iracdec,tmassJflux,tmasshflux,tmasskflux,irac1flux,irac2flux,irac3flux, irac4flux,tmassJmag,tmasshmag,tmasskmag,irac1mag,irac2mag,irac3mag, irac4mag,mips24, iracnndist, iracmatchdist,iracmatch, mips24xcenter, mips24ycenter, mips24ra, mips24dec, mips24flux, mips24fluxerr, mips24mag, mips24magerr, mips24bckgrnd, iter,sharp, mips24match, mips24nndist, mips24matchdist,mips70xcenter, mips70ycenter, mips70ra, mips70dec, mips70flux, mips70fluxerr, mips70mag, mips70magerr, mips70bckgrnd, mips70isoarea,mips70fwhm, mips70flags, mips70nndist, mips70matchdist, mips70match,zphot, chisq, prob, spt, age, av, Mabs, mass, plusmasserr, minusmasserr, masschi, massprob, massage, model,acsconcentration,acsgini,acsmag, acsmu,acsellip,acsflags,acsra,acsdec,acsisoarea,acstheta,acssnr,acsflux, acssky,acsskyrms,acscentralmu,acsassym,acshlightrad,acsm20,acsclassstar,acsmagerr,acsfwhm, acsfluxmax, acsnndist, acsmatchdist,acsmatch}

;----
;object= {ob, 0., 0., u:0.,0.,-99.,99.,99.,0.,0.,0.,-99., -99., -99. , 0., 0.,g: 0.,0.,-99.,99.,99.,0.,0.,0.,-99.,-99., -99., -99., 0., 0., r: 0.,0.,-99.,99.,99.,0.,0, 0.,-99., -99.,-99., 0., 0., i:0.,0.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0., z:0.,0.,-99.,99.,99.,-99.,99.,99.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0.,-99.,flamj:0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99. , wircj:0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99.,wirch:0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99.,wirck:0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99., irac:0.,0.,-99.,-99.,-99.,-99.,-99.,-99., -99.,99.,99.,99.,99.,99.,99., 99.,-99., 0., 0.,-99.,mips24: 0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0.,mips70:0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0., -99.,zphot:0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., "none",acs:0.,0.,99., 0.,0.,-99.,0.,0.,0.,0.,0.,-99., 0.,0.,0.,0.,0.,0.,-99.,99.,0., -99., 0., 0.,-99.}

;----------------------------------------------------------------
;u catalog has the smallest number of junk obs at the end of the catalog, so use it for size.
nu = 0.
openr, lung, "/Users/jkrick/palomar/LFC/catalog/SExtractor.u.cat", /get_lun
WHILE (NOT EOF(lung)) DO BEGIN
    READF, lung, o, xwcsg, ywcsg, xcenterg, ycenterg, fluxautog, magautog, magerrg, magaper, magapererr, fwhmg, isoareag, fluxmaxg, ellipu,gnum, class_star, flags

   if (magaper gt 21.0)  then begin
       magautog = magaper
       magerrg = magapererr
    endif

   if gnum ne 0 then begin
       uobject[gnum-1] ={uob, xwcsg, ywcsg,fluxautog, magautog, magerrg,fwhmg, isoareag, ellipu,fluxmaxg, class_star, flags,0.,0.}
       nu = nu + 1
    ENDIF

 endwhile
uobject = uobject[0:nu - 1]
print, "there are ",nu," u objects"
close, lung
free_lun, lung

;determine nndist
for nn = 0l, n_elements(uobject.ura) - 1 do begin
   a = where(sqrt( (uobject[nn].ura - uobject.ura)^2 + (uobject[nn].udec - uobject.udec)^2 ) le 0.0008 )
   uobject[nn].unndist = n_elements(a)
endfor


;----------------------------------------------------------------

lg = 0.
openr, lung, "/Users/jkrick/palomar/LFC/catalog/SExtractor.g.cat", /get_lun
WHILE (NOT EOF(lung)) DO BEGIN
    READF, lung, o, xwcsg, ywcsg, xcenterg, ycenterg, fluxautog, magautog, magerrg, magaper, magapererr,fwhmg, isoareag, fluxmaxg, ellipg, class_star, flags
   if (magaper gt 21.0)  then begin
       magautog = magaper
       magerrg = magapererr
    endif

       gobject[lg] ={gob, xwcsg, ywcsg,fluxautog, magautog, magerrg,fwhmg, isoareag, ellipg, fluxmaxg,-1, class_star, flags,0.,0.}
       lg = lg + 1

 endwhile
gobject = gobject[0:nu - 1]
print, "there are ",lg," g objects"
close, lung
free_lun, lung

 ;determine nndist
for nn = 0l, n_elements(gobject.gra) - 1 do begin
   a = where(sqrt( (gobject[nn].gra - gobject.gra)^2 + (gobject[nn].gdec - gobject.gdec)^2 ) le 0.0008 )
   gobject[nn].gnndist = n_elements(a)
endfor  
;----------------------------------------------------------------
i = 0.
openr, lunr, "/Users/jkrick/palomar/LFC/catalog/SExtractor.r.cat", /get_lun
WHILE (NOT EOF(lunr)) DO BEGIN
    READF, lunr, o, xwcsr, ywcsr, xcenterr, ycenterr, fluxautor, magautor, magautoerrr, magaper, magapererr,fwhmr, isoarear, fluxmaxr, ellipr,gnum, class_star, flags

    if (magaper gt 21.0)  then begin
       magautor = magaper
       magautoerrr = magapererr
    endif

    if gnum ne 0 then begin
       robject[gnum-1] ={rob, xwcsr, ywcsr,fluxautor, magautor, magautoerrr,fwhmr, isoarear, ellipr,fluxmaxr, class_star, flags,0.,0.}
       i = i + 1
    endif


 endwhile
robject = robject[0:nu-1]
print, "there are ",i," r objects"
close, lunr
free_lun, lunr

;determine nndist
for nn = 0l, n_elements(robject.rra) - 1 do begin
   a = where(sqrt( (robject[nn].rra - robject.rra)^2 + (robject[nn].rdec - robject.rdec)^2 ) le 0.0008 )
   robject[nn].rnndist = n_elements(a)
endfor
;--------------------------------------------------------
;am changing the zeropoints here, after the catalogs are made.

m = 0.
openr, lung, "/Users/jkrick/palomar/LFC/catalog/SExtractor.i.cat", /get_lun
WHILE (NOT EOF(lung)) DO BEGIN
    READF, lung, o, xwcsg, ywcsg, xcenterg, ycenterg, fluxautog, magautog, magerrg, magaper, magapererr, fwhmg, isoareag, fluxmaxg, ellipi,gnum, class_star, flags

   if (magaper gt 21.0)  then begin
       magautog = magaper
       magerrg = magapererr
    endif

   if gnum ne 0 then begin
       iobject[gnum-1] ={iob, xwcsg, ywcsg,fluxautog, magautog, magerrg,fwhmg, isoareag, ellipi,fluxmaxg, class_star, flags,0.,0.}
       m = m + 1
    ENDIF

 endwhile
iobject = iobject[0:nu - 1]
print, "there are ",m," i objects"
close, lung
free_lun, lung

;determine nndist
for nn = 0l, n_elements(iobject.ira) - 1 do begin
   a = where(sqrt( (iobject[nn].ira - iobject.ira)^2 + (iobject[nn].idec - iobject.idec)^2 ) le 0.0008 )
   iobject[nn].inndist = n_elements(a)
endfor
;----------------------------------------------------------------

p = 0.
openr, luni, "/Users/jkrick/mmt/SExtractor.cat", /get_lun

fits_read, '/Users/jkrick/mmt/IRACCF.z.coadd.cov.fits', zdata, zhead

WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, NUMBER,X_WORLD      ,Y_WORLD      ,X_IMAGE,Y_IMAGE,FLUX_AUTO,MAG_AUTO,MAGERR_AUTO,FLUX_APER,MAG_APER,MAGERR_APER,FLUX_BEST,MAG_BEST,MAGERR_BEST,FWHM_IMAGE,ISOAREA_IMAGE,FLUX_MAX,ELLIPTICITY,CLASS_STAR,FLAGS
   if zdata[x_image, y_image] gt 1 then begin
      zobject[p] ={zob, x_world, y_world, flux_auto, MAG_AUTO,MAGERR_AUTO,FLUX_APER,MAG_APER,MAGERR_APER,FLUX_BEST,MAG_BEST,MAGERR_BEST,FWHM_IMAGE,ISOAREA_IMAGE,ELLIPTICITY,flux_max,CLASS_STAR,FLAGS,0.,0.,-1.}
      p = p + 1
   endif

 endwhile
zobject =zobject[0:p - 1]
print, "there are ",p,"z objects"
close, luni
free_lun, luni

;determine nndist
for nn = 0l, n_elements(zobject.zra) - 1 do begin
   a = where(sqrt( (zobject[nn].zra - zobject.zra)^2 + (zobject[nn].zdec - zobject.zdec)^2 ) le 0.0008 )
   zobject[nn].znndist = n_elements(a)
endfor
;----------------------------------------------------------------
;ACS
mag = 0.D
p = 0l
openr, luni, "/Users/jkrick/hst/raw/wholeacs.cat", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   readf, luni, NUMBER,X_WORLD,Y_WORLD,X_IMAGE,Y_IMAGE,MAG_AUTO,FLUX_BEST,FLUXERR_BEST,MAG_BEST,MAGERR_BEST,BACKGROUND,FLUX_MAX,ISOAREA_IMAGE,ALPHA_J2000,DELTA_J2000,THETA_IMAGE,MU_THRESHOLD,MU_MAX,FLAGS,FWHM_IMAGE,CLASS_STAR,ELLIPTICITY

   if magerr_best eq 0. then magerr_best = 0.03
   acsobject[p]={acsob,0.,0., mag_best, 0., ellipticity, flags, x_world, y_world, isoarea_image, theta_image, 0., flux_best, background, 0., 0., 0., 0., 0., class_star, magerr_best,fwhm_image,flux_max,0.,0.,-1.}
   p = p + 1
   
 endwhile
acsobject = acsobject[0:p-1]
print, "there are ",p," acs objects"

close, luni
free_lun, luni

;determine nndist
for nn = 0l, n_elements(acsobject.acsra) - 1 do begin
   a = where(sqrt( (acsobject[nn].acsra - acsobject.acsra)^2 + (acsobject[nn].acsdec - acsobject.acsdec)^2 ) le 0.0008 )
   acsobject[nn].acsnndist = n_elements(a)
endfor
;----------------------------------------------------------------
;irac, 2mass, mips24(not best)
readcol, "/Users/jkrick/spitzer/irac/combined_catalog.prt", wcsra, wcsdec, xcenter, ycenter, tmassJa,tmassha,tmasska,irac1a,irac2a,irac3a, irac4a,mips24a

for p = 0, n_elements(wcsra) - 1 do begin

 if tmassja(p) le 0. then begin
      tmassja(p) = -99
      tmassjab = -99
   endif else begin
      tmassjab = -2.5*alog10(tmassja(p)) + 8.926 
   endelse

  if tmassha(p)  le 0. then begin
      tmassha(p) = -99
      tmasshab = -99
   endif else begin
      tmasshab = -2.5*alog10(tmassha(p)) + 8.926 
   endelse

  if tmasska(p) le 0. then begin
      tmasska(p) = -99
      tmasskab = -99
   endif else begin
      tmasskab = -2.5*alog10(tmasska(p)) + 8.926 
   endelse

   if irac1a(p) lt 0. then begin
      mag1 = -99. 
      irac1a(p)=-99.
   endif else begin
      mag1 = 8.926 - 2.5*alog10(1E-6*irac1a(p)) 
   endelse

   if irac2a(p) lt 0. then begin
      mag2 = -99. 
      irac2a(p)=-99.
   endif else begin
      mag2 = 8.926 - 2.5*alog10(1E-6*irac2a(p)) 
   endelse

   if irac3a(p) lt 0. then begin
      mag3 = -99. 
      irac3a(p)=-99.
   endif else begin
      mag3 = 8.926- 2.5*alog10(1E-6*irac3a(p))
   endelse

   if irac4a(p) lt 0. then begin
      mag4 = -99. 
      irac4a(p)=-99.
   endif else begin
      mag4 = 8.926 - 2.5*alog10(1E-6*irac4a(p))
   endelse


   iracobject[p] ={iracob, wcsra(p), wcsdec(p), tmassJa(p),tmassha(p),tmasska(p),irac1a(p),irac2a(p),irac3a(p), irac4a(p),tmassjab,tmasshab,tmasskab,mag1,mag2,mag3,mag4,mips24a(p),0.,0.,-1.}
   
endfor

print, "there are ",p," ir objects, plothist"
iracobject = iracobject[0:p-1]
;plothist, iracobject.tmasskflux, title = 'after reading file', xrange=[0,1], bin = 0.001, yrange=[0,100]

;determine nndist
for nn = 0l, n_elements(iracobject.iracra) - 1 do begin
   a = where(sqrt( (iracobject[nn].iracra - iracobject.iracra)^2 + (iracobject[nn].iracdec - iracobject.iracdec)^2 ) le 0.0008 )
   iracobject[nn].iracnndist = n_elements(a)
endfor
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

openr, lunw, "/Users/jkrick/Spitzer/mips/mips24/dao/mips24.ra.flux.phot", /get_lun
WHILE (NOT EOF(lunw)) DO BEGIN
;   print, w
   READF, lunw, num, ra, dec, xcen,ycen, mag, magerr, back, niter, sh, junk, junk, flux, fluxerr
;   print, xcen
   ;need to convert from magnitudes back to fluxes.
   mips24object[w] = {mips24ob, xcen, ycen, 0.0,0.0, 1.432*(10^((25-mag)/2.5)),  10^((25-(mag -magerr))/2.5) -(10^((25-mag)/2.5)) , mag,magerr,back, niter,sh,0.,0.,0.}
   w = w +1
endwhile
mips24object = mips24object[0:w - 1]
print, "there are ",w," mips24 objects"
close, lunw
free_lun, lunw

;need to go from x,y in pixels to ra, dec
FITS_READ, '/Users/jkrick/spitzer/mips/mips24/Combine/mosaic.wcs.fits',data, header
xyad, header, mips24object.mips24xcenter, mips24object.mips24ycenter, ra, dec
mips24object.mips24ra = ra
mips24object.mips24dec = dec

;determine nndist
for nn = 0l, n_elements(mips24object.mips24ra) - 1 do begin
   a = where(sqrt( (mips24object[nn].mips24ra - mips24object.mips24ra)^2 + (mips24object[nn].mips24dec - mips24object.mips24dec)^2 ) le 0.0008 )
   mips24object[nn].mips24nndist = n_elements(a)
endfor
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

   mips70object[p] ={mips70ob, x,y, ra, dec, fluxauto*23.5045, fluxerrauto*23.5045, fluxauto, fluxerrauto,back, isoarea, fwhm, flag,0.,0.,-1.}

   p = p + 1
   
 endwhile
mips70object = mips70object[0:p-1]
print, "there are ",p," mips70 objects"

close, luni
free_lun, luni

;determine nndist
for nn = 0l, n_elements(mips70object.mips70ra) - 1 do begin
   a = where(sqrt( (mips70object[nn].mips70ra - mips70object.mips70ra)^2 + (mips70object[nn].mips70dec - mips70object.mips70dec)^2 ) le 0.0008 )
   mips70object[nn].mips70nndist = n_elements(a)
endfor
;---------------------------------------------------------------

p = 0
openr, luni, "/Users/jkrick/palomar/wirc/jband_wirc_nohead.txt", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, junk,x_image, y_image,junk,junk,junk,junk,junk,junk,mag_best, magerr_best, junk,junk, junk,isoareaimage, ra,dec,junk,junk,junk, junk,fwhm_image, junk,ell, junk 

   wircjobject[p] ={wircjob, x_image, y_image, ra, dec, mag_best, magerr_best, fwhm_image, isoareaimage, ell,0.,0.,-1.}
   p = p + 1
   
 endwhile
wircjobject = wircjobject[0:p - 1]
print, "there are ",p," wirc j objects"
close, luni
free_lun, luni

;determine nndist
for nn = 0l, n_elements(wircjobject.wircjra) - 1 do begin
   a = where(sqrt( (wircjobject[nn].wircjra - wircjobject.wircjra)^2 + (wircjobject[nn].wircjdec - wircjobject.wircjdec)^2 ) le 0.0008 )
   wircjobject[nn].wircjnndist = n_elements(a)
endfor
;----------------------------------------------------------------

p = 0
openr, luni, "/Users/jkrick/palomar/wirc/hband_wirc_nohead.txt", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, junk,x_image, y_image,junk,junk,junk,junk,junk,junk,mag_best, magerr_best, junk,junk, junk,isoareaimage, ra,dec,junk,junk,junk, junk,fwhm_image, junk,ell, junk 

   wirchobject[p] ={wirchob, x_image, y_image, ra, dec, mag_best, magerr_best, fwhm_image, isoareaimage, ell,0.,0.,-1.}
   p = p + 1
   
 endwhile
wirchobject = wirchobject[0:p - 1]
print, "there are ",p," wirc h objects"
close, luni
free_lun, luni

;determine nndist
for nn = 0l, n_elements(wirchobject.wirchra) - 1 do begin
   a = where(sqrt( (wirchobject[nn].wirchra - wirchobject.wirchra)^2 + (wirchobject[nn].wirchdec - wirchobject.wirchdec)^2 ) le 0.0008 )
   wirchobject[nn].wirchnndist = n_elements(a)
endfor
 ;----------------------------------------------------------------

p = 0
openr, luni, "/Users/jkrick/palomar/wirc/kband_wirc_nohead.txt", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, junk,x_image, y_image,junk,junk,junk,junk,junk,junk,mag_best, magerr_best, junk,junk, junk,isoareaimage, ra,dec,junk,junk,junk, junk,fwhm_image, junk,ell, junk 

   wirckobject[p] ={wirckob, x_image, y_image, ra, dec, mag_best, magerr_best, fwhm_image, isoareaimage, ell,0.,0.,-1.}
   p = p + 1
   
 endwhile
wirckobject = wirckobject[0:p - 1]
print, "there are ",p," wirc k objects"
close, luni
free_lun, luni
;determine nndist
for nn = 0l, n_elements(wirckobject.wirckra) - 1 do begin
   a = where(sqrt( (wirckobject[nn].wirckra - wirckobject.wirckra)^2 + (wirckobject[nn].wirckdec - wirckobject.wirckdec)^2 ) le 0.0008 )
   wirckobject[nn].wircknndist = n_elements(a)
endfor
;----------------------------------------------------------------
p = 0.
openr, luni, "/Users/jkrick/nep/flamingos/jband_noheader.txt", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, junk,x_image, y_image,junk,junk,junk,junk,junk,junk,mag_best, magerr_best, junk,junk, junk,isoareaimage, ra,dec,junk,junk,junk, junk,fwhm_image, junk,ell, junk 

   flamjobject[p] ={flamjob, x_image, y_image, ra, dec, mag_best, magerr_best, fwhm_image, isoareaimage, ell,0.,0.,-1.}
   p = p + 1
   
 endwhile
flamjobject = flamjobject[0:p - 1]
print, "there are ",p,"flam  j objects"
close, luni
free_lun, luni
;determine nndist
for nn = 0l, n_elements(flamjobject.flamjra) - 1 do begin
   a = where(sqrt( (flamjobject[nn].flamjra - flamjobject.flamjra)^2 + (flamjobject[nn].flamjdec - flamjobject.flamjdec)^2 ) le 0.0008 )
   flamjobject[nn].flamjnndist = n_elements(a)
endfor
;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
;Matching--------------------------------------------------------------------
;-------------------------------------------------------------------------
;-------------------------------------------------------------------------



; match mips24 to irac 
ir=n_elements(mips24object.mips24ra)
irmatch=fltarr(ir)
irmatch[*]=-999
print,'Matching  mips24 to ir'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,ir-1 do begin

   dist=sphdist(mips24object[q].mips24ra,mips24object[q].mips24dec,iracobject.iracra,iracobject.iracdec, /degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008) then begin
      if iracobject[ind].iracmatch GE 0 then begin  ;if there was a previous match
         if sep lt iracobject[ind].iracmatchdist then begin  ;new match is a closer better match
           ;how do I put the old match into the unmatched list?  need to know its 'q'
            ;aka which number it is in mips24object
            irmatch[iracobject[ind].iracmatch] = -999
            ;put the new match into the matched pile
           irmatch[q]=ind
            iracobject[ind].iracmatch = q
            iracobject[ind].iracmatchdist = sep
         endif
      endif else begin    ;if there was no previous match
         irmatch[q]=ind
         iracobject[ind].iracmatch = q
         iracobject[ind].iracmatchdist = sep
      endelse

   endif 
endfor

print,"Finished at "+systime()
matched=where(irmatch GE 0)
nonmatched = where(irmatch lt 0)
print, n_elements(matched),"matched"
print, n_elements(nonmatched),"nonmatched"

;catalog  the matched and nonmatched objects 
for num=0, n_elements(matched) - 1 do begin
   
   object[num]= {ob, iracobject(irmatch[matched[num]]).iracra, iracobject(irmatch[matched[num]]).iracdec, 0.,0.,-99.,99.,99.,0.,0.,0.,-99., -99., -99. , 0., 0., 0.,0.,-99.,99.,99.,0.,0.,0.,-99.,-99., -99., -99., 0., 0.,  0.,0.,-99.,99.,99.,0.,0, 0.,-99., -99.,-99., 0., 0., 0.,0.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0., 0.,0.,-99.,99.,99.,-99.,99.,99.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99. , 0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99., iracobject(irmatch[matched[num]]).iracra, iracobject(irmatch[matched[num]]).iracdec,iracobject(irmatch[matched[num]]).tmassJflux, iracobject(irmatch[matched[num]]).tmasshflux, iracobject(irmatch[matched[num]]).tmasskflux, iracobject(irmatch[matched[num]]).irac1flux, iracobject(irmatch[matched[num]]).irac2flux,iracobject(irmatch[matched[num]]).irac3flux, iracobject(irmatch[matched[num]]).irac4flux, iracobject(irmatch[matched[num]]).tmassJmag, iracobject(irmatch[matched[num]]).tmasshmag, iracobject(irmatch[matched[num]]).tmasskmag, iracobject(irmatch[matched[num]]).irac1mag, iracobject(irmatch[matched[num]]).irac2mag,iracobject(irmatch[matched[num]]).irac3mag, iracobject(irmatch[matched[num]]).irac4mag,iracobject(irmatch[matched[num]]).mips24, iracobject(irmatch[matched[num]]).iracnndist, iracobject(irmatch[matched[num]]).iracmatchdist,iracobject(irmatch[matched[num]]).iracmatch,mips24object(matched[num]).mips24xcenter , mips24object(matched[num]).mips24ycenter , mips24object(matched[num]).mips24ra , mips24object(matched[num]).mips24dec , mips24object(matched[num]).mips24flux , mips24object(matched[num]).mips24fluxerr ,  mips24object(matched[num]).mips24mag , mips24object(matched[num]).mips24magerr ,mips24object(matched[num]).mips24bckgrnd , mips24object(matched[num]).iter ,mips24object(matched[num]).sharp, mips24object(matched[num]).mips24match, mips24object(matched[num]).mips24nndist, mips24object(matched[num]).mips24matchdist,0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0., -99.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., "none",0.,0.,99., 0.,0.,-99.,0.,0.,0.,0.,0.,-99., 0.,0.,0.,0.,0.,0.,-99.,99.,0., -99., 0., 0.,-99.}

endfor

count = n_elements(matched)
count = long(count)
print, 'count, n_elements(matched)', count, n_elements(matched)
;mips detected, irac non detected
for num= 0,  n_elements(nonmatched) -1. do begin
 


      object[count]= {ob, mips24object(nonmatched[num]).mips24ra , mips24object(nonmatched[num]).mips24dec , 0.,0.,-99.,99.,99.,0.,0.,0.,-99., -99., -99. , 0., 0., 0.,0.,-99.,99.,99.,0.,0.,0.,-99.,-99., -99., -99., 0., 0.,  0.,0.,-99.,99.,99.,0.,0, 0.,-99., -99.,-99., 0., 0., 0.,0.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0., 0.,0.,-99.,99.,99.,-99.,99.,99.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99. , 0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99., 0.,0.,-99.,-99.,-99.,-99.,-99.,-99., -99.,99.,99.,99.,99.,99.,99., 99.,-99., 0., 0.,-99.,mips24object(nonmatched[num]).mips24xcenter , mips24object(nonmatched[num]).mips24ycenter , mips24object(nonmatched[num]).mips24ra , mips24object(nonmatched[num]).mips24dec , mips24object(nonmatched[num]).mips24flux , mips24object(nonmatched[num]).mips24fluxerr ,  mips24object(nonmatched[num]).mips24mag , mips24object(nonmatched[num]).mips24magerr ,mips24object(nonmatched[num]).mips24bckgrnd , mips24object(nonmatched[num]).iter ,mips24object(nonmatched[num]).sharp, mips24object(nonmatched[num]).mips24match, mips24object(nonmatched[num]).mips24nndist, mips24object(nonmatched[num]).mips24matchdist,0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0., -99.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., "none",0.,0.,99., 0.,0.,-99.,0.,0.,0.,0.,0.,-99., 0.,0.,0.,0.,0.,0.,-99.,99.,0., -99., 0., 0.,-99.}

   count = count+1
endfor

print, 'number irac nonmatched', n_elements(where(iracobject.iracmatch lt 1))

;irac detected mips non detected
for num3 = 0l, n_elements(iracobject.iracmatch) -1 do begin
   if iracobject[num3].iracmatch lt 1 then begin
      ;this object did not match to a mips object
   object[count]= {ob, iracobject(num3).iracra, iracobject(num3).iracdec, 0.,0.,-99.,99.,99.,0.,0.,0.,-99., -99., -99. , 0., 0., 0.,0.,-99.,99.,99.,0.,0.,0.,-99.,-99., -99., -99., 0., 0.,  0.,0.,-99.,99.,99.,0.,0, 0.,-99., -99.,-99., 0., 0., 0.,0.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0., 0.,0.,-99.,99.,99.,-99.,99.,99.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99. , 0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99., iracobject(num3).iracra, iracobject(num3).iracdec,iracobject(num3).tmassJflux, iracobject(num3).tmasshflux, iracobject(num3).tmasskflux, iracobject(num3).irac1flux, iracobject(num3).irac2flux,iracobject(num3).irac3flux, iracobject(num3).irac4flux, iracobject(num3).tmassJmag, iracobject(num3).tmasshmag, iracobject(num3).tmasskmag, iracobject(num3).irac1mag, iracobject(num3).irac2mag,iracobject(num3).irac3mag, iracobject(num3).irac4mag,iracobject(num3).mips24, iracobject(num3).iracnndist, iracobject(num3).iracmatchdist,iracobject(num3).iracmatch, 0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0.,0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0., -99.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., "none",0.,0.,99., 0.,0.,-99.,0.,0.,0.,0.,0.,-99., 0.,0.,0.,0.,0.,0.,-99.,99.,0., -99., 0., 0.,-99.}

      count = count + 1
   endif
endfor

print, 'count going into 70micron ', count

;-------------------------------------------------------------------------
;-------------------------------------------------------------------------

;match mips70
m=n_elements(mips70object.mips70ra)
mmatch=fltarr(m)
mmatch[*]=-999

print,'Matching mips70 to object'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,m-1 do begin
   dist=sphdist( mips70object[q].mips70ra, mips70object[q].mips70dec,object.ra,object.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0011)  then begin
      if object[ind].mips70match GE 0 then begin           ;if there was a previous match
         if sep lt object[ind].mips70matchdist then begin  ;new match is a closer better match
           ;how do I put the old match into the unmatched list?  need to know its 'q'
            ;aka which number it is in mips24object
            mmatch[object[ind].mips70match] = -999
            ;put the new match into the matched pile
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
            object[ind].mips70nndist = mips70object[q].mips70nndist
            object[ind].mips70match = q
            object[ind].mips70matchdist = sep
         endif
      endif else begin          ;if there was no previous match
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
         object[ind].mips70nndist = mips70object[q].mips70nndist
         object[ind].mips70matchdist = sep
         object[ind].mips70match = q
      endelse
   endif 
endfor
matched=where(mmatch GE 0)
nonmatched = where(mmatch lt 0)

print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))
print,"Finished at "+systime()

;need to include 70 detected, object not detected.
for num= 0,  n_elements(nonmatched) -1. do begin
 
   object[count]= {ob, mips70object(nonmatched[num]).mips70ra, mips70object(nonmatched[num]).mips70dec, 0.,0.,-99.,99.,99.,0.,0.,0.,-99., -99., -99. , 0., 0., 0.,0.,-99.,99.,99.,0.,0.,0.,-99.,-99., -99., -99., 0., 0.,  0.,0.,-99.,99.,99.,0.,0, 0.,-99., -99.,-99., 0., 0., 0.,0.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0., 0.,0.,-99.,99.,99.,-99.,99.,99.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99. , 0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99., 0.,0.,-99.,-99.,-99.,-99.,-99.,-99., -99.,99.,99.,99.,99.,99.,99., 99.,-99., 0., 0.,-99., 0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0.,mips70object(nonmatched[num]).mips70xcenter, mips70object(nonmatched[num]).mips70ycenter, mips70object(nonmatched[num]).mips70ra, mips70object(nonmatched[num]).mips70dec, mips70object(nonmatched[num]).mips70flux, mips70object(nonmatched[num]).mips70fluxerr, mips70object(nonmatched[num]).mips70mag, mips70object(nonmatched[num]).mips70magerr, mips70object(nonmatched[num]).mips70bckgrnd, mips70object(nonmatched[num]).mips70isoarea,mips70object(nonmatched[num]).mips70fwhm, mips70object(nonmatched[num]).mips70flags, mips70object(nonmatched[num]).mips70nndist, mips70object(nonmatched[num]).mips70matchdist, mips70object(nonmatched[num]).mips70match,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., "none",0.,0.,99., 0.,0.,-99.,0.,0.,0.,0.,0.,-99., 0.,0.,0.,0.,0.,0.,-99.,99.,0., -99., 0., 0.,-99.}

   count = count+1
endfor


;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
print, count, ' going into flamj'
;match nearir
;match flamj
m=n_elements(flamjobject.flamjra)
mmatch=fltarr(m)
mmatch[*]=-999

print,'Matching flamj to object'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,m-1 do begin
   dist=sphdist( flamjobject[q].flamjra, flamjobject[q].flamjdec,object.ra,object.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin
      if object[ind].flamjmatch GE 0 then begin           ;if there was a previous match
         if sep lt object[ind].flamjmatchdist then begin  ;new match is a closer better match
           ;how do I put the old match into the unmatched list?  need to know its 'q'
            ;aka which number it is in mips24object
            mmatch[object[ind].flamjmatch] = -999
            ;put the new match into the matched pile
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
            object[ind].flamjnndist =flamjobject[q].flamjnndist
            object[ind].flamjmatch = q
            object[ind].flamjmatchdist = sep
         endif
      endif else begin          ;if there was no previous match
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
         object[ind].flamjnndist =flamjobject[q].flamjnndist
         object[ind].flamjmatchdist = sep
         object[ind].flamjmatch = q
      endelse
   endif 
endfor
matched=where(mmatch GE 0)
nonmatched = where(mmatch lt 0)

print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))
print,"Finished at "+systime()

;need to include 70 detected, object not detected.
for num= 0,  n_elements(nonmatched) -1. do begin
 

object[count]= {ob, flamjobject(nonmatched[num]).flamjra, flamjobject(nonmatched[num]).flamjdec, 0.,0.,-99.,99.,99.,0.,0.,0.,-99., -99., -99. , 0., 0., 0.,0.,-99.,99.,99.,0.,0.,0.,-99.,-99., -99., -99., 0., 0.,  0.,0.,-99.,99.,99.,0.,0, 0.,-99., -99.,-99., 0., 0., 0.,0.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0., 0.,0.,-99.,99.,99.,-99.,99.,99.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0.,-99.,flamjobject(nonmatched[num]).flamjxcenter, flamjobject(nonmatched[num]).flamjycenter, flamjobject(nonmatched[num]).flamjra, flamjobject(nonmatched[num]).flamjdec, flamjobject(nonmatched[num]).flamjmag, flamjobject(nonmatched[num]).flamjmagerr, flamjobject(nonmatched[num]).flamjfwhm, flamjobject(nonmatched[num]).flamjisoarea, flamjobject(nonmatched[num]).flamjellip, flamjobject(nonmatched[num]).flamjnndist, flamjobject(nonmatched[num]).flamjmatchdist,flamjobject(nonmatched[num]).flamjmatch, 0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99., 0.,0.,-99.,-99.,-99.,-99.,-99.,-99., -99.,99.,99.,99.,99.,99.,99., 99.,-99., 0., 0.,-99., 0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0.,0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0., -99.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., "none",0.,0.,99., 0.,0.,-99.,0.,0.,0.,0.,0.,-99., 0.,0.,0.,0.,0.,0.,-99.,99.,0., -99., 0., 0.,-99.}


   count = count+1
endfor


;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
print, count, ' going into wircj'
;match nearir
;match wircj
m=n_elements(wircjobject.wircjra)
mmatch=fltarr(m)
mmatch[*]=-999

print,'Matching wircj to object'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,m-1 do begin
   dist=sphdist( wircjobject[q].wircjra, wircjobject[q].wircjdec,object.ra,object.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin
      if object[ind].wircjmatch GE 0 then begin           ;if there was a previous match
         if sep lt object[ind].wircjmatchdist then begin  ;new match is a closer better match
           ;how do I put the old match into the unmatched list?  need to know its 'q'
            ;aka which number it is in mips24object
            mmatch[object[ind].wircjmatch] = -999
            ;put the new match into the matched pile
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
            object[ind].wircjnndist =wircjobject[q].wircjnndist
            object[ind].wircjmatch = q
            object[ind].wircjmatchdist = sep
         endif
      endif else begin          ;if there was no previous match
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
         object[ind].wircjnndist =wircjobject[q].wircjnndist
         object[ind].wircjmatchdist = sep
         object[ind].wircjmatch = q
      endelse
   endif 
endfor
matched=where(mmatch GE 0)
nonmatched = where(mmatch lt 0)

print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))
print,"Finished at "+systime()

;need to include 70 detected, object not detected.
for num= 0,  n_elements(nonmatched) -1. do begin
 
object[count]= {ob, wircjobject(nonmatched[num]).wircjra, wircjobject(nonmatched[num]).wircjdec, 0.,0.,-99.,99.,99.,0.,0.,0.,-99., -99., -99. , 0., 0., 0.,0.,-99.,99.,99.,0.,0.,0.,-99.,-99., -99., -99., 0., 0.,  0.,0.,-99.,99.,99.,0.,0, 0.,-99., -99.,-99., 0., 0., 0.,0.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0., 0.,0.,-99.,99.,99.,-99.,99.,99.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99. , wircjobject(nonmatched[num]).wircjxcenter, wircjobject(nonmatched[num]).wircjycenter, wircjobject(nonmatched[num]).wircjra, wircjobject(nonmatched[num]).wircjdec, wircjobject(nonmatched[num]).wircjmag, wircjobject(nonmatched[num]).wircjmagerr, wircjobject(nonmatched[num]).wircjfwhm, wircjobject(nonmatched[num]).wircjisoarea, wircjobject(nonmatched[num]).wircjellip, wircjobject(nonmatched[num]).wircjnndist, wircjobject(nonmatched[num]).wircjmatchdist,wircjobject(nonmatched[num]).wircjmatch,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99., 0.,0.,-99.,-99.,-99.,-99.,-99.,-99., -99.,99.,99.,99.,99.,99.,99., 99.,-99., 0., 0.,-99., 0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0.,0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0., -99.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., "none",0.,0.,99., 0.,0.,-99.,0.,0.,0.,0.,0.,-99., 0.,0.,0.,0.,0.,0.,-99.,99.,0., -99., 0., 0.,-99.}

   count = count+1
endfor


;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
print, count, ' going into wirch'
;match nearir
;match wirch
m=n_elements(wirchobject.wirchra)
mmatch=fltarr(m)
mmatch[*]=-999
print,'Matching wirch to object'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,m-1 do begin
   dist=sphdist( wirchobject[q].wirchra, wirchobject[q].wirchdec,object.ra,object.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin
      if object[ind].wirchmatch GE 0 then begin           ;if there was a previous match
         if sep lt object[ind].wirchmatchdist then begin  ;new match is a closer better match
           ;how do I put the old match into the unmatched list?  need to know its 'q'
            ;aka which number it is in mips24object
            mmatch[object[ind].wirchmatch] = -999
            ;put the new match into the matched pile
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
            object[ind].wirchnndist =wirchobject[q].wirchnndist
            object[ind].wirchmatch = q
            object[ind].wirchmatchdist = sep
         endif
      endif else begin          ;if there was no previous match
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
         object[ind].wirchnndist =wirchobject[q].wirchnndist
         object[ind].wirchmatchdist = sep
         object[ind].wirchmatch = q
      endelse
   endif 
endfor
matched=where(mmatch GE 0)
nonmatched = where(mmatch lt 0)

print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))
print,"Finished at "+systime()

;need to include 70 detected, object not detected.
for num= 0,  n_elements(nonmatched) -1. do begin

object[count]= {ob, wirchobject(nonmatched[num]).wirchra, wirchobject(nonmatched[num]).wirchdec, 0.,0.,-99.,99.,99.,0.,0.,0.,-99., -99., -99. , 0., 0., 0.,0.,-99.,99.,99.,0.,0.,0.,-99.,-99., -99., -99., 0., 0.,  0.,0.,-99.,99.,99.,0.,0, 0.,-99., -99.,-99., 0., 0., 0.,0.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0., 0.,0.,-99.,99.,99.,-99.,99.,99.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99. , 0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99.,wirchobject(nonmatched[num]).wirchxcenter, wirchobject(nonmatched[num]).wirchycenter, wirchobject(nonmatched[num]).wirchra, wirchobject(nonmatched[num]).wirchdec, wirchobject(nonmatched[num]).wirchmag, wirchobject(nonmatched[num]).wirchmagerr, wirchobject(nonmatched[num]).wirchfwhm, wirchobject(nonmatched[num]).wirchisoarea, wirchobject(nonmatched[num]).wirchellip, wirchobject(nonmatched[num]).wirchnndist, wirchobject(nonmatched[num]).wirchmatchdist, wirchobject(nonmatched[num]).wirchmatch,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99., 0.,0.,-99.,-99.,-99.,-99.,-99.,-99., -99.,99.,99.,99.,99.,99.,99., 99.,-99., 0., 0.,-99., 0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0.,0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0., -99.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., "none",0.,0.,99., 0.,0.,-99.,0.,0.,0.,0.,0.,-99., 0.,0.,0.,0.,0.,0.,-99.,99.,0., -99., 0., 0.,-99.}

   count = count+1
endfor


;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
print, count, ' going into wirck'
;match nearir
;match wirck
m=n_elements(wirckobject.wirckra)
mmatch=fltarr(m)
mmatch[*]=-999

print,'Matching wirck to object'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,m-1 do begin
   dist=sphdist( wirckobject[q].wirckra, wirckobject[q].wirckdec,object.ra,object.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin
      if object[ind].wirckmatch GE 0 then begin           ;if there was a previous match
         if sep lt object[ind].wirckmatchdist then begin  ;new match is a closer better match
           ;how do I put the old match into the unmatched list?  need to know its 'q'
            ;aka which number it is in mips24object
            mmatch[object[ind].wirckmatch] = -999
            ;put the new match into the matched pile
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
            object[ind].wircknndist =wirckobject[q].wircknndist
            object[ind].wirckmatch = q
            object[ind].wirckmatchdist = sep
         endif
      endif else begin          ;if there was no previous match
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
         object[ind].wircknndist =wirckobject[q].wircknndist
         object[ind].wirckmatchdist = sep
         object[ind].wirckmatch = q
      endelse
   endif 
endfor
matched=where(mmatch GE 0)
nonmatched = where(mmatch lt 0)

print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))
print,"Finished at "+systime()

;need to include 70 detected, object not detected.
for num= 0,  n_elements(nonmatched) -1. do begin
 

;----
object[count]= {ob, wirckobject(nonmatched[num]).wirckra, wirckobject(nonmatched[num]).wirckdec, 0.,0.,-99.,99.,99.,0.,0.,0.,-99., -99., -99. , 0., 0., 0.,0.,-99.,99.,99.,0.,0.,0.,-99.,-99., -99., -99., 0., 0.,  0.,0.,-99.,99.,99.,0.,0, 0.,-99., -99.,-99., 0., 0., 0.,0.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0., 0.,0.,-99.,99.,99.,-99.,99.,99.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99. , 0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99.,wirckobject(nonmatched[num]).wirckxcenter, wirckobject(nonmatched[num]).wirckycenter, wirckobject(nonmatched[num]).wirckra, wirckobject(nonmatched[num]).wirckdec, wirckobject(nonmatched[num]).wirckmag, wirckobject(nonmatched[num]).wirckmagerr, wirckobject(nonmatched[num]).wirckfwhm, wirckobject(nonmatched[num]).wirckisoarea, wirckobject(nonmatched[num]).wirckellip, wirckobject(nonmatched[num]).wircknndist, wirckobject(nonmatched[num]).wirckmatchdist, wirckobject(nonmatched[num]).wirckmatch, 0.,0.,-99.,-99.,-99.,-99.,-99.,-99., -99.,99.,99.,99.,99.,99.,99., 99.,-99., 0., 0.,-99., 0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0.,0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0., -99.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., "none",0.,0.,99., 0.,0.,-99.,0.,0.,0.,0.,0.,-99., 0.,0.,0.,0.,0.,0.,-99.,99.,0., -99., 0., 0.,-99.}

   count = count+1
endfor


;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
print, count, ' going into palomar optical'
;match g; and consequently u, r, i
m=n_elements(gobject.gra)
mmatch=lonarr(m)
mmatch[*]=-999

print,'Matching g to object'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0l,m-1 do begin
   if q eq fix(m/2) then print, 'halfway ', systime()
   dist=sphdist( gobject[q].gra, gobject[q].gdec,object.ra,object.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin
      if object[ind].gmatch GE 0 then begin           ;if there was a previous match
         if sep lt object[ind].gmatchdist then begin  ;new match is a closer better match
           ;how do I put the old match into the unmatched list?  need to know its 'q'
            ;aka which number it is in mips24object
            mmatch[object[ind].gmatch] = -999
            ;put the new match into the matched pile
            mmatch[q]=ind
            object[ind].gmatch = q
            object[ind].gra = gobject[q].gra
            object[ind].gdec = gobject[q].gdec
            object[ind].gflux = gobject[q].gflux
            object[ind].gfluxmax = gobject[q].gfluxmax
            object[ind].gmag = gobject[q].gmag
            object[ind].gmagerr = gobject[q].gmagerr
            object[ind].gfwhm = gobject[q].gfwhm
            object[ind].gisoarea = gobject[q].gisoarea
            object[ind].gellip = gobject[q].gellip
            object[ind].gclassstar = gobject[q].gclassstar
            object[ind].gflags = gobject[q].gflags
            object[ind].gnndist =gobject[q].gnndist
            object[ind].gmatchdist = sep
;--
            object[ind].ura = uobject[q].ura
            object[ind].udec = uobject[q].udec
            object[ind].uflux = uobject[q].uflux
            object[ind].ufluxmax = uobject[q].ufluxmax
            object[ind].umag = uobject[q].umag
            object[ind].umagerr = uobject[q].umagerr
            object[ind].ufwhm = uobject[q].ufwhm
            object[ind].uisoarea = uobject[q].uisoarea
            object[ind].uellip = uobject[q].uellip
            object[ind].uclassstar = uobject[q].uclassstar
            object[ind].uflags = uobject[q].uflags
            object[ind].unndist =uobject[q].unndist
            object[ind].umatchdist = sep
;--
            object[ind].rra = robject[q].rra
            object[ind].rdec = robject[q].rdec
            object[ind].rflux = robject[q].rflux
            object[ind].rfluxmax = robject[q].rfluxmax
            object[ind].rmag = robject[q].rmag
            object[ind].rmagerr = robject[q].rmagerr
            object[ind].rfwhm = robject[q].rfwhm
            object[ind].risoarea = robject[q].risoarea
            object[ind].rellip = robject[q].rellip
            object[ind].rclassstar = robject[q].rclassstar
            object[ind].rflags = robject[q].rflags
            object[ind].rnndist =robject[q].rnndist
            object[ind].rmatchdist = sep
;---
            object[ind].ira = iobject[q].ira
            object[ind].idec = iobject[q].idec
            object[ind].iflux = iobject[q].iflux
            object[ind].ifluxmax = iobject[q].ifluxmax
            object[ind].imag = iobject[q].imag
            object[ind].imagerr = iobject[q].imagerr
            object[ind].ifwhm = iobject[q].ifwhm
            object[ind].iisoarea = iobject[q].iisoarea
            object[ind].iellip = iobject[q].iellip
            object[ind].iclassstar = iobject[q].iclassstar
            object[ind].iflags = iobject[q].iflags
            object[ind].inndist =iobject[q].inndist
            object[ind].imatchdist = sep

         endif
      endif else begin          ;if there was no previous match
            mmatch[q]=ind
            object[ind].gmatch = q
            object[ind].gra = gobject[q].gra
            object[ind].gdec = gobject[q].gdec
            object[ind].gflux = gobject[q].gflux
            object[ind].gfluxmax = gobject[q].gfluxmax
            object[ind].gmag = gobject[q].gmag
            object[ind].gmagerr = gobject[q].gmagerr
            object[ind].gfwhm = gobject[q].gfwhm
            object[ind].gisoarea = gobject[q].gisoarea
            object[ind].gellip = gobject[q].gellip
            object[ind].gclassstar = gobject[q].gclassstar
            object[ind].gflags = gobject[q].gflags
            object[ind].gmatchdist = sep
            object[ind].gnndist =gobject[q].gnndist
;--
            object[ind].ura = uobject[q].ura
            object[ind].udec = uobject[q].udec
            object[ind].uflux = uobject[q].uflux
            object[ind].ufluxmax = uobject[q].ufluxmax
            object[ind].umag = uobject[q].umag
            object[ind].umagerr = uobject[q].umagerr
            object[ind].ufwhm = uobject[q].ufwhm
            object[ind].uisoarea = uobject[q].uisoarea
            object[ind].uellip = uobject[q].uellip
            object[ind].uclassstar = uobject[q].uclassstar
            object[ind].uflags = uobject[q].uflags
            object[ind].unndist =uobject[q].unndist
            object[ind].umatchdist = sep
;--
            object[ind].rra = robject[q].rra
            object[ind].rdec = robject[q].rdec
            object[ind].rflux = robject[q].rflux
            object[ind].rfluxmax = robject[q].rfluxmax
            object[ind].rmag = robject[q].rmag
            object[ind].rmagerr = robject[q].rmagerr
            object[ind].rfwhm = robject[q].rfwhm
            object[ind].risoarea = robject[q].risoarea
            object[ind].rellip = robject[q].rellip
            object[ind].rclassstar = robject[q].rclassstar
            object[ind].rflags = robject[q].rflags
            object[ind].rnndist =robject[q].rnndist
            object[ind].rmatchdist = sep
;---
            object[ind].ira = iobject[q].ira
            object[ind].idec = iobject[q].idec
            object[ind].iflux = iobject[q].iflux
            object[ind].ifluxmax = iobject[q].ifluxmax
            object[ind].imag = iobject[q].imag
            object[ind].imagerr = iobject[q].imagerr
            object[ind].ifwhm = iobject[q].ifwhm
            object[ind].iisoarea = iobject[q].iisoarea
            object[ind].iellip = iobject[q].iellip
            object[ind].iclassstar = iobject[q].iclassstar
            object[ind].iflags = iobject[q].iflags
            object[ind].inndist =iobject[q].inndist
            object[ind].imatchdist = sep
      endelse
   endif 
endfor
matched=where(mmatch GE 0)
nonmatched = where(mmatch lt 0)

print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))
print,"Finished at "+systime()

;need to include 70 detected, object not detected.
for num= 0l,  n_elements(nonmatched) -1. do begin
 
object[count]= {ob, gobject(nonmatched[num]).gra, gobject(nonmatched[num]).gdec, uobject(nonmatched[num]).ura,uobject(nonmatched[num]).udec,uobject(nonmatched[num]).uflux,uobject(nonmatched[num]).umag,uobject(nonmatched[num]).umagerr,uobject(nonmatched[num]).ufwhm,uobject(nonmatched[num]).uisoarea,uobject(nonmatched[num]).uellip,uobject(nonmatched[num]).ufluxmax, uobject(nonmatched[num]).uclassstar, uobject(nonmatched[num]).uflags , uobject(nonmatched[num]).unndist, uobject(nonmatched[num]).umatchdist, gobject(nonmatched[num]).gra,gobject(nonmatched[num]).gdec,gobject(nonmatched[num]).gflux,gobject(nonmatched[num]).gmag,gobject(nonmatched[num]).gmagerr,gobject(nonmatched[num]).gfwhm,gobject(nonmatched[num]).gisoarea,gobject(nonmatched[num]).gellip,gobject(nonmatched[num]).gfluxmax,gobject(nonmatched[num]).gmatch, gobject(nonmatched[num]).gclassstar, gobject(nonmatched[num]).gflags, gobject(nonmatched[num]).gnndist, gobject(nonmatched[num]).gmatchdist,  robject(nonmatched[num]).rra,robject(nonmatched[num]).rdec,robject(nonmatched[num]).rflux,robject(nonmatched[num]).rmag,robject(nonmatched[num]).rmagerr,robject(nonmatched[num]).rfwhm,robject(nonmatched[num]).risoarea, robject(nonmatched[num]).rellip,robject(nonmatched[num]).rfluxmax, robject(nonmatched[num]).rclassstar,robject(nonmatched[num]).rflags, robject(nonmatched[num]).rnndist, robject(nonmatched[num]).rmatchdist, iobject(nonmatched[num]).ira,iobject(nonmatched[num]).idec,iobject(nonmatched[num]).iflux,iobject(nonmatched[num]).imag,iobject(nonmatched[num]).imagerr,iobject(nonmatched[num]).ifwhm,iobject(nonmatched[num]).iisoarea, iobject(nonmatched[num]).iellip,iobject(nonmatched[num]).ifluxmax, iobject(nonmatched[num]).iclassstar, iobject(nonmatched[num]).iflags, iobject(nonmatched[num]).inndist, iobject(nonmatched[num]).imatchdist, 0.,0.,-99.,99.,99.,-99.,99.,99.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99. , 0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99., 0.,0.,-99.,-99.,-99.,-99.,-99.,-99., -99.,99.,99.,99.,99.,99.,99., 99.,-99., 0., 0.,-99., 0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0.,0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0., -99.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., "none",0.,0.,99., 0.,0.,-99.,0.,0.,0.,0.,0.,-99., 0.,0.,0.,0.,0.,0.,-99.,99.,0., -99., 0., 0.,-99.}

   count = count+1
endfor

;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
print, count, ' going into z'
;match nearir
;match wirck
m=n_elements(zobject.zra)
mmatch=lonarr(m)
mmatch[*]=-999

print,'Matching z to object'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0l,m-1 do begin
   dist=sphdist( zobject[q].zra, zobject[q].zdec,object.ra,object.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin
      if object[ind].zmatch GE 0 then begin           ;if there was a previous match
         if sep lt object[ind].zmatchdist then begin  ;new match is a closer better match
           ;how do I put the old match into the unmatched list?  need to know its 'q'
            ;aka which number it is in mips24object
            mmatch[object[ind].zmatch] = -999
            ;put the new match into the matched pile
            mmatch[q]=ind
            object[ind].zra = zobject[q].zra
            object[ind].zdec = zobject[q].zdec
            object[ind].zfluxauto = zobject[q].zfluxauto
            object[ind].zmagauto = zobject[q].zmagauto
            object[ind].zmagerrauto = zobject[q].zmagerrauto
            object[ind].zfluxaper = zobject[q].zfluxaper
            object[ind].zmagaper = zobject[q].zmagaper
            object[ind].zmagerraper = zobject[q].zmagerraper
            object[ind].zfluxbest = zobject[q].zfluxbest
            object[ind].zmagbest = zobject[q].zmagbest
            object[ind].zmagerrbest = zobject[q].zmagerrbest
            object[ind].zfwhm = zobject[q].zfwhm
            object[ind].zisoarea = zobject[q].zisoarea
            object[ind].zellip = zobject[q].zellip
            object[ind].zfluxmax = zobject[q].zfluxmax
            object[ind].zclassstar = zobject[q].zclassstar
            object[ind].zflags = zobject[q].zflags
            object[ind].znndist =zobject[q].znndist
           object[ind].zmatch = q
            object[ind].zmatchdist = sep
         endif
      endif else begin          ;if there was no previous match
         mmatch[q]=ind
         object[ind].zra = zobject[q].zra
         object[ind].zdec = zobject[q].zdec
         object[ind].zfluxauto = zobject[q].zfluxauto
         object[ind].zmagauto = zobject[q].zmagauto
         object[ind].zmagerrauto = zobject[q].zmagerrauto
         object[ind].zfluxaper = zobject[q].zfluxaper
         object[ind].zmagaper = zobject[q].zmagaper
         object[ind].zmagerraper = zobject[q].zmagerraper
         object[ind].zfluxbest = zobject[q].zfluxbest
         object[ind].zmagbest = zobject[q].zmagbest
         object[ind].zmagerrbest = zobject[q].zmagerrbest
         object[ind].zfwhm = zobject[q].zfwhm
         object[ind].zisoarea = zobject[q].zisoarea
         object[ind].zellip = zobject[q].zellip
         object[ind].zfluxmax = zobject[q].zfluxmax
         object[ind].zclassstar = zobject[q].zclassstar
         object[ind].zflags = zobject[q].zflags
         object[ind].znndist =zobject[q].znndist
         object[ind].zmatch = q
         object[ind].zmatchdist = sep
      endelse
   endif 
endfor
matched=where(mmatch GE 0)
nonmatched = where(mmatch lt 0)

print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))
print,"Finished at "+systime()

;need to include z detected, object not detected.
for num= 0,  n_elements(nonmatched) -1. do begin
 
object[count]= {ob,zobject(nonmatched[num]).zra, zobject(nonmatched[num]).zdec, 0.,0.,-99.,99.,99.,0.,0.,0.,-99., -99., -99. , 0., 0., 0.,0.,-99.,99.,99.,0.,0.,0.,-99.,-99., -99., -99., 0., 0.,  0.,0.,-99.,99.,99.,0.,0, 0.,-99., -99.,-99., 0., 0., 0.,0.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0.,  zobject(nonmatched[num]).zra,zobject(nonmatched[num]).zdec,zobject(nonmatched[num]).zfluxauto,zobject(nonmatched[num]).zmagauto,zobject(nonmatched[num]).zmagerrauto,zobject(nonmatched[num]).zfluxaper,zobject(nonmatched[num]).zmagaper,zobject(nonmatched[num]).zmagerraper,zobject(nonmatched[num]).zfluxbest,zobject(nonmatched[num]).zmagbest,zobject(nonmatched[num]).zmagerrbest,zobject(nonmatched[num]).zfwhm,zobject(nonmatched[num]).zisoarea, zobject(nonmatched[num]).zellip,zobject(nonmatched[num]).zfluxmax, zobject(nonmatched[num]).zclassstar, zobject(nonmatched[num]).zflags, zobject(nonmatched[num]).znndist, zobject(nonmatched[num]).zmatchdist,zobject(nonmatched[num]).zmatch,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99. , 0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99., 0.,0.,-99.,-99.,-99.,-99.,-99.,-99., -99.,99.,99.,99.,99.,99.,99., 99.,-99., 0., 0.,-99., 0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0.,0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0., -99.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., "none",0.,0.,99., 0.,0.,-99.,0.,0.,0.,0.,0.,-99., 0.,0.,0.,0.,0.,0.,-99.,99.,0., -99., 0., 0.,-99.}

   count = count+1
endfor

;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
print, count, ' going into acs'
;match nearir
;match wirck
m=n_elements(acsobject.acsra)
mmatch=lonarr(m)
mmatch[*]=-999

print,'Matching acs to object'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=long(0),m-1 do begin
   dist=sphdist( acsobject[q].acsra, acsobject[q].acsdec,object.ra,object.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin
      if object[ind].acsmatch GE 0 then begin           ;if there was a previous match
         if sep lt object[ind].acsmatchdist then begin  ;new match is a closer better match
           ;how do I put the old match into the unmatched list?  need to know its 'q'
            ;aka which number it is in mips24object
            mmatch[object[ind].acsmatch] = -999
            ;put the new match into the matched pile
            mmatch[q]=ind
            object[ind].acsconcentration = acsobject[q].acsconcentration
            object[ind].acsgini = acsobject[q].acsgini
            object[ind].acsmag = acsobject[q].acsmag
            object[ind].acsmu = acsobject[q].acsmu
            object[ind].acsellip = acsobject[q].acsellip
            object[ind].acsflags = acsobject[q].acsflags
            object[ind].acsra = acsobject[q].acsra
            object[ind].acsdec = acsobject[q].acsdec
            object[ind].acsisoarea = acsobject[q].acsisoarea
            object[ind].acstheta = acsobject[q].acstheta
            object[ind].acssnr = acsobject[q].acssnr
            object[ind].acsflux = acsobject[q].acsflux
            object[ind].acssky = acsobject[q].acssky
            object[ind].acsskyrms = acsobject[q].acsskyrms
            object[ind].acscentralmu = acsobject[q].acscentralmu
            object[ind].acsassym = acsobject[q].acsassym
            object[ind].acshlightrad= acsobject[q].acshlightrad
            object[ind].acsm20 = acsobject[q].acsm20
            object[ind].acsclassstar = acsobject[q].acsclassstar
            object[ind].acsmagerr = acsobject[q].acsmagerr
             object[ind].acsfwhm = acsobject[q].acsfwhm
            object[ind].acsfluxmax = acsobject[q].acsfluxmax
            object[ind].acsnndist =acsobject[q].acsnndist
            object[ind].acsmatch = q
            object[ind].acsmatchdist = sep
         endif
      endif else begin          ;if there was no previous match
         mmatch[q]=ind
            object[ind].acsconcentration = acsobject[q].acsconcentration
            object[ind].acsgini = acsobject[q].acsgini
            object[ind].acsmag = acsobject[q].acsmag
            object[ind].acsmu = acsobject[q].acsmu
            object[ind].acsellip = acsobject[q].acsellip
            object[ind].acsflags = acsobject[q].acsflags
            object[ind].acsra = acsobject[q].acsra
            object[ind].acsdec = acsobject[q].acsdec
            object[ind].acsisoarea = acsobject[q].acsisoarea
            object[ind].acstheta = acsobject[q].acstheta
            object[ind].acssnr = acsobject[q].acssnr
            object[ind].acsflux = acsobject[q].acsflux
            object[ind].acssky = acsobject[q].acssky
            object[ind].acsskyrms = acsobject[q].acsskyrms
            object[ind].acscentralmu = acsobject[q].acscentralmu
            object[ind].acsassym = acsobject[q].acsassym
            object[ind].acshlightrad= acsobject[q].acshlightrad
            object[ind].acsm20 = acsobject[q].acsm20
            object[ind].acsclassstar = acsobject[q].acsclassstar
            object[ind].acsmagerr = acsobject[q].acsmagerr
             object[ind].acsfwhm = acsobject[q].acsfwhm
            object[ind].acsfluxmax = acsobject[q].acsfluxmax
            object[ind].acsnndist =acsobject[q].acsnndist
            object[ind].acsmatch = q
            object[ind].acsmatchdist = sep
      endelse
   endif 
endfor
matched=where(mmatch GE 0)
nonmatched = where(mmatch lt 0)

print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))
print,"Finished at "+systime()

;need to include acs detected, object not detected.
for num=long(0),  n_elements(nonmatched) -1. do begin
 

;----
object[count]= {ob, acsobject(nonmatched[num]).acsra,acsobject(nonmatched[num]).acsdec,0.,0.,-99.,99.,99.,0.,0.,0.,-99., -99., -99. , 0., 0., 0.,0.,-99.,99.,99.,0.,0.,0.,-99.,-99., -99., -99., 0., 0.,  0.,0.,-99.,99.,99.,0.,0, 0.,-99., -99.,-99., 0., 0.,0.,0.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0., 0.,0.,-99.,99.,99.,-99.,99.,99.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99. , 0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99., 0.,0.,-99.,-99.,-99.,-99.,-99.,-99., -99.,99.,99.,99.,99.,99.,99., 99.,-99., 0., 0.,-99., 0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0.,0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0., -99.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., "none",acsobject(nonmatched[num]).acsconcentration,acsobject(nonmatched[num]).acsgini,acsobject(nonmatched[num]).acsmag, acsobject(nonmatched[num]).acsmu,acsobject(nonmatched[num]).acsellip,acsobject(nonmatched[num]).acsflags,acsobject(nonmatched[num]).acsra,acsobject(nonmatched[num]).acsdec,acsobject(nonmatched[num]).acsisoarea,acsobject(nonmatched[num]).acstheta,acsobject(nonmatched[num]).acssnr,acsobject(nonmatched[num]).acsflux, acsobject(nonmatched[num]).acssky,acsobject(nonmatched[num]).acsskyrms,acsobject(nonmatched[num]).acscentralmu,acsobject(nonmatched[num]).acsassym,acsobject(nonmatched[num]).acshlightrad,acsobject(nonmatched[num]).acsm20,acsobject(nonmatched[num]).acsclassstar,acsobject(nonmatched[num]).acsmagerr,acsobject(nonmatched[num]).acsfwhm, acsobject(nonmatched[num]).acsfluxmax, acsobject(nonmatched[num]).acsnndist, acsobject(nonmatched[num]).acsmatchdist,acsobject(nonmatched[num]).acsmatch}


   count = count+1
endfor


;-------------------------------------------------------------------------
;-------------------------------------------------------------------------

print, 'There area ', count, ' objects in the catalog'
object = object[0:count-1]
save, object, filename='/Users/jkrick/idlbin/object.spitzer.sav'
;save, infraredobject, filename='/Users/jkrick/idlbin/infraredobject.sav'

;close, outlun
;free_lun, outlun
;close, outlun2
;free_lun, outlun2
;close, outlun3
;free_lun, outlun3
;close, outlun4
;free_lun, outlun4


undefine, robject
undefine, gobject
undefine, object
undefine, iobject
undefine, uobject
undefine, zobject
undefine, acsobject
undefine, infraredobject
undefine, mips24object
undefine, flamjobject
undefine, wircjobject
undefine, wirchobject
undefine, wirckobject
end



;----------------------------------------------------------------
; nearir combined data
; have used topcat to match flamj, wircj, wirch, wirck
;have not matched on flux, but could add that to the topcat match
;actually I think I do want to match all the individual catalogs so that I can still have match distances....

;openr, luni, "/Users/jkrick/nep/nearir.topcat.txt", /get_lun, error=err
;if (err ne 0.) then print, "file didn't open", ﻿!ERROR_STATE.MSG  
;p = 0l
;WHILE (NOT EOF(luni)) DO BEGIN

;   READF, luni, junk,x_image, y_image,junk,junk,junk,junk,junk,junk,mag_best, magerr_best, junk,junk, junk,isoareaimage, ra,dec,junk,junk,junk, junk,fwhm_image, junk,ell, junk, junk,x_image2, y_image2,junk,junk,junk,junk,junk,junk,mag_best2, magerr_best2, junk,junk, junk,isoareaimage2, ra2,dec2,junk,junk,junk, junk,fwhm_image2, junk,ell2, junk , junk,x_image3, y_image3,junk,junk,junk,junk,junk,junk,mag_best3, magerr_best3, junk,junk, junk,isoareaimage3, ra3,dec3,junk,junk,junk, junk,fwhm_image3, junk,ell3, junk ,  junk,x_image4, y_image4,junk,junk,junk,junk,junk,junk,mag_best4, magerr_best4, junk,junk, junk,isoareaimage4, ra4,dec4,junk,junk,junk, junk,fwhm_image4, junk,ell4, junk 

;if ra4 gt 200 then begin
;   ora = ra4
;   odec = dec4
;endif
;if ra3 gt 200 then begin
;   ora = ra3
;   odec = dec3
;endif
;if ra2 gt 200 then begin
;   ora = ra2
;   odec = dec2
;endif
;if ra gt 200 then begin
;   ora = ra
;   odec = dec
;endif 

;   nearirobject[p] ={nearirob, ora, odec, x_image, y_image, ra, dec, mag_best, magerr_best, fwhm_image, isoareaimage, ell,0.,0.,-1., x_image2, y_image2, ra2, dec2, mag_best2, magerr_best2, fwhm_image2, isoareaimage2, ell2,0.,0.,-1.,x_image3, y_image3, ra3, dec3, mag_best3, magerr_best3, fwhm_image3, isoareaimage3, ell3,0.,0.,-1., x_image4, y_image4, ra4, dec4, mag_best4, magerr_best4, fwhm_image4, isoareaimage4, ell4,0.,0.,-1.}


;   p = p + 1
   
; endwhile
;nearirobject = nearirobject[0:p-1]
;print, "there are ",p," nearir objects"

;close, luni
;free_lun, luni
;----------------------------------------------------------------
