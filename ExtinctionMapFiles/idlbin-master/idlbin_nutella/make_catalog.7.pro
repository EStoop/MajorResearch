pro make_catalog


close, /all

!P.multi = [0,2,2]

numoobjects = 200000
colorarr = fltarr(numoobjects)
magarr = fltarr(numoobjects)
FITS_READ, '/Users/jkrick/palomar/lfc/coadd_g.fits',gdata, gheader
FITS_READ, '/Users/jkrick/spitzer/irac/ch1/mosaic.fits',irdata, irheader

uobject = replicate({uob, uxcenter:0D,uycenter:0D,ufluxa:0D,umaga:0D,umagerra:0D,ufwhm:0D,uisoarea:0D,uellip:0D,ufluxmax:0D, uclassstar:0D, uflags:0D, unndist:0D, umatchdist:0D},numoobjects)
gobject = replicate({gob, gxcenter:0D,gycenter:0D,gfluxa:0D,gmaga:0D,gmagerra:0D,gfwhm:0D,gisoarea:0D,gellip:0D,gfluxmax:0D,gmatch:0D, gclassstar:0D, gflags:0D, gnndist:0D, gmatchdist:0D},numoobjects)
robject = replicate({rob, rxcenter:0D,rycenter:0D,rfluxa:0D,rmaga:0D,rmagerra:0D,rfwhm:0D,risoarea:0D, rellip:0D,rfluxmax:0D, rclassstar:0D,rflags:0D, rnndist:0D, rmatchdist:0D},numoobjects)
iobject = replicate({iob, ixcenter:0D,iycenter:0D,ifluxa:0D,imaga:0D,imagerra:0D,ifwhm:0D,iisoarea:0D, iellip:0D,ifluxmax:0D, iclassstar:0D, iflags:0D, inndist:0D, imatchdist:0D},numoobjects)

zobject = replicate({zob, zra:0D,zdec:0D,zfluxauto:0D,zmagauto:0D,zmagerrauto:0D,zfluxaper:0D,zmagaper:0D,zmagerraper:0D,zfluxbest:0D,zmagbest:0D,zmagerrbest:0D,zfwhm:0D,zisoarea:0D, zellip:0D,zfluxmax:0D, zclassstar:0D, zflags:0D, znndist:0D, zmatchdist:0D,zmatch:0D},numoobjects)

acsobject =replicate({acsob,acsconcentration:0D,acsgini:0D,acsmag:0D, acsmu:0D,acsellip:0D,acsflag:0D,acsra:0D,acsdec:0D,acsisoarea:0D,acstheta:0D,acssnr:0D,acsflux:0D, acssky:0D,acsskyrms:0D,acscentralmu:0D,acsassym:0D,acshlightrad:0D,acsm20:0D,acsclassstar:0D,acsmagerr:0D,acsfwhm:0D, acsfluxmax:0D, acsnndist:0D, acsmatchdist:0D, acsmatch:0D},80000)
;-------
flamjobject = replicate({flamjob,flamjxcenter:0D, flamjycenter:0D, flamjra:0D, flamjdec:0D, flamjmag:0D, flamjmagerr:0D, flamjfwhm:0D, flamjisoarea:0D, flamjellip:0D, flamjnndist:0D, flamjmatchdist:0D, flamjmatch:0D },numoobjects)
wircjobject = replicate({wircjob, wircjxcenter:0D, wircjycenter:0D, wircjra:0D, wircjdec:0D, wircjmag:0D, wircjmagerr:0D, wircjfwhm:0D, wircjisoarea:0D, wircjellip:0D, wircjnndist:0D, wircjmatchdist:0D, wircjmatch:0D },numoobjects)
wirchobject = replicate({wirchob, wirchxcenter:0D, wirchycenter:0D, wirchra:0D, wirchdec:0D, wirchmag:0D, wirchmagerr:0D, wirchfwhm:0D, wirchisoarea:0D, wirchellip:0D , wirchnndist:0D, wirchmatchdist:0D, wirchmatch:0D},numoobjects)
wirckobject = replicate({wirckob, wirckxcenter:0D, wirckycenter:0D, wirckra:0D, wirckdec:0D, wirckmag:0D, wirckmagerr:0D, wirckfwhm:0D, wirckisoarea:0D, wirckellip:0D, wircknndist:0D, wirckmatchdist:0D , wirckmatch:0D},numoobjects)

mips24object = replicate({mips24ob, mips24xcenter:0D, mips24ycenter:0D, mips24ra:0D, mips24dec:0D, mips24flux:0D, mips24fluxerr:0D, mips24mag:0D, mips24magerr:0D, mips24bckgrnd:0D, iter:0D,sharp:0D, mips24match:0D, mips24nndist:0D, mips24matchdist:0D},numoobjects)

mips70object = replicate({mips70ob, mips70xcenter:0D, mips70ycenter:0D, mips70ra:0D, mips70dec:0D, mips70flux:0D, mips70fluxerr:0D, mips70mag:0D, mips70magerr:0D, mips70bckgrnd:0D, mips70isoarea:0D,mips70fwhm:0D, mips70flags:0D, mips70nndist:0D, mips70matchdist:0D, mips70match:0D},numoobjects)

iracobject = replicate({iracob, irxcenter:0D,irycenter:0D,tmassJflux:0D,tmasshflux:0D,tmasskflux:0D,irac1flux:0D,irac2flux:0D,irac3flux:0D, irac4flux:0D,tmassJmag:0D,tmasshmag:0D,tmasskmag:0D,irac1mag:0D,irac2mag:0D,irac3mag:0D, irac4mag:0D,mips24:0D, iracnndist:0D, iracmatchdist:0D, iracmatch:0D},numoobjects)

infraredobject = replicate({spitzob, INHERITS iracob,INHERITS mips24ob,INHERITS mips70ob},numoobjects)
;=======
object= replicate({ob, ra:0D, dec:0D, INHERITS uob, INHERITS gob, INHERITS rob, INHERITS iob, INHERITS zob, INHERITS flamjob, INHERITS wircjob, INHERITS wirchob, INHERITS wirckob, INHERITS spitzob, zphot:0D, chisq:0D, prob:0D, spt:0D, age:0D, av:0D, Mabs:0D, mass:0D, plusmasserr:0D, minusmasserr:0D, masschi:0D, massprob:0D, massage:0D, model:' ',INHERITS acsob}, numoobjects)

;----------------------------------------------------------------
;template objects
;object= {ob, ra, dec, uxcenter,uycenter,ufluxa,umaga,umagerra,ufwhm,uisoarea,uellip,ufluxmax, uclassstar, uflags , unndist, umatchdist, gxcenter,gycenter,gfluxa,gmaga,gmagerra,gfwhm,gisoarea,gellip,gfluxmax,gmatch, gclassstar, gflags, gnndist, gmatchdist,  rxcenter,rycenter,rfluxa,rmaga,rmagerra,rfwhm,risoarea, rellip,rfluxmax, rclassstar,rflags, rnndist, rmatchdist, ixcenter,iycenter,ifluxa,imaga,imagerra,ifwhm,iisoarea, iellip,ifluxmax, iclassstar, iflags, inndist, imatchdist, zra,zdec,zfluxauto,zmagauto,zmagerrauto,zfluxaper,zmagaper,zmagerraper,zfluxbest,zmagbest,zmagerrbest,zfwhm,zisoarea, zellip,zfluxmax, zclassstar, zflags, znndist, zmatchdist,zmatch,flamjxcenter, flamjycenter, flamjra, flamjdec, flamjmag, flamjmagerr, flamjfwhm, flamjisoarea, flamjellip, flamjnndist, flamjmatchdist,flamjmatch , wircjxcenter, wircjycenter, wircjra, wircjdec, wircjmag, wircjmagerr, wircjfwhm, wircjisoarea, wircjellip, wircjnndist, wircjmatchdist,wircjmatch,wirchxcenter, wirchycenter, wirchra, wirchdec, wirchmag, wirchmagerr, wirchfwhm, wirchisoarea, wirchellip, wirchnndist, wirchmatchdist, wirchmatch,wirckxcenter, wirckycenter, wirckra, wirckdec, wirckmag, wirckmagerr, wirckfwhm, wirckisoarea, wirckellip, wircknndist, wirckmatchdist, wirckmatch, irxcenter,irycenter,tmassJflux,tmasshflux,tmasskflux,irac1flux,irac2flux,irac3flux, irac4flux,tmassJmag,tmasshmag,tmasskmag,irac1mag,irac2mag,irac3mag, irac4mag,mips24, iracnndist, iracmatchdist,iracmatch, mips24xcenter, mips24ycenter, mips24ra, mips24dec, mips24flux, mips24fluxerr, mips24mag, mips24magerr, mips24bckgrnd, iter,sharp, mips24match, mips24nndist, mips24matchdist,mips70xcenter, mips70ycenter, mips70ra, mips70dec, mips70flux, mips70fluxerr, mips70mag, mips70magerr, mips70bckgrnd, mips70isoarea,mips70fwhm, mips70flags, mips70nndist, mips70matchdist, mips70match,zphot, chisq, prob, spt, age, av, Mabs, mass, plusmasserr, minusmasserr, masschi, massprob, massage, model,acsconcentration,acsgini,acsmag, acsmu,acsellip,acsflag,acsra,acsdec,acsisoarea,acstheta,acssnr,acsflux, acssky,acsskyrms,acscentralmu,acsassym,acshlightrad,acsm20,acsclassstar,acsmagerr,acsfwhm, acsfluxmax, acsnndist, acsmatchdist,acsmatch}

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
;----------------------------------------------------------------


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


;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
;Matching--------------------------------------------------------------------
;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
; match mips24 to irac 
ir=n_elements(iracobject.irxcenter)
irmatch=fltarr(ir)
irmatch[*]=-999
print,'Matching  mips24 to ir'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,ir-1 do begin

   dist=sphdist(iracobject[q].irxcenter,iracobject[q].irycenter, mips24object.mips24ra,mips24object.mips24dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008) then begin
      if mips24object[ind].mips24match eq 1 then begin
         print, 'already matched', iracobject[q].irxcenter, iracobject[q].irycenter, mips24object[ind].mips24ra, mips24object[ind].mips24dec, mips24object[ind].mips24matchdist, sep
      endif

      irmatch[q]=ind
      mips24object[ind].mips24match = 1
      mips24object[ind].mips24matchdist = sep
   endif 
endfor

print,"Finished at "+systime()
matched=where(irmatch GT 0)
nonmatched = where(irmatch lt 0)
print, n_elements(matched),"matched"
print, n_elements(nonmatched),"nonmatched"

;catalog  the matched and nonmatched objects 
for num=0, n_elements(matched) - 1 do begin
   
   object[num]= {ob, iracobject(matched[num]).irxcenter, iracobject(matched[num]).irycenter, 0.,0.,-99.,99.,99.,0.,0.,0.,-99., -99., -99. , 0., 0., 0.,0.,-99.,99.,99.,0.,0.,0.,-99.,-99., -99., -99., 0., 0.,  0.,0.,-99.,99.,99.,0.,0, 0.,-99., -99.,-99., 0., 0., 0.,0.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0., 0.,0.,-99.,99.,99.,-99.,99.,99.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99. , 0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99., iracobject(matched[num]).irxcenter, iracobject(matched[num]).irycenter,iracobject(matched[num]).tmassJflux, iracobject(matched[num]).tmasshflux, iracobject(matched[num]).tmasskflux, iracobject(matched[num]).irac1flux, iracobject(matched[num]).irac2flux,iracobject(matched[num]).irac3flux, iracobject(matched[num]).irac4flux, iracobject(matched[num]).tmassJmag, iracobject(matched[num]).tmasshmag, iracobject(matched[num]).tmasskmag, iracobject(matched[num]).irac1mag, iracobject(matched[num]).irac2mag,iracobject(matched[num]).irac3mag, iracobject(matched[num]).irac4mag,iracobject(matched[num]).mips24, iracobject(matched[num]).iracnndist, iracobject(matched[num]).iracmatchdist,iracobject(matched[num]).iracmatch,mips24object(irmatch[matched[num]]).mips24xcenter , mips24object(irmatch[matched[num]]).mips24ycenter , mips24object(irmatch[matched[num]]).mips24ra , mips24object(irmatch[matched[num]]).mips24dec , mips24object(irmatch[matched[num]]).mips24flux , mips24object(irmatch[matched[num]]).mips24fluxerr ,  mips24object(irmatch[matched[num]]).mips24mag , mips24object(irmatch[matched[num]]).mips24magerr ,mips24object(irmatch[matched[num]]).mips24bckgrnd , mips24object(irmatch[matched[num]]).iter ,mips24object(irmatch[matched[num]]).sharp, mips24object(irmatch[matched[num]]).mips24match, mips24object(irmatch[matched[num]]).mips24nndist, mips24object(irmatch[matched[num]]).mips24matchdist,0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0., -99.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., "none",0.,0.,99., 0.,0.,-99.,0.,0.,0.,0.,0.,-99., 0.,0.,0.,0.,0.,0.,-99.,99.,0., -99., 0., 0.,-99.}

endfor

count = n_elements(matched)
for num= 0,  n_elements(nonmatched) -1. do begin
 
   object[count]= {ob, iracobject(nonmatched[num]).irxcenter, iracobject(nonmatched[num]).irycenter, 0.,0.,-99.,99.,99.,0.,0.,0.,-99., -99., -99. , 0., 0., 0.,0.,-99.,99.,99.,0.,0.,0.,-99.,-99., -99., -99., 0., 0.,  0.,0.,-99.,99.,99.,0.,0, 0.,-99., -99.,-99., 0., 0., 0.,0.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0., 0.,0.,-99.,99.,99.,-99.,99.,99.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99. , 0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99., iracobject(nonmatched[num]).irxcenter, iracobject(nonmatched[num]).irycenter,iracobject(nonmatched[num]).tmassJflux, iracobject(nonmatched[num]).tmasshflux, iracobject(nonmatched[num]).tmasskflux, iracobject(nonmatched[num]).irac1flux, iracobject(nonmatched[num]).irac2flux,iracobject(nonmatched[num]).irac3flux, iracobject(nonmatched[num]).irac4flux, iracobject(nonmatched[num]).tmassJmag, iracobject(nonmatched[num]).tmasshmag, iracobject(nonmatched[num]).tmasskmag, iracobject(nonmatched[num]).irac1mag, iracobject(nonmatched[num]).irac2mag,iracobject(nonmatched[num]).irac3mag, iracobject(nonmatched[num]).irac4mag,iracobject(nonmatched[num]).mips24, iracobject(nonmatched[num]).iracnndist, iracobject(nonmatched[num]).iracmatchdist,iracobject(nonmatched[num]).iracmatch, 0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0.,0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0., -99.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., "none",0.,0.,99., 0.,0.,-99.,0.,0.,0.,0.,0.,-99., 0.,0.,0.,0.,0.,0.,-99.,99.,0., -99., 0., 0.,-99.}
   count = count+1
endfor

print, 'number 24 nonmatched', n_elements(where(mips24object.mips24match lt 1))
for num3 = 0l, n_elements(mips24object.mips24match) -1 do begin
   if mips24object[num3].mips24match lt 1 then begin
      ;this object did not match to an irac object
      ;is it a real object?


      object[count]= {ob, mips24object(num3).mips24ra , mips24object(num3).mips24dec , 0.,0.,-99.,99.,99.,0.,0.,0.,-99., -99., -99. , 0., 0., 0.,0.,-99.,99.,99.,0.,0.,0.,-99.,-99., -99., -99., 0., 0.,  0.,0.,-99.,99.,99.,0.,0, 0.,-99., -99.,-99., 0., 0., 0.,0.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0., 0.,0.,-99.,99.,99.,-99.,99.,99.,-99.,99.,99.,0.,0., 0.,-99., -99., -99., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99. , 0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0.,-99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99.,0., 0., 0., 0., 99., 99., 0., 0., 0., 0., 0., -99., 0.,0.,-99.,-99.,-99.,-99.,-99.,-99., -99.,99.,99.,99.,99.,99.,99., 99.,-99., 0., 0.,-99.,mips24object(num3).mips24xcenter , mips24object(num3).mips24ycenter , mips24object(num3).mips24ra , mips24object(num3).mips24dec , mips24object(num3).mips24flux , mips24object(num3).mips24fluxerr ,  mips24object(num3).mips24mag , mips24object(num3).mips24magerr ,mips24object(num3).mips24bckgrnd , mips24object(num3).iter ,mips24object(num3).sharp, mips24object(num3).mips24match, mips24object(num3).mips24nndist, mips24object(num3).mips24matchdist,0., 0., 0., 0., -99., -99., 99., 99., 0., 0.,0., -99., 0., 0., -99.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., "none",0.,0.,99., 0.,0.,-99.,0.,0.,0.,0.,0.,-99., 0.,0.,0.,0.,0.,0.,-99.,99.,0., -99., 0., 0.,-99.}

      count = count + 1
   endif
endfor


;-------------------------------------------------------------------------
;-------------------------------------------------------------------------

;match mips70
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
   if (sep LT 0.0011)  then begin
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
      object[ind].mips70matchdist = sep
      mips70object[q].mips70match = 1

   endif 
endfor
matched=where(mmatch GT 0)
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

