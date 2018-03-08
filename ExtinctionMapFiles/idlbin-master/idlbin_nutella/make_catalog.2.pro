pro make_catalog


close, /all

numoobjects = 34000
colorarr = fltarr(numoobjects)
magarr = fltarr(numoobjects)
FITS_READ, '/Users/jkrick/palomar/lfc/coadd_g.fits',gdata, gheader
FITS_READ, '/Users/jkrick/spitzer/irac/ch1/mosaic.fits',irdata, irheader

robject = replicate({rob, rxcenter:0D,rycenter:0D,rfluxa:0D,rmaga:0D,rmagerra:0D,rfwhm:0D,risoarea:0D, rellip:0D},numoobjects)
gobject = replicate({gob, gxcenter:0D,gycenter:0D,gfluxa:0D,gmaga:0D,gmagerra:0D,gfwhm:0D,gisoarea:0D,gellip:0D,match:0D},numoobjects)
iobject = replicate({iob, ixcenter:0D,iycenter:0D,ifluxa:0D,imaga:0D,imagerra:0D,ifwhm:0D,iisoarea:0D, iellip:0D},numoobjects)
uobject = replicate({uob, uxcenter:0D,uycenter:0D,ufluxa:0D,umaga:0D,umagerra:0D,ufwhm:0D,uisoarea:0D,uellip:0D},numoobjects)
;-------
jobject = replicate({job, jxcenter:0D, jycenter:0D, jra:0D, jdec:0D, jmag:0D, jmagerr:0D, jfwhm:0D, jisoarea:0D, jellip:0D },numoobjects)

;------
mips24object = replicate({mips24ob, xcenter:0D, ycenter:0D, ra:0D, dec:0D, flux:0D, fluxerr:0D, bckgrnd:0D, iter:0D,sharp:0D},numoobjects)

iracobject = replicate({iracob, irxcenter:0D,irycenter:0D,tmassJ:0D,tmassh:0D,tmassk:0D,irac1:0D,irac2:0D,irac3:0D, irac4:0D,mips24:0D},numoobjects)

infraredobject = replicate({spitzob, INHERITS job, INHERITS iracob,INHERITS mips24ob},numoobjects)
;=======
object= replicate({ob, INHERITS uob, INHERITS gob, INHERITS rob, INHERITS iob, INHERITS spitzob, zphot:0D, chisq:0D, prob:0D, spt:0D, age:0D, av:0D, Mabs:0D}, numoobjects)


i = 0.
openr, lunr, "/Users/jkrick/palomar/LFC/catalog/SExtractor.r.cat", /get_lun
WHILE (NOT EOF(lunr)) DO BEGIN
    READF, lunr, o, xwcsr, ywcsr, xcenterr, ycenterr, fluxautor, magautor, magerrr, fwhmr, isoarear, fluxmaxr, ellipr

 ;   if (magautor LT 50)  then begin
       robject[i] ={rob, xwcsr, ywcsr,fluxautor, magautor, magerrr,fwhmr, isoarear, ellipr}
       i = i + 1
 ;   ENDIF

 endwhile
robject = robject[0:i - 1]
print, "there are ",i," r objects"
close, lunr
free_lun, lunr
;----------------------------------------------------------------

l = 0.
openr, lung, "/Users/jkrick/palomar/LFC/catalog/SExtractor.g.cat", /get_lun
WHILE (NOT EOF(lung)) DO BEGIN
    READF, lung, o, xwcsg, ywcsg, xcenterg, ycenterg, fluxautog, magautog, magerrg, fwhmg, isoareag, fluxmaxg, ellipg
;    if o eq 20389 then print, xwcsg, ywcsg, fluxmaxg
  ;  if (fluxautog gt 0 )   then begin
       gobject[l] ={gob, xwcsg, ywcsg,fluxautog, magautog, magerrg,fwhmg, isoareag, ellipg, -1}
;       if o eq 20389 then print, "found g object", l, gobject[l].gxcenter, gobject[l].gycenter
       l = l + 1
;    ENDIF

 endwhile
gobject = gobject[0:l - 1]
print, "there are ",l," g objects"
close, lung
free_lun, lung

;set them all as unmatched objects
;gobject.match = -1
   
;----------------------------------------------------------------

m = 0.
openr, lung, "/Users/jkrick/palomar/LFC/catalog/SExtractor.i.cat", /get_lun
WHILE (NOT EOF(lung)) DO BEGIN
    READF, lung, o, xwcsg, ywcsg, xcenterg, ycenterg, fluxautog, magautog, magerrg, fwhmg, isoareag, fluxmaxg, ellipi

;    if (fluxautog gt 0 )   then begin
       iobject[m] ={iob, xwcsg, ywcsg,fluxautog, magautog, magerrg,fwhmg, isoareag, ellipi}
       m = m + 1
 ;   ENDIF

 endwhile
iobject = iobject[0:m - 1]
print, "there are ",m," i objects"
close, lung
free_lun, lung
;----------------------------------------------------------------

n = 0.
openr, lung, "/Users/jkrick/palomar/LFC/catalog/SExtractor.u.cat", /get_lun
WHILE (NOT EOF(lung)) DO BEGIN
    READF, lung, o, xwcsg, ywcsg, xcenterg, ycenterg, fluxautog, magautog, magerrg, fwhmg, isoareag, fluxmaxg, ellipu

;    if (fluxautog gt 0 )   then begin
       uobject[n] ={uob, xwcsg, ywcsg,fluxautog, magautog, magerrg,fwhmg, isoareag, ellipu}
       n = n + 1
;    ENDIF

 endwhile
uobject = uobject[0:n - 1]
print, "there are ",n," u objects"
close, lung
free_lun, lung
;----------------------------------------------------------------

p = 0
openr, luni, "/Users/jkrick/palomar/LFC/catalog/jband_noheader.txt", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, junk,x_image, y_image,junk,junk,junk,junk,junk,junk,mag_best, magerr_best, junk,junk, junk,isoareaimage, ra,dec,junk,junk,junk, junk,fwhm_image, junk,ell, junk 

   jobject[p] ={job, x_image, y_image, ra, dec, mag_best, magerr_best, fwhm_image, isoareaimage, ell}
   p = p + 1
   
 endwhile
jobject = jobject[0:p - 1]
print, "there are ",p," jobjects"
close, luni
free_lun, luni
;----------------------------------------------------------------

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
w = 0

openr, luni, "/Users/jkrick/Spitzer/mips/mips24/dao/mips24.phot", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, num, xcen,ycen, mag, magerr, back, niter, sh,junk,junk
   ;need to convert from magnitudes back to fluxes.
   mips24object[w] = {mips24ob, xcen, ycen, 0.0,0.0, 1.432*(10^((25-mag)/2.5)),  10^((25-(mag -magerr))/2.5) -(10^((25-mag)/2.5)) , back, niter,sh}
   w = w + 1

endwhile
mips24object = mips24object[0:w - 1]
print, "there are ",w," mips24 objects"
close, luni
free_lun, luni

;need to go from x,y in pixels to ra, dec
FITS_READ, '/Users/jkrick/spitzer/mips/mips24/daophot/mosaic.fits',data, header
xyad, header, mips24object.xcenter, mips24object.ycenter, ra, dec
mips24object.ra = ra
mips24object.dec = dec

;now match mips to irac before I even deal with the optical.
;there shouldn't be any mips detections which are not detected in irac.
; create initial arrays
m=n_elements(mips24object.ra)
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

   dist=sphdist(iracobject[q].irxcenter,iracobject[q].irycenter, mips24object.ra,mips24object.dec,/degrees)
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
   
   infraredobject[num] = {spitzob,  iracobject(matched[num]).irxcenter, iracobject(matched[num]).irycenter,iracobject(matched[num]).tmassJ, iracobject(matched[num]).tmassh, iracobject(matched[num]).tmassk, iracobject(matched[num]).irac1, iracobject(matched[num]).irac2,iracobject(matched[num]).irac3, iracobject(matched[num]).irac4,iracobject(matched[num]).mips24, mips24object(irmatch[matched[num]]).xcenter , mips24object(irmatch[matched[num]]).ycenter , mips24object(irmatch[matched[num]]).ra , mips24object(irmatch[matched[num]]).dec , mips24object(irmatch[matched[num]]).flux , mips24object(irmatch[matched[num]]).fluxerr , mips24object(irmatch[matched[num]]).bckgrnd , mips24object(irmatch[matched[num]]).iter ,mips24object(irmatch[matched[num]]).sharp}

;print, "matched", num, infraredobject[num].irxcenter, infraredobject[num].irycenter, infraredobject[num].ra, infraredobject[num].dec

endfor
count = n_elements(matched)
for num= 0 , n_elements(nonmatched) -1 do begin
   
   infraredobject[count] = {spitzob,  iracobject(nonmatched[num]).irxcenter, iracobject(nonmatched[num]).irycenter,iracobject(nonmatched[num]).tmassJ,iracobject(nonmatched[num]).tmassh,iracobject(nonmatched[num]).tmassk,iracobject(nonmatched[num]).irac1,iracobject(nonmatched[num]).irac2,iracobject(nonmatched[num]).irac3, iracobject(nonmatched[num]).irac4,iracobject(nonmatched[num]).mips24, 0.,0.,0.,0.,0.,0.,0.,0.,0.}
count = count + 1

endfor

infraredobject = infraredobject[0:n_elements(matched) + n_elements(nonmatched)]
;test
;print, "testing infraredob, irac vs. mips24"
;print, infraredobject[100].irxcenter, infraredobject[17].irycenter, infraredobject[17].ra, infraredobject[17].dec
;print,  infraredobject[10000].irxcenter, infraredobject[10000].irycenter, infraredobject[10000].ra, infraredobject[10000].dec
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
   if (sep LT 0.0005) then begin
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
   
   object[num] = {ob,  uobject(irmatch[matched[num]]).uxcenter,uobject(irmatch[matched[num]]).uycenter,uobject(irmatch[matched[num]]).ufluxa,uobject(irmatch[matched[num]]).umaga,uobject(irmatch[matched[num]]).umagerra,uobject(irmatch[matched[num]]).ufwhm,uobject(irmatch[matched[num]]).uisoarea,uobject(irmatch[matched[num]]).uellip, gobject(irmatch[matched[num]]).gxcenter,gobject(irmatch[matched[num]]).gycenter,gobject(irmatch[matched[num]]).gfluxa,gobject(irmatch[matched[num]]).gmaga,gobject(irmatch[matched[num]]).gmagerra,gobject(irmatch[matched[num]]).gfwhm,gobject(irmatch[matched[num]]).gisoarea,gobject(irmatch[matched[num]]).gellip,gobject(irmatch[matched[num]]).match , robject(irmatch[matched[num]]).rxcenter,robject(irmatch[matched[num]]).rycenter,robject(irmatch[matched[num]]).rfluxa,robject(irmatch[matched[num]]).rmaga,robject(irmatch[matched[num]]).rmagerra,robject(irmatch[matched[num]]).rfwhm,robject(irmatch[matched[num]]).risoarea, robject(irmatch[matched[num]]).rellip, iobject(irmatch[matched[num]]).ixcenter,iobject(irmatch[matched[num]]).iycenter,iobject(irmatch[matched[num]]).ifluxa,iobject(irmatch[matched[num]]).imaga,iobject(irmatch[matched[num]]).imagerra,iobject(irmatch[matched[num]]).ifwhm,iobject(irmatch[matched[num]]).iisoarea, iobject(irmatch[matched[num]]).iellip, infraredobject(matched[num]).irxcenter,infraredobject(matched[num]).irycenter,infraredobject(matched[num]).tmassJ,infraredobject(matched[num]).tmassh,infraredobject(matched[num]).tmassk,infraredobject(matched[num]).irac1,infraredobject(matched[num]).irac2,infraredobject(matched[num]).irac3, infraredobject(matched[num]).irac4,infraredobject(matched[num]).mips24, infraredobject(matched[num]).xcenter , infraredobject(matched[num]).ycenter , infraredobject(matched[num]).ra , infraredobject(matched[num]).dec , infraredobject(matched[num]).flux , infraredobject(matched[num]).fluxerr , infraredobject(matched[num]).bckgrnd , infraredobject(matched[num]).iter ,infraredobject(matched[num]).sharp,0.,0.,0.,0.,0.,0.,0.}

 

endfor

;print the ir objects which don't have an optical match 
for num2=0, n_elements(nonmatched) - 1 do begin
   
   adxy, gheader, infraredobject(nonmatched[num2]).irxcenter,infraredobject(nonmatched[num2]).irycenter, xcen, ycen
   printf,outlun3, xcen,ycen
   adxy, irheader, infraredobject(nonmatched[num2]).irxcenter,infraredobject(nonmatched[num2]).irycenter, xcen2, ycen2
   printf,outlun4, xcen2,ycen2


   object[num]={ob, -1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,infraredobject(nonmatched[num2]).irxcenter ,infraredobject(nonmatched[num2]).irycenter ,infraredobject(nonmatched[num2]).tmassJ ,infraredobject(nonmatched[num2]).tmassh ,infraredobject(nonmatched[num2]).tmassk ,infraredobject(nonmatched[num2]).irac1 ,infraredobject(nonmatched[num2]).irac2 ,infraredobject(nonmatched[num2]).irac3 , infraredobject(nonmatched[num2]).irac4,infraredobject(nonmatched[num2]).mips24, infraredobject(nonmatched[num2]).xcenter , infraredobject(nonmatched[num2]).ycenter , infraredobject(nonmatched[num2]).ra , infraredobject(nonmatched[num2]).dec , infraredobject(nonmatched[num2]).flux , infraredobject(nonmatched[num2]).fluxerr , infraredobject(nonmatched[num2]).bckgrnd , infraredobject(nonmatched[num2]).iter ,infraredobject(nonmatched[num2]).sharp,0.,0.,0.,0.,0.,0.,0.}
 
   num = num+1

endfor


;;print the optical objects which are detected in 2 optical bands and not in the IR catalog

for num3 = 0.D, n_elements(gobject.match) -1 do begin
   if gobject[num3].match lt 0 then begin
      ;this object did not match to an irac object
      ;is it a real object?

      if (gobject[num3].gfluxa GT 0) and (gobject[num3].gmaga lt 30)$  
         and (gobject[num3].gmagerra lt 2) $
         and (robject[num3].rfluxa GT 0) and (robject[num3].rmaga lt 30)$  
         and (robject[num3].rmagerra lt 2) $ and (iobject[num3].ifluxa GT 0) and (iobject[num3].imaga lt 30)$  
         and (iobject[num3].imagerra lt 2) $
      then begin
         ;it is both in r and g catalogs and apparently not noise
         ;so add it to the catalog
          adxy, gheader, gobject[num3].gxcenter,gobject[num3].gycenter, xcen, ycen
 ;       printf,outlun3, xcen,ycen
          adxy, irheader, gobject[num3].gxcenter,gobject[num3].gycenter, xcen2, ycen2
 ;       printf,outlun4, xcen2,ycen2

         object[num]={ob, uobject[num3].uxcenter,uobject[num3].uycenter,uobject[num3].ufluxa ,uobject[num3].umaga ,uobject[num3].umagerra ,uobject[num3].ufwhm ,uobject[num3].uisoarea ,uobject[num3].uellip ,gobject[num3].gxcenter,gobject[num3].gycenter,gobject[num3].gfluxa ,gobject[num3].gmaga ,gobject[num3].gmagerra ,gobject[num3].gfwhm ,gobject[num3].gisoarea ,gobject[num3].gellip, gobject[num3].match,robject[num3].rxcenter,robject[num3].rycenter,robject[num3].rfluxa ,robject[num3].rmaga ,robject[num3].rmagerra ,robject[num3].rfwhm ,robject[num3].risoarea ,robject[num3].rellip, iobject[num3].ixcenter,iobject[num3].iycenter,iobject[num3].ifluxa ,iobject[num3].imaga ,iobject[num3].imagerra ,iobject[num3].ifwhm ,iobject[num3].iisoarea ,iobject[num3].iellip, -1, -1, -1, -1, -1, -1, -1, -1, -1,-1,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,0.,0.,0.,0.,0.,0.,0.}


         num = num + 1
      endif
   endif
endfor



;----------------------------------------------------------
object = object[0:num - 1]



save, object, filename='object.sav'


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
end


