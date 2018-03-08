pro plotsed, good, psname

restore, '/Users/jkrick/idlbin/objectnew.sav'

ps_open, filename=psname,/portrait,/square, xsize=4, ysize=4;,/color

!P.multi=[0,1,1]
;-----------------------------------------------------------------
;what do these stars look like in ACS, K band, and Irac channel 1?

;
acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits');
adxy, acshead, objectnew[good].ra,objectnew[good].dec , xcenter, ycenter

fits_read, '/Users/jkrick/palomar/wirc/coadd_k.fits', kdata, khead
adxy, khead, objectnew[good].ra,objectnew[good].dec , kxcenter, kycenter

fits_read, '/Users/jkrick/spitzer/irac/ch1/mosaic.fits', irac1data, irac1head
adxy, irac1head, objectnew[good].ra,objectnew[good].dec , irac1xcenter, irac1ycenter



size = 10
x = [.3540,.4660,.6255,.7680,.8330,1.196, 1.25,1.24,1.635,1.65,2.15,2.16,3.6000,4.5000,5.8,8.]    ;wavelengths of our photometry in microns
y = x
yerr = x
for n = 0,n_elements(good) - 1 do begin

;output the acs image of the objectnew
   print, "working on n ", n

   if objectnew[good(n)].acsmag gt 0 and objectnew[good(n)].acsmag lt 90 then begin
;     print,    xcenter(n) - size/0.05, xcenter(n)+ size/0.05  ,ycenter(n) -size/0.05, ycenter(n)+size/0.05
      acsdata = mrdfits('/Users/jkrick/hst/raw/wholeacs.fits', range=[ycenter(n) -size/0.05, ycenter(n)+size/0.05])
      plotimage, xrange=[xcenter(n) - size/0.05, xcenter(n)+ size/0.05],$
;                 yrange=[ycenter(n) -size/0.05, ycenter(n)+size/0.05], $
                 bytscl(acsdata, min = -0.01, max = 0.1),$
                 /preserve_aspect, /noaxes, ncolors=60
;      xyouts, xcenter(n)- 0.6*size/0.05, -10., strcompress(string(good(n))+ '    ACS' + string(objectnew[good(n)].acsmag) + string( (objectnew[good(n)].wirckmag) - (objectnew[good(n)].irac2mag - 3.25))),charthick = 3
;      xyouts, xcenter(n)- 0.6*size/0.05, -10., string(good(n)),charthick = 3
   endif

;output the k band image

   if objectnew[good(n)].wirckmag gt 0 and objectnew[good(n)].wirckmag lt 90 then begin
;     print,    kxcenter(n) - size/0.25, kxcenter(n)+ size/0.25  ,kycenter(n) -size/0.25, kycenter(n)+size/0.25
      plotimage, xrange=[kxcenter(n) - size/0.25, kxcenter(n)+ size/0.25],$
                 yrange=[kycenter(n) -size/0.25, kycenter(n)+size/0.25], $
                 bytscl(kdata, min = -21, max = 21),$
                 /preserve_aspect, /noaxes, ncolors=60
   endif

;output the ch1 image

   if objectnew[good(n)].wirckmag gt 0 and objectnew[good(n)].wirckmag lt 90 then begin
;     print,    irac1xcenter(n) - size/0.6, irac1xcenter(n)+ size/0.6  ,irac1ycenter(n) -size/0.6, irac1ycenter(n)+size/0.6
      plotimage, xrange=[irac1xcenter(n) - size/0.6, irac1xcenter(n)+ size/0.6],$
                 yrange=[irac1ycenter(n) -size/0.6, irac1ycenter(n)+size/0.6], $
                 bytscl(irac1data, min = 0., max = 0.3),$
                 /preserve_aspect, /noaxes, ncolors=60
   endif


;-----------------------------------------------------------------
;now lets  look at some SED's.
;want all to be in fluxes in microjanskies
   
;flamingos and wirc data is in vega in the catalog, so first convert that to AB, even though that means converting through flux twice
  if objectnew[good(n)].flamjmag gt 0 and objectnew[good(n)].flamjmag ne 99 then begin
      fab = 1594.*10^(objectnew[good(n)].flamjmag/(-2.5))
      jmagab = -2.5*alog10(fab) +8.926
   endif else begin
      jmagab = objectnew[good(n)].flamjmag
   endelse

   if objectnew[good(n)].wircjmag gt 0 and objectnew[good(n)].wircjmag ne 99 then begin
      wircjab = 1594.*10^(objectnew[good(n)].wircjmag/(-2.5))
      wircjmagab = -2.5*alog10(wircjab) +8.926
   endif else begin
      wircjmagab = objectnew[good(n)].wircjmag
   endelse

   if objectnew[good(n)].wirchmag gt 0 and objectnew[good(n)].wirchmag ne 99 then begin
      wirchab = 1024.*10^(objectnew[good(n)].wirchmag/(-2.5))
      wirchmagab = -2.5*alog10(wirchab) +8.926
   endif else begin
      wirchmagab = objectnew[good(n)].wirchmag
   endelse

   if objectnew[good(n)].wirckmag gt 0 and objectnew[good(n)].wirckmag ne 99 then begin
      wirckab = 666.8*10^(objectnew[good(n)].wirckmag/(-2.5))
      wirckmagab = -2.5*alog10(wirckab) +8.926
   endif else begin
      wirckmagab = objectnew[good(n)].wirckmag
   endelse


   yab = [objectnew[good(n)].umag, objectnew[good(n)].gmag,objectnew[good(n)].rmag,objectnew[good(n)].imag,objectnew[good(n)].acsmag,jmagab,wircjmagab,objectnew[good(n)].tmassjmag,wirchmagab,objectnew[good(n)].tmasshmag,wirckmagab,objectnew[good(n)].tmasskmag,objectnew[good(n)].irac1mag,objectnew[good(n)].irac2mag,objectnew[good(n)].irac3mag,objectnew[good(n)].irac4mag,objectnew[good(n)].mips24mag,objectnew[good(n)].mips70mag ]

   yflux = 1E6*10^((yab - 8.926)/(-2.5))

   abtovega = [0.94,-0.08,0.17,0.4,0.43,-0.92,-0.92,-0.92,-1.4,-1.4,-1.87,-1.87,-2.7,-3.25,-3.73,-4.37,0,0]
   yvega = yab + abtovega

   yaberr = [objectnew[good(n)].umagerr, objectnew[good(n)].gmagerr,objectnew[good(n)].rmagerr,objectnew[good(n)].imagerr,objectnew[good(n)].acsmagerr,objectnew[good(n)].flamjmagerr,objectnew[good(n)].wircjmagerr,0.02,objectnew[good(n)].wirchmagerr,0.02,objectnew[good(n)].wirckmagerr,0.02,0.05,0.05,0.05,0.05, 0.05,0.05]

   ypluserr = 1E6*10^(((yab + yaberr) - 8.926)/(-2.5))
   yminuserr = 1E6*10^(((yab - yaberr) - 8.926)/(-2.5))

   for count = 0, n_elements(yab) - 1 do begin
      if yab(count) gt 90. or yab(count) lt 0 then begin
         ypluserr(count) = 0
         yminuserr(count) = 0
      endif
   endfor


   plot,(x), alog10(yflux), psym = 5, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
        xtitle = "microns", ytitle = "log(flux(microjansky))", title=strcompress("objectnew " + string(good(n)) + string( (objectnew[good(n)].wirckmag) - (objectnew[good(n)].irac2mag - 3.25))), xrange=[0,10], yrange=[-1,2.5], ystyle = 1
   
   errplot, x, alog10(yminuserr), alog10(ypluserr)
   
;   for count = 0, 5 do xyouts, 0.8 + 1.4*(count), -0.7, yvega(count)
;   for count = 6, 11 do xyouts, 0.8 + 1.4*(count - 6), -0.8, yvega(count)
;   for count = 12, 17 do xyouts, 0.8 + 1.4*(count - 12), -0.9, yflux(count)


endfor

ps_close, /noprint,/noid


end
