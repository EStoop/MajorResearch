pro browndwarfs
close, /all
help, /memory
;device, true_color=24
;device, decomposed=0


;colors = GetColor(/load, Start=1)

;AB - Vega
;ch 1 = 2.7
;ch2 = 3.25
;ch3 = 3.73
;ch4 = 4.37

ps_open, filename='/Users/jkrick/nep/bdwarfs3.ps',/portrait,/square;,/color

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

!P.thick = 3
!P.charthick = 3

restore, '/Users/jkrick/idlbin/object.sav'

notstars = [898,922,1181,1658,5990,6082,6239,6881,7519,7636,8164]
object[notstars].acsclassstar = 0.
object[notstars].uclassstar = 0.
object[notstars].gclassstar = 0.
object[notstars].rclassstar = 0.
object[notstars].iclassstar = 0.

;-----------------------------------------------------------------
; first try to recreate plots of Patten et al with only the stars.  
;be wary of flagged items
point8 = where(object.acsclassstar eq 1 or object.uclassstar eq 1 or object.gclassstar eq 1 or object.rclassstar eq 1 or object.iclassstar eq 1 and object.irac1mag gt 0 and object.irac1mag lt 90 and object.irac2mag gt 0 and object.irac2mag lt 90 and object.irac3mag gt 0 and object.irac3mag lt 90 and object.irac4mag gt 0 and object.irac4mag lt 90  and object.acsflag eq 0 and object.uflags eq 0 and object.gflags eq 0 and object.rflags eq 0 and object.iflags eq 0)

print, "point8", n_elements(point8)

;-----------------------------------------------------------------
;then loosen the criterion to be anything with some ellip and fwhm (let go of classstar)

relax = where(object.irac1mag gt 0 and object.irac1mag lt 90 and object.irac2mag gt 0 and object.irac2mag lt 90 and object.irac3mag gt 0 and object.irac3mag lt 90 and object.acsflag eq 0 and object.uflags eq 0 and object.gflags eq 0 and object.rflags eq 0 and object.iflags eq 0)
print, "relax", n_elements(relax)

;-----------------------------------------------------------------

;bring in near-IR K band

withk = where(object.acsmag gt 0 and object.acsmag lt 90 and object.acsclassstar eq 1 and object.irac1mag gt 0 and object.irac1mag lt 90 and object.irac2mag gt 0 and object.irac2mag lt 90 and object.wirckmag gt 0. and object.wirckmag lt 90. )

print, "withk", n_elements(withk)

good = where(  (object[withk].wirckmag) - (object[withk].irac2mag-3.25) gt 0.8  )
print, "good", n_elements(good)


;-----------------------------------------------------------------

plot, (object[withk].irac1mag-2.7) - (object[withk].irac2mag-3.25), (object[withk].wirckmag) - (object[withk].irac2mag - 3.25), psym = 2, xrange=[-0.5, 2.0], yrange=[0, 4.0], xtitle = '[3.6]-[4.5]', ytitle='K-[4.5]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1

oplot, findgen(20)/ 10, 1.05*(findgen(20)/10) + 1.2

oplot, (object[withk[good]].irac1mag-2.7) - (object[withk[good]].irac2mag-3.25), (object[withk[good]].wirckmag) - (object[withk[good]].irac2mag - 3.25), psym = 2,color = redcolor


plot, (object[withk].irac1mag-2.7) - (object[withk].irac2mag-3.25), (object[withk].wirckmag) - (object[withk].irac1mag - 2.7), psym = 2, xrange=[-0.5, 2.0], yrange=[0, 1.8], xtitle = '[3.6]-[4.5]', ytitle='K-[3.6]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1

oplot, (object[withk[good]].irac1mag-2.7) - (object[withk[good]].irac2mag-3.25), (object[withk[good]].wirckmag) - (object[withk[good]].irac1mag - 2.7), psym = 2, color = redcolor


;-----------------------------------------------------------------

;now if these are our candidates, where do they fall in the other color color plots?
;fortunately they are all detected in ch3 and ch 4

plot, (object[withk[good]].irac2mag-3.25) - (object[withk[good]].irac3mag-3.73), (object[withk[good]].irac1mag - 2.7) - (object[withk[good]].irac2mag - 3.25), psym = 2, xrange=[-0.8, 0.6], yrange=[-0.5, 2.5], xtitle = '[4.5]-[5.8]', ytitle='[3.6]-[4.5]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1

plot, (object[withk[good]].irac3mag-3.73) - (object[withk[good]].irac4mag-4.37), (object[withk[good]].irac1mag - 2.7) - (object[withk[good]].irac2mag - 3.25), psym = 2, xrange=[-0.2, 1.0], yrange=[-0.5, 2.5], xtitle = '[5.8]-[8.0]', ytitle='[3.6]-[4.5]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1

plot, (object[withk[good]].irac1mag-2.7) - (object[withk[good]].irac4mag-4.37), (object[withk[good]].irac2mag - 3.25) - (object[withk[good]].irac3mag - 3.73), psym = 2, xrange=[-0.2, 2.3], yrange=[-0.7, 0.6], xtitle = '[3.6]-[8.0]', ytitle='[4.5]-[5.8]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1





;-----------------------------------------------------------------
;what do these stars look like in ACS, K band, and Irac channel 1?

;
acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits');
adxy, acshead, object[withk[good]].ra,object[withk[good]].dec , xcenter, ycenter

fits_read, '/Users/jkrick/palomar/wirc/coadd_k.fits', kdata, khead
adxy, khead, object[withk[good]].ra,object[withk[good]].dec , kxcenter, kycenter

fits_read, '/Users/jkrick/spitzer/irac/ch1/mosaic.fits', irac1data, irac1head
adxy, irac1head, object[withk[good]].ra,object[withk[good]].dec , irac1xcenter, irac1ycenter



size = 10
x = [.3540,.4660,.6255,.7680,.8330,1.196, 1.25,1.24,1.635,1.65,2.15,2.16,3.6000,4.5000,5.8,8.]    ;wavelengths of our photometry in microns
y = x
yerr = x
openw, outlun, '/Users/jkrick/nep/bdwarfs.reg', /get_lun
for n = 0,n_elements(good) - 1 do begin

;output the acs image of the object
   print, "working on n ", n

   if object[withk[good(n)]].acsmag gt 0 and object[withk[good(n)]].acsmag lt 90 then begin
;     print,    xcenter(n) - size/0.05, xcenter(n)+ size/0.05  ,ycenter(n) -size/0.05, ycenter(n)+size/0.05
      acsdata = mrdfits('/Users/jkrick/hst/raw/wholeacs.fits', range=[ycenter(n) -size/0.05, ycenter(n)+size/0.05])
      plotimage, xrange=[xcenter(n) - size/0.05, xcenter(n)+ size/0.05],$
;                 yrange=[ycenter(n) -size/0.05, ycenter(n)+size/0.05], $
                 bytscl(acsdata, min = -0.01, max = 0.1),$
                 /preserve_aspect, /noaxes, ncolors=60
;      xyouts, xcenter(n)- 0.6*size/0.05, -10., strcompress(string(withk[good(n)])+ '    ACS' + string(object[withk[good(n)]].acsmag) + string( (object[withk[good(n)]].wirckmag) - (object[withk[good(n)]].irac2mag - 3.25))),charthick = 3
      xyouts, xcenter(n)- 0.6*size/0.05, -10., string(withk[good(n)]),charthick = 3
   endif
   printf, outlun, 'circle(',xcenter(n), ycenter(n),' 300)'
;output the k band image

   if object[withk[good(n)]].wirckmag gt 0 and object[withk[good(n)]].wirckmag lt 90 then begin
;     print,    kxcenter(n) - size/0.25, kxcenter(n)+ size/0.25  ,kycenter(n) -size/0.25, kycenter(n)+size/0.25
      plotimage, xrange=[kxcenter(n) - size/0.25, kxcenter(n)+ size/0.25],$
                 yrange=[kycenter(n) -size/0.25, kycenter(n)+size/0.25], $
                 bytscl(kdata, min = -21, max = 21),$
                 /preserve_aspect, /noaxes, ncolors=60
   endif

;output the ch1 image

   if object[withk[good(n)]].wirckmag gt 0 and object[withk[good(n)]].wirckmag lt 90 then begin
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
  if object[withk[good(n)]].flamjmag gt 0 and object[withk[good(n)]].flamjmag ne 99 then begin
      fab = 1594.*10^(object[withk[good(n)]].flamjmag/(-2.5))
      jmagab = -2.5*alog10(fab) +8.926
   endif else begin
      jmagab = object[withk[good(n)]].flamjmag
   endelse

   if object[withk[good(n)]].wircjmag gt 0 and object[withk[good(n)]].wircjmag ne 99 then begin
      wircjab = 1594.*10^(object[withk[good(n)]].wircjmag/(-2.5))
      wircjmagab = -2.5*alog10(wircjab) +8.926
   endif else begin
      wircjmagab = object[withk[good(n)]].wircjmag
   endelse

   if object[withk[good(n)]].wirchmag gt 0 and object[withk[good(n)]].wirchmag ne 99 then begin
      wirchab = 1024.*10^(object[withk[good(n)]].wirchmag/(-2.5))
      wirchmagab = -2.5*alog10(wirchab) +8.926
   endif else begin
      wirchmagab = object[withk[good(n)]].wirchmag
   endelse

   if object[withk[good(n)]].wirckmag gt 0 and object[withk[good(n)]].wirckmag ne 99 then begin
      wirckab = 666.8*10^(object[withk[good(n)]].wirckmag/(-2.5))
      wirckmagab = -2.5*alog10(wirckab) +8.926
   endif else begin
      wirckmagab = object[withk[good(n)]].wirckmag
   endelse


   yab = [object[withk[good(n)]].umaga, object[withk[good(n)]].gmaga,object[withk[good(n)]].rmaga,object[withk[good(n)]].imaga,object[withk[good(n)]].acsmag,jmagab,wircjmagab,object[withk[good(n)]].tmassjmag,wirchmagab,object[withk[good(n)]].tmasshmag,wirckmagab,object[withk[good(n)]].tmasskmag,object[withk[good(n)]].irac1mag,object[withk[good(n)]].irac2mag,object[withk[good(n)]].irac3mag,object[withk[good(n)]].irac4mag,object[withk[good(n)]].mips24mag,object[withk[good(n)]].mips70mag ]

   yflux = 1E6*10^((yab - 8.926)/(-2.5))

   yaberr = [object[withk[good(n)]].umagerra, object[withk[good(n)]].gmagerra,object[withk[good(n)]].rmagerra,object[withk[good(n)]].imagerra,object[withk[good(n)]].acsmagerr,object[withk[good(n)]].flamjmagerr,object[withk[good(n)]].wircjmagerr,0.02,object[withk[good(n)]].wirchmagerr,0.02,object[withk[good(n)]].wirckmagerr,0.02,0.05,0.05,0.05,0.05, 0.05,0.05]

   ypluserr = 1E6*10^(((yab + yaberr) - 8.926)/(-2.5))
   yminuserr = 1E6*10^(((yab - yaberr) - 8.926)/(-2.5))

   for count = 0, n_elements(yab) - 1 do begin
      if yab(count) gt 90. or yab(count) lt 0 then begin
         ypluserr(count) = 0
         yminuserr(count) = 0
      endif
   endfor


   plot,(x), alog10(yflux), psym = 5, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
        xtitle = "microns", ytitle = "log(flux(microjansky))", title=strcompress("object " + string(withk[good(n)]) + string( (object[withk[good(n)]].wirckmag) - (object[withk[good(n)]].irac2mag - 3.25))), xrange=[0,10], yrange=[-1,2.5], ystyle = 1
   
   errplot, x, alog10(yminuserr), alog10(ypluserr)
   
   for count = 0, 5 do xyouts, 0.8 + 1.4*(count), -0.7, yab(count)
   for count = 6, 11 do xyouts, 0.8 + 1.4*(count - 6), -0.8, yab(count)
   for count = 12, 17 do xyouts, 0.8 + 1.4*(count - 12), -0.9, yab(count)


endfor

close, outlun
free_lun, outlun

;-----------------------------------------------------------------

;try selecting only based on acs detection, acs starlike, and ch1 and ch2 detection

withoutk = where( object.acsmag gt 0 and object.acsmag lt 90 and object.acsfwhm lt 3.2 and object.acsellip lt 0.19 and object.irac1mag gt 0 and object.irac1mag lt 90 and object.irac2mag gt 0 and object.irac2mag lt 90 )

print, "withoutk", n_elements(withoutk)

good2 = where(  (object[withoutk].irac1mag - 2.7 ) - (object[withoutk].irac2mag-3.25) gt 0.2  )
print, "good2", n_elements(good2)

;-----------------------------------------------------------------
;how do they fit on those plots again
plot, (object[withoutk].irac2mag-3.25) - (object[withoutk].irac3mag-3.73), (object[withoutk].irac1mag - 2.7) - (object[withoutk].irac2mag - 3.25), psym = 2, xrange=[-0.8, 0.6], yrange=[-0.5, 2.5], xtitle = '[4.5]-[5.8]', ytitle='[3.6]-[4.5]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1

oplot, (object[withoutk[good2]].irac2mag-3.25) - (object[withoutk[good2]].irac3mag-3.73), (object[withoutk[good2]].irac1mag - 2.7) - (object[withoutk[good2]].irac2mag - 3.25), psym = 2, color = redcolor

;where are they on the acs image

adxy, acshead, object[withoutk[good2]].ra, object[withoutk[good2]].dec, xcenter3, ycenter3

openw, outlun, '/users/jkrick/nep/bdwarfs3.reg', /get_lun
for count = 0, n_elements(good2) - 1 do printf, outlun, 'circle( ', xcenter3(count), ycenter3(count), ' 100)'
close, outlun
free_lun, outlun 


;----------------------------------------------------------------
; try selecting objects detected in ch2, acs, and are acs starlike but not detected in ch1

good3 = where(object.irac2mag gt 0 and object.irac2mag lt 90 and object.irac1mag lt 0 and object.acsmag gt 0 and object.acsmag lt 90 and object.acsfwhm lt 3.2 and object.acsellip lt 0.19)

good4 = where(object.irac2mag gt 0 and object.irac2mag lt 90 and object.irac1mag gt 90 and object.acsmag gt 0 and object.acsmag lt 90 and object.acsfwhm lt 3.2 and object.acsellip lt 0.19)
print, good3, good4
print, 'good3', n_elements(good3) + n_elements(good4)

;adxy, acshead, object[good3].ra, object[good3].dec, xcenter3, ycenter3

;openw, outlun, '/users/jkrick/nep/bdwarfs4.reg', /get_lun
;for count = 0, n_elements(good3) - 1 do printf, outlun, 'circle( ', xcenter3(count), ycenter3(count), ' 100)'
;close, outlun
;free_lun, outlun 

a = where(object.irac2mag lt object.irac1mag and object.irac1mag gt 0 and object.irac1mag lt 90 and object.irac2mag gt 0 and object.irac2mag lt 90 and object.wirckmag gt 0 and object.wirckmag lt 90 )
print, n_elements(a)

b = where(object[a].wirckmag - object[a].irac2mag + 3.25 gt 1.0)
print, n_elements(b)

ps_close, /noprint,/noid

end
