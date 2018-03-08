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

ps_open, filename='/Users/jkrick/nep/bdwarfs.ps',/portrait,/square;,/color

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

print, object[withk[good]].irac3mag
print, object[withk[good]].irac4mag

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
;now lets  look at some SED's.
;want all to be in fluxes in microjanskies

uab = 10^((umaga - 8.926)/(-2.5))


;-----------------------------------------------------------------

openw, outlunh, "/Users/jkrick/ZPHOT/hyperz_cat_bdwarf.txt",/get_lun
for num = 0, n_elements(good) - 1 do begin

   if object[withk[good(num)]].flamjmag gt 0 and object[withk[good(num)]].flamjmag ne 99 then begin
      fab = 1594.*10^(object[withk[good(num)]].flamjmag/(-2.5))
      jmagab = -2.5*alog10(fab) +8.926
   endif else begin
      jmagab = object[withk[good(num)]].flamjmag
   endelse

   if object[withk[good(num)]].wircjmag gt 0 and object[withk[good(num)]].wircjmag ne 99 then begin
      wircjab = 1594.*10^(object[withk[good(num)]].wircjmag/(-2.5))
      wircjmagab = -2.5*alog10(wircjab) +8.926
   endif else begin
      wircjmagab = object[withk[good(num)]].wircjmag
   endelse

   if object[withk[good(num)]].wirchmag gt 0 and object[withk[good(num)]].wirchmag ne 99 then begin
      wirchab = 1024.*10^(object[withk[good(num)]].wirchmag/(-2.5))
      wirchmagab = -2.5*alog10(wirchab) +8.926
   endif else begin
      wirchmagab = object[withk[good(num)]].wirchmag
   endelse

   if object[withk[good(num)]].wirckmag gt 0 and object[withk[good(num)]].wirckmag ne 99 then begin
      wirckab = 666.8*10^(object[withk[good(num)]].wirckmag/(-2.5))
      wirckmagab = -2.5*alog10(wirckab) +8.926
   endif else begin
      wirckmagab = object[withk[good(num)]].wirckmag
   endelse

   if object[withk[good(num)]].irac1flux lt 0 then err1 = -1. else err1 = 0.05
   if object[withk[good(num)]].irac2flux lt 0 then err2 = -1. else err2 =  0.05;8.926 - 2.5*alog10(0.05E-6*object[withk[good(num)]].irac2)
   if object[withk[good(num)]].irac3flux lt 0 then err3 = -1. else err3 =  0.05;8.926 - 2.5*alog10(0.05E-6*object[withk[good(num)]].irac3)
   if object[withk[good(num)]].irac4flux lt 0 then err4 = -1. else err4 =  0.05;8.926 - 2.5*alog10(0.05E-6*object[withk[good(num)]].irac4)
   
   if object[withk[good(num)]].imagerra gt 1000. then object[withk[good(num)]].imagerra = 1000.
   if object[withk[good(num)]].gmagerra gt 1000. then object[withk[good(num)]].gmagerra = 1000.
   if object[withk[good(num)]].rmagerra gt 1000. then object[withk[good(num)]].rmagerra = 1000.
   if object[withk[good(num)]].umagerra gt 1000. then object[withk[good(num)]].umagerra = 1000.
   

    printf, outlunh, format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',$
                 withk[good(num)], object[withk[good(num)]].umaga, object[withk[good(num)]].gmaga, object[withk[good(num)]].rmaga, $
                 object[withk[good(num)]].imaga,  object[withk[good(num)]].acsmag,  jmagab, wircjmagab, $
                 object[withk[good(num)]].tmassjmag,  wirchmagab, object[withk[good(num)]].tmasshmag, wirckmagab, $
                 object[withk[good(num)]].tmasskmag, object[withk[good(num)]].irac1mag,object[withk[good(num)]].irac2mag,$
                 object[withk[good(num)]].irac3mag,object[withk[good(num)]].irac4mag,   $
                 object[withk[good(num)]].umagerra, object[withk[good(num)]].gmagerra, $
                 object[withk[good(num)]].rmagerra, object[withk[good(num)]].imagerra, object[withk[good(num)]].acsmagerr, $
            object[withk[good(num)]].flamjmagerr, $
                 object[withk[good(num)]].wircjmagerr,0.02, object[withk[good(num)]].wirchmagerr, 0.02, object[withk[good(num)]].wirckmagerr,$
                 0.02,err1,err2,err3,err4
 endfor
close, outlunh
free_lun, outlunh

;run hyperz

readcol,'/Users/jkrick/ZPHOT/hyperz_swire_bdwarf.z_phot',idz, zphota, chia, proba, specta, nagea, agea, $
        ava, ba, zinf99a,zsup99a,zinf90a,zsup90a,zinf68a,zsup68a,zwma,probwma,Mabsa,$
        zphot2a,prob2a,format="A"
readlargecol, '/Users/jkrick/ZPHOT/hyperz_swire_bdwarf.obs_sed',idz, u,g,r,i,acs,j,wj,tj,wh,th,wk,tk,one,two,three, four,$
         uerr, gerr, rerr, ierr,acserr, jerr, wjerr,tjerr,wherr,therr,wkerr,tkerr,oneerr, twoerr, threeerr, fourerr,format="A",/debug

;spectral type key
sptkey=[" Arp220 " ," Ell13  " ," Ell2  " ," I19254 " ," M82  " ," Mrk231" ," QSO1 " ," QSO2 " ," S0   " ," Sa  " ," Sb  " ," Sc  " ," Sd  " ," Sey2  " ," Torus" ]     

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

interest =indgen(n_elements(idz)) 
for n = 0,n_elements(idz) - 1 do begin

;output the acs image of the object
   print, "working on n ", n

   if object[withk[good(n)]].acsmag gt 0 and object[withk[good(n)]].acsmag lt 90 then begin
     print,    xcenter(n) - size/0.05, xcenter(n)+ size/0.05  ,ycenter(n) -size/0.05, ycenter(n)+size/0.05
      acsdata = mrdfits('/Users/jkrick/hst/raw/wholeacs.fits', range=[ycenter(n) -size/0.05, ycenter(n)+size/0.05])
      plotimage, xrange=[xcenter(n) - size/0.05, xcenter(n)+ size/0.05],$
;                 yrange=[ycenter(n) -size/0.05, ycenter(n)+size/0.05], $
                 bytscl(acsdata, min = -0.01, max = 0.1),$
                 /preserve_aspect, /noaxes, ncolors=60
      xyouts, xcenter(n)- 0.6*size/0.05, -10., strcompress(string(idz(n))+ '    ACS' + string(object[idz(n)].acsmag) + string( (object[withk[good(n)]].wirckmag) - (object[withk[good(n)]].irac2mag - 3.25))),charthick = 3
   endif

;output the k band image

   if object[withk[good(n)]].wirckmag gt 0 and object[withk[good(n)]].wirckmag lt 90 then begin
     print,    kxcenter(n) - size/0.25, kxcenter(n)+ size/0.25  ,kycenter(n) -size/0.25, kycenter(n)+size/0.25
      plotimage, xrange=[kxcenter(n) - size/0.25, kxcenter(n)+ size/0.25],$
                 yrange=[kycenter(n) -size/0.25, kycenter(n)+size/0.25], $
                 bytscl(kdata, min = -21, max = 21),$
                 /preserve_aspect, /noaxes, ncolors=60
      xyouts, kxcenter(n)- 0.6*size/0.25, -10., strcompress(string(idz(n))+ '    K' + string(object[idz(n)].wirckmag) + string( (object[withk[good(n)]].wirckmag) - (object[withk[good(n)]].irac2mag - 3.25))),charthick = 3
   endif

;output the ch1 image

   if object[withk[good(n)]].wirckmag gt 0 and object[withk[good(n)]].wirckmag lt 90 then begin
     print,    irac1xcenter(n) - size/0.6, irac1xcenter(n)+ size/0.6  ,irac1ycenter(n) -size/0.6, irac1ycenter(n)+size/0.6
      plotimage, xrange=[irac1xcenter(n) - size/0.6, irac1xcenter(n)+ size/0.6],$
                 yrange=[irac1ycenter(n) -size/0.6, irac1ycenter(n)+size/0.6], $
                 bytscl(irac1data, min = 0., max = 0.3),$
                 /preserve_aspect, /noaxes, ncolors=60
      xyouts, irac1xcenter(n)- 0.6*size/0.6, -10., strcompress(string(idz(n))+ '    IRAC1' + string(object[idz(n)].irac1mag) + string( (object[withk[good(n)]].wirckmag) - (object[withk[good(n)]].irac2mag - 3.25))),charthick = 3
   endif


;plot the fitted SED along with photometry of an individual galaxy

   y = [u(interest(n)),g(interest(n)),r(interest(n)),i(interest(n)),acs(interest(n)), j(interest(n)),wj(interest(n)),tj(interest(n)),wh(interest(n)),th(interest(n)),wk(interest(n)),tk(interest(n)),one(interest(n)),two(interest(n)),three(interest(n)), four(interest(n))] ;photometry
   yerr = [uerr(interest(n)),gerr(interest(n)),rerr(interest(n)),ierr(interest(n)),acserr(interest(n)),jerr(interest(n)),wjerr(interest(n)),tjerr(interest(n)), wherr(interest(n)),therr(interest(n)),wkerr(interest(n)), tkerr(interest(n)),oneerr(interest(n)),twoerr(interest(n)),threeerr(interest(n)), fourerr(interest(n))]

   plot,(x), alog10(y), psym = 5, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
        xtitle = "microns", ytitle = "log(flux(microjansky))", title=strcompress("object " + string(idz(n)) +string(object[idz(n)].wirckmag) ), xrange=[0,10]
   
   errplot, x, alog10(y - yerr), alog10(y + yerr)
;
   readcol,strcompress('/Users/jkrick/ZPHOT/' + string(idz(n)) + '.spe', /remove_all),x2, y2,format="A",/silent
;   oplot, alog10((x2/1E4)), alog10(y2), thick = 3
   
endfor


ps_close, /noprint,/noid

end
