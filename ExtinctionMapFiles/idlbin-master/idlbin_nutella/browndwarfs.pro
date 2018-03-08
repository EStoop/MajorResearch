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

withk = where(object.acsclassstar eq 1 or object.uclassstar eq 1 or object.gclassstar eq 1 or object.rclassstar eq 1 or object.iclassstar eq 1 and object.irac1mag gt 0 and object.irac1mag lt 90 and object.irac2mag gt 0 and object.irac2mag lt 90  and object.acsflag eq 0 and object.uflags eq 0 and object.gflags eq 0 and object.rflags eq 0 and object.iflags eq 0 and object.wirckmag gt 0. and object.wirckmag lt 90. )

print, "withk", n_elements(withk)

good = where(  (object[withk].wirckmag) - (object[withk].irac2mag-3.25) gt 0.8  )
print, "good", n_elements(good)
print, object[withk[good]].wirckmag

relaxk = where(object.acsmag gt 0 and object.acsmag lt 90 and object.acsclassstar eq 1 and object.irac1mag gt 0 and object.irac1mag lt 90 and object.irac2mag gt 0 and object.irac2mag lt 90  and object.wirckmag gt 0. and object.wirckmag lt 90. )

;-----------------------------------------------------------------
plot, (object[point8].irac2mag-3.25) - (object[point8].irac3mag-3.73), (object[point8].irac1mag - 2.7) - (object[point8].irac2mag - 3.25), psym = 2, xrange=[-0.8, 0.6], yrange=[-0.5, 2.5], xtitle = '[4.5]-[5.8]', ytitle='[3.6]-[4.5]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1



oplot, (object[withk[good]].irac2mag-3.25) - (object[withk[good]].irac3mag-3.73), (object[withk[good]].irac1mag - 2.7) - (object[withk[good]].irac2mag - 3.25), psym = 2,color = redcolor

;plot, (object[relax].irac2mag-3.25) - (object[relax].irac3mag-3.73), (object[relax].irac1mag - 2.7) - (object[relax].irac2mag - 3.25), psym = 2, xrange=[-0.8, 0.6], yrange=[-0.5, 2.5], xtitle = '[4.5]-[5.8]', ytitle='[3.6]-[4.5]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1


;im_hessplot, (object[relax].irac2mag-3.25) - (object[relax].irac3mag-3.73), (object[relax].irac1mag - 2.7) - (object[relax].irac2mag - 3.25),  xrange=[-0.8, 0.6], yrange=[-0.5, 2.5], xtitle = '[4.5]-[5.8]', ytitle='[3.6]-[4.5]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1,/negative, nbin2d = 45
;-----------------------------------------------------------------
plot, (object[point8].irac3mag-3.73) - (object[point8].irac4mag-4.37), (object[point8].irac1mag - 2.7) - (object[point8].irac2mag - 3.25), psym = 2, xrange=[-0.2, 1.0], yrange=[-0.5, 2.5], xtitle = '[5.8]-[8.0]', ytitle='[3.6]-[4.5]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1

oplot, (object[withk[good]].irac3mag-3.73) - (object[withk[good]].irac4mag-4.37), (object[withk[good]].irac1mag - 2.7) - (object[withk[good]].irac2mag - 3.25), psym = 2,color=redcolor


;-----------------------------------------------------------------

plot, (object[point8].irac1mag-2.7) - (object[point8].irac4mag-4.37), (object[point8].irac2mag - 3.25) - (object[point8].irac3mag - 3.73), psym = 2, xrange=[-0.2, 2.3], yrange=[-0.7, 0.6], xtitle = '[3.6]-[8.0]', ytitle='[4.5]-[5.8]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1

oplot, (object[withk[good]].irac1mag-2.7) - (object[withk[good]].irac4mag-4.37), (object[withk[good]].irac2mag - 3.25) - (object[withk[good]].irac3mag - 3.73), psym = 2,color = redcolor
;-----------------------------------------------------------------

plot, (object[withk].irac1mag-2.7) - (object[withk].irac2mag-3.25), (object[withk].wirckmag) - (object[withk].irac1mag - 2.7), psym = 2, xrange=[-0.5, 2.0], yrange=[0, 1.8], xtitle = '[3.6]-[4.5]', ytitle='K-[3.6]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1

oplot, (object[withk[good]].irac1mag-2.7) - (object[withk[good]].irac2mag-3.25), (object[withk[good]].wirckmag) - (object[withk[good]].irac1mag - 2.7), psym = 2, color = redcolor

;-----------------------------------------------------------------

;plot, (object[withk].irac1mag-2.7) - (object[withk].irac2mag-3.25), (object[withk].wirckmag) - (object[withk].irac2mag - 3.25), psym = 2, xrange=[-0.5, 2.0], yrange=[0, 4.0], xtitle = '[3.6]-[4.5]', ytitle='K-[4.5]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1

;oplot, (object[withk[good]].irac1mag-2.7) - (object[withk[good]].irac2mag-3.25), (object[withk[good]].wirckmag) - (object[withk[good]].irac2mag - 3.25), psym = 2,color = redcolor

;really only need k and ch1 and ch2 for this plot


good2 = where(  (object[relaxk].wirckmag) - (object[relaxk].irac2mag-3.25) gt 0.9  )
print, "good2", n_elements(good2)
plot, (object[relaxk].irac1mag-2.7) - (object[relaxk].irac2mag-3.25), (object[relaxk].wirckmag) - (object[relaxk].irac2mag - 3.25), psym = 2, xrange=[-0.5, 2.0], yrange=[0, 4.0], xtitle = '[3.6]-[4.5]', ytitle='K-[4.5]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1

oplot, findgen(20)/ 10, 1.05*(findgen(20)/10) + 1.2

oplot, (object[relaxk[good2]].irac1mag-2.7) - (object[relaxk[good2]].irac2mag-3.25), (object[relaxk[good2]].wirckmag) - (object[relaxk[good2]].irac2mag - 3.25), psym = 2,color = bluecolor

;-----------------------------------------------------------------
;now lets  look at some SED's.

openw, outlunh, "/Users/jkrick/ZPHOT/hyperz_cat_bdwarf.txt",/get_lun
for num = 0, n_elements(good2) - 1 do begin

   if object[relaxk[good2(num)]].flamjmag gt 0 and object[relaxk[good2(num)]].flamjmag ne 99 then begin
      fab = 1594.*10^(object[relaxk[good2(num)]].flamjmag/(-2.5))
      jmagab = -2.5*alog10(fab) +8.926
   endif else begin
      jmagab = object[relaxk[good2(num)]].flamjmag
   endelse

   if object[relaxk[good2(num)]].wircjmag gt 0 and object[relaxk[good2(num)]].wircjmag ne 99 then begin
      wircjab = 1594.*10^(object[relaxk[good2(num)]].wircjmag/(-2.5))
      wircjmagab = -2.5*alog10(wircjab) +8.926
   endif else begin
      wircjmagab = object[relaxk[good2(num)]].wircjmag
   endelse

   if object[relaxk[good2(num)]].wirchmag gt 0 and object[relaxk[good2(num)]].wirchmag ne 99 then begin
      wirchab = 1024.*10^(object[relaxk[good2(num)]].wirchmag/(-2.5))
      wirchmagab = -2.5*alog10(wirchab) +8.926
   endif else begin
      wirchmagab = object[relaxk[good2(num)]].wirchmag
   endelse

   if object[relaxk[good2(num)]].wirckmag gt 0 and object[relaxk[good2(num)]].wirckmag ne 99 then begin
      wirckab = 666.8*10^(object[relaxk[good2(num)]].wirckmag/(-2.5))
      wirckmagab = -2.5*alog10(wirckab) +8.926
   endif else begin
      wirckmagab = object[relaxk[good2(num)]].wirckmag
   endelse

   if object[relaxk[good2(num)]].irac1flux lt 0 then err1 = -1. else err1 = 0.05
   if object[relaxk[good2(num)]].irac2flux lt 0 then err2 = -1. else err2 =  0.05;8.926 - 2.5*alog10(0.05E-6*object[relaxk[good2(num)]].irac2)
   if object[relaxk[good2(num)]].irac3flux lt 0 then err3 = -1. else err3 =  0.05;8.926 - 2.5*alog10(0.05E-6*object[relaxk[good2(num)]].irac3)
   if object[relaxk[good2(num)]].irac4flux lt 0 then err4 = -1. else err4 =  0.05;8.926 - 2.5*alog10(0.05E-6*object[relaxk[good2(num)]].irac4)
   
   if object[relaxk[good2(num)]].imagerra gt 1000. then object[relaxk[good2(num)]].imagerra = 1000.
   if object[relaxk[good2(num)]].gmagerra gt 1000. then object[relaxk[good2(num)]].gmagerra = 1000.
   if object[relaxk[good2(num)]].rmagerra gt 1000. then object[relaxk[good2(num)]].rmagerra = 1000.
   if object[relaxk[good2(num)]].umagerra gt 1000. then object[relaxk[good2(num)]].umagerra = 1000.
   

    printf, outlunh, format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',$
                 relaxk[good2(num)], object[relaxk[good2(num)]].umaga, object[relaxk[good2(num)]].gmaga, object[relaxk[good2(num)]].rmaga, $
                 object[relaxk[good2(num)]].imaga,  object[relaxk[good2(num)]].acsmag,  jmagab, wircjmagab, $
                 object[relaxk[good2(num)]].tmassjmag,  wirchmagab, object[relaxk[good2(num)]].tmasshmag, wirckmagab, $
                 object[relaxk[good2(num)]].tmasskmag, object[relaxk[good2(num)]].irac1mag,object[relaxk[good2(num)]].irac2mag,$
                 object[relaxk[good2(num)]].irac3mag,object[relaxk[good2(num)]].irac4mag,   $
                 object[relaxk[good2(num)]].umagerra, object[relaxk[good2(num)]].gmagerra, $
                 object[relaxk[good2(num)]].rmagerra, object[relaxk[good2(num)]].imagerra, object[relaxk[good2(num)]].acsmagerr, object[relaxk[good2(num)]].flamjmagerr, $
                 object[relaxk[good2(num)]].wircjmagerr,0.02, object[relaxk[good2(num)]].wirchmagerr, 0.02, object[relaxk[good2(num)]].wirckmagerr,$
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

help, /memory
;fits_read, '/Users/jkrick/hst/raw/wholeacs.fits', acsdata, acshead
;fits_read, '/Users/jkrick/palomar/lfc/coadd_i.fits', idata, ihead
help, /memory

acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits');
adxy, acshead, object[relaxk[good2]].ra,object[relaxk[good2]].dec , xcenter, ycenter
;headeri = headfits('/Users/jkrick/palomar/lfc/coadd_i.fits');
;adxy, headeri, object[relaxk[good2]].ra,object[relaxk[good2]].dec , xcenteri, ycenteri



size = 10
x = [.3540,.4660,.6255,.7680,.8330,1.196, 1.25,1.24,1.635,1.65,2.15,2.16,3.6000,4.5000,5.8,8.]    ;wavelengths of our photometry in microns
y = x
yerr = x

interest =indgen(n_elements(idz)) 
for n = 0,n_elements(idz) - 1 do begin

;output the acs image of the object if it exists
   print, "working on n ", n

   if object[relaxk[good2(n)]].acsmag gt 0 and object[relaxk[good2(n)]].acsmag lt 90 then begin
;      if object[relaxk[good2(n)]].imaga gt 0 and object[relaxk[good2(n)]].imaga lt 90 then begin
;      print,    xcenteri(n) - size/0.18, xcenteri(n)+ size/0.18  ,ycenteri(n) -size/0.18, ycenteri(n)+size/0.18, "i"
;         plotimage, xrange=[xcenteri(n) - size/0.18, xcenteri(n)+ size/0.18],$
;                 yrange=[ycenteri(n) -size/0.18, ycenteri(n)+size/0.18], $
;                 bytscl(idata, min = -0.01, max = 0.1),$
;                 /preserve_aspect, /noaxes, ncolors=60
;      xyouts, xcenteri(n)- 0.6*size/0.18, ycenteri(n)-1.2*size/0.18, strcompress(string(relaxk[good2(n)])+ 'i' + string(object[idz(n)].imaga)),charthick = 3
;      endif
;   endif else begin

      print,    xcenter(n) - size/0.05, xcenter(n)+ size/0.05  ,ycenter(n) -size/0.05, ycenter(n)+size/0.05
      acsdata = mrdfits('/Users/jkrick/hst/raw/wholeacs.fits', range=[ycenter(n) -size/0.05, ycenter(n)+size/0.05])
      plotimage, xrange=[xcenter(n) - size/0.05, xcenter(n)+ size/0.05],$
;                 yrange=[ycenter(n) -size/0.05, ycenter(n)+size/0.05], $
                 bytscl(acsdata, min = -0.01, max = 0.1),$
                 /preserve_aspect, /noaxes, ncolors=60
      xyouts, xcenter(n)- 0.6*size/0.05, -10., strcompress(string(idz(n))+ '    ACS' + string(object[idz(n)].acsmag) + string( (object[relaxk[good2(n)]].wirckmag) - (object[relaxk[good2(n)]].irac2mag - 3.25))),charthick = 3
;      xyouts, xcenter(n)- 0.6*size/0.05, 100, strcompress(string(idz(n))+ 'ACS' + string(object[idz(n)].acsmag)),charthick = 3
;   endelse
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
