pro zcat

close, /all
;device, true=24
;device, decomposed=0
colors = GetColor(/load, Start=1)

!p.multi = [0, 0, 1]
ps_open, file = "/Users/jkrick/ZPHOT/zhist.swire.ps", /portrait, xsize = 4, ysize = 4,/color


restore, '/Users/jkrick/idlbin/object.sav'

readcol,'/Users/jkrick/ZPHOT/hyperz_swire.z_phot',idz, zphota, chia, proba, specta, nagea, agea, $
        ava, ba, zinf99a,zsup99a,zinf90a,zsup90a,zinf68a,zsup68a,zwma,probwma,Mabsa,$
        zphot2a,prob2a,format="A"

plot, chia, proba, psym = 3, xtitle = "chia", ytitle = "proba", xrange = [0,10]

;readlargecol, '/Users/jkrick/ZPHOT/hyperz_swire.obs_sed',idz, u,g,r,i,j,wj,tj,wh,th,wk,tk,one,two,three, four,$
;         uerr, gerr, rerr, ierr, jerr, wjerr,tjerr,wherr,therr,wkerr,tkerr,oneerr, twoerr, threeerr, fourerr,format="A",/debug
readlargecol, '/Users/jkrick/ZPHOT/hyperz_swire.obs_sed',idz, u,g,r,i,acs,z,j,wj,tj,wh,th,wk,tk,one,two,three, four,$
         uerr, gerr, rerr, ierr,zerr,acserr, jerr, wjerr,tjerr,wherr,therr,wkerr,tkerr,oneerr, twoerr, threeerr, fourerr,format="A",/debug

print, n_elements(zphota), n_elements(u), n_elements(object.rmaga)
for count = 0l, n_elements(zphota) - 1 do begin
;   print, count
   object[idz].zphot = zphota(count)
   object[idz].chisq = chia(count)
   object[idz].prob = proba(count)
   object[idz].spt = specta(count)
   object[idz].age = nagea(count)
   object[idz].av = ava(count)
   object[idz].Mabs = Mabsa(count) 
;   object[idz(count)].zphot = zphota(count)
;   object[idz(count)].chisq = chia(count)
;   object[idz(count)].prob = proba(count)
;   object[idz(count)].spt = specta(count)
;   object[idz(count)].age = nagea(count)
;   object[idz(count)].av = ava(count)
;   object[idz(count)].Mabs = Mabsa(count) 
endfor

safety = fltarr(n_elements(zphota))
rarr = fltarr(n_elements(object.rmaga))
garr = fltarr(n_elements(object.gmaga))
absmag = fltarr(n_elements(object.gmaga))
z = fltarr(n_elements(object.gmaga))
irac1m = fltarr(n_elements(object.gmaga))
irac2m = fltarr(n_elements(object.gmaga))
numarr = intarr(n_elements(object.gmaga))
iarr = fltarr(n_elements(object.gmaga))
acsarr = fltarr(n_elements(object.gmaga))

;spectral type key
sptkey=[" Arp220 " ," Ell13  " ," Ell2  " ," I19254 " ," M82  " ," Mrk231" ," QSO1 " ," QSO2 " ," S0   " ," Sa  " ," Sb  " ," Sc  " ," Sd  " ," Sey2  " ," Torus" ]     
;---------------------------------------------------------------------
;get good info out of catalog
newz = object.zphot
count = 0l
is = 0l
ys = 0l
c = 0
zs=0
for num=0l, n_elements(object.zphot) - 1 do begin
;   if newz(num) LE 0.1 OR newz(num) gt 5.9 then newz(num) = -10.
   if object[num].prob lt 50.0 then newz(num) = -10
   if object[num].acsclassstar gt 0.7  then newz(num) = -10

    nflux = 0
   if object[num].umaga gt 0 and object[num].umaga ne 99 then nflux = nflux + 1
   if object[num].gmaga gt 0 and object[num].gmaga ne 99 then nflux = nflux + 1
   if object[num].rmaga gt 0 and object[num].rmaga ne 99 then nflux = nflux + 1
   if object[num].imaga gt 0 and object[num].imaga ne 99 then nflux = nflux + 1
   if object[num].acsmag gt 0 and object[num].acsmag ne 99 then nflux = nflux + 1
   if object[num].zmagbest gt 0 and object[num].zmagbest ne 99 then nflux = nflux + 1
   if object[num].flamjmag gt 0 and object[num].flamjmag ne 99 then nflux = nflux + 1
   if object[num].wircjmag gt 0 and object[num].wircjmag ne 99 then nflux = nflux + 1
   if object[num].wirchmag gt 0 and object[num].wirchmag ne 99 then nflux = nflux + 1
   if object[num].wirckmag gt 0 and object[num].wirckmag ne 99 then nflux = nflux + 1
   if object[num].tmassjmag gt 0 and object[num].tmassjmag ne 99 then nflux = nflux + 1
   if object[num].tmasshmag gt 0 and object[num].tmasshmag ne 99 then nflux = nflux + 1
   if object[num].tmasskmag gt 0 and object[num].tmasskmag ne 99 then nflux = nflux + 1
   if object[num].irac1mag gt 0 and object[num].irac1mag ne 99 then nflux = nflux + 1
   if object[num].irac2mag gt 0 and object[num].irac2mag ne 99 then nflux = nflux + 1
   if object[num].irac3mag gt 0 and object[num].irac3mag ne 99 then nflux = nflux + 1
   if object[num].irac4mag gt 0 and object[num].irac4mag ne 99 then nflux = nflux + 1
   if object[num].mips24mag gt 0 and object[num].mips24mag ne 99 then nflux = nflux + 1

   if nflux gt 4 and object[num].acsclassstar lt 0.9 and object[num].prob gt 50. then begin
      safety(is) = object[num].zphot
      is = is + 1
   endif

;don't include the stars
   if object[num].acsclassstar gt 0.9 and object[num].rmaga gt 0. and object[num].gmaga gt 0. and object[num].irac1mag gt 0. and object[num].irac2mag gt 0. then begin ;and object[num].prob gt 0. then begin
      rarr(ys) = object[num].rmaga
      garr(ys) = object[num].gmaga
      absmag(ys) = object[num].Mabs
      z(ys) = object[num].zphot
      irac1m(ys) = object[num].irac1mag
      irac2m(ys) = object[num].irac2mag
      ys = ys + 1
   endif

;test zphot
   if object[num].rmaga lt 22 and num lt 5900 then begin
      numarr(c) = num
      c=c+1
   endif

;test acs
if object[num].acsclassstar gt 0.9 and object[num].acsmag gt 0 and object[num].acsmag lt 90 and object[num].imaga gt 0 and object[num].imaga lt 90 then begin
   iarr(zs) = object[num].imaga
;   flux = 10^((25.937-object[num].acsmag)/4)
;   acsarr(zs) = 25.937 - 2.5*alog10(2*flux)
   acsarr(zs) = object[num].acsmag
   zs = zs + 1
endif

endfor

safety = safety[0:is-1]
rarr=rarr[0:ys-1]
garr = garr[0:ys-1]
absmag = absmag[0:ys-1]
z = z[0:ys-1]
irac1m = irac1m[0:ys-1]
irac2m = irac2m[0:ys-1]
numarr = numarr[0:c-1]
print, "n_elements(numarr)",n_elements(numarr)
;-----------------------------------------------------------------------
;plot the redshift distribution of the sample with decent redshift determinations without stars

plothist, newz, xhist, yhist, bin=0.05, /noprint,/noplot

plot, xhist, yhist,title = "z distribution ",$
thick = 3, charthick = 3, xthick = 3, ythick = 3, xtitle = "z",$
xstyle = 1, xrange=[0.0,6.0];, yrange=[0,50]

;-----------------------------------------------------------------------
;plot the redshift distribution of only the sources with at least 5 bands and a good z determination without stars

;plothist, safety, xhist, yhist, bin=0.05, /noprint,/noplot;,xrange=[0,480]
;print, "n_elenments(safety)", n_elements(safety)

;plot, xhist, yhist,title = "redshift distribution ",$
;thick = 3, charthick = 3, xthick = 3, ythick = 3, xtitle = "z",$
;xstyle = 1, xrange=[0,3], yrange=[0,50], ytitle = "Number"


;-----------------------------------------------------------------------
;plot the distribution of r magnitudes

;plothist, rarr, xhist, yhist, bin=0.10, /noprint,/noplot;,xrange=[0,480]

;plot, xhist, yhist,title = "rmag distribution ",$
;thick = 3, charthick = 3, xthick = 3, ythick = 3, xtitle = "r mag",$
;xstyle = 1, xrange=[10.0,30.0], yrange=[0,200]

;------------------------------------
;CMD
;plot, rarr,garr - rarr,  psym = 2, xrange=[9,30],xstyle = 1,$
;      yrange=[-4,4], xtitle="r magnitude", ytitle = "g-r", title="CMD"

;absolute magnitude vs. z
;plot, z, absmag, thick = 3, psym = 3, xtitle = "redshift", ytitle="r band absolute mag", $
 ;     yrange=[0,-40], xrange=[0,6]

;color color plot
;plot, garr-rarr, irac1m-irac2m, thick = 3, psym = 3, xtitle="g - r", $
;      ytitle = "irac1 - irac2", title="color color", xrange=[-2,3],yrange=[-2,4]

; colors versus redshift
;plot, z, rarr - irac1m, psym = 2, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
;      xtitle = "phot-z", ytitle = "r-3.6", xrange=[0,6], ystyle = 1, xstyle = 1,yrange=[-5,5]

;plot, z, garr-rarr, psym = 2, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
;      xtitle = "phot-z", ytitle = "g-r", xrange=[0,6],yrange=[0,4],ystyle = 1, xstyle = 1

;fukugita g-r
x = [0.,0.2,0.5,0.8]
y = [0.77,1.31,1.78,1.76]
;oplot, x, y, color = 255

;acs vs palomar
;plot, acsarr, acsarr - iarr, psym  = 2, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
;      xtitle = "acs mag 814", ytitle = "acs814 - i", xrange=[10,30], ystyle = 1, xstyle = 1,yrange=[-5,5]
;plot, iarr, acsarr - iarr, psym  = 2, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
 ;     xtitle = "palomar i", ytitle = "acs814 - i", xrange=[10,30], ystyle = 1, xstyle = 1,yrange=[-5,5]

;-----------------------------------------------------------------------
;plot the fitted SED along with photometry of an individual galaxy

;x = [3540.,4660.,6255.,7680.,36000.,45000.,58000.,80000.]    ;wavelengths of our photometry in angs.
;x = [.3540,.4660,.6255,.7680,.8140,1.196, 1.25,1.24,1.635,1.65,2.15,2.16,3.6000,4.5000,5.8,8.,24.]    ;wavelengths of our photometry in microns
x = [.3540,.4660,.6255,.7680,.8330,.9250,1.196, 1.25,1.24,1.635,1.65,2.15,2.16,3.6000,4.5000,5.8,8.]    ;wavelengths of our photometry in microns
;x = [.3540,.4660,.6255,.7680,1.196, 1.25,1.24,1.635,1.65,2.15,2.16,3.6000,4.5000,5.8,8.]    ;wavelengths of our photometry in microns
y = x
yerr = x
;ashby
;interest = [1602,1367,132,1643]

;jason
;interest=[798,5000,848,633,2900,963,3017,1935,1057]
interest = numarr
interest =indgen(1000)
interest = 0
for n = 0,n_elements(interest) - 1 do begin
  y = [u(interest(n)),g(interest(n)),r(interest(n)),i(interest(n)),acs(interest(n)), z(interest(n)),j(interest(n)),wj(interest(n)),tj(interest(n)),wh(interest(n)),th(interest(n)),wk(interest(n)),tk(interest(n)),one(interest(n)),two(interest(n)),three(interest(n)), four(interest(n))] ;photometry
  yerr = [uerr(interest(n)),gerr(interest(n)),rerr(interest(n)),ierr(interest(n)),acserr(interest(n)),zerr(interest(n)),jerr(interest(n)),wjerr(interest(n)),tjerr(interest(n)), wherr(interest(n)),therr(interest(n)),wkerr(interest(n)), tkerr(interest(n)),oneerr(interest(n)),twoerr(interest(n)),threeerr(interest(n)), fourerr(interest(n))]
;   y = [u(interest(n)),g(interest(n)),r(interest(n)),i(interest(n)), j(interest(n)),wj(interest(n)),tj(interest(n)),wh(interest(n)),th(interest(n)),wk(interest(n)),tk(interest(n)),one(interest(n)),two(interest(n)),three(interest(n)), four(interest(n))] ;photometry
;   yerr = [uerr(interest(n)),gerr(interest(n)),rerr(interest(n)),ierr(interest(n)),jerr(interest(n)),wjerr(interest(n)),tjerr(interest(n)), wherr(interest(n)),therr(interest(n)),wkerr(interest(n)), tkerr(interest(n)),oneerr(interest(n)),twoerr(interest(n)),threeerr(interest(n)), fourerr(interest(n))]

   plot,alog10(x), alog10(y), psym = 5, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
        xtitle = "log(microns)", ytitle = "log(flux(microjansky))", title=strcompress("object " + string(interest(n)) +string(zphota(interest(n))) + string(proba(interest(n))) + string(sptkey(specta(interest(n))-1)))

   errplot, alog10(x), alog10(y - yerr), alog10(y + yerr)
;
   readcol,strcompress('/Users/jkrick/ZPHOT/' + string(interest(n)) + '.spe', /remove_all),x2, y2,format="A",/silent
   oplot, alog10((x2/1E4)), alog10(y2), thick = 3

endfor


;-----------------------------------------------------------------------

ps_close, /noprint, /noid


;now it has the photo-z info in it.
save, object, filename='/Users/jkrick/idlbin/object.test.sav'

END

;y = [object[n-1].umaga,object[n-1].gmaga,object[n-1].rmaga,object[n-1].imaga,$
;     6.12 - 2.5*alog10(1E-6*object[n-1].irac1) , 5.64 - 2.5*alog10(1E-6*object[n-1].irac2) ,$
;     5.16 - 2.5*alog10(1E-6*object[n-1].irac3), 4.52 - 2.5*alog10(1E-6*object[n-1].irac4)]

;yerr = [object[n-1].umagerra,object[n-1].gmagerra,object[n-1].rmagerra,object[n-1].imagerra,$
;        0.05*(6.12 - 2.5*alog10(1E-6*object[n-1].irac1)), 0.05*(5.64 - 2.5*alog10(1E-6*object[n-1].irac2)),$
;        0.05*(5.16 - 2.5*alog10(1E-6*object[n-1].irac3)),0.05*(4.52 - 2.5*alog10(1E-6*object[n-1].irac4))]

;print, idz

;readcol, '/Users/jkrick/palomar/lfc/catalog/junk',idc, umag,gmag,rmag,imag,jmag,tjmag,thmag,tkmag,onemag,twomag,threemag, fourmag,umagerr, gmagerr, $
;         rmagerr, imagerr, jmagerr, tjmagerr,thmagerr,tkmagerr,onemagerr, twomagerr, threemagerr, fourmagerr,format="A"
;readcol, '/Users/jkrick/ZPHOT/ugrijirac.temp_sed',id, u,g,r,i,j,tj,th,tk,one,two,three, four,format="A"
;y = [u(n-1),g(n-1),r(n-1),i(n-1),j(n-1),tj(n-1),th(n-1),tk(n-1),one(n-1),two(n-1),three(n-1), four(n-1)]
;oplot,alog10(x), alog10(y), psym = 2, thick = 10
;n = 1027
;   y = [u(n-1),g(n-1),r(n-1),i(n-1),one(n-1),two(n-1),three(n-1), four(n-1)] ;photometry
;    yerr = [uerr(n-1),gerr(n-1),rerr(n-1),ierr(n-1),oneerr(n-1),twoerr(n-1),threeerr(n-1), fourerr(n-1)]

;   oplot,alog10(x), alog10(y), psym = 6, thick = 3
;   errplot, alog10(x), alog10(y - yerr), alog10(y + yerr)

;   readcol,strcompress('/Users/jkrick/ZPHOT/' + string(n) + '.spe', /remove_all),x2, y2,format="A"
;   oplot, alog10((x2/1E4)), alog10(y2), thick = 3, color = colors.blue

;for n = 1, 48 do begin    ;which galaxy do I want to look at

;   y = [u(n-1),g(n-1),r(n-1),i(n-1),one(n-1),two(n-1),three(n-1), four(n-1)] ;photometry
;   yerr = [uerr(n-1),gerr(n-1),rerr(n-1),ierr(n-1),oneerr(n-1),twoerr(n-1),threeerr(n-1), fourerr(n-1)]

;   plot,alog10(x), alog10(y), psym = 2, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
;;        title=strcompress("spectrum of object " + string(n)),
; xtitle = "log(microns)", ytitle = "log(flux)"
;   errplot, alog10(x), alog10(y - yerr), alog10(y + yerr)

;   readcol,strcompress('/Users/jkrick/ZPHOT/' + string(n) + '.spe', /remove_all),x2, y2,format="A"
;   oplot, alog10((x2/1E4)), alog10(y2), thick = 3

;readcol, '/Users/jkrick/ZPHOT/ugrijirac.temp_sed',id, u,g,r,i,one,two,three, four,format="A"
;y = [u(n-1),g(n-1),r(n-1),i(n-1),one(n-1),two(n-1),three(n-1), four(n-1)]
;oplot,alog10(x), alog10(y), psym = 2, thick = 10

;endfor
;-----------------------------------------------------------------------
;test make_catalog within zphot

;readcol,'/Users/jkrick/ZPHOT/ugrijirac_test.z_phot',idz, zphota, chia, proba, specta, nagea, agea, $
;        ava, ba, zinf99a,zsup99a,zinf90a,zsup90a,zinf68a,zsup68a,zwma,probwma,Mabsa,$
;        zphot2a,prob2a,format="A"

;readcol,'/Users/jkrick/ZPHOT/models/ugriirac.model',id, zmodel, spt, junk,mu,mg,mr,mi,m1,m2,m3,m4,mj,av,format="A"

;;plot, zmodel, zphota, psym = 3, thick = 3, xtitle = "z model", ytitle = "z phot",xrange=[0,2],yrange=[0,2]

;newzmodel = fltarr(n_elements(zmodel))
;newzphota =  fltarr(n_elements(zmodel))
;c=0
;for i = 0, n_elements(zmodel) -1 do begin

;   if zmodel(i) lt 2.0 then begin
;      newzmodel(c) = zmodel(i)
;      newzphota(c) = zphota(i)
;      c = c+1
;   endif

;endfor
;newzmodel = newzmodel[0:c-1]
;newzphota=newzphota[0:c-1]


;coeff1 = ROBUST_LINEFIT(newzmodel,newzphota, yfit, sig, coeff_sig)
;print, "sig" , sig
;oplot,newzmodel, yfit, thick = 3;, color = colors.yellow
;oplot, newzmodel, yfit + sig;, thick = 3, color = colors.yellow
;oplot, newzmodel, yfit - sig;, thick = 3, color = colors.yellow

;;plot, zphota, Mabsa, psym = 3, thick = 3, xtitle = "redshift test", ytitle = "r band abs mag test", yrange=[0,-50]

;;plot, zmodel, mg-mr, psym = 3, thick =0, xtitle = "z model", ytitle = "g-r model",yrange=[0,4];, color = colors.white
;;oplot, zphota, mg-mr, psym = 3, color = colors.orange
