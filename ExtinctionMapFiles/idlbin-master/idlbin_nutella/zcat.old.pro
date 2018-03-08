pro zcat

close, /all
;device, true=24
;device, decomposed=0
colors = GetColor(/load, Start=1)

!p.multi = [0, 0, 1]
ps_open, file = "/Users/jkrick/palomar/LFC/catalog/zhist.swire.ps", /portrait, xsize = 4, ysize = 4,/color


restore, '/Users/jkrick/idlbin/object.sav'



readcol,'/Users/jkrick/ZPHOT/cat.swire.z_phot',idz, zphota, chia, proba, specta, nagea, agea, $
        ava, ba, zinf99a,zsup99a,zinf90a,zsup90a,zinf68a,zsup68a,zwma,probwma,Mabsa,$
        zphot2a,prob2a,format="A"
safety = zphota

readlargecol, '/Users/jkrick/ZPHOT/cat.swire.obs_sed',idz, u,g,r,i,j,wj,tj,wh,th,wk,tk,one,two,three, four,$
         uerr, gerr, rerr, ierr, jerr, wjerr,tjerr,wherr,therr,wkerr,tkerr,oneerr, twoerr, threeerr, fourerr,format="A",/debug

;print, idz

;readcol, '/Users/jkrick/palomar/lfc/catalog/junk',idc, umag,gmag,rmag,imag,jmag,tjmag,thmag,tkmag,onemag,twomag,threemag, fourmag,umagerr, gmagerr, $
;         rmagerr, imagerr, jmagerr, tjmagerr,thmagerr,tkmagerr,onemagerr, twomagerr, threemagerr, fourmagerr,format="A"



;really only need to do this if I re-run hyperz, and want to keep the new info
for count = 0, n_elements(object.rmaga) - 1 do begin
   object[idz[count]].prob = proba(count)
   object[idz[count]].spt = specta(count)
   object[idz[count]].age = nagea(count)
   object[idz[count]].av = ava(count)
   object[idz[count]].Mabs = Mabsa(count)
endfor

;-----------------------------------------------------------------------
;plot the redshift distribution of the sample with decent redshift determinations
newz = object.zphot
count = 0
for num=0, n_elements(object.zphot) - 1 do begin
   if newz(num) LE 0.1 OR newz(num) gt 5.9 then newz(num) = -10.
   if object[num].prob lt 90.0 then begin
      newz(num) = -10
   endif
;   if object[num].rmaga gt 22.0 or object[num].rmaga lt 10 then begin
;      newz(num) = -10
;   endif

   if object[num].rmaga le 21.5 and object[num].rmaga gt 10 and object[num].irac1 gt 0  and object[num].mips24flux gt 0 then begin
      count = count + 1
   endif

endfor
 
print, "there are ",count, " galaxies with rmag le 21"

plothist, newz, xhist, yhist, bin=0.10, /noprint,xrange=[0,3]

plot, xhist, yhist,title = "z distribution ",$
thick = 3, charthick = 3, xthick = 3, ythick = 3, xtitle = "z",$
xstyle = 1, xrange=[0.0,3.0];, yrange=[0,50]


plothist, object.rmaga, xhist, yhist, bin=0.10, /noprint,/noplot;,xrange=[0,480]

plot, xhist, yhist,title = "rmag distribution ",$
thick = 3, charthick = 3, xthick = 3, ythick = 3, xtitle = "r mag",$
xstyle = 1, xrange=[10.0,29.0], yrange=[0,150]

;-----------------------------------------------------------------------
;plot the redshift distribution of only the sources with all 8 bands
for jc = 0, n_elements(safety) - 1 do begin
   if proba(jc) lt 90. then begin
      safety(jc) = -10 
   endif
endfor

plothist, safety, xhist, yhist, bin=0.1, /noprint,/noplot;,xrange=[0,480]
print, "n_elenments(safety)", n_elements(safety)
plot, xhist, yhist,$
;title = "redshift distribution ",$
thick = 3, charthick = 3, xthick = 3, ythick = 3, xtitle = "z",$
xstyle = 1, xrange=[0,2], yrange=[0,1000], ytitle = "Number"

;-----------------------------------------------------------------------
;plot the fitted SED along with photometry of an individual galaxy
;x = fltarr(8)
;x = [3540.,4660.,6255.,7680.,36000.,45000.,58000.,80000.]    ;wavelengths of our photometry in angs.
x = [.3540,.4660,.6255,.7680,1.196, 1.24,1.65,2.16,3.6000,4.5000,5.8,8.]    ;wavelengths of our photometry in microns
y = x
yerr = x

for n = 1,1650 do begin
   y = [u(n),g(n),r(n),i(n),j(n),tj(n),th(n),tk(n),one(n),two(n),three(n), four(n)] ;photometry
   yerr = [uerr(n),gerr(n),rerr(n),ierr(n),jerr(n),tjerr(n), therr(n),tkerr(n),oneerr(n),twoerr(n),threeerr(n), fourerr(n)]

   plot,alog10(x), alog10(y), psym = 5, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
        xtitle = "log(microns)", ytitle = "log(flux)", title=strcompress("spectrum of object " + string(idz(n-1)) +string(zphota(n-1)) + string(proba(n-1)))

   errplot, alog10(x), alog10(y - yerr), alog10(y + yerr)

   readcol,strcompress('/Users/jkrick/ZPHOT/' + string(n+1) + '.spe', /remove_all),x2, y2,format="A",/silent
   oplot, alog10((x2/1E4)), alog10(y2), thick = 3

;readcol, '/Users/jkrick/ZPHOT/ugrijirac.temp_sed',id, u,g,r,i,j,tj,th,tk,one,two,three, four,format="A"
;y = [u(n-1),g(n-1),r(n-1),i(n-1),j(n-1),tj(n-1),th(n-1),tk(n-1),one(n-1),two(n-1),three(n-1), four(n-1)]
;oplot,alog10(x), alog10(y), psym = 2, thick = 10
endfor

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
;plot colors versus redshift

znew = object.zphot

for count = 0, n_elements(object.zphot) - 1 do begin
   if object[count].gmaga gt 20 or object[count].gmaga lt 8 then znew(count) = -10
   if object[count].rmaga gt 20 or object[count].rmaga lt 8 then znew(count) = -10
;   if object[count].umaga gt 20 or object[count].umaga lt 8 then znew(count) = -10
;   if object[count].prob lt 90. then znew(count) = -10
endfor

;plot, znew, object.rmaga - (6.12 - 2.5*alog10(1E-6*object.irac1)), psym = 2, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
;      xtitle = "phot-z", ytitle = "r-3.6", xrange=[0,1], ystyle = 1, xstyle = 1

plot, znew, object.gmaga - object.rmaga, psym = 2, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
      xtitle = "phot-z", ytitle = "g-r", xrange=[0,6],yrange=[0,4],ystyle = 1, xstyle = 1


;fukugita g-r
x = [0.,0.2,0.5,0.8]
y = [0.77,1.31,1.78,1.76]
;oplot, x, y

;oplot, object.zphot, object.gmaga- object.rmaga, thick =3, psym = 2, color = colors.red
;------------------------------------
;CMD
plot, object.rmaga,object.gmaga - object.rmaga,  psym = 2, xrange=[9,30],xstyle = 1,$
      yrange=[-4,4], xtitle="r magnitude", ytitle = "g-r", title="CMD"

;absolute magnitude vs. z
plot, object.zphot, object.Mabs, thick = 3, psym = 3, xtitle = "redshift", ytitle="r band absolute mag", yrange=[0,-50]

;color color plot
plot, object.gmaga - object.rmaga, object.irac1 - object.irac2, thick = 3, psym = 3, xtitle="g - r", $
      ytitle = "irac1 - irac2", title="color color", xrange=[-2,3],yrange=[-2,4]
;-----------------------------------------------------------------------
;test make_catalog within zphot

readcol,'/Users/jkrick/ZPHOT/ugrijirac_test.z_phot',idz, zphota, chia, proba, specta, nagea, agea, $
        ava, ba, zinf99a,zsup99a,zinf90a,zsup90a,zinf68a,zsup68a,zwma,probwma,Mabsa,$
        zphot2a,prob2a,format="A"

readcol,'/Users/jkrick/ZPHOT/models/ugriirac.model',id, zmodel, spt, junk,mu,mg,mr,mi,m1,m2,m3,m4,mj,av,format="A"

plot, zmodel, zphota, psym = 3, thick = 3, xtitle = "z model", ytitle = "z phot",xrange=[0,2],yrange=[0,2]

newzmodel = fltarr(n_elements(zmodel))
newzphota =  fltarr(n_elements(zmodel))
c=0
for i = 0, n_elements(zmodel) -1 do begin

   if zmodel(i) lt 2.0 then begin
      newzmodel(c) = zmodel(i)
      newzphota(c) = zphota(i)
      c = c+1
   endif

endfor
newzmodel = newzmodel[0:c-1]
newzphota=newzphota[0:c-1]


coeff1 = ROBUST_LINEFIT(newzmodel,newzphota, yfit, sig, coeff_sig)
print, "sig" , sig
;oplot,newzmodel, yfit, thick = 3;, color = colors.yellow
;oplot, newzmodel, yfit + sig;, thick = 3, color = colors.yellow
;oplot, newzmodel, yfit - sig;, thick = 3, color = colors.yellow

plot, zphota, Mabsa, psym = 3, thick = 3, xtitle = "redshift test", ytitle = "r band abs mag test", yrange=[0,-50]

plot, zmodel, mg-mr, psym = 3, thick =0, xtitle = "z model", ytitle = "g-r model",yrange=[0,4];, color = colors.white
oplot, zphota, mg-mr, psym = 3, color = colors.orange
;-----------------------------------------------------------------------

ps_close, /noprint, /noid


;now it has the photo-z info in it.
save, object, filename='/Users/jkrick/idlbin/object.sav'

END

;y = [object[n-1].umaga,object[n-1].gmaga,object[n-1].rmaga,object[n-1].imaga,$
;     6.12 - 2.5*alog10(1E-6*object[n-1].irac1) , 5.64 - 2.5*alog10(1E-6*object[n-1].irac2) ,$
;     5.16 - 2.5*alog10(1E-6*object[n-1].irac3), 4.52 - 2.5*alog10(1E-6*object[n-1].irac4)]

;yerr = [object[n-1].umagerra,object[n-1].gmagerra,object[n-1].rmagerra,object[n-1].imagerra,$
;        0.05*(6.12 - 2.5*alog10(1E-6*object[n-1].irac1)), 0.05*(5.64 - 2.5*alog10(1E-6*object[n-1].irac2)),$
;        0.05*(5.16 - 2.5*alog10(1E-6*object[n-1].irac3)),0.05*(4.52 - 2.5*alog10(1E-6*object[n-1].irac4))]
