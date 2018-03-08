pro zcat

close, /all
;device, true=24
;device, decomposed=0
colors = GetColor(/load, Start=1)

!p.multi = [0, 0, 1]
ps_open, file = "/Users/jkrick/ZPHOT/sed.small.ps", /portrait, xsize = 4, ysize = 4,/color


;restore, '/Users/jkrick/idlbin/object.sav'

readcol,'/Users/jkrick/ZPHOT/hyperz_swire_small.z_phot',idz, zphota, chia, proba, specta, nagea, agea, $
        ava, ba, zinf99a,zsup99a,zinf90a,zsup90a,zinf68a,zsup68a,zwma,probwma,Mabsa,$
        zphot2a,prob2a,format="A"

;plot, chia, proba, psym = 3, xtitle = "chia", ytitle = "proba", xrange = [0,10]

;readlargecol, '/Users/jkrick/ZPHOT/hyperz_swire_small.obs_sed',idz, u,g,r,i,j,wj,tj,wh,th,wk,tk,one,two,three, four,$
;         uerr, gerr, rerr, ierr, jerr, wjerr,tjerr,wherr,therr,wkerr,tkerr,oneerr, twoerr, threeerr, fourerr,format="A",/debug
readlargecol, '/Users/jkrick/ZPHOT/hyperz_swire_small.obs_sed',idz, u,g,r,i,acs,j,wj,tj,wh,th,wk,tk,one,two,three, four,$
         uerr, gerr, rerr, ierr,acserr, jerr, wjerr,tjerr,wherr,therr,wkerr,tkerr,oneerr, twoerr, threeerr, fourerr,format="A",/debug



;spectral type key
sptkey=[" Arp220 " ," Ell13  " ," Ell2  " ," I19254 " ," M82  " ," Mrk231" ," QSO1 " ," QSO2 " ," S0   " ," Sa  " ," Sb  " ," Sc  " ," Sd  " ," Sey2  " ," Torus" ]     


;-----------------------------------------------------------------------
;plot the fitted SED along with photometry of an individual galaxy

;x = [3540.,4660.,6255.,7680.,36000.,45000.,58000.,80000.]    ;wavelengths of our photometry in angs.
x = [.3540,.4660,.6255,.7680,.8330,1.196, 1.25,1.24,1.635,1.65,2.15,2.16,3.6000,4.5000,5.8,8.,24.]    ;wavelengths of our photometry in microns
;x = [.3540,.4660,.6255,.7680,1.196, 1.25,1.24,1.635,1.65,2.15,2.16,3.6000,4.5000,5.8,8.]    ;wavelengths of our photometry in microns
y = x
yerr = x
;ashby
;interest = [1602,1367,132,1643]

;jason
;interest=[798,5000,848,633,2900,963,3017,1935,1057]

interest =indgen(3000) 
for n = 0,n_elements(interest) - 1 do begin
   y = [u(interest(n)),g(interest(n)),r(interest(n)),i(interest(n)),acs(interest(n)), j(interest(n)),wj(interest(n)),tj(interest(n)),wh(interest(n)),th(interest(n)),wk(interest(n)),tk(interest(n)),one(interest(n)),two(interest(n)),three(interest(n)), four(interest(n))] ;photometry
   yerr = [uerr(interest(n)),gerr(interest(n)),rerr(interest(n)),ierr(interest(n)),acserr(interest(n)),jerr(interest(n)),wjerr(interest(n)),tjerr(interest(n)), wherr(interest(n)),therr(interest(n)),wkerr(interest(n)), tkerr(interest(n)),oneerr(interest(n)),twoerr(interest(n)),threeerr(interest(n)), fourerr(interest(n))]
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
;save, object, filename='/Users/jkrick/idlbin/object.sav'

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
