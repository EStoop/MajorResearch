pro zcat

close, /all
;device, true=24
;device, decomposed=0
colors = GetColor(/load, Start=1)

!p.multi = [0, 0, 1]
ps_open, file = "/Users/jkrick/ZPHOT/sed.jason.ps", /portrait, xsize = 4, ysize = 4,/color


restore, '/Users/jkrick/idlbin/object.sav'

readcol,'/Users/jkrick/ZPHOT/hyperz_swire_jason.z_phot',idz, zphota, chia, proba, specta, nagea, agea, $
        ava, ba, zinf99a,zsup99a,zinf90a,zsup90a,zinf68a,zsup68a,zwma,probwma,Mabsa,$
        zphot2a,prob2a,format="A"

plot, chia, proba, psym = 3, xtitle = "chia", ytitle = "proba", xrange = [0,10]

readlargecol, '/Users/jkrick/ZPHOT/hyperz_swire_jason.obs_sed',idz, u,g,r,i,j,wj,tj,wh,th,wk,tk,one,two,three, four,$
         uerr, gerr, rerr, ierr, jerr, wjerr,tjerr,wherr,therr,wkerr,tkerr,oneerr, twoerr, threeerr, fourerr,format="A",/debug
;readlargecol, '/Users/jkrick/ZPHOT/hyperz_swire_jason.obs_sed',idz, u,g,r,i,acs,j,wj,tj,wh,th,wk,tk,one,two,three, four,$
;         uerr, gerr, rerr, ierr,acserr, jerr, wjerr,tjerr,wherr,therr,wkerr,tkerr,oneerr, twoerr, threeerr, fourerr,format="A",/debug

print, n_elements(zphota), n_elements(u), n_elements(object.rmaga)
;for count = 0l, n_elements(object.rmaga) - 5 do begin
 ;  object[count].zphot = zphota(count)
;   object[count].chisq = chia(count)
;   object[count].prob = proba(count)
;   object[count].spt = specta(count)
;   object[count].age = nagea(count)
;   object[count].av = ava(count)
;   object[count].Mabs = Mabsa(count) ;;
;endfor

;safety = fltarr(n_elements(zphota))
;rarr = fltarr(n_elements(object.rmaga))
;garr = fltarr(n_elements(object.gmaga))
;absmag = fltarr(n_elements(object.gmaga))
;z = fltarr(n_elements(object.gmaga))
;irac1m = fltarr(n_elements(object.gmaga))
;irac2m = fltarr(n_elements(object.gmaga))
;numarr = intarr(n_elements(object.gmaga))
;iarr = fltarr(n_elements(object.gmaga))
;acsarr = fltarr(n_elements(object.gmaga))

;spectral type key
sptkey=[" Arp220 " ," Ell13  " ," Ell2  " ," I19254 " ," M82  " ," Mrk231" ," QSO1 " ," QSO2 " ," S0   " ," Sa  " ," Sb  " ," Sc  " ," Sd  " ," Sey2  " ," Torus" ]     
;---------------------------------------------------------------------
c = 0l
delta = fltarr(n_elements(object.acsmag))
mag = delta
for j2 = 0l, n_elements(object.acsmag) - 1 do begin
   if object[j2].acsmag gt 0 and object[j2].acsmag lt 90 and object[j2].imaga gt 0 and object[j2].imaga lt 90 then begin
      delta[c] = object[j2].acsmag - object[j2].imaga
      mag[c] = object[j2].acsmag
      c = c + 1
   endif
endfor
delta = delta[0:c-1]
mag = mag[0:c-1]

 plot, mag, delta, psym = 2, xrange=[10,30], yrange=[-5,5]

;-----------------------------------------------------------------------

;-----------------------------------------------------------------------
;plot the fitted SED along with photometry of an individual galaxy

;x = [3540.,4660.,6255.,7680.,36000.,45000.,58000.,80000.]    ;wavelengths of our photometry in angs.
;x = [.3540,.4660,.6255,.7680,.8140,1.196, 1.25,1.24,1.635,1.65,2.15,2.16,3.6000,4.5000,5.8,8.,24.]    ;wavelengths of our photometry in microns
;x = [.3540,.4660,.6255,.7680,.8330,1.196, 1.25,1.24,1.635,1.65,2.15,2.16,3.6000,4.5000,5.8,8.]    ;wavelengths of our photometry in microns
x = [.3540,.4660,.6255,.7680,1.196, 1.25,1.24,1.635,1.65,2.15,2.16,3.6000,4.5000,5.8,8.]    ;wavelengths of our photometry in microns
y = x
yerr = x
;ashby
;interest = [1602,1367,132,1643]

;jason
;interest=[798,5000,848,633,2900,963,3017,1935,1057]
;interest = numarr
;interest =indgen(276)
for n = 0,n_elements(idz) - 1 do begin
;  y = [u(interest(n)),g(interest(n)),r(interest(n)),i(interest(n)),acs(interest(n)), j(interest(n)),wj(interest(n)),tj(interest(n)),wh(interest(n)),th(interest(n)),wk(interest(n)),tk(interest(n)),one(interest(n)),two(interest(n)),three(interest(n)), four(interest(n))] ;photometry
;  yerr = [uerr(interest(n)),gerr(interest(n)),rerr(interest(n)),ierr(interest(n)),acserr(interest(n)),jerr(interest(n)),wjerr(interest(n)),tjerr(interest(n)), wherr(interest(n)),therr(interest(n)),wkerr(interest(n)), tkerr(interest(n)),oneerr(interest(n)),twoerr(interest(n)),threeerr(interest(n)), fourerr(interest(n))]
   y = [u(n),g(n),r(n),i(n), j(n),wj(n),tj(n),wh(n),th(n),wk(n),tk(n),one(n),two(n),three(n), four(n)] ;photometry
   yerr = [uerr(n),gerr(n),rerr(n),ierr(n),jerr(n),wjerr(n),tjerr(n), wherr(n),therr(n),wkerr(n), tkerr(n),oneerr(n),twoerr(n),threeerr(n), fourerr(n)]

;   print, n, x, y

   plot,alog10(x), alog10(y), psym = 5, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
        xtitle = "log(microns)", ytitle = "log(flux(microjansky))", title=strcompress("object " + string(idz(n)) +string(zphota(n)) + string(proba(n)) + string(sptkey(specta(n)-1)))

   errplot, alog10(x), alog10(y - yerr), alog10(y + yerr)
;
   readcol,strcompress('/Users/jkrick/ZPHOT/' + string(idz(n)) + '.spe', /remove_all),x2, y2,format="A",/silent
   oplot, alog10((x2/1E4)), alog10(y2), thick = 3

endfor


;-----------------------------------------------------------------------

ps_close, /noprint, /noid


;
END

