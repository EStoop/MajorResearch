pro plothyperz
filenameplot ='/Users/jkrick/palomar/cosmic/mask3.thumb.ps'
apnum = [2,3,5,6,7,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,36,37,38,39,40,42,43]

cosmicnum=[5065,5066,5064,1,2081,334,2076,2078,333,2082,1,2080,3081,1,1,23,332,20,3085,2017,17,1102,1018,1096,1100,303,335,336,1014,1093,1008,2012,301,350,302,26,1017]

objnum = [19544,19613,19373,60100,4970,4941,4313,4697,4645,5027,13626,4871,3617,3608,4131,4413,3621,4062,4103,484,3919,5150,473,3766,4532,401,5064,5413,373,3636,280,221,183,271,367,5191,465]

print, n_elements(apnum), n_elements(cosmicnum), n_elements(objnum)

keeper = objnum

print, "keeper", keeper
restore, '/Users/jkrick/idlbin/object.sav'
;----------------------------------------------------------------------
;make hyperz input file
;----------------------------------------------------------------------

openw, outlunh, "/Users/jkrick/ZPHOT/hyperz_cat_target.txt",/get_lun
for num = 0, n_elements(keeper) - 1 do begin
   print, "working on", num, keeper(num)
   if object[keeper(num)].flamjmag gt 0 and object[keeper(num)].flamjmag ne 99 then begin
      fab = 1594.*10^(object[keeper(num)].flamjmag/(-2.5))
      jmagab = -2.5*alog10(fab) +8.926
   endif else begin
      jmagab = object[keeper(num)].flamjmag
   endelse

   if object[keeper(num)].wircjmag gt 0 and object[keeper(num)].wircjmag ne 99 then begin
      wircjab = 1594.*10^(object[keeper(num)].wircjmag/(-2.5))
      wircjmagab = -2.5*alog10(wircjab) +8.926
   endif else begin
      wircjmagab = object[keeper(num)].wircjmag
   endelse

   if object[keeper(num)].wirchmag gt 0 and object[keeper(num)].wirchmag ne 99 then begin
      wirchab = 1024.*10^(object[keeper(num)].wirchmag/(-2.5))
      wirchmagab = -2.5*alog10(wirchab) +8.926
   endif else begin
      wirchmagab = object[keeper(num)].wirchmag
   endelse

   if object[keeper(num)].wirckmag gt 0 and object[keeper(num)].wirckmag ne 99 then begin
      wirckab = 666.8*10^(object[keeper(num)].wirckmag/(-2.5))
      wirckmagab = -2.5*alog10(wirckab) +8.926
   endif else begin
      wirckmagab = object[keeper(num)].wirckmag
   endelse

   if object[keeper(num)].irac1flux lt 0 then err1 = -1. else err1 = 0.05
   if object[keeper(num)].irac2flux lt 0 then err2 = -1. else err2 =  0.05;8.926 - 2.5*alog10(0.05E-6*object[keeper(num)].irac2)
   if object[keeper(num)].irac3flux lt 0 then err3 = -1. else err3 =  0.05;8.926 - 2.5*alog10(0.05E-6*object[keeper(num)].irac3)
   if object[keeper(num)].irac4flux lt 0 then err4 = -1. else err4 =  0.05;8.926 - 2.5*alog10(0.05E-6*object[keeper(num)].irac4)
   
   if object[keeper(num)].imagerra gt 1000. then object[keeper(num)].imagerra = 1000.
   if object[keeper(num)].gmagerra gt 1000. then object[keeper(num)].gmagerra = 1000.
   if object[keeper(num)].rmagerra gt 1000. then object[keeper(num)].rmagerra = 1000.
   if object[keeper(num)].umagerra gt 1000. then object[keeper(num)].umagerra = 1000.
   
   if object[keeper(num)].tmassjmag lt 0 then tmassjerr = 99 else tmassjerr = 0.02
   if object[keeper(num)].tmasshmag lt 0 then tmassherr = 99 else tmassherr = 0.02
   if object[keeper(num)].tmasskmag lt 0 then tmasskerr = 99 else tmasskerr = 0.02

    printf, outlunh, format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',$
                 keeper(num), object[keeper(num)].umaga, object[keeper(num)].gmaga, object[keeper(num)].rmaga, $
                 object[keeper(num)].imaga,  object[keeper(num)].acsmag, object[keeper(num)].zmagbest, jmagab, wircjmagab, $
                 object[keeper(num)].tmassjmag,  wirchmagab, object[keeper(num)].tmasshmag, wirckmagab, $
                 object[keeper(num)].tmasskmag, object[keeper(num)].irac1mag,object[keeper(num)].irac2mag,$
                 object[keeper(num)].irac3mag,object[keeper(num)].irac4mag,   $
                 object[keeper(num)].umagerra, object[keeper(num)].gmagerra, $
                 object[keeper(num)].rmagerra, object[keeper(num)].imagerra, object[keeper(num)].acsmagerr, object[keeper(num)].zmagerrbest, object[keeper(num)].flamjmagerr, $
                 object[keeper(num)].wircjmagerr,tmassjerr, object[keeper(num)].wirchmagerr, tmassherr, object[keeper(num)].wirckmagerr,$
                 tmasskerr,err1,err2,err3,err4
 endfor
close, outlunh
free_lun, outlunh

;----------------------------------------------------------------------
;run hyperz
;----------------------------------------------------------------------

commandline = '~/bin/hyperz /Users/jkrick/ZPHOT/zphot.target.param'
;spawn, commandline
;spawn, 'mv *.spe /Users/jkrick/ZPHOT/'

;----------------------------------------------------------------------
;plot SED's of keeper objects
;----------------------------------------------------------------------
!p.multi = [0, 0, 1]

ps_open, file =filenameplot, /portrait, xsize = 6, ysize =6,/color
;ps_open, file = "/Users/jkrick/nep/target.ps", /portrait, xsize = 6, ysize =6,/color

readcol,'/Users/jkrick/ZPHOT/hyperz_swire_target.z_phot',idz, zphota, chia, proba, specta, nagea, agea, $
        ava, ba, zinf99a,zsup99a,zinf90a,zsup90a,zinf68a,zsup68a,zwma,probwma,Mabsa,$
        zphot2a,prob2a,format="A"

readlargecol, '/Users/jkrick/ZPHOT/hyperz_swire_target.obs_sed',idz, u,g,r,i,acs,z,j,wj,tj,wh,th,wk,tk,one,two,three, four,$
         uerr, gerr, rerr, ierr,zerr,acserr, jerr, wjerr,tjerr,wherr,therr,wkerr,tkerr,oneerr, twoerr, threeerr, fourerr,format="A",/debug

print, n_elements(zphota), n_elements(u)

;spectral type key
sptkey=[" Arp220 " ," Ell13  " ," Ell2  " ," I19254 " ," M82  " ," Mrk231" ," QSO1 " ," QSO2 " ," S0   " ," Sa  " ," Sb  " ," Sc  " ," Sd  " ," Sey2  " ," Torus" ]     
;---------------------------------------------------------------------

acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits');
adxy, acshead, object[keeper].ra,object[keeper].dec , xcenter, ycenter

size = 10

x = [.3540,.4660,.6255,.7680,.8330,0.86,1.196, 1.25,1.24,1.635,1.65,2.15,2.16,3.6000,4.5000,5.8,8.]    ;wavelengths of our photometry in microns
y = x
yerr = x

fits_read, '/Users/jkrick/spitzer/mips/mips24/Combine/mosaic.fits',mipsdata, mipsheader
adxy, mipsheader, object[keeper].ra,object[keeper].dec , mipsxcenter, mipsycenter

fits_read, '/Users/jkrick/spitzer/IRAC/ch1/mosaic.fits',iracdata, iracheader
adxy, iracheader, object[keeper].ra,object[keeper].dec , iracxcenter, iracycenter
fits_read, '/Users/jkrick/mmt/IRACCF.z.coadd.fits',zdata,zheader
adxy,zheader, object[keeper].ra,object[keeper].dec , zxcenter,zycenter

for n = 0,n_elements(keeper) - 1 do begin

;output the acs image of the object
   print, "working on n ", n
   if object[keeper(n)].acsmag gt 0  and object[keeper(n)].acsmag lt 90  and ycenter(n) - size/0.05 gt 0 and ycenter(n) + size/0.05 lt 20300 then begin;and object[keeper(n)].acsmag lt 90 
      print, 'object[keeper(n)].acsmag',object[keeper(n)].acsmag , xcenter(n) , ycenter(n), size/0.05
      acsdata = mrdfits('/Users/jkrick/hst/raw/wholeacs.fits', range=[ycenter(n) -size/0.05, ycenter(n)+size/0.05])
      plotimage, xrange=[xcenter(n) - size/0.05, xcenter(n)+ size/0.05],$
  ;               yrange=[ycenter(n) -size/0.05, ycenter(n)+size/0.05], $
                 bytscl(acsdata, min = -0.01, max = 0.1),$
                 /preserve_aspect, /noaxes, ncolors=60
;      xyouts, xcenter(n)- 0.6*size/0.05, -10., string(keeper(n)),charthick = 3
   endif else begin

      if object[keeper(n)].zmagbest gt 0 then begin ;and object[keeper(n)].acsmag lt 90 
         plotimage, xrange=[zxcenter(n) - size/0.31,zxcenter(n)+ size/.31],$
                    yrange=[zycenter(n) -size/.31,zycenter(n)+size/.31], $
                    bytscl(zdata, min =-100, max = 500),$
                    /preserve_aspect, /noaxes, ncolors=60
;      xyouts, xcenter(n)- 0.6*size/0.05, -10., string(keeper(n)),charthick = 3
      endif
   endelse
   
   if object[keeper(n)].irac1mag gt 0 then begin;and object[keeper(n)].acsmag lt 90 
      plotimage, xrange=[iracxcenter(n) - size/1.2, iracxcenter(n)+ size/1.2],$
                 yrange=[iracycenter(n) -size/1.2, iracycenter(n)+size/1.2], $
                 bytscl(iracdata, min =0.04, max = 0.06),$
                 /preserve_aspect, /noaxes, ncolors=60
;      xyouts, xcenter(n)- 0.6*size/0.05, -10., string(keeper(n)),charthick = 3
   endif
;   if object[keeper(n)].mips24mag gt 0 then begin ;and object[keeper(n)].acsmag lt 90 
;      plotimage, xrange=[mipsxcenter(n) - size/2.5, mipsxcenter(n)+ size/2.5],$
;                 yrange=[mipsycenter(n) -size/2.5, mipsycenter(n)+size/2.5], $
;                 bytscl(mipsdata, min = 15.3, max = 15.9),$
;                 /preserve_aspect, /noaxes, ncolors=60
;;      xyouts, xcenter(n)- 0.6*size/0.05, -10., string(keeper(n)),charthick = 3
;   endif

 ;plot the fitted SED along with photometry of an individual galaxy

  y = [u(n),g(n),r(n),i(n),acs(n), z(n),j(n),wj(n),tj(n),wh(n),th(n),wk(n),tk(n),one(n),two(n),three(n), four(n)] ;photometry
  yerr = [uerr(n),gerr(n),rerr(n),ierr(n),acserr(n),zerr(n), jerr(n),wjerr(n),tjerr(n), wherr(n),therr(n),wkerr(n), tkerr(n),oneerr(n),twoerr(n),threeerr(n), fourerr(n)]

   plot,alog10(x), alog10(y), psym = 5, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
        xtitle = "log(microns)", ytitle = "log(flux(microjansky))", title=strcompress("ap" + string(apnum(n)) +string(zphota(n)) + string(proba(n)) + string(sptkey(specta(n)-1))), xrange=[-0.5,1.5]
;   plot,(x)*10000., (y), psym = 5, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
;        xtitle = "angstroms", ytitle = "flux(microjansky))", title=strcompress("object " + string(idz(n)) +string(zphota(n)) + string(proba(n)) + string(sptkey(specta(n)-1))), xrange=[4000,8000]

   xyouts, alog10(24), alog10(object[keeper[n]].mips24flux), 'x', charthick = 3
   errplot, alog10(x), alog10(y - yerr), alog10(y + yerr)
;
   readcol,strcompress('/Users/jkrick/ZPHOT/' + string(idz(n)) + '.spe', /remove_all),x2, y2,format="A",/silent
   oplot, alog10((x2/1E4)), alog10(y2), thick = 3
 
   
endfor


;-----------------------------------------------------------------------

ps_close, /noprint, /noid

help, /memory


end
