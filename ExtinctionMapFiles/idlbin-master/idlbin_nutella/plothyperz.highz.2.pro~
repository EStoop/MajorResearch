pro plothyperz, keeper, filenameplot
;/Users/jkrick/zphot/zphot.target.param
;print, "keeper", keeper
print, '/Users/jkrick/zphot/zphot.target.param'
restore, '/Users/jkrick/idlbin/objectnew.sav'
;restore, '/Users/jkrick/idlbin/object.
;----------------------------------------------------------------------
;make hyperz input file
;----------------------------------------------------------------------

openw, outlunh, "/Users/jkrick/ZPHOT/hyperz_cat_target.txt",/get_lun
for num = long(0), n_elements(keeper) - 1 do begin
 ;  print, "working on", num, keeper(num)
   if objectnew[keeper(num)].flamjmag gt 0 and objectnew[keeper(num)].flamjmag ne 99 then begin
      fab = 1594.*10^(objectnew[keeper(num)].flamjmag/(-2.5))
      jmagab = -2.5*alog10(fab) +8.926
   endif else begin
      jmagab = objectnew[keeper(num)].flamjmag
   endelse

   if objectnew[keeper(num)].wircjmag gt 0 and objectnew[keeper(num)].wircjmag ne 99 then begin
      wircjab = 1594.*10^(objectnew[keeper(num)].wircjmag/(-2.5))
      wircjmagab = -2.5*alog10(wircjab) +8.926
   endif else begin
      wircjmagab = objectnew[keeper(num)].wircjmag
   endelse

   if objectnew[keeper(num)].wirchmag gt 0 and objectnew[keeper(num)].wirchmag ne 99 then begin
      wirchab = 1024.*10^(objectnew[keeper(num)].wirchmag/(-2.5))
      wirchmagab = -2.5*alog10(wirchab) +8.926
   endif else begin
      wirchmagab = objectnew[keeper(num)].wirchmag
   endelse

   if objectnew[keeper(num)].wirckmag gt 0 and objectnew[keeper(num)].wirckmag ne 99 then begin
      wirckab = 666.8*10^(objectnew[keeper(num)].wirckmag/(-2.5))
      wirckmagab = -2.5*alog10(wirckab) +8.926
   endif else begin
      wirckmagab = objectnew[keeper(num)].wirckmag
   endelse

   if objectnew[keeper(num)].irac1flux lt 0 then err1 = -1. else err1 = 0.05
   if objectnew[keeper(num)].irac2flux lt 0 then err2 = -1. else err2 =  0.05;8.926 - 2.5*alog10(0.05E-6*objectnew[keeper(num)].irac2)
   if objectnew[keeper(num)].irac3flux lt 0 then err3 = -1. else err3 =  0.05;8.926 - 2.5*alog10(0.05E-6*objectnew[keeper(num)].irac3)
   if objectnew[keeper(num)].irac4flux lt 0 then err4 = -1. else err4 =  0.05;8.926 - 2.5*alog10(0.05E-6*objectnew[keeper(num)].irac4)
   
   if objectnew[keeper(num)].imagerra gt 1000. then objectnew[keeper(num)].imagerra = 1000.
   if objectnew[keeper(num)].gmagerra gt 1000. then objectnew[keeper(num)].gmagerra = 1000.
   if objectnew[keeper(num)].rmagerra gt 1000. then objectnew[keeper(num)].rmagerra = 1000.
   if objectnew[keeper(num)].umagerra gt 1000. then objectnew[keeper(num)].umagerra = 1000.
   
   if objectnew[keeper(num)].tmassjmag lt 0 then tmassjerr = 99 else tmassjerr = 0.02
   if objectnew[keeper(num)].tmasshmag lt 0 then tmassherr = 99 else tmassherr = 0.02
   if objectnew[keeper(num)].tmasskmag lt 0 then tmasskerr = 99 else tmasskerr = 0.02

    printf, outlunh, format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',$
                 keeper(num), objectnew[keeper(num)].umaga, objectnew[keeper(num)].gmaga, objectnew[keeper(num)].rmaga, $
                 objectnew[keeper(num)].imaga,  objectnew[keeper(num)].acsmag, objectnew[keeper(num)].zmagbest, jmagab, wircjmagab, $
                 objectnew[keeper(num)].tmassjmag,  wirchmagab, objectnew[keeper(num)].tmasshmag, wirckmagab, $
                 objectnew[keeper(num)].tmasskmag, objectnew[keeper(num)].irac1mag,objectnew[keeper(num)].irac2mag,$
                 objectnew[keeper(num)].irac3mag,objectnew[keeper(num)].irac4mag,   $
                 objectnew[keeper(num)].umagerra, objectnew[keeper(num)].gmagerra, $
                 objectnew[keeper(num)].rmagerra, objectnew[keeper(num)].imagerra, objectnew[keeper(num)].acsmagerr, objectnew[keeper(num)].zmagerrbest, objectnew[keeper(num)].flamjmagerr, $
                 objectnew[keeper(num)].wircjmagerr,tmassjerr, objectnew[keeper(num)].wirchmagerr, tmassherr, objectnew[keeper(num)].wirckmagerr,$
                 tmasskerr,objectnew[keeper(num)].irac1magerr,objectnew[keeper(num)].irac2magerr,$
                 objectnew[keeper(num)].irac3magerr, objectnew[keeper(num)].irac4magerr;err1,err2,err3,err4;
 endfor
close, outlunh
free_lun, outlunh

;----------------------------------------------------------------------
;run hyperz
;----------------------------------------------------------------------

commandline = '~/bin/hyperz /Users/jkrick/ZPHOT/zphot.target.param'
spawn, commandline
spawn, 'mv *.spe /Users/jkrick/ZPHOT/'

;----------------------------------------------------------------------
;plot SED's of keeper objectnews
;----------------------------------------------------------------------
;!p.multi = [0, 2, 1]
;!p.multi = [0,3, 3]
!p.multi=[0,3,3]
;ps_open, file =filenameplot, /portrait, xsize = 7, ysize =3,/color
;ps_open, file =filenameplot, /portrait;, xsize = 6, ysize =6,/color
;ps_open, file = "/Users/jkrick/nep/target.ps", /portrait, xsize = 6, ysize =6,/color

readcol,'/Users/jkrick/ZPHOT/hyperz_swire_target.z_phot',idz, zphota, chia, proba, specta, nagea, agea, $
        ava, ba, zinf99a,zsup99a,zinf90a,zsup90a,zinf68a,zsup68a,zwma,probwma,Mabsa,$
        zphot2a,prob2a,format="A"

readlargecol, '/Users/jkrick/ZPHOT/hyperz_swire_target.obs_sed',idz, u,g,r,i,acs,z,j,wj,tj,wh,th,wk,tk,one,two,three, four,$
         uerr, gerr, rerr, ierr,zerr,acserr, jerr, wjerr,tjerr,wherr,therr,wkerr,tkerr,oneerr, twoerr, threeerr, fourerr,format="A",/debug

print, n_elements(zphota), n_elements(u)

;spectral type key
sptkey=[" Arp220 " ," Ell13  " ," Ell2  " ," I19254 " ," M82  " ," Mrk231" ," QSO1 " ," QSO2 " ," S0   " ," Sa  " ," Sb  " ," Sc  " ," Sd  " ," Sey2  " ," Torus" ]     

newsptkey=[" Arp220 " ," Ell13  " ," Ell2  "  ," M82  " ," Mrk231" ," QSO1 " ," QSO2 " ," S0   " ," Sa  " ," Sb  " ," Sc  " ," Sd  " ,"Sey18", " Sey2  " ," Torus" ]  
;---------------------------------------------------------------------

acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits');
adxy, acshead, objectnew[keeper].ra,objectnew[keeper].dec , xcenter, ycenter

size = 10

x = [.3540,.4660,.6255,.7680,.8330,0.86,1.196, 1.25,1.24,1.635,1.65,2.15,2.16,3.6000,4.5000,5.8,8.]    ;wavelengths of our photometry in microns
y = x
yerr = x

fits_read, '/Users/jkrick/spitzer/mips/mips24/Combine/mosaic.fits',mipsdata, mipsheader
adxy, mipsheader, objectnew[keeper].ra,objectnew[keeper].dec , mipsxcenter, mipsycenter

fits_read, '/Users/jkrick/spitzer/IRAC/ch1/mosaic.fits',iracdata, iracheader
adxy, iracheader, objectnew[keeper].ra,objectnew[keeper].dec , iracxcenter, iracycenter
fits_read, '/Users/jkrick/spitzer/IRAC/ch2/mosaic.fits',irac2data, irac2header
adxy, irac2header, objectnew[keeper].ra,objectnew[keeper].dec , irac2xcenter, irac2ycenter
fits_read, '/Users/jkrick/spitzer/IRAC/ch3/mosaic.fits',irac3data, irac3header
adxy, irac3header, objectnew[keeper].ra,objectnew[keeper].dec , irac3xcenter, irac3ycenter
fits_read, '/Users/jkrick/spitzer/IRAC/ch4/mosaic.fits',irac4data, irac4header
adxy, irac4header, objectnew[keeper].ra,objectnew[keeper].dec , irac4xcenter, irac4ycenter
fits_read, '/Users/jkrick/mmt/IRACCF.z.coadd.fits',zdata,zheader
adxy,zheader, objectnew[keeper].ra,objectnew[keeper].dec , zxcenter,zycenter
fits_read, '/Users/jkrick/palomar/lfc/coadd_r.fits',rdata,rheader
adxy,rheader, objectnew[keeper].ra,objectnew[keeper].dec , rxcenter,rycenter
fits_read, '/Users/jkrick/mmt/jband/NEPswmos.fits',jdata,jheader
adxy,jheader, objectnew[keeper].ra,objectnew[keeper].dec , jxcenter,jycenter

for n = 0,n_elements(keeper) - 1 do begin

;output the acs image of the objectnew
   print, "working on n ", n

   if objectnew[keeper(n)].rmaga gt 0 then begin;and objectnew[keeper(n)].acsmag lt 90 
      plotimage, xrange=[rxcenter(n) - size/0.18,rxcenter(n)+ size/.18],$
                yrange=[rycenter(n) -size/.18,rycenter(n)+size/.18], $
                 bytscl(rdata, min =-0.2, max = 0.2),title='r',$
                 /preserve_aspect, /noaxes, ncolors=60
 ;     xyouts, xcenter(n)- 0.6*size/0.05, -10., string(keeper(n)),charthick = 3
   endif else begin
      ;need a placeholder
      plot, findgen(10), findgen(10), thick=0, xthick=0,ythick=0, charthick=0, charsize=0, xstyle=4, ystyle=4, psym=3
      xyouts, 4, 4, 'r', charthick=3
   endelse

   if objectnew[keeper(n)].acsmag gt 0  and ycenter(n) - size/0.05 gt 0 and ycenter(n) + size/0.05 lt 20300 then begin;and objectnew[keeper(n)].acsmag lt 90 
      acsdata = mrdfits('/Users/jkrick/hst/raw/wholeacs.fits', range=[ycenter(n) -size/0.05, ycenter(n)+size/0.05])
      plotimage, xrange=[xcenter(n) - size/0.05, xcenter(n)+ size/0.05],$
;;                 yrange=[ycenter(n) -size/0.05, ycenter(n)+size/0.05], $
                 bytscl(acsdata, min = -0.01, max = 0.1),title='acs',$
                 /preserve_aspect, /noaxes, ncolors=60;, position=[0,0.5,0.5,1.0]
      xyouts, xcenter(n)- 0.6*size/0.05, -16., strcompress( string(objectnew[keeper(n)].ra) +'   '+ string(objectnew[keeper(n)].dec)),charthick = 3
   endif else begin
      ;need a placeholder
      plot, findgen(10), findgen(10), thick=0, xthick=0,ythick=0, charthick=0, charsize=0, xstyle=4, ystyle=4, psym=3
      xyouts, 4, 4, 'acs', charthick=3
   endelse

   if objectnew[keeper(n)].zmagbest gt 0 then begin;and objectnew[keeper(n)].acsmag lt 90 
      plotimage, xrange=[zxcenter(n) - size/0.31,zxcenter(n)+ size/.31],$
                 yrange=[zycenter(n) -size/.31,zycenter(n)+size/.31], $
                 bytscl(zdata, min =-100, max = 198),title='z',$
                 /preserve_aspect, /noaxes, ncolors=60
;      xyouts, 0, -10., strcompress( string(objectnew[keeper(n)].ra) +'   '+ string(objectnew[keeper(n)].dec)),charthick = 3
   endif else begin
      ;need a placeholder
      plot, findgen(10), findgen(10), thick=0, xthick=0,ythick=0, charthick=0, charsize=0, xstyle=4, ystyle=4, psym=3
      xyouts, 4, 4, 'z', charthick=3
   endelse

   combfilename=strcompress("/Users/jkrick/nep/highz/f"+string(keeper(n)) + "irac1.ps", /remove_all)
   ps_open, filename=combfilename, /portrait
   if objectnew[keeper(n)].irac1mag gt 0 then begin;and objectnew[keeper(n)].acsmag lt 90 
      plotimage, xrange=[iracxcenter(n) - size/0.6, iracxcenter(n)+ size/0.6],$
                 yrange=[iracycenter(n) -size/0.6, iracycenter(n)+size/0.6], $
                 bytscl(iracdata, min =0.04, max = 0.06),title='irac1',$
                 /preserve_aspect, /noaxes, ncolors=60;,position=[0.5,1.0,0.5,1.0]
 ;     xyouts, xcenter(n)- 0.6*size/0.05, -10., string(keeper(n)),charthick = 3
   endif else begin
      ;need a placeholder
      plot, findgen(10), findgen(10), thick=0, xthick=0,ythick=0, charthick=0, charsize=0, xstyle=4, ystyle=4, psym=3
      xyouts, 4, 4, 'irac1', charthick=3
   endelse
   ps_close, /noprint, /noid

   combfilename=strcompress("/Users/jkrick/nep/highz/f"+string(keeper(n)) + "irac2.ps", /remove_all)
   ps_open, filename=combfilename, /portrait
   if objectnew[keeper(n)].irac2mag gt 0 then begin;and objectnew[keeper(n)].acsmag lt 90 
      plotimage, xrange=[irac2xcenter(n) - size/0.6, irac2xcenter(n)+ size/0.6],$
                 yrange=[irac2ycenter(n) -size/0.6, irac2ycenter(n)+size/0.6], $
                 bytscl(irac2data, min =0.2, max = 0.32),title='irac2',$
                 /preserve_aspect, /noaxes, ncolors=60;,position=[0.5,1.0,0.5,1.0]
 ;     xyouts, xcenter(n)- 0.6*size/0.05, -10., string(keeper(n)),charthick = 3
   endif else begin
      ;need a placeholder
      plot, findgen(10), findgen(10), thick=0, xthick=0,ythick=0, charthick=0, charsize=0, xstyle=4, ystyle=4, psym=3
      xyouts, 4, 4, 'irac2', charthick=3
   endelse
   ps_close, /noprint, /noid

   combfilename=strcompress("/Users/jkrick/nep/highz/f"+string(keeper(n)) + "irac3.ps", /remove_all)
   ps_open, filename=combfilename, /portrait
   if objectnew[keeper(n)].irac3mag gt 0 then begin;and objectnew[keeper(n)].acsmag lt 90 
      plotimage, xrange=[irac3xcenter(n) - size/0.6, irac3xcenter(n)+ size/0.6],$
                 yrange=[irac3ycenter(n) -size/0.6, irac3ycenter(n)+size/0.6], $
                 bytscl(irac3data, min =2.0, max = 2.08),title='irac3',$
                 /preserve_aspect, /noaxes, ncolors=60;,position=[0.5,1.0,0.5,1.0]
 ;     xyouts, xcenter(n)- 0.6*size/0.05, -10., string(keeper(n)),charthick = 3
   endif else begin
      ;need a placeholder
      plot, findgen(10), findgen(10), thick=0, xthick=0,ythick=0, charthick=0, charsize=0, xstyle=4, ystyle=4, psym=3
      xyouts, 4, 4, 'irac3', charthick=3
   endelse
   ps_close, /noprint, /noid

   combfilename=strcompress("/Users/jkrick/nep/highz/f"+string(keeper(n)) + "irac4.ps", /remove_all)
   ps_open, filename=combfilename, /portrait
   if objectnew[keeper(n)].irac4mag gt 0 then begin;and objectnew[keeper(n)].acsmag lt 90 
      plotimage, xrange=[irac4xcenter(n) - size/0.6, irac4xcenter(n)+ size/0.6],$
                 yrange=[irac4ycenter(n) -size/0.6, irac4ycenter(n)+size/0.6], $
                 bytscl(irac4data, min =5.29, max = 5.39),title='irac4',$
                 /preserve_aspect, /noaxes, ncolors=60;,position=[0.5,1.0,0.5,1.0]
 ;     xyouts, xcenter(n)- 0.6*size/0.05, -10., string(keeper(n)),charthick = 3
   endif else begin
      ;need a placeholder
      plot, findgen(10), findgen(10), thick=0, xthick=0,ythick=0, charthick=0, charsize=0, xstyle=4, ystyle=4, psym=3
      xyouts, 4, 4, 'irac4', charthick=3
   endelse
   ps_close, /noprint, /noid

   if objectnew[keeper(n)].mips24mag gt 0 then begin ;and objectnew[keeper(n)].acsmag lt 90 
      plotimage, xrange=[mipsxcenter(n) - size/2.5, mipsxcenter(n)+ size/2.5],$
                 yrange=[mipsycenter(n) -size/2.5, mipsycenter(n)+size/2.5], $
                 bytscl(mipsdata, min = 15.45, max = 15.6),title='mips24',$
                 /preserve_aspect, /noaxes, ncolors=60;,position=[0,0.5,0,0.5]
;      xyouts, xcenter(n)- 0.6*size/0.05, -10., string(keeper(n)),charthick = 3
   endif else begin
      ;need a placeholder
      plot, findgen(10), findgen(10), thick=0, xthick=0,ythick=0, charthick=0, charsize=0, xstyle=4, ystyle=4, psym=3
      xyouts, 4, 4, 'mips24', charthick=3
   endelse


 ;plot the fitted SED along with photometry of an individual galaxy

  y = [u(n),g(n),r(n),i(n),acs(n), z(n),j(n),wj(n),tj(n),wh(n),th(n),wk(n),tk(n),one(n),two(n),three(n), four(n)] ;photometry
  yerr = [uerr(n),gerr(n),rerr(n),ierr(n),acserr(n),zerr(n), jerr(n),wjerr(n),tjerr(n), wherr(n),therr(n),wkerr(n), tkerr(n),oneerr(n),twoerr(n),threeerr(n), fourerr(n)]

   plot,alog10(x), alog10(y), psym = 5, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
        xtitle = "log(microns)", ytitle = "log(flux(microjansky))", xrange=[-0.5,1.5], title=strcompress("objectnew " + string(idz(n)) +string(zphota(n)) + string(proba(n)) + string(sptkey(specta(n)-1)));,position=[0.5,1.0,0,0.5]
;   plot,(x)*10000., (y), psym = 5, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
;        xtitle = "angstroms", ytitle = "flux(microjansky))", title=strcompress("objectnew " + string(idz(n)) +string(zphota(n)) + string(proba(n)) + string(sptkey(specta(n)-1))), xrange=[4000,8000]

   xyouts, alog10(24), alog10(objectnew[keeper[n]].mips24flux), 'x', charthick = 3
   errplot, alog10(x), alog10(y - yerr), alog10(y + yerr)

;
   readcol,strcompress('/Users/jkrick/ZPHOT/' + string(idz(n)) + '.spe', /remove_all),x2, y2,format="A",/silent
   oplot, alog10((x2/1E4)), alog10(y2), thick = 3
 
   
endfor


;-----------------------------------------------------------------------

;ps_close, /noprint, /noid

;---------
;make .reg file for viewing objectnews in ds9
;openw, outlunred, '/Users/jkrick/nep/mips24noacs.reg', /get_lun
;printf, outlunred, 'fk5'
;for rc=0, n_elements(keeper) -1 do  printf, outlunred, 'circle( ', objectnew[keeper[rc]].ra, objectnew[keeper[rc]].dec, ' 3")'
;close, outlunred
;free_lun, outlunred



help, /memory


end
