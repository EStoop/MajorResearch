pro cluster_hyperz
restore, '/Users/jkrick/idlbin/objectnew.sav'

keeper=[2160,1250,2007,1087]
plothyperz, keeper, 'junk.ps'

!p.multi = [0,0, 1]

ps_open, file ='/Users/jkrick/noao/gemini/2008b/fig3.ps', /portrait,/square, xsize = 4, ysize =4;,/color
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
adxy, acshead, objectnew[keeper].ra,objectnew[keeper].dec , xcenter, ycenter

size = 10

x = [.3540,.4660,.6255,.7680,.8330,0.86,1.196, 1.25,1.24,1.635,1.65,2.15,2.16,3.6000,4.5000,5.8,8.]    ;wavelengths of our photometry in microns
y = x
yerr = x

fits_read, '/Users/jkrick/spitzer/mips/mips24/Combine/mosaic.fits',mipsdata, mipsheader
adxy, mipsheader, objectnew[keeper].ra,objectnew[keeper].dec , mipsxcenter, mipsycenter

mipsxcenter = mipsxcenter + 0.5
mipsycenter = mipsycenter + 0.5

fits_read, '/Users/jkrick/spitzer/IRAC/ch1/mosaic.fits',iracdata, iracheader
adxy, iracheader, objectnew[keeper].ra,objectnew[keeper].dec , iracxcenter, iracycenter
;fits_read, '/Users/jkrick/spitzer/IRAC/ch2/mosaic.fits',irac2data, irac2header
;adxy, irac2header, objectnew[keeper].ra,objectnew[keeper].dec , irac2xcenter, irac2ycenter
;fits_read, '/Users/jkrick/spitzer/IRAC/ch3/mosaic.fits',irac3data, irac3header
;adxy, irac3header, objectnew[keeper].ra,objectnew[keeper].dec , irac3xcenter, irac3ycenter
;fits_read, '/Users/jkrick/spitzer/IRAC/ch4/mosaic.fits',irac4data, irac4header
;adxy, irac4header, objectnew[keeper].ra,objectnew[keeper].dec , irac4xcenter, irac4ycenter
;fits_read, '/Users/jkrick/mmt/IRACCF.z.coadd.fits',zdata,zheader
;adxy,zheader, objectnew[keeper].ra,objectnew[keeper].dec , zxcenter,zycenter
;fits_read, '/Users/jkrick/palomar/lfc/coadd_r.fits',rdata,rheader
;adxy,rheader, objectnew[keeper].ra,objectnew[keeper].dec , rxcenter,rycenter

for n = 0,n_elements(keeper) - 1 do begin

;output the acs image of the objectnew
   print, "working on n ", n

 
   if objectnew[keeper(n)].acsmag gt 0  and ycenter(n) - size/0.05 gt 0 and ycenter(n) + size/0.05 lt 20300 then begin;and objectnew[keeper(n)].acsmag lt 90 
      acsdata = mrdfits('/Users/jkrick/hst/raw/wholeacs.fits', range=[ycenter(n) -size/0.05, ycenter(n)+size/0.05])
      plotimage, xrange=[xcenter(n) - size/0.05, xcenter(n)+ size/0.05],$
;;                 yrange=[ycenter(n) -size/0.05, ycenter(n)+size/0.05], $
                 bytscl(acsdata, min = -0.01, max = 0.1),title='acs',$
                 /preserve_aspect, /noaxes, ncolors=60, position=[0.1,0.1,0.9,0.9]
;      xyouts, xcenter(n)- 0.6*size/0.05, -16., strcompress( string(objectnew[keeper(n)].ra) +'   '+ string(objectnew[keeper(n)].dec)),charthick = 3
   endif else begin
      ;need a placeholder
      plot, findgen(10), findgen(10), thick=0, xthick=0,ythick=0, charthick=0, charsize=0, xstyle=4, ystyle=4, psym=3
      xyouts, 4, 4, 'acs', charthick=3
   endelse





 ;plot the fitted SED along with photometry of an individual galaxy

  y = [u(n),g(n),r(n),i(n),acs(n), z(n),j(n),wj(n),tj(n),wh(n),th(n),wk(n),tk(n),one(n),two(n),three(n), four(n)] ;photometry
  yerr = [uerr(n),gerr(n),rerr(n),ierr(n),acserr(n),zerr(n), jerr(n),wjerr(n),tjerr(n), wherr(n),therr(n),wkerr(n), tkerr(n),oneerr(n),twoerr(n),threeerr(n), fourerr(n)]

   plot,alog10(x), alog10(y), psym = 5, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
        xtitle = "log(microns)", ytitle = "log(flux(microjansky))", xrange=[-0.5,1.5],position=[0.1,0.1,0.9,0.9];, title=strcompress("objectnew " + string(idz(n)) +string(zphota(n)) + string(proba(n)) + string(sptkey(specta(n)-1)));,position=[0.5,1.0,0,0.5]
;   plot,(x)*10000., (y), psym = 5, thick = 3, charthick = 3, xthick = 3, ythick = 3, $
;        xtitle = "angstroms", ytitle = "flux(microjansky))", title=strcompress("objectnew " + string(idz(n)) +string(zphota(n)) + string(proba(n)) + string(sptkey(specta(n)-1))), xrange=[4000,8000]

   xyouts, alog10(24), alog10(objectnew[keeper[n]].mips24flux), 'x', charthick = 3
   errplot, alog10(x), alog10(y - yerr), alog10(y + yerr)

;
   readcol,strcompress('/Users/jkrick/ZPHOT/' + string(idz(n)) + '.spe', /remove_all),x2, y2,format="A",/silent
   oplot, alog10((x2/1E4)), alog10(y2), thick = 3
   
   if n gt 1 then begin
      readcol,strcompress('/Users/jkrick/ZPHOT/' + string(idz(n)) + '_highz.spe', /remove_all),x2, y2,format="A",/silent
      oplot, alog10((x2/1E4)), alog10(y2), thick = 3, linestyle = 1
   endif

endfor


;-----------------------------------------------------------------------

ps_close, /noprint, /noid

end
