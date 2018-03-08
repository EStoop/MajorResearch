pro plotthumb
;this code plots thumbnails based on location;
;not catalog number like plothyperz.pro


restore, '/Users/jkrick/idlbin/objectnew.sav'
;restore, '/Users/jkrick/idlbin/object.

filenameplot='/Users/jkrick/nep/mips24nomatch.ps'

readcol, '/Users/jkrick/spitzer/mips/mips24/nomatch.txt', ra, dec,  format="A"

print, ra[64], dec[64]
;----------------------------------------------------------------------
;plot SED's of keeper objectnews
;----------------------------------------------------------------------
;!p.multi = [0, 2, 1]
!p.multi = [0,3, 3]

;ps_open, file =filenameplot, /portrait, xsize = 7, ysize =3,/color
ps_open, file =filenameplot, /portrait;, xsize = 6, ysize =6,/color
;ps_open, file = "/Users/jkrick/nep/target.ps", /portrait, xsize = 6, ysize =6,/color


acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits');
adxy, acshead, ra, dec , xcenter, ycenter
naxis1acs=sxpar(acshead,'NAXIS1')
naxis2acs=sxpar(acshead,'NAXIS2')

size = 10


fits_read, '/Users/jkrick/spitzer/mips/mips24/Combine/mosaic.fits',mipsdata, mipsheader
adxy, mipsheader, ra,dec , mipsxcenter, mipsycenter
naxis1mips=sxpar(mipsheader,'NAXIS1')
naxis2mips=sxpar(mipsheader,'NAXIS2')
fits_read, '/Users/jkrick/spitzer/IRAC/ch1/mosaic.fits',iracdata, iracheader
adxy, iracheader, ra,dec , iracxcenter, iracycenter
naxis1irac=sxpar(iracheader,'NAXIS1')
naxis2irac=sxpar(iracheader,'NAXIS2')
fits_read, '/Users/jkrick/spitzer/IRAC/ch2/mosaic.fits',irac2data, irac2header
adxy, irac2header, ra,dec , irac2xcenter, irac2ycenter
naxis1irac2=sxpar(irac2header,'NAXIS1')
naxis2irac2=sxpar(irac2header,'NAXIS2')
fits_read, '/Users/jkrick/spitzer/IRAC/ch3/mosaic.fits',irac3data, irac3header
adxy, irac3header, ra,dec , irac3xcenter, irac3ycenter
naxis1irac3=sxpar(irac3header,'NAXIS1')
naxis2irac3=sxpar(irac3header,'NAXIS2')
fits_read, '/Users/jkrick/spitzer/IRAC/ch4/mosaic.fits',irac4data, irac4header
adxy, irac4header, ra,dec , irac4xcenter, irac4ycenter
naxis1irac4=sxpar(irac4header,'NAXIS1')
naxis2irac4=sxpar(irac4header,'NAXIS2')
fits_read, '/Users/jkrick/mmt/IRACCF.z.coadd.fits',zdata,zheader
adxy,zheader, ra,dec , zxcenter,zycenter
naxis1z=sxpar(zheader,'NAXIS1')
naxis2z=sxpar(zheader,'NAXIS2')
fits_read, '/Users/jkrick/palomar/lfc/coadd_r.fits',rdata,rheader
adxy,rheader, ra,dec , rxcenter,rycenter
naxis1r=sxpar(rheader,'NAXIS1')
naxis2r=sxpar(rheader,'NAXIS2')
fits_read, '/Users/jkrick/palomar/lfc/coadd_u.fits',udata,uheader
adxy,uheader, ra,dec , uxcenter,uycenter
naxis1u=sxpar(uheader,'NAXIS1')
naxis2u=sxpar(uheader,'NAXIS2')

for n = 0,n_elements(ra) - 1 do begin

;output the acs image of the objectnew
   print, "working on n ", n

   if uxcenter(n) - size/.18 gt 0 and uxcenter(n) - size/.18 lt naxis1u and  uycenter(n) - size/.18 gt 0 and uycenter(n) - size/.18 lt naxis2u then begin
      plotimage, xrange=[uxcenter(n) - size/0.18,uxcenter(n)+ size/.18],$
                yrange=[uycenter(n) -size/.18,uycenter(n)+size/.18], $
                 bytscl(udata, min =-0.01, max = 0.01),title='r',$
                 /preserve_aspect, /noaxes, ncolors=60
   endif else begin
      ;need a placeholder
      plot, findgen(10), findgen(10), thick=0, xthick=0,ythick=0, charthick=0, charsize=0, xstyle=4, ystyle=4, psym=3
      xyouts, 4, 4, 'u', charthick=3
   endelse

   if rxcenter(n) - size/.18 gt 0 and rxcenter(n) - size/.18 lt naxis1r and  rycenter(n) - size/.18 gt 0 and rycenter(n) - size/.18 lt naxis2r then begin
      plotimage, xrange=[rxcenter(n) - size/0.18,rxcenter(n)+ size/.18],$
                yrange=[rycenter(n) -size/.18,rycenter(n)+size/.18], $
                 bytscl(rdata, min =-0.2, max = 0.2),title='r',$
                 /preserve_aspect, /noaxes, ncolors=60
   endif else begin
      ;need a placeholder
      plot, findgen(10), findgen(10), thick=0, xthick=0,ythick=0, charthick=0, charsize=0, xstyle=4, ystyle=4, psym=3
      xyouts, 4, 4, 'r', charthick=3
   endelse

   if xcenter(n) - size/.05 gt 0 and xcenter(n) - size/.05 lt naxis1acs and  ycenter(n) - size/.05 gt 0 and ycenter(n) - size/.05 lt naxis2acs  then begin
      acsdata = mrdfits('/Users/jkrick/hst/raw/wholeacs.fits', range=[ycenter(n) -size/0.05, ycenter(n)+size/0.05])
      plotimage, xrange=[xcenter(n) - size/0.05, xcenter(n)+ size/0.05],$
;;                 yrange=[ycenter(n) -size/0.05, ycenter(n)+size/0.05], $
                 bytscl(acsdata, min = -0.01, max = 0.1),title='acs',$
                 /preserve_aspect, /noaxes, ncolors=60;, position=[0,0.5,0.5,1.0]
      xyouts, xcenter(n)- 0.6*size/0.05, -16., strcompress( string(ra(n)) +'   '+ string(dec(n))),charthick = 3
   endif else begin
      ;need a placeholder
      plot, findgen(10), findgen(10), thick=0, xthick=0,ythick=0, charthick=0, charsize=0, xstyle=4, ystyle=4, psym=3
      xyouts, 4, 4, 'acs', charthick=3
   endelse


   if zxcenter(n) - size/.31 gt 0 and zxcenter(n) - size/.31 lt naxis1z and  zycenter(n) - size/.31 gt 0 and zycenter(n) - size/.31 lt naxis2z  then begin
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


   if iracxcenter(n) - size/0.6 gt 0 and iracxcenter(n) - size/0.6 lt naxis1irac and  iracycenter(n) - size/0.6 gt 0 and iracycenter(n) - size/0.6 lt naxis2irac then begin
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

   if  irac2xcenter(n) - size/0.6 gt 0 and irac2xcenter(n) - size/0.6 lt naxis1irac2 and  irac2ycenter(n) - size/0.6 gt 0 and irac2ycenter(n) - size/0.6 lt naxis2irac2 then begin
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

   if  irac3xcenter(n) - size/0.6 gt 0 and irac3xcenter(n) - size/0.6 lt naxis1irac3 and  irac3ycenter(n) - size/0.6 gt 0 and irac3ycenter(n) - size/0.6 lt naxis2irac3 then begin
      plotimage, xrange=[irac3xcenter(n) - size/0.6, irac3xcenter(n)+ size/0.6],$
                 yrange=[irac3ycenter(n) -size/0.6, irac3ycenter(n)+size/0.6], $
                 bytscl(irac3data, min =1.95, max = 2.2),title='irac3',$
                 /preserve_aspect, /noaxes, ncolors=60;,position=[0.5,1.0,0.5,1.0]
 ;     xyouts, xcenter(n)- 0.6*size/0.05, -10., string(keeper(n)),charthick = 3
   endif else begin
      ;need a placeholder
      plot, findgen(10), findgen(10), thick=0, xthick=0,ythick=0, charthick=0, charsize=0, xstyle=4, ystyle=4, psym=3
      xyouts, 4, 4, 'irac3', charthick=3
   endelse

   if  irac4xcenter(n) - size/0.6 gt 0 and irac4xcenter(n) - size/0.6 lt naxis1irac4 and  irac4ycenter(n) - size/0.6 gt 0 and irac4ycenter(n) - size/0.6 lt naxis2irac4 then begin
      plotimage, xrange=[irac4xcenter(n) - size/0.6, irac4xcenter(n)+ size/0.6],$
                 yrange=[irac4ycenter(n) -size/0.6, irac4ycenter(n)+size/0.6], $
                 bytscl(irac4data, min =5.25, max = 5.4),title='irac4',$
                 /preserve_aspect, /noaxes, ncolors=60;,position=[0.5,1.0,0.5,1.0]
 ;     xyouts, xcenter(n)- 0.6*size/0.05, -10., string(keeper(n)),charthick = 3
   endif else begin
      ;need a placeholder
      plot, findgen(10), findgen(10), thick=0, xthick=0,ythick=0, charthick=0, charsize=0, xstyle=4, ystyle=4, psym=3
      xyouts, 4, 4, 'irac4', charthick=3
   endelse

   if  mipsxcenter(n) - size/2.5 gt 0 and mipsxcenter(n) - size/2.5 lt naxis1mips and  mipsycenter(n) - size/2.5 gt 0 and mipsycenter(n) - size/2.5 lt naxis2mips then begin
      plotimage, xrange=[mipsxcenter(n) - size/2.5, mipsxcenter(n)+ size/2.5],$
                 yrange=[mipsycenter(n) -size/2.5, mipsycenter(n)+size/2.5], $
                 bytscl(mipsdata, min = 15.45, max = 15.6),title='mips24',$
                 /preserve_aspect, /noaxes, ncolors=60 ;,position=[0,0.5,0,0.5]
      xyouts, 0, -16., strcompress( string(ra(n)) +'   '+ string(dec(n))),charthick = 3

   endif else begin
      ;need a placeholder
      plot, findgen(10), findgen(10), thick=0, xthick=0,ythick=0, charthick=0, charsize=0, xstyle=4, ystyle=4, psym=3
      xyouts, 4, 4, 'mips24', charthick=3
   endelse


 
endfor


;-----------------------------------------------------------------------

ps_close, /noprint, /noid

;---------
;make .reg file for viewing objectnews in ds9
;openw, outlunred, '/Users/jkrick/nep/mips24noacs.reg', /get_lun
;printf, outlunred, 'fk5'
;for rc=0, n_elements(keeper) -1 do  printf, outlunred, 'circle( ', objectnew[keeper[rc]].ra, objectnew[keeper[rc]].dec, ' 3")'
;close, outlunred
;free_lun, outlunred



help, /memory


end
