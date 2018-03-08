pro thumbnail, ra, dec

mydevice = !D.NAME
!p.multi = [0, 4,3]
SET_PLOT, 'ps'

device, filename = '/Users/jkrick/nep/thumbail.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

uimage = '/Users/jkrick/palomar/lfc/coadd_u.fits'
gimage = '/Users/jkrick/palomar/lfc/coadd_g.fits'
rimage = '/Users/jkrick/palomar/lfc/coadd_r.fits'
iimage = '/Users/jkrick/palomar/lfc/coadd_i.fits'
jimage = '/Users/jkrick/palomar/wirc/jband_swarp_median_gaiawcs.fits'
himage = '/Users/jkrick/palomar/wirc/hband_swarp_median_gaiawcs.fits'
kimage = '/Users/jkrick/palomar/wirc/coadd_K.fits'
irac1image = '/Users/jkrick/spitzer/irac/ch1/mosaic.fits'
irac2image = '/Users/jkrick/spitzer/irac/ch2/mosaic.fits'
irac3image = '/Users/jkrick/spitzer/irac/ch3/mosaic.fits'
irac4image = '/Users/jkrick/spitzer/irac/ch4/mosaic.fits'
mips24image= '/Users/jkrick/spitzer/mips/mips24/Combine/mosaic.fits'

;first find in the catalog and print relevant info
 deltara= 5.0E-4
deltadec=5.0E-4

restore, '/Users/jkrick/idlbin/object.sav'

print, "num, ra, dec, phot z, umag, umag err, gmag, gmagerr, rmag, rmagerr, imag, imagerr, acsmag, acsmagerr, flamingo J mag, J err, wirc J mag, J err, wirc H mag, H err, wirc K mag, K err, irac1, irac2, irac3, irac4, mips24, mips24 err"

for i =0l, n_elements(object.ra) - 1 do begin

   if object[i].ra gt ra-deltara and  object[i].ra lt ra+ deltara then begin
      if object[i].dec GT dec-deltadec and  object[i].dec lt dec+deltadec then begin
            
          print, i, object[i].ra, object[i].dec, object[i].zphot, object[i].umaga, object[i].umagerra, object[i].gmaga, object[i].gmagerra,object[i].rmaga, object[i].rmagerra,object[i].imaga, object[i].imagerra, object[i].acsmag, object[i].acsmagerr, object[i].flamjmag, object[i].flamjmagerr, object[i].wircjmag, object[i].wircjmagerr, object[i].wirchmag, object[i].wirchmagerr, object[i].wirckmag, object[i].wirckmagerr, object[i].irac1mag,object[i].irac2mag,object[i].irac3mag,object[i].irac4mag, object[i].mips24mag, object[i].mips24magerr,format = '(I10,F10.5,F10.6,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)'

          num = i
      endif
   endif

endfor

;how big should the thumbnails be? from least to most desirable
size = 5 ;in arcseconds
if object[num].gmaga gt 0 and object[num].gmaga ne 99 then size = 0.18*0.25*object[num].gisoarea
if object[num].acsmag gt 0 and object[num].acsmag ne 99 then  size = 0.05*0.25*object[num].acsisoarea 

;now print out thumbnails where the object is on the frame
;--------------
if object[num].umaga gt 0 then begin
   fits_read, uimage, udata, uheader
   adxy,uheader,ra,dec,x,y
   plotimage, xrange=[x - size/0.18, x+size/0.18],$
              yrange=[y-size/0.18, y+size/0.18], $
              bytscl(udata, min = -.01, max = object[num].ufluxmax) ,/preserve_aspect, $
              /noaxes, ncolors=60
   xyouts, x, y-1.2*size/0.18, 'u',charthick = 3
endif
;--------------
if object[num].gmaga gt 0 then begin
   fits_read, gimage, gdata, gheader
   adxy,gheader,ra,dec,x,y
   plotimage, xrange=[x - size/0.18, x+size/0.18],$
              yrange=[y-size/0.18, y+size/0.18], $
              bytscl(gdata, min = -.01, max = object[num].gfluxmax) ,/preserve_aspect, $
              /noaxes, ncolors=60
   xyouts, x, y-1.2*size/0.18, 'g',charthick = 3
endif
;--------------
if object[num].rmaga gt 0 then begin
   fits_read, rimage, rdata, rheader
   adxy,rheader,ra,dec,x,y
   plotimage, xrange=[x - size/0.18, x+size/0.18],$
              yrange=[y-size/0.18, y+size/0.18], $
              bytscl(rdata, min = -.01, max = object[num].rfluxmax) ,/preserve_aspect, $
              /noaxes, ncolors=60
   xyouts, x, y-1.2*size/0.18, 'r',charthick = 3
endif

;--------------

if object[num].acsmag gt 0 then begin
   ;first find which frame the object is on

   images = ['/Users/jkrick/hst/a_drz.fits','/Users/jkrick/hst/b_drz.fits','/Users/jkrick/hst/c_drz.fits',$
             '/Users/jkrick/hst/d_drz.fits','/Users/jkrick/hst/e_drz.fits','/Users/jkrick/hst/f_drz.fits',$
             '/Users/jkrick/hst/g_drz.fits']

   nfiles=n_elements(images)
   
   for i=0,nfiles-1 do begin

      bcd_hd=headfits(images[i],exten=1,/silent)

      naxis1=sxpar(bcd_hd,'NAXIS1')
      naxis2=sxpar(bcd_hd,'NAXIS2')
      crval1=sxpar(bcd_hd,'CRVAL1')
      crval2=sxpar(bcd_hd,'CRVAL2')

      hix=0.95*naxis1
      hiy=0.95*naxis2
      lox=0.05*naxis1
      loy=0.05*naxis2


      adxy,bcd_hd,ra,dec,x,y
      d=sphdist(ra,dec,crval1,crval2,/DEGREES)

      if ((x[0] LT hix) AND (x[0] GT lox) AND (y[0] LT hiy) AND (y[0] GT loy)) then begin
         fits_read, images[i], acsdata, acsheader, exten_no=1
         if acsdata[x,y] gt 0 and acsdata[x+20,y+20] gt 0 and acsdata[x-20,y-20] gt 0 then  begin
            plotimage, xrange=[x - size/0.05, x+ size/0.05],$
                       yrange=[y -size/0.05, y+size/0.05], $
                       bytscl(acsdata, min = -0.01, max = object[num].acsflux/object[num].acsisoarea) ,$
                       /preserve_aspect, /noaxes, ncolors=60
            xyouts, x, y-1.2*size/0.05, 'ACS I',charthick = 3
         endif


      endif

   endfor

endif
;--------------
if object[num].wircjmag gt 0 then begin
   fits_read, jimage, jdata, jheader
   adxy,jheader,ra,dec,x,y
   plotimage, xrange=[x - size/0.25, x+size/0.25],$
              yrange=[y-size/0.25, y+size/0.25], $
              bytscl(jdata, min = -75., max = 70.) ,/preserve_aspect, $
              /noaxes, ncolors=60
   xyouts, x, y-1.2*size/0.25, 'J',charthick = 3
endif
;--------------
if object[num].wirchmag gt 0 then begin
   fits_read, himage, hdata, hheader
   adxy,hheader,ra,dec,x,y
   plotimage, xrange=[x - size/0.25, x+size/0.25],$
              yrange=[y-size/0.25, y+size/0.25], $
              bytscl(hdata, min = -55., max = 50.) ,/preserve_aspect, $
              /noaxes, ncolors=60
   xyouts, x, y-1.2*size/0.25, 'H',charthick = 3
endif
;--------------
if object[num].wirckmag gt 0 then begin
   fits_read, kimage, kdata, kheader
   adxy,kheader,ra,dec,x,y
   plotimage, xrange=[x -size/0.25, x+size/0.25],$
              yrange=[y-size/0.25, y+size/0.25], $
              bytscl(kdata, min = -35, max = 30.) ,/preserve_aspect, $
              /noaxes, ncolors=60
   xyouts, x,y-1.2*size/0.25, 'K',charthick = 3
endif
;--------------
if object[num].irac1mag gt 0 then begin
   fits_read, irac1image, irac1data, irac1header
   adxy,irac1header,ra,dec,x,y
   plotimage, xrange=[x -size/1.2, x+size/1.2],$
              yrange=[y-size/1.2, y+size/1.2], $
              bytscl(irac1data, min = 0., max = 0.2) ,/preserve_aspect, $
              /noaxes, ncolors=60
   xyouts, x,y-1.2*size/1.2, 'IRAC1',charthick = 3
endif
;--------------
if object[num].irac2mag gt 0 then begin
   fits_read, irac2image, irac2data, irac2header
   adxy,irac2header,ra,dec,x,y
   print, x, y
   plotimage, xrange=[x -size/1.2, x+size/1.2],$
              yrange=[y-size/1.2, y+size/1.2], $
              bytscl(irac2data, min = 0., max = 0.2) ,/preserve_aspect, $
              /noaxes, ncolors=60
   xyouts, x,y-1.2*size/1.2, 'IRAC2',charthick = 3
endif
;--------------
if object[num].irac3mag gt 0 then begin
   fits_read, irac3image, irac3data, irac3header
   adxy,irac3header,ra,dec,x,y
   plotimage, xrange=[x -size/1.2, x+size/1.2],$
              yrange=[y-size/1.2, y+size/1.2], $
              bytscl(irac3data, min = 0., max = 0.2) ,/preserve_aspect, $
              /noaxes, ncolors=60
   xyouts, x,y-1.2*size/1.2, 'IRAC3',charthick = 3
endif
;--------------
if object[num].irac4mag gt 0 then begin
   fits_read, irac4image, irac4data, irac4header
   adxy,irac4header,ra,dec,x,y
   plotimage, xrange=[x -size/1.2, x+size/1.2],$
              yrange=[y-size/1.2, y+size/1.2], $
              bytscl(irac4data, min = 0., max = 0.2) ,/preserve_aspect, $
              /noaxes, ncolors=60
   xyouts, x,y-1.2*size/1.2, 'IRAC4',charthick = 3
endif



device, /close
set_plot, mydevice

end
