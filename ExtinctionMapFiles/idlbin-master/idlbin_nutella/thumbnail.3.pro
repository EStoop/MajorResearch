pro thumbnail, ra, dec
;takes ra and dec and outputs a ps file with thumbnail images of the object where available
close, /all
mydevice = !D.NAME
!p.multi = [0, 4,3]
SET_PLOT, 'ps'
device, filename = '/Users/jkrick/nep/thumbail.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

imagelist = [ '/Users/jkrick/palomar/lfc/coadd_u.fits', '/Users/jkrick/palomar/lfc/coadd_g.fits', $
'/Users/jkrick/palomar/lfc/coadd_r.fits', '/Users/jkrick/palomar/lfc/coadd_i.fits', $
'doesnt matter there is no image here','/Users/jkrick/palomar/wirc/jband_swarp_median_gaiawcs.fits', $
'/Users/jkrick/palomar/wirc/hband_swarp_median_gaiawcs.fits', '/Users/jkrick/palomar/wirc/coadd_K.fits', $
'/Users/jkrick/spitzer/irac/ch1/mosaic.fits','/Users/jkrick/spitzer/irac/ch2/mosaic.fits', $
'/Users/jkrick/spitzer/irac/ch3/mosaic.fits', '/Users/jkrick/spitzer/irac/ch4/mosaic.fits', $
'/Users/jkrick/spitzer/mips/mips24/Combine/mosaic.fits']

images = ['/Users/jkrick/hst/a_drz.fits','/Users/jkrick/hst/b_drz.fits','/Users/jkrick/hst/c_drz.fits',$
          '/Users/jkrick/hst/d_drz.fits','/Users/jkrick/hst/e_drz.fits','/Users/jkrick/hst/f_drz.fits',$
          '/Users/jkrick/hst/g_drz.fits']

labellist=['u','g','r','i','ACS i','J','H','K','irac1','irac2','irac3','irac4','mips24']
sizelist=[0.18,0.18,0.18,0.18,0.05,0.25,0.25,0.25,1.2,1.2,1.2,1.2,2.5]

restore, '/Users/jkrick/idlbin/object.sav'


;first find in the catalog and print relevant info
deltara= 5.0E-4
deltadec=5.0E-4

true = 0

for i =0, n_elements(object.ra) - 1 do begin

   if object[i].ra gt ra-deltara and  object[i].ra lt ra+ deltara then begin
      if object[i].dec GT dec-deltadec and  object[i].dec lt dec+deltadec then begin
         
         print, "num, ra, dec, phot z, umag, umag err, gmag, gmagerr, rmag, rmagerr, imag, imagerr, acsmag, acsmagerr, flamingo J mag, J err, wirc J mag, J err, wirc H mag, H err, wirc K mag, K err, irac1, irac2, irac3, irac4, mips24, mips24 err"
         
         print, i, object[i].ra, object[i].dec, object[i].zphot, object[i].umaga, object[i].umagerra, object[i].gmaga, object[i].gmagerra,object[i].rmaga, object[i].rmagerra,object[i].imaga, object[i].imagerra, object[i].acsmag, object[i].acsmagerr, object[i].flamjmag, object[i].flamjmagerr, object[i].wircjmag, object[i].wircjmagerr, object[i].wirchmag, object[i].wirchmagerr, object[i].wirckmag, object[i].wirckmagerr, object[i].irac1mag,object[i].irac2mag,object[i].irac3mag,object[i].irac4mag, object[i].mips24mag, object[i].mips24magerr,format = '(I10,F10.5,F10.6,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)'
         
         num = i
         true = 1
      endif
   endif
   
endfor
if true lt 1 then begin
   print, "There is no object in the catalog at that location"
endif else begin
;---------------------------------------------------------
   minlist = [-0.01,-0.1,-0.2,-0.2,-0.05,-75,-55,-35,0.04,0.21,2.0,5.3,15.46]
   maxlist = [0.01,0.09,0.2,0.22,0.1,70,50,30,0.06,0.23,2.05,5.35,15.56]
   maglist=[object[num].umaga, object[num].gmaga, object[num].rmaga, object[num].imaga, object[num].acsmag, object[num].wircjmag, object[num].wirchmag, object[num].wirckmag, object[num].irac1mag, object[num].irac2mag, object[num].irac3mag, object[num].irac4mag, object[num].mips24mag]
;---------------------------------------------------------
;how big should the thumbnails be? 
   size = 30                    ;in arcseconds
;   if object[num].gmaga gt 0 and object[num].gmaga ne 99 then size = 0.18*sqrt(object[num].gisoarea)
;   if object[num].acsmag gt 0 and object[num].acsmag ne 99 then  size = 0.05*sqrt(object[num].acsisoarea)
;use the isoarea sizes to get zoomed in, but in general too small for irac/mips
;---------------------------------------------------------
;now print out thumbnails where the object is on the frame
   
   for i = 0, n_elements(imagelist) - 1 do begin
      
      if i eq 4 then begin
       ;acs tricky
       ;first find which frame the object is on
         
         for j=0,n_elements(images)-1 do begin
            
            bcd_hd=headfits(images[j],exten=1,/silent)
            
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
               fits_read, images[j], acsdata, acsheader, exten_no=1
               if acsdata[x,y] gt 0 and acsdata[x+50,y] ne 0 and acsdata[x-50,y] ne 0 then  begin
                  plotimage, xrange=[x - size/0.05, x+ size/0.05],$
                             yrange=[y -size/0.05, y+size/0.05], $
                             bytscl(acsdata, min = minlist[i], max = maxlist[i]),$
                             /preserve_aspect, /noaxes, ncolors=60
                  xyouts, x, y-1.2*size/sizelist[i], 'ACS I',charthick = 3
               endif
               
            endif
            
         endfor
         
      endif else begin
         ;non acs, if the object is on the frame then print out an image
         if maglist[i] gt 0 then begin
            fits_read, imagelist[i], data, header
            adxy,header,ra,dec,x,y
            plotimage, xrange=[x - size/sizelist[i], x+size/sizelist[i]],$
                       yrange=[y-size/sizelist[i], y+size/sizelist[i]], $
                       bytscl(data, min = minlist[i], max = maxlist[i]) ,/preserve_aspect, $
                       /noaxes, ncolors=60
;            plotimage, xrange=[x - size/sizelist[i], x+size/sizelist[i]],$
;                       yrange=[y-size/sizelist[i], y+size/sizelist[i]], $
;                       bytscl(data, min = 0.95*median(data[x - size/sizelist[i]:x+size/sizelist[i],y-size/sizelist[i]:y+size/sizelist[i]]), max = 1.05*median(data[x - size/sizelist[i]:x+size/sizelist[i],y-size/sizelist[i]:y+size/sizelist[i]])) ,/preserve_aspect, $
;                       /noaxes, ncolors=60
           
            xyouts, x,y-1.2*size/sizelist[i], labellist[i],charthick = 3
         endif
         
      endelse
      
   endfor
   
endelse


device, /close
set_plot, mydevice

end
