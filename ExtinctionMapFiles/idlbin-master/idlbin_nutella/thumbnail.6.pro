pro thumbnail
;takes ra and dec and outputs a ps file with thumbnail images of the object where available
close, /all
mydevice = !D.NAME
!p.multi = [0, 4,3]
SET_PLOT, 'ps'
device, filename = '/Users/jkrick/nep/thumbail.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

imagelist = [ '/Users/jkrick/palomar/lfc/coadd_u.fits', '/Users/jkrick/palomar/lfc/coadd_g.fits', $
'/Users/jkrick/palomar/lfc/coadd_r.fits', '/Users/jkrick/palomar/lfc/coadd_i.fits', $
'/Users/jkrick/hst/raw/wholeacs.fits','/Users/jkrick/palomar/wirc/jband_swarp_median_gaiawcs.fits', $
'/Users/jkrick/palomar/wirc/hband_swarp_median_gaiawcs.fits', '/Users/jkrick/palomar/wirc/coadd_K.fits', $
'/Users/jkrick/spitzer/irac/ch1/mosaic.fits','/Users/jkrick/spitzer/irac/ch2/mosaic.fits', $
'/Users/jkrick/spitzer/irac/ch3/mosaic.fits', '/Users/jkrick/spitzer/irac/ch4/mosaic.fits', $
'/Users/jkrick/spitzer/mips/mips24/Combine/mosaic.fits']

images = ['/Users/jkrick/hst/raw/a.fits','/Users/jkrick/hst/raw/b.fits','/Users/jkrick/hst/raw/c.fits',$
          '/Users/jkrick/hst/raw/d.fits']

;fits_read, '/Users/jkrick/hst/raw/a.fits', dataa, headera
;fits_read, '/Users/jkrick/hst/raw/b.fits', datab, headerb
;fits_read, '/Users/jkrick/hst/raw/c.fits', datac, headerc
;fits_read,  '/Users/jkrick/hst/raw/d.fits', datad, headerd


labellist=['u','g','r','i','ACS i','J','H','K','irac1','irac2','irac3','irac4','mips24']
sizelist=[0.18,0.18,0.18,0.18,0.05,0.25,0.25,0.25,1.2,1.2,1.2,1.2,2.5]

restore, '/Users/jkrick/idlbin/object.sav'
maglist = fltarr(n_elements(labellist))
c = 0

;readcol, '/Users/jkrick/hst/good.radec.txt', raarr, decarr, format="A"

raarr = [264.68143,264.68782,264.68663,264.69259,264.68869,264.67843]
decarr = [69.04381,69.04228,69.04543,69.04577,69.04399,69.04753]

raarr=[264.98380]
decarr=[69.002106]


readcol, '/Users/jkrick/nep/jasontarget.list', numarr, raarr, decarr

;first find in the catalog and print relevant info
deltara= 5.0E-4
deltadec=5.0E-4

true = 0

for index = 0, n_elements(raarr) - 1 do begin
   ra = raarr(index)
   dec = decarr(index)

   dist=sphdist( ra,dec, object.ra,object.dec,/degrees)
   sep = min(dist, ind)
   if (sep lt 0.006) then begin
;   if object[i].ra gt ra-deltara and  object[i].ra lt ra+ deltara then begin
;      if object[i].dec GT dec-deltadec and  object[i].dec lt dec+deltadec then begin
         
         print, "num, ra, dec, phot z, umag, umag err, gmag, gmagerr, rmag, rmagerr, imag, imagerr, acsmag, acsmagerr, flamingo J mag, J err, wirc J mag, J err, wirc H mag, H err, wirc K mag, K err, irac1, irac2, irac3, irac4, mips24, mips24 err"
         
         print, ind, object[ind].ra, object[ind].dec, object[ind].zphot, object[ind].umaga, object[ind].umagerra, object[ind].gmaga, object[ind].gmagerra,object[ind].rmaga, object[ind].rmagerra,object[ind].imaga, object[ind].imagerra, object[ind].acsmag, object[ind].acsmagerr, object[ind].flamjmag, object[ind].flamjmagerr, object[ind].wircjmag, object[ind].wircjmagerr, object[ind].wirchmag, object[ind].wirchmagerr, object[ind].wirckmag, object[ind].wirckmagerr, object[ind].irac1mag,object[ind].irac2mag,object[ind].irac3mag,object[ind].irac4mag, object[ind].mips24mag, object[ind].mips24magerr,format = '(I10,F10.5,F10.6,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)'
         print, object[ind].acsclassstar

         maglist= [object[ind].umaga, object[ind].gmaga, object[ind].rmaga, object[ind].imaga,  object[ind].acsmag,  object[ind].flamjmag, object[ind].wircjmag,  object[ind].wirchmag,  object[ind].wirckmag,  object[ind].irac1mag,object[ind].irac2mag,object[ind].irac3mag,object[ind].irac4mag, object[ind].mips24mag]
 
         ralist = [object[ind].uxcenter, object[ind].gxcenter, object[ind].rxcenter, object[ind].ixcenter,  object[ind].acsra,  object[ind].flamjxcenter, object[ind].wircjxcenter,  object[ind].wirchxcenter,  object[ind].wirckxcenter,  object[ind].irxcenter,object[ind].irxcenter,object[ind].irxcenter,object[ind].irxcenter, object[ind].mips24ra]
         declist = [object[ind].uycenter, object[ind].gycenter, object[ind].rycenter, object[ind].iycenter,  object[ind].acsdec,  object[ind].flamjycenter, object[ind].wircjycenter,  object[ind].wirchycenter,  object[ind].wirckycenter,  object[ind].irycenter,object[ind].irycenter,object[ind].irycenter,object[ind].irycenter, object[ind].mips24dec]
         num = ind
         true = 1
      endif
   

true = 1
if true lt 1 then begin
   print, "There is no object in the catalog at that location"
endif else begin
;---------------------------------------------------------
   minlist = [-0.01,-0.1,-0.2,-0.2,-0.05,-75,-55,-35,0.04,0.21,2.0,5.3,15.46]
   maxlist = [0.01,0.09,0.2,0.22,0.08,70,50,30,0.06,0.23,2.05,5.35,15.56]
;   maglist=[object[num].umaga, object[num].gmaga, object[num].rmaga, object[num].imaga, object[num].acsmag, object[num].wircjmag, object[num].wirchmag, object[num].wirckmag, object[num].irac1mag, object[num].irac2mag, object[num].irac3mag, object[num].irac4mag, object[num].mips24mag]
;---------------------------------------------------------
;how big should the thumbnails be? 
   size = 30                    ;in arcseconds
   size = 10
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
            
            bcd_hd=headfits(images[j]);,exten=1,/silent)
            
            naxis1=sxpar(bcd_hd,'NAXIS1')
            naxis2=sxpar(bcd_hd,'NAXIS2')
            crval1=sxpar(bcd_hd,'CRVAL1')
            crval2=sxpar(bcd_hd,'CRVAL2')
            
            hix=0.95*naxis1
            hiy=0.95*naxis2
            lox=0.05*naxis1
            loy=0.05*naxis2
            
            
            adxy,bcd_hd,ralist[i],declist[i],x,y

            if ((x[0] LT hix) AND (x[0] GT lox) AND (y[0] LT hiy) AND (y[0] GT loy)) then begin
               print, images[j], x, y
               fits_read, images[j], acsdata, acsheader
               if acsdata[x,y] gt 0 and acsdata[x+50,y] ne 0 and acsdata[x-50,y] ne 0 then  begin
                  plotimage, xrange=[x - size/0.05, x+ size/0.05],$
                             yrange=[y -size/0.05, y+size/0.05], $
                             bytscl(acsdata, min = minlist[i], max = maxlist[i]),$
                             /preserve_aspect, /noaxes, ncolors=60
                  xyouts, x- 0.6*size/sizelist[i], y-1.2*size/sizelist[i], strcompress('ACS I' + string(maglist[i])),charthick = 3
                  
                                ;if I have plotted the image, then don't need to see the same thing on another frame
                  j = 4
                ;  xyouts, x, y, '+', charthick = 3,alignment = 0.5
               endif
               
            endif
            
         endfor
         
      endif else begin
         ;non acs, if the object is on the frame then print out an image

  ;          fits_read, imagelist[i], data, header
  ;          if i lt 4 or i gt 7 then begin
  ;              adxy,header,ralist[i],declist[i],x,y 
;              print, "changing ra to x", imagelist[i], ralist[i], declist[i], x, y
 ;           endif else begin
 ;              x = ralist[i]
 ;              y = declist[i]
 ;           endelse
 ;           print, imagelist[i], x, y
;            plotimage, xrange=[x - size/sizelist[i], x+size/sizelist[i]],$
;                       yrange=[y-size/sizelist[i], y+size/sizelist[i]], $
;                       bytscl(data, min = minlist[i], max = maxlist[i]) ,/preserve_aspect, $
;                       /noaxes, ncolors=60
        
;            xyouts, x- 0.6*size/sizelist[i],y-1.2*size/sizelist[i], strcompress(labellist[i] + string(maglist[i])), charthick = 3
;            xyouts, x, y, '+', charthick = 3,alignment = 0.5

         ner = 2
      endelse
      
   endfor
   
endelse

endfor

device, /close
set_plot, mydevice

end
