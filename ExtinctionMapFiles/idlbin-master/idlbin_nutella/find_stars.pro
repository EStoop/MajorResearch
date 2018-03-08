pro find_stars
close, /all
restore, '/Users/jkrick/idlbin/object.sav'
!P.multi = [0,1,1]

;ps_open, filename='/Users/jkrick/nep/pop3/stars.ps',/portrait,/square,/color

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)


;plothist, object.acsclassstar, xhist, yhist, bin=0.01, thick = 3, /noprint

;plot, object.acsclassstar, object.acsfwhm, thick = 3, psym = 2, yrange = [0,30] , xtitle='acs class star', ytitle='acs fwhm'
;fits_read, '/Users/jkrick/spitzer/irac/ch1/mosaic.fits', data, header
;fits_read, '/Users/jkrick/hst/raw/wholeacs.fits', dataacs, headeracs
;fits_read, '/Users/jkrick/palomar/lfc/coadd_i.fits', dataacs, headeracs

targ0 = where(object.acsclassstar gt 0.98 or object.uclassstar gt 0.98 or object.gclassstar gt 0.98 or object.rclassstar gt 0.98 or object.iclassstar gt 0.98)
targ1 = where(object.iclassstar gt 0.90 and object.acsfwhm lt 3.2 and object.acsellip lt 0.19)
targ2 = where(object.iclassstar gt 0.90  and object.ifwhm lt 6.5 and object.iellip lt 0.20)
;not good targ = where(object.rclassstar gt 0.90 );and object.rfwhm lt 6.5 and object.rellip lt 0.18)
targ3 = where(object.ufluxmax gt 20 or object.gfluxmax gt 78 or object.rfluxmax gt 200. or object.ifluxmax gt 180 or object.acsfluxmax gt 124)

print, n_elements(targ0), n_elements(targ1), n_elements(targ2), n_elements(targ3)
targ = [targ0,targ1,targ2,targ3]


;adxy, headeracs, object[targ].ra, object[targ].dec, x, y

;plot, object[targ].ifwhm, object[targ].iellip, psym = 2, thick = 3, xrange=[0,15]
;im_hessplot, object[targ].rfwhm, object[targ].rellip, xrange=[4,10]

;----------------------------------------------------
;from jason:
; define what is a star and what is not 
; based on the bifurcation in the mag vs isoarea diagram 
; similar to the swire estimator 
; stars have smaller isoarea than galaxies at the same mag 
; 



;window,xsize=640,ysize=480 
plot,object.acsmag,alog10(object.acsisoarea),psym=3,xrange=[12,28] , xtitle="acs mag (AB)", ytitle="acs isoarea", charthick = 1


; define the bounding polygonal area for stars 
px=[12.5,25, 25, 10, 12.5] ; acsmag 
py=[6e5 ,60, 35, 3e5, 6e5] ; ascisoarea 
py=alog10(py) ; switch to log space 


oplot,px,py 

roi = Obj_New('IDLanROI', px, py) 
draw_roi,roi,/line_fill,color='444444'x 


; test to see who is inside the polygon 
star=roi->ContainsPoints(object.acsmag,alog10(object.acsisoarea)) 


; at this point, inside the region is a star, i.e. star>1 

temp=where(star GT 0) 
oplot,object[temp].acsmag,alog10(object[temp].acsisoarea),psym=3,color=bluecolor
;maketake,'star.take',object[temp].acsra,object[temp].acsdec 

;plot, object.acsfwhm, object.acsellip, xrange=[0,10], psym = 3, xtitle='acs fwhm', ytitle='acs ellip'
;oplot,object[temp].acsfwhm,object[temp].acsellip,psym=3,color='FF0000'x 

headeracs = headfits('/Users/jkrick/HST/raw/wholeacs.fits')
adxy, headeracs, object[temp].ra, object[temp].dec, x, y

openw, outlun, '/Users/jkrick/nep/stars.isoarea.reg', /get_lun
for j = 0, n_elements(temp) -1 do begin
   printf, outlun, 'circle(', x(j),y(j), 20,') #color=blue'
endfor
close, outlun
free_lun, outlun


a = where(object.acsfwhm lt 2.9 and object.acsellip lt 0.19 and object.acsfwhm gt 0 and object.acsellip gt 0)
adxy, headeracs, object[a].ra, object[a].dec, x, y

openw, outlun, '/Users/jkrick/nep/stars.fwhm.reg', /get_lun
for j = 0, n_elements(a) -1 do begin
   printf, outlun, 'circle(', x(j),y(j), 15,')'
endfor
close, outlun
free_lun, outlun

oplot, object[a].acsmag,alog10(object[a].acsisoarea),psym=3,color=greencolor

openw, outlun, '/Users/jkrick/nep/stars.dim.reg', /get_lun
printf, outlun, 'fk5'
b = where(object[a].acsmag gt 25 and object[a].acsmag lt 26)

for c = 0, n_elements(b) - 1 do printf, outlun, 'circle (', object[a[b[c]]].ra, object[a[b[c]]].dec, ' 3.5") #color=red'
close, outlun
free_lun, outlun

;ps_close, /noprint,/noid

end


;xarr = fltarr(n_elements(object.ra))
;yarr = xarr
;num = xarr
;j = 0l
;for i = 0l,n_elements(object.ra) -1 do begin
 ;  if object[i].acsmag gt 0 and object[i].acsmag ne 99. then begin
      ;make an array of all the x, y values.
      ; could just read in the irac catalog here, not the whole thing
;      xarr[j] = object[i].acsra
;      yarr[j] = object[i].acsdec
;      num[j] = i
;      j = j + 1
;   endif
;endfor
;xarr = xarr[0:j-1]
;yarr = yarr[0:j-1]

 ;do aperture photometry at the location of the objects in the catalog at two radii, 3", 12"

;adxy, headeracs, xarr, yarr, x, y
;aper, dataacs, x, y, mag, magerr, sky, skyerr, 1, [4.0,15], [25,35], [0,0], /NAN, /exact,/silent

;print, "mag", mag

;print, "mag[0,*]", mag[0,*];
;print, "mag[1,*]", mag[1,*]
;print, "mag[0,*] - mag[1,*],", mag[0,*] - mag[1,*]
;plot 3"aper - 12"aper vs 12"aper

;plot, mag[1,*] , mag[0,*] - mag[1,*], psym = 2, yrange=[-1,5], xrange=[8,26]

;for c= 0, j-1 do begin

;   if mag[1,c] lt 18.0 and mag[0,c] - mag[1,c] lt 0.65 then print, object[num(c)].acsmag,  object[num(c)].acsclassstar
;endfor

; at this point, inside the region is a star, i.e. star>1 
