PRO catalog_test

close, /all
colors = GetColor(/load, Start=1)
mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/Users/jkrick/palomar/wirc/comparej.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

restore, '/Users/jkrick/idlbin/object.sav'

plot, object.wircjmag, object.flamjmag, psym = 2, xrange=[15,21], xtitle = "wirc J mag", ytitle = "flamingos J mag",$
      yrange=[15,21], thick = 3, xthick = 3, ythick=3, charthick = 3, xstyle = 1, ystyle = 1
oplot, findgen(30), findgen(30)

;print, object(sort(object.flamjdec)).flamjdec




;look at wirc k band distribution to find rough detection thresholds
plothist, object.wirckmag, xhist, yhist,/noprint,/noplot, bin = 0.1, thick = 3, xthick=3,ythick=3, xtitle = "wirc k mag"

plot, xhist, yhist, xrange=[15,20], thick = 3

;look at overlay of wirc and flamingo J band to test matching of catalogs.
fits_read, "/Users/jkrick/palomar/LFC/coadd_r.fits", data, header
plotimage, xrange=[2000,6000],yrange=[2000,6000], bytscl(data, min = -0.2, max = 0.2) ,$; bytscl(data, min=1,max=3),$
 /preserve_aspect, /noaxes, ncolors=8

for i = 0, n_elements(object.ra) - 1 do begin

   if (object[i].wircjmag gt 0) and (object[i].wircjmag lt 99) then begin
      adxy, header, object[i].wircjra, object[i].wircjdec, x,y

      xyouts,  x,y, 'w', alignment=0.5
      
      if (object[i].flamjmag lt 10) then begin
         adxy, header, object[i].wircjra, object[i].wircjdec, x,y
         xyouts,  x,y, 'w', alignment=0.5, color = colors.red

         


      endif
   endif

   if (object[i].flamjmag gt 0) and (object[i].flamjmag lt 99) then begin
      adxy, header, object[i].flamjra, object[i].flamjdec, x,y

      xyouts,  x,y, 'f', alignment=0.5
      
   endif



endfor




device, /close
set_plot, mydevice

END
