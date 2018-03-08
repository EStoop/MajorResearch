pro make_filter_figure
!p.thick=3
!p.charthick = 3

ps_open, filename="/Users/jkrick/nep/datapaper/filters_new.ps", /portrait, xsize = 6, ysize = 6
!P.multi=[0,1,2]

filename="hello"
readcol, '/Users/jkrick/ZPHOT/filters/myfilters.list', filename,format="A"

plot,findgen(5), findgen(5), xrange=[2E-1,3], yrange=[0,1.2], /xlog, /nodata, xstyle = 1, xthick =3, ythick = 3, xtitle='Wavelength (microns)', ytitle='Transmission'

for i = 0, n_elements(filename) - 1 do begin
   print, 'working on', i, filename(i)
   readcol, filename(i), num, lambda, transmission,format="A",/silent

   ;work in micron
   lambda = lambda / 1E4

   ;normalize the array
   transmission =  transmission / max(transmission)

    ;plot
   oplot, lambda, transmission
 
endfor

centrallambda = [0.35,0.47,0.63,0.77,0.97,1.20,1.55,2.10]
centrallambda = centrallambda 
filtername = ['u','g','r','i',' z',' J',' H',' K']
xyouts, 0.73, 0.18, 'F814W'
for filter = 0,n_elements(centrallambda) - 1 do begin
   xyouts, centrallambda(filter), 0.1, filtername(filter), alignment=0.5
endfor
;---------------------------------------------------------
;do it again to cover the mid-far IR
plot,findgen(5), findgen(5), xrange=[3, 1E2], yrange=[0,1.2], /xlog, /nodata, xstyle = 1, xthick =3, ythick = 3, xtitle='Wavelength (microns)', ytitle='Transmission'


for i = 0, n_elements(filename) - 1 do begin
   print, 'working on', i, filename(i)
   readcol, filename(i), num, lambda, transmission,format="A",/silent

  ;work in micron
   lambda = lambda / 1E4

   ;normalize the array
   transmission =  transmission / max(transmission)

    ;plot
   oplot, lambda, transmission
 
endfor

centrallambda = [3.5,4.4,5.6,8.00,24.00,70.00]
centrallambda = centrallambda 
filtername = [' ch1',' ch2',' ch3',' ch4',' mips24',' mips70']
print, 'filtername', filtername(5)
for filter = 0,n_elements(centrallambda) - 1 do begin
   xyouts, centrallambda(filter), 0.1, filtername(filter), alignment=0.5
endfor


;---------------------------------------------------------
;add akari curves
readcol, '/Users/jkrick/ZPHOT/filters/akarifilters.list', filename,format="A"

for i = 0, n_elements(filename) - 1 do begin
   print, 'working on', i, filename(i)
   readcol, filename(i), lambda, transmission,format="A",/silent

     ;plot
   oplot, lambda, transmission, linestyle = 2
 
endfor
centrallambda = [4.2,10.5,15.0,19]
centrallambda = centrallambda 
filtername = [' a4',' a11',' a15',' a18']
;print, 'filtername', filtername(5)
for filter = 0,n_elements(centrallambda) - 1 do begin
   xyouts, centrallambda(filter), 0.25, filtername(filter), alignment=0.5
endfor
;---------------------------------------------------------


ps_close, /noprint, /noid

end
