pro test
ps_open, filename='/Users/jkrick/IRAC/iwic210/noise.ps',/portrait,/color

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)
graycolor = FSC_COLOR("gray", !D.Table_Size-9)
browncolor = FSC_COLOR("brown", !D.Table_Size-10)

;colorname = ['redcolor','bluecolor','greencolor','yellowcolor','cyancolor','orangecolor','purplecolor','graycolor','browncolor']
colorname = [ graycolor, redcolor , bluecolor , greencolor , yellowcolor , cyancolor , orangecolor , purplecolor ,  browncolor ]

dir_name = '/Users/jkrick/IRAC/iwic210/r7665664'

command1 = ' find '+dir_name+' -name "*.fits" > ' + dir_name + '/files_test.list'
command2 = 'grep ch1 < '+dir_name+'/files_test.list |  grep _bcd.fits > '+dir_name+'/ch1_bcd.list'
command3 = 'grep ch2 < '+dir_name+'/files_test.list |  grep _bcd.fits > '+dir_name+'/ch2_bcd.list'

a = [command1, command2, command3]
for i = 0, n_elements(a) -1 do spawn, a(i)

readcol, dir_name+'/ch1_bcd.list', bcdname_ch1, format="A"
readcol, dir_name+'/ch2_bcd.list', bcdname_ch2, format="A"

;two stars
ra = [271.015833, 270.93958333]
dec = [66.928333, 66.93416667 ]

gain = [3.3,3.7]



nimages = n_elements(bcdname_ch1)
;two stars to keep track of
xcenarr_1 = (fltarr(nimages))
ycenarr_1 = fltarr(nimages)
xcenarr_2 = (fltarr(nimages))
ycenarr_2 = fltarr(nimages)

;for each channel
for k = 0, 1 do begin
print, 'working on k', k
   if k eq 0 then bcdname = bcdname_ch1
   if k eq 1 then bcdname = bcdname_ch2
print, 'bcdname', bcdname
;find the locations of the stars on all frames
for i = 0, nimages - 1 do begin
   header = headfits(bcdname(i))    ;read in the header of the image on which we are working
   adxy, header, ra, dec, xcen, ycen
   ;keep track of where the two stars are
   xcenarr_1(i) = xcen(0)
   ycenarr_1(i) = ycen(0)
   xcenarr_2(i) = xcen(1)
   ycenarr_2(i) = ycen(1)

endfor

fluxarr_1 = fltarr(nimages,nimages)
fluxarr_2 = fltarr(nimages,nimages)

;measure the flux at each of those positions
for i = 0, nimages - 1 do begin
   fits_read, bcdname(i), data, header
   ;star 1
   aper,data, xcenarr_1, ycenarr_1, flux , fluxerr, sky, $
       skyerr, 3.3, [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent
;   flux(where(finite(flux) lt 1)) = 1
   fluxarr_1(i,*) = flux
  ;star 2
   aper,data, xcenarr_2, ycenarr_2, flux2 , fluxerr2, sky2, $
       skyerr2, 3.3, [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent
   nanflux = where(finite(flux2) lt 1)

   fluxarr_2(i,*) = flux2
 
endfor
print, ' fluxarr_2', fluxarr_2
print, max(fluxarr_2), 'max'

plot, findgen(nimages), fluxarr_1(0,*), yrange= [1,max(fluxarr_1)],/ylog, thick = 3, charthick = 3, $
      xthick = 3, ythick = 3, title = 'STAR 1'
oplot, findgen(nimages), fluxarr_1(0,*), psym = 2, thick = 3

for i =1, nimages - 2 do begin
   oplot, findgen(nimages), fluxarr_1(i,*), color = colorname[i], thick = 3
   oplot, findgen(nimages), fluxarr_1(i,*), psym = 2, color = colorname[i], thick = 3
endfor



plot, findgen(nimages), fluxarr_2(0,*), yrange= [1,max(fluxarr_2)],/ylog, thick = 3, charthick = 3, $
      xthick = 3, ythick = 3, title = 'STAR 2'
oplot, findgen(nimages), fluxarr_2(0,*), psym = 2, thick = 3

for i =1, nimages - 2 do begin
   oplot, findgen(nimages), fluxarr_2(i,*), color = colorname[i], thick = 3
   oplot, findgen(nimages), fluxarr_2(i,*), psym = 2, color = colorname[i], thick = 3
endfor
endfor

ps_close, /noprint,/noid

end

