pro size_sigma_single
nan = alog10(-1)
ch = 0

if ch eq 0 then begin
   fits_read, '/Users/jkrick/Virgo/IRAC/s18p14/r35320320/ch1/bcd/SPITZER_I1_35320320_0010_0000_ndark_ff_match.fits', data, header
endif

if ch eq 1 then begin
   print, '****ch2*****'
endif

;need to mask.
CD, '/Users/jkrick/Virgo/IRAC/s18p14/r35320320/ch1/bcd'
command = 'sex SPITZER_I1_35320320_0010_0000_ndark_ff_match.fits -c ./default.sex'
spawn, command

fits_read,  '/Users/jkrick/Virgo/IRAC/s18p14/r35320320/ch1/bcd/segmentation.fits', segdata, segheader

;try converting the data image into electrons.
;sbtoe = gain*exptime/flux_conv
;Mjypersr * sbtoe = electrons 
;sbtoe = (covdata*100) * 3.7 / .1198
;data = data * sbtoe

;first use the segmentation image as a mask
   segmask = where(segdata gt 0)
   data(segmask) = nan
   fits_write, strcompress('/Users/jkrick/Virgo/IRAC/s18p14/r35320320/ch1/bcd/ch1_10_seg.fits',/remove_all), data, header

;xcen = [996,668,1427,1338,591,1391,2369,3363,1149]
;ycen = [5399,4427,4323,3560,3034,2938,1852,1836,949]

xcen = randomu(12, 9)*100. + 50.
ycen = randomu(120, 9)*100 + 50.
print, 'xcen', xcen, 'ycen', ycen

e = 50
sigmaarr = fltarr(e)
sizearr = fltarr(e)

for n = 1, e do begin 
   size = 1*n
   sumarr = fltarr(n_elements(xcen))
   for i = 0, n_elements(xcen) - 1 do begin 
      goodarea = where(finite(data[xcen[i] - size/2.:xcen[i] + size/2., ycen[i] -size/2.:ycen[i] + size/2.]) gt 0, goodcount)
      sumarr(i) = total(data[xcen[i] - size/2.:xcen[i] + size/2., ycen[i] -size/2.:ycen[i] + size/2.],/nan) / goodcount
   endfor
   
   sigmaarr(n-1) = stddev(sumarr)
   sizearr(n-1) = (sqrt(goodcount))*0.6
endfor

ps_start, filename= '/Users/jkrick/Virgo/IRAC/S18p14/size_sigma_single.ps'

!P.multi = [0,1,1]
plot, sizearr, sigmaarr, thick = 3, charthick = 3,/ylog, ytitle = 'Noise in Mjy/sr/pixel', xtitle = 'Binning Length (arcsec on bottom, native pixels on top)',  psym = 2, /xlog, xstyle = 9, xrange = [1,100], xthick = 3, ythick =3

;axis, xaxis = 1, xrange = [1/1.2,1000/1.2], xstyle = 1,/xlog, xthick = 3

y = 0.006/sizearr
oplot, sizearr, y, thick = 3
print, 'sizearr', sizearr
print, '0.006/ sizearr', 0.006/ sizearr
print, 'sigmaarr',sigmaarr
ps_end
end
