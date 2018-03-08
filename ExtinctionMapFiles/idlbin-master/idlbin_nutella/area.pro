pro area

fits_read, '/Users/jkrick/hst/raw/wholeacs.fits', data, header

adxy, header, 264.68160, 69.04481, clusterx, clustery
print, clusterx, clustery
x = lindgen(3000) + 18945
y = lindgen(3000) + 12216

;a = x(20440)
;b = y(13700)
;print, sqrt( (clustery - b)^2 + (clusterx - a)^2 )
if  data[21431,13475] gt 0. then print, 'zero gt 0 = problem'
if  data[20363,13799] gt 0. then print, 'nonzero gt 0 = good'

hasdata = long(0)
allarea = long(0)
for i = 0, n_elements(x) -1  do begin
   for j = 0, n_elements(y) - 1 do begin
      if sqrt( (clustery - y(j))^2 + (clusterx - x(i))^2 ) lt 1224. and data[x(i),y(j)] eq 0 then hasdata = hasdata + 1
      if sqrt( (clustery - y(j))^2 + (clusterx - x(i))^2 ) lt 1224. then allarea = allarea + 1
   endfor
endfor

;hasdata = where(sqrt( (clustery - y)^2 + (clusterx - x)^2 ) lt 1224. and data[x,y] gt 0.)
;allarea = where( sqrt( (clustery - y)^2 + (clusterx - x)^2 ) lt 1224.)

print,hasdata, allarea
end
