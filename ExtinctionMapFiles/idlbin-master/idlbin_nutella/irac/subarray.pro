pro subarray, data


;take a subarray raw frame 256x256 and covert that into a 32x32x64
;image.

data3d = fltarr(32, 32, 64)
;data = findgen(256,256)

for j = 0, 255 do begin
   for i = 0,255 do begin
;         print, 'i,j', i, j
;         print, 'data(i,j)', data(i,j)
;         print, 'I,J,K' , i mod 32, 8*(j mod 4) + floor(i/32), floor(j/4)
         data3d[i mod 32, 8*(j mod 4.) + floor(i/32.), floor(j/4.)] = data(i,j)

   endfor
endfor


;now I want to measure the std dev within each pixel that gets imaged
;64 times
c = long(0)
stddevarr = fltarr(1024)
print, n_elements(stddevarr)
for i = 0, 31 do begin
   for j = 0, 31 do begin
      print, 'begin', c, i, j, stddev(data3d[i,j,*])
      stddevarr(c) =  stddev(data3d[i,j,*])
      c = c +1
   endfor
endfor


;plot the histogram of standard deviations
plothist, stddevarr

end
