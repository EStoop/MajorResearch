pro cosmic_removerepeats

readcol,'/Users/jkrick/palomar/cosmic/slitmasks/june08/masklist.txt' ,masternum,h,m,s,dd,dm,ds,rmag, format='A'

for i = 0, n_elements(masternum) - 1 do begin
   j = where(s eq s(i))
   a = float(masternum(j))
   print, masternum(i), masternum(j)
   masternum(i) = long(min(float(masternum(j))))
   print, masternum(i)
endfor


openw, outlun, '/Users/jkrick/palomar/cosmic/slitmasks/june08/masklist_norepeat.txt',/get_lun
for k = 0, n_elements(masternum) - 1 do begin
   printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',masternum(k),h(k),m(k),s(k),dd(k),dm(k),ds(k),rmag(k)
endfor

close, outlun
free_lun, outlun



end
