pro readbc03
close, /all


model = fltarr(7,1086,30)

readcol,'/Users/jkrick/bin/bc03/filelist',filename,format="A",/silent
filename = strcompress(filename + ",0.5,1,3,5,8,12")

for i = 0, n_elements(filename) -1 do begin
   readcol, filename(i), lambda, age0, age1, age2,age3,age4,age5,age6,format="A",/silent

   model[0,*,i] = age0
   model[1,*,i] = age1
   model[2,*,i] = age2
   model[3,*,i] = age3
   model[4,*,i] = age4
   model[5,*,i] = age5
   model[6,*,i] = age6
endfor




end
