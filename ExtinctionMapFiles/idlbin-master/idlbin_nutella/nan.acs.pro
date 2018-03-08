pro nan_acs
filelist = '/Users/jkrick/hst/drzlist'
NaN = alog10(-1)
readcol,filelist,filename,format="A"

print, filename[0]

for j=0,n_elements(filename)-1 do begin

i = 1
print, j, i, "j, i"
data1 =readfits(filename(j), head1)

data1(where(data1 eq 0)) = NaN

print, data1[1000,1]
writefits,filename(j),data1, head1

endfor

end
