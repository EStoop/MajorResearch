pro cosmic_slitdone

filenames =['/Users/jkrick/palomar/cosmic/slitmasks/june08/mask1.txt', '/Users/jkrick/palomar/cosmic/slitmasks/june08/mask2.txt', '/Users/jkrick/palomar/cosmic/slitmasks/june08/mask3.txt', '/Users/jkrick/palomar/cosmic/slitmasks/june08/mask4.txt', '/Users/jkrick/palomar/cosmic/slitmasks/june08/mask5.txt', '/Users/jkrick/palomar/cosmic/slitmasks/june08/mask6.txt', '/Users/jkrick/palomar/cosmic/slitmasks/june08/mask7.txt']
donenum = [0,1]

for i = 0, n_elements(filenames) - 1 do begin
   readcol, filenames(i) ,donenum0, junk, junk,junk, junk,junk, junk, junk,junk, junk,junk ,format="A"
   donenum = [donenum, donenum0]
endfor

donenum = donenum(where(donenum gt 99))

print, donenum

readcol,'/Users/jkrick/palomar/cosmic/slitmasks/june08/masklist_norepeat.txt' ,masternum,h,m,s,dd,dm,ds,rmag ,format="A"
print, n_elements(masternum), n_elements(donenum)

for j = 0, n_elements(donenum) - 1 do begin
   h = h(where(masternum ne donenum(j)) )
   m = m(where(masternum ne donenum(j)) );
   s = s(where(masternum ne donenum(j)) )
   dd = dd(where(masternum ne donenum(j)) )
   dm = dm(where(masternum ne donenum(j)) )
   ds = ds(where(masternum ne donenum(j)) )
   rmag = rmag(where(masternum ne donenum(j)) )
   masternum = masternum(where(masternum ne donenum(j)) )

endfor

print, n_elements(masternum)
openw, outlun, '/Users/jkrick/palomar/cosmic/slitmasks/june08/masklist8.txt', /get_lun
for k = 0, n_elements(masternum) - 1 do begin
   printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',masternum(k),h(k),m(k),s(k),dd(k),dm(k),ds(k),rmag(k)
endfor

close, outlun
free_lun, outlun
end
