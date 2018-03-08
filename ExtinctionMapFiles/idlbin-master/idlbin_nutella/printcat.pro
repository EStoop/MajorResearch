pro printcat
;device, true=24
;device, decomposed=0

close, /all


restore, '/Users/jkrick/idlbin/object.sav'
count = 0
openw, outlun, '/Users/jkrick/nep/catalog.txt', /get_lun
for i =0, n_elements(object.rmaga) - 1 do begin
   
         printf,outlun, count, object[i].ra, object[i].dec, object[i].umaga, object[i].umagerra, object[i].gmaga, object[i].gmagerra,object[i].rmaga, object[i].rmagerra,object[i].imaga, object[i].imagerra, object[i].flamjmag, object[i].flamjmagerr, object[i].wircjmag, object[i].wircjmagerr, object[i].wirchmag, object[i].wirchmagerr, object[i].wirckmag, object[i].wirckmagerr, object[i].irac1,object[i].irac2,object[i].irac3,object[i].irac4, object[i].mips24flux, object[i].mips24fluxerr,format = '(I10,F10.5,F10.5,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)'

count = count + 1 


endfor
print, count
print, object[22941].imaga
close, outlun
free_lun, outlun
END

