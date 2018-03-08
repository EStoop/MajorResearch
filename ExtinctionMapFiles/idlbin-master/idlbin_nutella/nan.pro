pro nan, filelist, bpmlist, outlist
NaN = alog10(-1)

readcol,filelist,mosaics,format="A"
readcol,bpmlist,masks,FORMAT="A"
readcol,outlist,procs,FORMAT="A"


nfiles=n_elements(mosaics)
for j=0,nfiles-1 do begin

filename = mosaics[j]           ;"/Users/jkrick/palomar/LFC/22jul04/FFdata/ccd.050.fits"
bpmname = masks[j]            ;"/Users/jkrick/palomar/LFC/CRmasks/bpm22jul04.050.fits"
outname = procs[j]              ;"test.fits"
i = 1

bpmname = strmid(bpmname, 0,50) + strcompress("_"+string(i) + ".fits",/remove_all)
data1 =readfits(filename, head1, exten=1)
mask1 = readfits(bpmname, hd)
data1(where(mask1 ne 0)) = NaN

bpmname = strmid(bpmname, 0,48) + strcompress("_"+string(i+1) + ".fits",/remove_all)
data2 =readfits(filename, head2, exten=2)
mask1 = readfits(bpmname, hd)
data2(where(mask1 ne 0)) = NaN

bpmname = strmid(bpmname, 0,48) + strcompress("_"+string(i+2) + ".fits",/remove_all)
print, bpmname
data3 =readfits(filename, head3, exten=3)
mask1 = readfits(bpmname, hd)
data3(where(mask1 ne 0)) = NaN

bpmname = strmid(bpmname, 0,48) + strcompress("_"+string(i+3) + ".fits",/remove_all)
data4 =readfits(filename, head4, exten=4)
mask1 = readfits(bpmname, hd)
data4(where(mask1 ne 0)) = NaN

bpmname = strmid(bpmname, 0,48) + strcompress("_"+string(i+4) + ".fits",/remove_all)
data5 =readfits(filename, head5, exten=5)
mask1 = readfits(bpmname, hd)
data5(where(mask1 ne 0)) = NaN

bpmname = strmid(bpmname, 0,48) + strcompress("_"+string(i+5) + ".fits",/remove_all)
data6 =readfits(filename, head6, exten=6)
mask1 = readfits(bpmname, hd)
data6(where(mask1 ne 0)) = NaN


print, data1(26,951)


mkhdr, head0,",/exten

writefits,outname, ",head0
writefits,outname,data1, head1,/app
writefits,outname,data2, head2,/app
writefits,outname,data3, head3,/app
writefits,outname,data4, head4,/app
writefits,outname,data5, head5,/app
writefits,outname,data6, head6,/app

endfor

end
