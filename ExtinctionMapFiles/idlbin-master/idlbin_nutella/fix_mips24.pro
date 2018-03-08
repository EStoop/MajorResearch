pro fix_mips24

restore, '/Users/jkrick/idlbin/objectnew.sav'

newflux = [1251,1253,1187,1868,1469,1472,917,939,952,1765,993,1768,1005,1006,1007,1042,1501,1064,1101,1092,79395,2405,14648,1825,1202,1227,1215,1228,2368,1324,2257,1530,1847,1279,2662,2629,2037]

oldflux = [33904,34349,44719,45263,73365,73840,74008,74878,74897,76126,76278,76338,76400,76680,76859,77502,77526,78242,78729,79276,79340,79460,79705,81250,81697,82093,82146,82556,82908,83340,83963,84418,84437,84897,85269,85574,85945]

print, 'equal?' , n_elements(newflux), n_elements(oldflux)

for i = 0, n_elements(newflux) -1 do begin
   objectnew[newflux(i)].mips24ra = objectnew[oldflux(i)].mips24ra
   objectnew[newflux(i)].mips24dec = objectnew[oldflux(i)].mips24dec
   objectnew[newflux(i)].mips24xcenter = objectnew[oldflux(i)].mips24xcenter 
   objectnew[newflux(i)].mips24ycenter = objectnew[oldflux(i)].mips24ycenter 
   objectnew[newflux(i)].mips24flux = objectnew[oldflux(i)].mips24flux 
   objectnew[newflux(i)].mips24fluxerr = objectnew[oldflux(i)].mips24fluxerr 
   objectnew[newflux(i)].mips24mag = objectnew[oldflux(i)].mips24mag 
   objectnew[newflux(i)].mips24magerr = objectnew[oldflux(i)].mips24magerr
   objectnew[newflux(i)].mips24bckgrnd = objectnew[oldflux(i)].mips24bckgrnd 
   objectnew[newflux(i)].iter = objectnew[oldflux(i)].iter
   objectnew[newflux(i)].sharp = objectnew[oldflux(i)].sharp 
   objectnew[oldflux(i)].mips24ra =0
   objectnew[oldflux(i)].mips24dec = 0
   objectnew[oldflux(i)].mips24xcenter = 0
   objectnew[oldflux(i)].mips24ycenter = 0
   objectnew[oldflux(i)].mips24flux = -99
   objectnew[oldflux(i)].mips24fluxerr = -99
   objectnew[oldflux(i)].mips24mag = 99
   objectnew[oldflux(i)].mips24magerr =99
   objectnew[oldflux(i)].mips24bckgrnd = 0
   objectnew[oldflux(i)].iter = 0
   objectnew[oldflux(i)].sharp = 0
endfor

print, 'test', objectnew[14648].mips24flux, objectnew[79705].mips24flux

;save, objectnew, filename='/Users/jkrick/idlbin/objectnew.sav'

end
