function cluster_mips24_sample

common sharevariables

cluster24 = fltarr(400)
clusterall = fltarr(8000)
dist24 = fltarr(400)
distall = fltarr(8000)

openw, outlunred, '/Users/jkrick/nep/clusters/mips24/all.reg', /get_lun
printf, outlunred, 'fk5'


j = 0
l = 0
;for all clusters
for i = 0, n_elements(candra) - 1 do begin

   ;those objects within cluster area that have mips detections
   good = where(sphdist(objectnew.ra, objectnew.dec, candra[i], canddec[i], /degrees) lt rvir(i) $
                and objectnew.mips24flux ge 17.3 and objectnew.mips24mag lt 90  )
;   print, 'good', good

   ;sort based on rband magnitude for comparison purposes.
   good = good(sort(objectnew[good].mips24flux))
   if i eq 0 then begin
      filename = '/Users/jkrick/nep/clusters/mips24/newastrometry/mips24_0.ps'
      textname = '/Users/jkrick/nep/clusters/mips24/newastrometry/cluster0.txt'
   endif
   if i eq 1 then begin
      filename = '/Users/jkrick/nep/clusters/mips24/newastrometry/mips24_1.ps'
      textname = '/Users/jkrick/nep/clusters/mips24/newastrometry/cluster1.txt'
   endif
   if i eq 2 then begin
      filename = '/Users/jkrick/nep/clusters/mips24/newastrometry/mips24_2.ps'
      textname = '/Users/jkrick/nep/clusters/mips24/newastrometry/cluster2.txt'
   endif

;   plothyperz, good, filename

   openw, textlun,textname, /get_lun
   for j = 0, n_elements(good) -1 do begin
      printf, textlun,  format='(I10,F10.5,F10.5,F10.2,F10.2,F10.2,F10.2,F10.2)', good[j],  objectnew[good[j]].ra, objectnew[good[j]].dec, objectnew[good[j]].zphot, objectnew[good[j]].rmag,objectnew[good[j]].acsmag, objectnew[good[j]].irac1mag, objectnew[good[j]].mips24flux
      printf, outlunred, 'circle( ', objectnew[good[j]].ra, objectnew[good[j]].dec, ' 3")'

   endfor
   close, textlun
   free_lun, textlun

   ;--------------------------------------------------------
   ;those objects within cluster area that have irac detections and NO mips detections
   goodall = where(sphdist(objectnew.ra, objectnew.dec, candra[i], canddec[i], /degrees) lt sep $
                   and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.mips24mag gt 90)

   ;make combined arrays for all three clusters, and distance arrays
   for k = 0, n_elements(good) -1 do begin
      cluster24[j] = good[k]
      dist24[j] = sphdist(objectnew[good[k]].ra, objectnew[good[k]].dec, candra(i), canddec(i), /degrees) 
      j = j + 1
   endfor
   for m = 0, n_elements(goodall) -1 do begin
      clusterall[l] = goodall[m]
      distall[l] = sphdist(objectnew[goodall[m]].ra, objectnew[goodall[m]].dec, candra(i), canddec(i), /degrees) 
      l = l + 1
   endfor


endfor

cluster24 = cluster24[0:j-1]
dist24 = dist24[0:j-1]
clusterall = clusterall[0:l-1]
distall = distall[0:l-1]


;sort good on rmag
cluster24 = cluster24[sort(objectnew[cluster24].rmag)]
clusterall = clusterall[sort(objectnew[clusterall].rmag)]

;plothyperz, cluster24, '/Users/jkrick/nep/clusters/mips24/mips24_cluster.ps'

;after looking at the output of all of this, the final lists are
;see /Users/jkrick/nep/clusters/mips24/newastrometry/cluster*.xls
;allmips = all objects except obvious stars
;allgals = all objects except obvious stars and sey or QSO and noise and impossible matches
; allcluster= all objects except obvious stars and sey or QSO  and noise and impossible matches and they have phot-z's between 0.8 and 1.3
;this is a broad generous range.


allmips = [1039,2314,2394,984,1983,958,999,1753,942,2393,1984,1985,74320,2395,1752,985,12401,12611,940,998,1075,912,13342,2228,2133,1060,1076,1026,2978,1759,75343,1490,880,921,879,1079,1998,1986,13901,1505,2145,13913,1040,1496,1987,1742,1751,986,1493,1023,2360,1497,1474,1491,1790,1495,2750,74900,73709,922, 2007,3124,1073,2160,1776,2010,2034,2536,1818,1533,79340,1226,2251,1258,2246,2036,81250,1010,1500,1087,1144,1090,1525,1507,1817,2941,2035,32717,2030,1145,2006,1244,83580,2252,2161,1141,1519,1830,2029,1043,1092,2008,1531,1088,1178,2026,84553,15680,2025,32364,1823,1265,1821,15404,14740,2162,15426,1523,2253,2028,2027,1176,1256,1091,2881,15688,2009,1044,1179,1177,1222,32283,2005,1826,1175,1168,2323,15647,1796,78504,81755,78984,1257,1995,1181,1055,2152,1235,32544,1250,2326,1524,1840,1196,1094]


allgals = [2007,1073,2160,1776,2010,2034,1818,1533,1226,2251,1258,2036,1010,1500,1087,1144,1090,1525,1507,2941,2030,1145,1244,83580,2252,2161,1141,1519,1830,1043,1092,2008,1531,1088,1178,84553,15680,2025,1265,1821,15404,14740,2162,1523,2253,2027,1176,1256,1091,15688,1044,1179,1177,1222,1175,2323,1257,1995,1181,1055,1235,32544,1250,1524,1840,1094,1039,2314,1983,958,942,2393,1984,1985,74320,985,12401,12611,940,1075,912,13342,2228,2133,1026,1759,1490,921,879,1998,1986,1505,2145,13913,1040,1496,1987,1742,1751,986,1023,2360,1497,1474,1491,1790,1495,922]


allcluster = [1983,958,942,2393,1984,12611,1075,912,13342,1026,1490,1505,13913,1040,1496,1987,1742,1751,986,2360,1491,1495,922,2160,1776,1533,1226,2251,1010,1500,1090,1525,2030,1244,1830,2008,1531,1178,2162,1523,2027,1176,1091,15688,1044,1055,32544,1250,1840]

print, 'allmips', n_elements(allmips)
print, 'allgals', n_elements(allgals)
print, 'allcluster', n_elements(allcluster)

;figure out what the distances to the cluster center are for each
;object in the three different samples.  This is more complex since I
;have not so far kept track of which cluster they are in, so first
;test to see which one(or two) they belong to.
distallmips = fltarr(n_elements(allmips) * 2)
j = 0
for k = 0, n_elements(allmips) -1 do begin
   zero =   sphdist(objectnew[allmips[k]].ra, objectnew[allmips[k]].dec, candra(0), canddec(0), /degrees) 
   one =   sphdist(objectnew[allmips[k]].ra, objectnew[allmips[k]].dec, candra(1), canddec(1), /degrees) 
   two =   sphdist(objectnew[allmips[k]].ra, objectnew[allmips[k]].dec, candra(2), canddec(2), /degrees) 

   if zero le sep then begin
      distallmips[j] = zero
      j = j + 1
   endif
   if one le sep then begin
      distallmips[j] = one
      j = j + 1
   endif
   if two le sep then begin
      distallmips[j] = two
      j = j + 1
   endif
endfor
distallmips = distallmips[0:j - 1]

distallgals = fltarr(n_elements(allgals) * 2)
j = 0
for k = 0, n_elements(allgals) -1 do begin
   zero =   sphdist(objectnew[allgals[k]].ra, objectnew[allgals[k]].dec, candra(0), canddec(0), /degrees) 
   one =   sphdist(objectnew[allgals[k]].ra, objectnew[allgals[k]].dec, candra(1), canddec(1), /degrees) 
   two =   sphdist(objectnew[allgals[k]].ra, objectnew[allgals[k]].dec, candra(2), canddec(2), /degrees) 
   if zero le sep then begin
      distallgals[j] = zero
      j = j + 1
   endif
   if one le sep then begin
      distallgals[j] = one
      j = j + 1
   endif
   if two le sep then begin
      distallgals[j] = two
      j = j + 1
   endif
endfor
distallgals = distallgals[0:j - 1]

distallcluster = fltarr(n_elements(allcluster) * 2)
j = 0
for k = 0, n_elements(allcluster) -1 do begin
   zero =   sphdist(objectnew[allcluster[k]].ra, objectnew[allcluster[k]].dec, candra(0), canddec(0), /degrees) 
   one =   sphdist(objectnew[allcluster[k]].ra, objectnew[allcluster[k]].dec, candra(1), canddec(1), /degrees) 
   two =   sphdist(objectnew[allcluster[k]].ra, objectnew[allcluster[k]].dec, candra(2), canddec(2), /degrees) 
   if zero le sep then begin
      distallcluster[j] = zero
      j = j + 1
   endif
   if one le sep then begin
      distallcluster[j] = one
      j = j + 1
   endif
   if two le sep then begin
      distallcluster[j] = two
      j = j + 1
   endif
endfor
distallcluster = distallcluster[0:j - 1]

close, outlunred
free_lun, outlunred
return, 0
end

;  if i eq 0 then save, good, filename='/Users/jkrick/nep/clusters/mips24/good_0.sav'
;   if i eq 1 then save, good, filename='/Users/jkrick/nep/clusters/mips24/good_1.sav'
;   if i eq 2 then save, good, filename='/Users/jkrick/nep/clusters/mips24/good_2.sav'

