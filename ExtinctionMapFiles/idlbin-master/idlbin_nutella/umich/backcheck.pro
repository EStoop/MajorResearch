PRO backcheck
close, /all
;device, true=24
;device, decomposed=0

sig = fltarr(27)
yaxis = fltarr(30)
xaxis = fltarr(30)
fits_read, "/n/Godiva1/jkrick/A3888/A3888V.div.fits", data, header
;fits_read, "/n/Godiva1/jkrick/A3888/largeV.block.fits", data, header
;find average stddev of 25 "blank" regions with no smoothing
;sig(0) = stddev(data[2326:2345,655:674])
;sig(1) = stddev(data[2256:2275,444:463])
;sig(2) = stddev(data[1869:1888,440:459])
;sig(3) = stddev(data[1841:1860,780:799])
;sig(4) = stddev(data[2115:2134,1009:1028])
;sig(5) = stddev(data[1816:1835,794:813])
;sig(6) = stddev(data[2170:2189,1195:1214])
;sig(7) = stddev(data[1681:1700,2032:2051])
;sig(8) = stddev(data[2025:2044,1974:1993])
;sig(9) = stddev(data[2217:2236,2226:2245])
;sig(10) = stddev(data[2377:2396,2173:2192])
;sig(11) = stddev(data[2285:2304,1870:1889])
;sig(12) = stddev(data[2340:2359,1458:1477])
;sig(13) = stddev(data[2107:2126,1381:1400])
;sig(14) = stddev(data[2408:2427,1366:1385])
;sig(15) = stddev(data[1182:1201,2210:2229])
;sig(16) = stddev(data[796:815,2030:2049])
;sig(17) = stddev(data[1004:1023,1997:2016])
;sig(18) = stddev(data[599:618,1956:1975])
;sig(19) = stddev(data[525:544,1894:1913])
;sig(20) = stddev(data[613:632,1649:1668])
;sig(21) = stddev(data[562:581,1505:1524])
;sig(22) = stddev(data[826:845,1659:1678])
;sig(23) = stddev(data[865:884,1179:1198])
;sig(24) = stddev(data[892:911,868:887])
;sig(25) = stddev(data[1088:1107,593:612])
;sig(26) = stddev(data[951:970,497:516])
sig(0) = stddev(data[2826:2845,1355:1374])  
sig(1) = stddev(data[2756:2775,1144:1163])  
sig(2) = stddev(data[2369:2388,1140:1159])  
sig(3) = stddev(data[2341:2360,1480:1499])  
sig(4) = stddev(data[2615:2634,1709:1728])  
sig(5) = stddev(data[2316:2335,1494:1513])  
sig(6) = stddev(data[2670:2689,1895:1914])  
sig(7) = stddev(data[2181:2200,2732:2751])  
sig(8) = stddev(data[2525:2544,2674:2693])  
sig(9) = stddev(data[2717:2736,2926:2945])  
sig(10) = stddev(data[2877:2896,2873:2892])  
sig(11) = stddev(data[2785:2804,2570:2589])  
sig(12) = stddev(data[2840:2859,2158:2177])  
sig(13) = stddev(data[2607:2626,2081:2100])  
sig(14) = stddev(data[2908:2927,2066:2085])  
sig(15) = stddev(data[1682:1701,2910:2929])  
sig(16) = stddev(data[1296:1315,2730:2749])  
sig(17) = stddev(data[1504:1523,2697:2716])  
sig(18) = stddev(data[1099:1118,2656:2675])  
sig(19) = stddev(data[1025:1044,2594:2613])  
sig(20) = stddev(data[1113:1132,2349:2368])  
sig(21) = stddev(data[1062:1081,2205:2224])  
sig(22) = stddev(data[1326:1345,2359:2378])  
sig(23) = stddev(data[1365:1384,1879:1898])  
sig(24) = stddev(data[1392:1411,1568:1587])  
sig(25) = stddev(data[1588:1607,1293:1312])  
sig(26) = stddev(data[1451:1470,1197:1216])  



yaxis(0) = mean(sig)
xaxis(0) = 1

;now smooth on varying scales and recalculate sigma every time
c = 1
FOR i = 3,59,2 DO BEGIN

    data2 = median(data, i)
    fits_write, "/n/Godiva1/jkrick/A3888/background/testsmooth.fits", data2, header
 ;   sig(0) = stddev(data2[2326:2345,655:674])
 ;   sig(1) = stddev(data2[2256:2275,444:463])
 ;   sig(2) = stddev(data2[1869:1888,440:459])
 ;   sig(3) = stddev(data2[1841:1860,780:799])
 ;   sig(4) = stddev(data2[2115:2134,1009:1028])
 ;   sig(5) = stddev(data2[1816:1835,794:813])
 ;   sig(6) = stddev(data2[2170:2189,1195:1214])
 ;   sig(7) = stddev(data2[1681:1700,2032:2051])
 ;   sig(8) = stddev(data2[2025:2044,1974:1993])
 ;   sig(9) = stddev(data2[2217:2236,2226:2245])
 ;   sig(10) = stddev(data2[2377:2396,2173:2192])
 ;   sig(11) = stddev(data2[2285:2304,1870:1889])
 ;   sig(12) = stddev(data2[2340:2359,1458:1477])
 ;   sig(13) = stddev(data2[2107:2126,1381:1400])
 ;   sig(14) = stddev(data2[2408:2427,1366:1385])
 ;   sig(15) = stddev(data2[1182:1201,2210:2229])
 ;   sig(16) = stddev(data2[796:815,2030:2049])
 ;   sig(17) = stddev(data2[1004:1023,1997:2016])
 ;   sig(18) = stddev(data2[599:618,1956:1975])
 ;   sig(19) = stddev(data2[525:544,1894:1913])
 ;   sig(20) = stddev(data2[613:632,1649:1668])
 ;   sig(21) = stddev(data2[562:581,1505:1524])
 ;   sig(22) = stddev(data2[826:845,1659:1678])
 ;   sig(23) = stddev(data2[865:884,1179:1198])
 ;   sig(24) = stddev(data2[892:911,868:887])
 ;   sig(25) = stddev(data2[1088:1107,593:612])
 ;   sig(26) = stddev(data2[951:970,497:516])

    sig(0) = stddev(data2[2826:2845,1355:1374])  
    sig(1) = stddev(data2[2756:2775,1144:1163])  
    sig(2) = stddev(data2[2369:2388,1140:1159])  
    sig(3) = stddev(data2[2341:2360,1480:1499])  
    sig(4) = stddev(data2[2615:2634,1709:1728])  
    sig(5) = stddev(data2[2316:2335,1494:1513])  
    sig(6) = stddev(data2[2670:2689,1895:1914])  
    sig(7) = stddev(data2[2181:2200,2732:2751])  
    sig(8) = stddev(data2[2525:2544,2674:2693])  
    sig(9) = stddev(data2[2717:2736,2926:2945])  
    sig(10) = stddev(data2[2877:2896,2873:2892])  
    sig(11) = stddev(data2[2785:2804,2570:2589])  
    sig(12) = stddev(data2[2840:2859,2158:2177])  
    sig(13) = stddev(data2[2607:2626,2081:2100])  
    sig(14) = stddev(data2[2908:2927,2066:2085])  
    sig(15) = stddev(data2[1682:1701,2910:2929])  
    sig(16) = stddev(data2[1296:1315,2730:2749])  
    sig(17) = stddev(data2[1504:1523,2697:2716])  
    sig(18) = stddev(data2[1099:1118,2656:2675])  
    sig(19) = stddev(data2[1025:1044,2594:2613])  
    sig(20) = stddev(data2[1113:1132,2349:2368])  
    sig(21) = stddev(data2[1062:1081,2205:2224])  
    sig(22) = stddev(data2[1326:1345,2359:2378])  
    sig(23) = stddev(data2[1365:1384,1879:1898])  
    sig(24) = stddev(data2[1392:1411,1568:1587])  
    sig(25) = stddev(data2[1588:1607,1293:1312])  
    sig(26) = stddev(data2[1451:1470,1197:1216])  
    
    
    yaxis(c) = mean(sig)
    xaxis(c) = i
    c = c + 1
ENDFOR

;ps_open, file = "/n/Godiva1/jkrick/A3888/backcheck.ps", /portrait, xsize = 6, ysize = 6
mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A3888/background/backcheckVtest.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plot, xaxis* 0.259, 24.3- 2.5*alog10(yaxis/(0.259^2)), linestyle = 0,title = "Combined V image",$
thick = 3, charthick = 3, xthick = 3, ythick = 3, xtitle = "smoothing length (arcsec)", $
ytitle = "stddev(mag/arcsec^2)", yrange =[24,34],/xlog
print, yaxis

oplot, xaxis*0.259, 24.3 - 2.5*alog10(1/(22* xaxis* (0.259^2))), linestyle = 2
;ps_close, /noprint, /noid
device, /close
set_plot, mydevice

END

;sig(0) = stddev(data[1393:1397,1898:1902])  
;sig(1) = stddev(data[1565:1569,1522:1526])  
;sig(2) = stddev(data[1851:1855,1231:1235])  
;sig(3) = stddev(data[1114:1118,2341:2345])  
;sig(4) = stddev(data[1188:1192,2542:2546])  
;sig(5) = stddev(data[1296:1300,2747:2751])  
;sig(6) = stddev(data[1560:1564,2612:2616])  
;sig(7) = stddev(data[1698:1702,2925:2929])  
;sig(8) = stddev(data[1724:1728,3106:3110])  
;sig(9) = stddev(data[1831:1835,3591:3595])  
;sig(10) = stddev(data[2888:2892,2913:2917])  
;sig(11) = stddev(data[2956:2960,2773:2777])  
;sig(12) = stddev(data[2804:2808,2569:2573])  
;s;ig(13) = stddev(data[2528:2532,2743:2747])  
;s;ig(14) = stddev(data[2786:2790,2587:2591])  
;sig(15) = stddev(data[2960:2964,2274:2278])  
;sig(16) = stddev(data[2620:2624,2077:2081])  
;sig(17) = stddev(data[2921:2925,2069:2073])  
;sig(18) = stddev(data[2765:2769,1818:1822])  
;sig(19) = stddev(data[2362:2366,1429:1433])  
;sig(20) = stddev(data[2210:2214,1257:1261])  
;sig(21) = stddev(data[2350:2354,1098:1102])    
;sig(22) = stddev(data[2800:2804,1247:1251])    
;sig(23) = stddev(data[2779:2783,1419:1423])    
;sig(24) = stddev(data[2945:2949,1489:1493])    
