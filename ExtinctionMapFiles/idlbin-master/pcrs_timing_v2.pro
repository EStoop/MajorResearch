pro pcrs_timing_v2
!p.multi = [0, 1, 2]
ps_open, filename='/Users/jkrick/irac_warm/pcrsstare/fancy_timing_y.ps',/portrait,/square,/color


;CH1 SUBARRAY

ch1_sub_sclk = [999405565.D,999444230.D,999693661.D,999732318.D,999751864.D,999756044.D,999760507.D,999769922.D,999771702.D,999781845.D]
ch1_sub_x=[14.875,14.851,14.909,14.842,14.924,14.843,14.869,14.864,14.874,14.922]
ch1_sub_y=[15.005,15.097,15.027,15.115,15.097,15.098,15.128,15.128,15.096,15.210]
ch1_sub_XMODE=      14.8813
ch1_sub_YMODE =      15.1101
ch1_sub_dayarr = ch1_sub_sclk
for i = 0, n_elements(ch1_sub_sclk) - 1 do begin
   year = (ch1_sub_sclk[i] / 60./60./24./365.) + 1980.
   day = year - fix(year)
   day = day *365.
   ch1_sub_dayarr[i]  = day
endfor

;plot, ch1_sub_dayarr, ch1_sub_x- ch1_sub_xmode, psym = 1, xtitle = 'Time (day of year)', xrange = [250, 280], yrange = [-0.06, 0.06], ytitle = 'X distance from mode', title = 'Ch1'
plot, ch1_sub_dayarr, ch1_sub_y- ch1_sub_ymode, psym = 1, xtitle = 'Time (day of year)', xrange = [245, 280], yrange = [-0.06, 0.06], ytitle = 'Y distance from mode', title = 'Ch1'

;CH1 FULL ARRAY

ch1_full_sclk=[  999546815.D,  999553526.D,  999577309.D,  999644194.D,  999686640.D,  999689118.D,  999782053.D,  999904668.D, 1000072460.D,  1000440302.D]
ch1_full_x=[ 129.879    , 129.901    , 129.897    , 129.941    , 129.937    , 129.947    , 129.889    , 129.770 , 129.865 ,  129.892]
ch1_full_y=[ 126.908    , 126.852    , 126.940    , 126.859    , 126.926    , 126.852    , 126.966    , 126.914 , 126.860 ,  126.889]
ch1_full_XMODE=       129.916
ch1_full_YMODE=      126.880
ch1_full_dayarr = ch1_full_sclk
for i = 0, n_elements(ch1_full_sclk) - 1 do begin
   year = (ch1_full_sclk[i] / 60./60./24./365.) + 1980.
   day = year - fix(year)
   day = day *365.
   ch1_full_dayarr[i]  = day
endfor
;oplot, ch1_full_dayarr, ch1_full_x- ch1_full_xmode, psym = 1
oplot, ch1_full_dayarr, ch1_full_y- ch1_full_ymode, psym = 1

;CH1 FULL ARRAY OBS to SUBARRAY FOV
 
ch1_fullsub_sclk=[  999380117.D,  999630508.D,  999685707.D,  999755111.D,  999759574.D,  999770131.D,  1000257191.D, 1000205883.D,  1000587195.D,  1000590410.D]
ch1_fullsub_x=[ 23.094 , 23.114 , 23.084 , 23.029 , 23.135 , 23.153 ,  23.113, 23.118 ,  22.937,  22.961]
ch1_fullsub_y=[ 230.958 , 231.035 , 231.058 , 231.036 , 231.086 , 231.058 ,  231.043, 231.026 ,  230.949,  230.986]
ch1_fullsub_XMODE=       23.1070
ch1_fullsub_YMODE=       231.048
    
ch1_fullsub_dayarr = ch1_fullsub_sclk
for i = 0, n_elements(ch1_fullsub_sclk) - 1 do begin
   year = (ch1_fullsub_sclk[i] / 60./60./24./365.) + 1980.
   day = year - fix(year)
   day = day *365.
   ch1_fullsub_dayarr[i]  = day
endfor
;oplot, ch1_fullsub_dayarr, ch1_fullsub_x- ch1_fullsub_xmode, psym = 1
oplot, ch1_fullsub_dayarr, ch1_fullsub_y- ch1_fullsub_ymode, psym = 1

yflat = fltarr(300) + 0.
oplot,findgen(300), yflat
yflat = fltarr(300) + 0.02
oplot,findgen(300), yflat, linestyle = 2
yflat = fltarr(300) -0.02
oplot,findgen(300), yflat, linestyle = 2

;-------------------------------------
;here are the gyro updates
t = ['245:17:24:05.500','246:17:48:05.300','247:18:37:05.200','249:18:23:05.000','251:02:39:04.800','253:01:21:04.500','254:00:59:04.300','255:00:49:04.200','257:00:34:03.900','258:01:59:03.700','259:01:00:03.500','259:23:50:03.300','261:09:01:03.100','263:01:55:02.700','264:01:19:02.500','266:01:24:02.100','267:00:00:01.900','268:03:44:01.700','269:04:24:01.400','270:06:29:01.200','271:00:09:01.000']

yflat = findgen(10) -5.
dayarr = fltarr(n_elements(t))
for j = 0, n_elements(t) - 1 do begin
   time = strsplit(t(j), ':',/extract)
   day = time(0) + (time(1) + time(2)/60.) /24.
   dayarr[j] = day
   oplot, fltarr(10) + day, yflat
endfor


;----------------------------------------------------------------------------------------------------
;CH2 stats from Jim
;----------------------------------------------------------------------------------------
;CH2 Subarray (INC_POINT_A)

ch2_sub_novw_sclk=[998852063.D,998905935.D,999018147.D,999085967.D,999134696.D,999139922.D,999141794.D,999184171.D,999185144.D,999217431.D,999265374.D,999273037.D,999382085.D,999528273.D]
ch2_sub_novw_x=[14.720,14.671,14.773,14.597,14.512,14.440,14.568,14.540,14.506,14.761,14.813,14.721,14.764,14.735]
ch2_sub_novw_y=[15.083  , 15.130 , 14.985 , 15.176 , 15.112 , 15.171 , 15.111 , 15.096 , 15.172 , 15.029 , 14.947 , 15.084 , 15.086 , 15.041 ]
ch2_sub_novw_XMODE = 14.7556 
ch2_sub_novw_YMODE=15.0946 

ch2_sub_novw_dayarr = ch2_sub_novw_sclk
for i = 0, n_elements(ch2_sub_novw_sclk) - 1 do begin
   year = (ch2_sub_novw_sclk[i] / 60./60./24./365.) + 1980.
   day = year - fix(year)
   day = day *365.
   ch2_sub_novw_dayarr[i]  = day
endfor

;plot, ch2_sub_novw_dayarr, ch2_sub_novw_x- ch2_sub_novw_xmode, psym = 1, xtitle = 'Time (day of year)', xrange = [250, 280], yrange = [-0.06, 0.06], ytitle = 'X distance from mode', title = 'Ch2'

plot, ch2_sub_novw_dayarr, ch2_sub_novw_y- ch2_sub_novw_ymode, psym = 1, xtitle = 'Time (day of year)', xrange = [245, 280], yrange = [-0.06, 0.06], ytitle = 'Y distance from mode', title = 'Ch2'

ch2_sub_vw_sclk=[998852449.D,998905550.D,999018534.D,999085582.D,999135082.D,999139537.D,999142180.D,999183786.D,999185531.D,999217045.D,999265760.D,999272652.D,999382471.D,999527888.D]
ch2_sub_vw_x=[14.991,14.986,15.114,14.927,14.969,14.918,14.905,14.922,14.908,15.090,15.119,15.083,15.029,15.034]
ch2_sub_vw_y=[15.182,15.162,15.062,15.131,15.150,15.101,15.118,15.117,15.205,15.093,15.019,15.094,15.072,15.030]
ch2_sub_vw_XMODE =    15.1200 
ch2_sub_vw_YMODE=15.0850 

ch2_sub_vw_dayarr = ch2_sub_vw_sclk
for i = 0, n_elements(ch2_sub_vw_sclk) - 1 do begin
   year = (ch2_sub_vw_sclk[i] / 60./60./24./365.) + 1980.
   day = year - fix(year)
   day = day *365.
   ch2_sub_vw_dayarr[i]  = day
endfor

;oplot, ch2_sub_vw_dayarr, ch2_sub_vw_x- ch2_sub_vw_xmode, psym = 1
oplot, ch2_sub_vw_dayarr, ch2_sub_vw_y- ch2_sub_vw_ymode, psym = 1


;CH2 SUBARRAY (INC_POINT_C)
ch2_sub_incc_sclk=[999379515.D,999391143.D,999395786.D,999400218.D,999405722.D,999444386.D,999452723.D,999529793.D,999540774.D,999542758.D]
ch2_sub_incc_x=[14.806,14.754,14.811,14.906,14.742,14.642,14.801,14.693,14.808,14.774]
ch2_sub_incc_y=[15.118,15.051,14.979,14.916,15.013,15.160,14.987,15.018,14.977,15.063]
ch2_sub_incc_XMODE=       14.7955
ch2_sub_incc_YMODE=       15.0196

ch2_sub_incc_dayarr = ch2_sub_incc_sclk
for i = 0, n_elements(ch2_sub_incc_sclk) - 1 do begin
   year = (ch2_sub_incc_sclk[i] / 60./60./24./365.) + 1980.
   day = year - fix(year)
   day = day *365.
   ch2_sub_incc_dayarr[i]  = day
endfor

;oplot, ch2_sub_incc_dayarr, ch2_sub_incc_x- ch2_sub_incc_xmode, psym = 1
oplot, ch2_sub_incc_dayarr, ch2_sub_incc_y- ch2_sub_incc_ymode, psym = 1

;CH2 FULL ARRAY

ch2_full_sclk=[  999444987.D,  999541376.D,  999552600.D,  999651532.D,  999688191.D,  999696224.D, 1000204956.D, 1000258124.D,  1000552456.D,  1000589483.D]
ch2_full_x=[ 126.939 , 127.055 , 127.007 , 127.131 , 127.040 , 126.993 , 127.094 , 127.026 ,  126.911,  126.841]
ch2_full_y=[ 127.933, 127.899, 127.949, 127.756, 127.750, 127.905, 127.760, 127.770,  127.88,  127.69]
ch2_full_XMODE=       127.027
ch2_full_YMODE=       127.758

ch2_full_dayarr = ch2_full_sclk
for i = 0, n_elements(ch2_full_sclk) - 1 do begin
   year = (ch2_full_sclk[i] / 60./60./24./365.) + 1980.
   day = year - fix(year)
   day = day *365.
   ch2_full_dayarr[i]  = day
endfor

;oplot, ch2_full_dayarr, ch2_full_x- ch2_full_xmode, psym = 1
oplot, ch2_full_dayarr, ch2_full_y- ch2_full_ymode, psym = 1

;CH2 FULL ARRAY, WITH OFFSET TO SUBARRAY FOV

ch2_fullsub_sclk=[   999543360.D,   999545886.D,   999582183.D,   999802515.D,   999905598.D,  1000159590.D,  1000207404.D,   1000439372.D,   1000489782.D,   1000588128.D]
ch2_fullsub_x=[ 23.788 , 23.880 , 23.750 , 23.828 , 23.760 , 23.851 , 23.770 ,  23.730,  23.790,  23.584]
ch2_fullsub_y=[ 231.368 , 231.402 , 231.440 , 231.474 , 231.472 , 231.349 , 231.325 ,  231.434,  231.513,  231.341]
ch2_fullsub_XMODE=       23.7984
ch2_fullsub_YMODE=       231.400

ch2_fullsub_dayarr = ch2_fullsub_sclk
for i = 0, n_elements(ch2_fullsub_sclk) - 1 do begin
   year = (ch2_fullsub_sclk[i] / 60./60./24./365.) + 1980.
   day = year - fix(year)
   day = day *365.
   ch2_fullsub_dayarr[i]  = day
endfor

;oplot, ch2_fullsub_dayarr, ch2_fullsub_x- ch2_fullsub_xmode, psym = 1
oplot, ch2_fullsub_dayarr, ch2_fullsub_y- ch2_fullsub_ymode, psym = 1

;ch2_sub_sweet spot
 
ch2_sub_sweet_sclk=[1001257540.D,1001270945.D,1001276019.D,1001329342.D,1001367863.D,1001373884.D,1001374819.D,1001375824.D,1001379917.D,1001576151.D]
ch2_sub_sweet_x=[15.090,15.127,15.104,15.078,15.105,15.115,15.097,15.085,15.079,15.096]
ch2_sub_sweet_y=[15.138,15.113,15.171,15.068,15.033,14.966,15.011,15.057,15.043,15.080]
ch2_sub_sweet_XMODE=       15.0983
ch2_sub_sweet_YMODE=       15.0542
   print, 'sweet', ch2_sub_sweet_y - ch2_sub_sweet_ymode


ch2_sub_sweet_dayarr = ch2_sub_sweet_sclk
for i = 0, n_elements( ch2_sub_sweet_sclk) - 1 do begin
   year = ( ch2_sub_sweet_sclk[i] / 60./60./24./365.) + 1980.
   day = year - fix(year)
   day = day *365.
   ch2_sub_sweet_dayarr[i]  = day
endfor

;oplot, ch2_sub_sweet_dayarr, ch2_sub_sweet_x- ch2_sub_sweet_xmode, psym = 1
oplot, ch2_sub_sweet_dayarr, ch2_sub_sweet_y- ch2_sub_sweet_ymode, psym = 1


yflat = fltarr(300) + 0.
oplot,findgen(300), yflat
yflat = fltarr(300) + 0.02
oplot,findgen(300), yflat, linestyle = 2
yflat = fltarr(300) -0.02
oplot,findgen(300), yflat, linestyle = 2

;-------------------------------------
;here are the gyro updates
t = ['245:17:24:05.500','246:17:48:05.300','247:18:37:05.200','249:18:23:05.000','251:02:39:04.800','253:01:21:04.500','254:00:59:04.300','255:00:49:04.200','257:00:34:03.900','258:01:59:03.700','259:01:00:03.500','259:23:50:03.300','261:09:01:03.100','263:01:55:02.700','264:01:19:02.500','266:01:24:02.100','267:00:00:01.900','268:03:44:01.700','269:04:24:01.400','270:06:29:01.200','271:00:09:01.000']

yflat = findgen(10) -5.
dayarr = fltarr(n_elements(t))
for j = 0, n_elements(t) - 1 do begin
   time = strsplit(t(j), ':',/extract)
   day = time(0) + (time(1) + time(2)/60.) /24.
   dayarr[j] = day
   oplot, fltarr(10) + day, yflat
endfor

ps_close, /noprint,/noid

end



