pro iwic850
;36 hr stability test
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
;now for the mean level per amplifier from the skydarks
ps_open, filename='/Users/jkrick/iwic/iwic850_2.ps',/portrait,/square,/color
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
orangecolor = FSC_COLOR("Orange", !D.Table_Size-6)


;read in all the darks
dir_name = '/Users/jkrick/iwic/iwicts390/IRAC025600'
aor_list = 'skydark_all.txt'

readcol, dir_name + '/'+ aor_list, aorname, format="A",/silent

for ch = 0,  1 do begin
   for n = 0, 0 do begin        ;the dark and the uncertainty
      c4=0
      c12=0
      c100 =0
      
      med_4 = fltarr(4, 400)
      med_12 = fltarr(4, 400)
      med_100 = fltarr(4, 400)
      t_4 = fltarr(400)
      t_12 = fltarr(400)
      t_100 = fltarr(400)
      zodiarr = fltarr(400)

      for pos = 0,    n_elements(aorname) - 1  do begin ; 
         cd, dir_name+'/cal/' +aorname(pos)
;         print, 'working on ', aorname(pos)
         if ch eq 0 then command1 =  ' find . -name "IRAC.1*.fits" >  list.txt'
         if ch eq 1 then command1 =  ' find . -name "IRAC.2*.fits" >  list.txt'
         spawn, command1
         readcol, dir_name + '/cal/' + aorname(pos) + '/list.txt', fitsname , format="A",/silent ;will be both 0.4 and 12s and 100s
         
         for f = 0, n_elements(fitsname) - 1 do begin
            if strmid(fitsname(f), 13, 1,/reverse) eq '4' or strmid(fitsname(f), 13, 1,/reverse) eq '3' then begin
               fits_read, fitsname(f), data, header
               
                                ;which exptime is this
               exp = fxpar(header, 'EXPTIME')
               time  = fxpar (header, 'SCLK_CAL')
               zodi = fxpar(header, 'ZODY_EST')
;               print, 'zodi', zodi
               if ch eq 0 then factor = 17100
               if ch eq 1 then factor = 1130

                                ;divide into 4 readouts
               cnvrt_stdim_to_64_256_4, data[*,*,n], data_4
               
               for read = 0, 3 do begin
                  mmm, data_4[*,*,read], skymode, sigma, skew,/silent
                 ; median_read = median(data_4[*,*,read])
                  median_read = skymode
 ;                 print, 'read, median', read, median_read

                                ;need some way of keeping track of these
                  if exp eq 0.4 then begin
                     med_4(read,c4) = median_read
                     t_4(c4) = time
                     c4 = c4 + 1
                  endif
                  if exp eq 10.4 then begin
                     med_12(read,c12) = median_read
                     t_12(c12) = time
                     c12 = c12 + 1
                  endif
                  if exp gt 92 then begin
                     med_100(read,c100) = median_read
                     t_100(c100) = time
                     zodiarr(c100) = zodi *factor
                     c100 = c100 + 1
                     
                  endif
                  
               endfor
            endif
            
         endfor
      endfor
        
      t_4 = t_4[0:c4-1]
      t_12 = t_12[0:c12-1]
      t_100 = t_100[0:c100-1]
      med_4 = med_4[*,0:c4-1]
      med_12 = med_12[*,0:c12-1]
      med_100= med_100[*,0:c100-1]
      zodiarr = zodiarr[0:c100-1]

      !P.multi =[0,1,3]
      !P.charsize = 1.5
      if n eq 0 and ch eq 0 then yr = [140,155]
      if n eq 0 and ch eq 1 then yr = [6,10]
      if n eq 1 and ch eq 0 then yr = [1.5,2.0]
      if n eq 1 and ch eq 1 then yr = [1.0,1.5]
      if n eq 0 then yname = 'mode of background regions (DN)'
      if n eq 1 then yname = 'mode of uncertainty'
      plot, (t_4 - min(t_4))/3600., med_4(0,*),  psym = 1, xtitle = 'Time (hrs)', $
            ytitle =yname,$
            title = "channel "+string( ch+1)+ "    exptime 0.4", thick=3, $
            xthick=3,ythick=3,charthick=3,$
            yrange = yr
      oplot,  (t_4 - min(t_4))/3600., med_4(1,*),  psym = 2, thick = 3
      oplot,  (t_4 - min(t_4))/3600., med_4(2,*),  psym = 4, thick = 3
      oplot,  (t_4 - min(t_4))/3600., med_4(3,*),  psym = 5, thick = 3
      
      if n eq 0 and ch eq 0 then yr = [430,460]
      if n eq 0 and ch eq 1 then yr = [46,52]
      if n eq 1 and ch eq 0 then yr = [1.4,1.9]
      if n eq 1 and ch eq 1 then yr = [.75,1.25]
      plot, (t_12 - min(t_12))/3600., med_12(0,*),  psym = 1, xtitle = 'Time (hrs)', $
            ytitle = yname,$
            title = "channel "+string( ch+1)+ "    exptime 12", thick=3, $
            xthick=3,ythick=3,charthick=3,$
            yrange = yr
      oplot,  (t_12 - min(t_12))/3600., med_12(1,*),  psym = 2, thick = 3
      oplot,  (t_12 - min(t_12))/3600., med_12(2,*),  psym = 4, thick = 3
      oplot,  (t_12 - min(t_12))/3600., med_12(3,*),  psym = 5, thick = 3
      
      if n eq 0 and ch eq 0 then yr = [700,740]
      if n eq 0 and ch eq 1 then yr = [235,250]
      if n eq 1 and ch eq 0 then yr = [2.0,2.5]
      if n eq 1 and ch eq 1 then yr = [1.75,2.25]
 plot, (t_100 - min(t_100))/3600., med_100(0,*),  psym = 1, xtitle = 'Time (hrs)', $
            ytitle =yname,$
            title = "channel "+string( ch+1)+ "    exptime 100", thick=3, $
            xthick=3,ythick=3,charthick=3,$
            yrange = yr
      oplot,  (t_100 - min(t_100))/3600., med_100(1,*),  psym = 2, thick = 3
      oplot,  (t_100 - min(t_100))/3600., med_100(2,*),  psym = 4, thick = 3
      oplot,  (t_100 - min(t_100))/3600., med_100(3,*),  psym = 5, thick = 3
      oplot, (t_100 - min(t_100))/3600., zodiarr,   thick = 3
   endfor
endfor

;ps_close, /noprint,/noid
   
   
   
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------

!P.multi =[0,1,1]

channel = [1, 2]

gain = [3.3,3.7]
gain = [3.82, 3.82]
flux_conv = [.1277,.1461] ;from Mark's most recent post on MODA

fits_read, '/Users/jkrick/iwic/ch1_photcorr_rj.fits', ch1_arraycorr, ch1head
fits_read, '/Users/jkrick/iwic/ch2_photcorr_rj.fits', ch2_arraycorr, ch2head
fits_read, '/Users/jkrick/iwic/test_large.fits', ch2_arraycorr, ch2head

;start with NPM1p67.0536
dir_name = '/Users/jkrick/iwic/iwicts390/IRAC025600'
aor_list1 = 'NPMcalstar.txt'
aor_list2 = 'KFcalstar.txt'
cal_list = ['NPMcalstar.txt','KFcalstar.txt']
;filenames, dir_name, aor_list


;ps_open, filename='/Users/jkrick/iwic/phot_stabil_NPM.ps',/portrait,/square,/color
;ps_open, filename='/Users/jkrick/iwic/phot_stabil_KF.ps',/portrait,/square,/color

for cal = 0, n_elements(cal_list) - 1 do begin
   aor_list = cal_list[cal]
   readcol, dir_name + '/'+ aor_list, aorname, format="A",/silent
   
   if cal eq 0 then ra = 269.72731
   if cal eq 0 then dec = 67.79371
   
   if cal eq 1 then ra = 269.49386
   if cal eq 1 then dec = 66.87483


   for ch = 0,  1 do begin
      ;restore Jim's pixel phase files
      if ch eq 0 then restore, '/Users/jkrick/iwic/pixel_phase_img_ch1.sav'   
      if ch eq 1 then restore, '/Users/jkrick/iwic/pixel_phase_img_ch2.sav'   

      aperarr = fltarr(900)
      eaperarr = fltarr(900)
      timearr  = lonarr(900)
      c  = 0
      
      for pos = 0,   n_elements(aorname) - 1  do begin ; 
;         print, 'working on position ',pos, aorname(pos)
         
         
         cd, dir_name+'/bcd/' +aorname(pos)
         if cal eq 0 then readcol,strcompress('bcdlist_ch'+string(channel(ch)) + '_0pt4s.txt',/remove_all), fitsname, format = 'A',/silent
         if cal eq 1 then  readcol,strcompress('bcdlist_ch'+string(channel(ch)) + '_12s.txt',/remove_all), fitsname, format = 'A',/silent     
         for f = 0, n_elements(fitsname) - 1 do begin
            if strmid(fitsname(f), 12, 1,/reverse) eq '4' or strmid(fitsname(f), 12, 1,/reverse) eq '3' then begin
;               print, 'working on fitsname', fitsname(f)
               fits_read, fitsname(f), data, header
                                ;do aperture photometry at the location of the calstar
               
                                ;put the data in units of electrons so
                                ;that noise calculation is correct
               data = data * gain(ch)
               data = data * fxpar(header, 'EXPTIME')
               data = data / flux_conv(ch)
               
            ;add the zodi back in for error purposes
               data = data + fxpar(header, 'ZODY_EST')

            ;apply array-location-dependent correction from cryo mission
               if ch eq 0 then data = data * ch1_arraycorr
               if ch eq 1 then data = data * ch2_arraycorr
               
               
               adxy, header, ra, dec, x, y
               cntrd, data, x, y, xcen, ycen, 2.0
               aper, data, xcen,ycen, flux, eflux, sky, skyerr, gain(ch), 10, [12,20], [-32767,80000],/exact, /flux,/silent
               

                                ;don't keep the fluxes that are NaN's
                                ;because the star was being observed
                                ;in the other channel
               if finite(flux) gt 0 then begin

                                ;apply pixel phase correction from cryo mission
                                ;only ch1
                  ;p = sqrt((xcen - fix(xcen))^2 + (ycen - fix(ycen))^2)
                  ;fphase = 1 + 0.0535*(1/sqrt(2*!Pi) - p)
                  ;if ch eq 0 then flux = flux / fphase
               
                                ;make a correction for pixel phase 
                  xphase = (xcen mod 1) - NINT(xcen mod 1)
                  yphase = (ycen mod 1) - NINT(ycen mod 1)
                  interpolated_relative_flux = interp2d(relative_flux,x_center,y_center,xphase,yphase,/regular)
                  flux = flux / interpolated_relative_flux

               
                                ;also keep track of the time of the observation
                  time  = fxpar (header, 'SCLK_OBS')
                  
                  aperarr(c) = flux
                  eaperarr(c) = eflux
                  timearr(c) = time
                  c = c + 1
               endif
            endif
            
         endfor
      endfor
      
      timearr = timearr[0:c-1]
      aperarr = aperarr[0:c-1]
      eaperarr = eaperarr[0:c-1]
      
      if cal eq 0 and ch eq 0 then yr = [1E5,1.4E5]
      if cal eq 0 and ch eq 1 then yr = [6E4, 7.5E4]
      if cal eq 1 and ch eq 0 then yr = [8E4,14E4]
      if cal eq 1 and ch eq 1 then yr = [5E4, 7.5E4]
      plot, (timearr - min(timearr) )/ 3600., aperarr,  yrange=yr, psym = 2, xtitle = 'time (hrs)', ytitle = 'flux (e-?)', charsize = 1.5$
      , title = strcompress('ch ' + string(ch + 1)) + aor_list; ' KF06T1' ; ' NPM1p67.0536'
      oploterror, (timearr - min(timearr)) / 3600., aperarr, eaperarr, nskip = 5, psym = 2
      
      ;want first one of every observation suite
      ;obs 0, 5, 10, ....
      x = (timearr - min(timearr) )/ 3600.
      c1 = 0
      c2 = 0
      c3 = 0
      c4 = 0
      c5=0
      aper1 = fltarr(n_elements(aperarr))
      t1 = fltarr(n_elements(aperarr))
      aper2 = fltarr(n_elements(aperarr))
      t2 = fltarr(n_elements(aperarr))
      aper3 = fltarr(n_elements(aperarr))
      t3 = fltarr(n_elements(aperarr))
      aper4 = fltarr(n_elements(aperarr))
      t4 = fltarr(n_elements(aperarr))
     aper5 = fltarr(n_elements(aperarr))
      t5 = fltarr(n_elements(aperarr))

       for cn = 0, n_elements(aperarr) - 1, 5 do begin
         aper1[c1] = aperarr[cn]
         t1[c1] = x[cn]
         c1 = c1 + 1
      endfor
      t1 = t1[0:c1-1]
      aper1 = aper1[0:c1-1]
      
      for cn = 1, n_elements(aperarr) - 1, 5 do begin
         aper2[c2] = aperarr[cn]
         t2[c2] = x[cn]
         c2 = c2 + 1
      endfor
      t2 = t1[0:c2-1]
      aper2 = aper2[0:c2-1]

      for cn = 2, n_elements(aperarr) - 1, 5 do begin
         aper3[c3] = aperarr[cn]
         t3[c3] = x[cn]
         c3 = c3 + 1
      endfor
      t3 = t3[0:c3-1]
      aper3 = aper3[0:c3-1]

      for cn = 3, n_elements(aperarr) - 1, 5 do begin
         aper4[c4] = aperarr[cn]
         t4[c4] = x[cn]
         c4 = c4 + 1
      endfor
      t4 = t4[0:c4-1]
      aper4 = aper4[0:c4-1]

      for cn = 4, n_elements(aperarr) - 1, 5 do begin
         aper5[c5] = aperarr[cn]
         t5[c5] = x[cn]
         c5 = c5 + 1
      endfor
      t5 = t5[0:c5-1]
      aper5 = aper5[0:c5-1]
      ;test this out
 
      oplot, t1, aper1, psym = 2, color = redcolor
      oplot, t2, aper2, psym = 2, color = orangecolor
      oplot, t3, aper3, psym = 2, color =yellowcolor
      oplot, t4, aper4, psym = 2, color = greencolor
      oplot, t5, aper5, psym = 2, color = bluecolor



                                ; robust linear regression
      x = (timearr - min(timearr) )/ 3600.
      y =  aperarr
      
      a = where(finite(y) gt 0)
      x = x(a)
      y = y(a)
      
      save, filename = '/Users/jkrick/iwic/iwicts390/IRAC025600/variables.sav', x, y
      coeff = robust_linefit(x, y, yfit, sig, coef_sig)
      oplot, x, coeff(0) + coeff(1)*x, color = redcolor
      flaty = fltarr(n_elements(y)) + coeff(0)
      oplot, x, flaty , color = bluecolor
      print, 'coeff', coeff
   endfor
endfor


ps_close, /noprint,/noid





end
