pro iwic430_latent
  ;ps_open, filename='/Users/jkrick/iwic/longer_latent.ps',/color
  ps_start,filename='/Users/jkrick/iwic/longer_latent_servs.ps'
  !P.multi =[0,0,1]
  !P.charthick = 3
  !P.thick = 3
  !X.thick = 3
  !Y.thick = 3

  redcolor = FSC_COLOR("Red", !D.Table_Size-2)
  bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
  greencolor = FSC_COLOR("Green", !D.Table_Size-4)
  yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
  cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
  orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
  purplecolor = FSC_COLOR("purple", !D.Table_Size-8)
  blackcolor = FSC_COLOR("black", !D.Table_Size-9)

;look at latent test 3.5 and 4.5 mag stars
  dir_name =[ '/Users/jkrick/iwic/iwicDAW/IRAC022000/bcd/0037381632','/Users/jkrick/iwic/iwicDAW/IRAC022000/bcd/0037379840', '/Users/jkrick/iwic/iwicDAW/IRAC022000/bcd/0037381888','/Users/jkrick/iwic/iwicDAW/IRAC022000/bcd/0037380864','/Users/jkrick/iwic/iwicDAW/IRAC022000/bcd/0037382144','/Users/jkrick/iwic/iwicDAY/IRAC025000/bcd/0037399552']
  bcdch1list = ' '
  bcdch2list = ' '

  for j = 0, n_elements(dir_name) - 1 do begin
     CD,dir_name(j)
     
     command1 =  ' find '+ dir_name(j) +' -name "*bcd_fp.fits" >  bcd_list.txt'
     command2 = 'cat bcd_list.txt | grep IRAC.1. > bcd_ch1.txt'
     command3 = 'cat bcd_list.txt | grep IRAC.2. | grep -v IRAC.1 > bcd_ch2.txt'
     
     commands = [command1, command2, command3]
     for i = 0, n_elements(commands) -1 do spawn, commands(i)
     
     readcol,dir_name(j) + '/bcd_ch1.txt', bcdch1list_a, format="A", /silent
     readcol,dir_name(j) + '/bcd_ch2.txt', bcdch2list_a, format="A", /silent

     bcdch1list = [bcdch1list , bcdch1list_a]
     bcdch2list = [bcdch2list , bcdch2list_a]

  endfor
  
  ;print, bcdch1list

  xcen = [165,166]
  ycen = [228,203]
  flux1arr_0p6 =fltarr(n_elements(bcdch1list))
  flux2arr_0p6 =fltarr(n_elements(bcdch1list))
  time_0p6 = fltarr(n_elements(bcdch1list))
  count_0p6 = 0
  flux1arr_12 =fltarr(n_elements(bcdch1list))
  flux2arr_12 =fltarr(n_elements(bcdch1list))
  time_12 = fltarr(n_elements(bcdch1list))
  count_12 = 0
  flux1arr_100 =fltarr(n_elements(bcdch1list))
  flux2arr_100 =fltarr(n_elements(bcdch1list))
  time_100 = fltarr(n_elements(bcdch1list))
  count_100 = 0
  
  for i = 30, n_elements(bcdch1list) - 1  do begin
     fits_read, bcdch1list(i), data, imheader
     
                                ;do aperture photometry on star position 
     object = data[xcen(0)-2:xcen(0)+2,ycen(0)-2:ycen(0)+2]
     nobject = n_elements(object)
     sky = total(data[xcen(0)-23:xcen(0)+23, ycen(0)-23:ycen(0)+23])- total(data[xcen(0)-15:xcen(0)+15,ycen(0)-15:ycen(0)+15])  
     nannulus = n_elements(data[xcen(0)-23:xcen(0)+23, ycen(0)-23:ycen(0)+23]) - n_elements(data[xcen(0)-15:xcen(0)+15,ycen(0)-15:ycen(0)+15])    
     sky = sky / nannulus       ;mean sky flux
     flux1 =  total( object)  - sky*nobject
     
     object = data[xcen(1)-2:xcen(1)+2,ycen(1)-2:ycen(1)+2]
     nobject = n_elements(object)
     sky = total(data[xcen(1)-23:xcen(1)+23, ycen(1)-23:ycen(1)+23])- total(data[xcen(1)-15:xcen(1)+15,ycen(1)-15:ycen(1)+15])  
     nannulus = n_elements(data[xcen(1)-23:xcen(1)+23, ycen(1)-23:ycen(1)+23]) - n_elements(data[xcen(1)-15:xcen(1)+15,ycen(1)-15:ycen(1)+15])    
     sky = sky / nannulus       ;mean sky flux
     flux2 =  total( object)  - sky*nobject
     
                                ;keep track of time of observation
     time= fxpar(imheader, 'SCLK_OBS')
     framtime = fxpar(imheader, 'FRAMTIME')
     
     if framtime eq 0.6 then  begin
        flux1arr_0p6(count_0p6) = flux1
        flux2arr_0p6(count_0p6) = flux2
        time_0p6(count_0p6) = time
        count_0p6 = count_0p6+ 1
     endif
     if framtime eq 12 then  begin
        flux1arr_12(count_12) = flux1
        flux2arr_12(count_12) = flux2
        time_12(count_12) = time
        count_12 = count_12+ 1
        
     endif
     if framtime eq 100 then  begin
        flux1arr_100(count_100) = flux1
        flux2arr_100(count_100) = flux2
        time_100(count_100) = time
        count_100 = count_100+ 1
     endif
     
     
  endfor
  
  flux1arr_0p6 = flux1arr_0p6[0:count_0p6 - 1]
  flux2arr_0p6 = flux2arr_0p6[0:count_0p6 - 1]
  time_0p6 = time_0p6[0:count_0p6 - 1]
  flux1arr_12 = flux1arr_12[0:count_12 - 1]
  flux2arr_12 = flux2arr_12[0:count_12 - 1]
  time_12 = time_12[0:count_12 - 1]
  flux1arr_100 = flux1arr_100[0:count_100 - 1]
  flux2arr_100 = flux2arr_100[0:count_100 - 1]
  time_100 = time_100[0:count_100 - 1]

plot, (time_100 - time_100(0))/3600., flux1arr_100/ flux1arr_100(0), psym = 2, xrange=[0,3], /ylog, yrange=[1E-4, 1],$
        xtitle = 'hours', ytitle = 'fractional flux', title = 'Latents vs. time',/nodata
  oplot, (time_12 - time_12(0))/3600., flux1arr_12/ flux1arr_12(0), psym = 2, color = bluecolor
   
;  oplot, (time_100 - time_100(0))/3600., flux1arr_100/ flux1arr_100(0)
  oplot, (time_12 - time_12(0))/3600., flux1arr_12/ flux1arr_12(0),color = bluecolor
   
;  oplot, (time_100 - time_100(0))/3600., flux2arr_100/ flux2arr_100(0), psym = 6
  oplot, (time_12 - time_12(0))/3600., flux2arr_12/ flux2arr_12(0), psym = 6, color = bluecolor
   
;  oplot, (time_100 - time_100(0))/3600., flux2arr_100/ flux2arr_100(0)
  oplot, (time_12 - time_12(0))/3600., flux2arr_12/ flux2arr_12(0),  color = bluecolor

;-------------------------------------------------------------
;look at star in Mark's SERVS data

  dir_name =['/Users/jkrick/iwic/iwicts390/IRAC025600/bcd/0037713920'];, '/Users/jkrick/iwic/iwicts390/IRAC025600/bcd/0037714432'];,'/Users/jkrick/iwic/iwicts390/IRAC025600/bcd/0037714176']
  bcdch1list = ' '
  bcdch2list = ' '

  for j = 0, n_elements(dir_name) - 1 do begin
     CD,dir_name(j)
     
     command1 =  ' find '+ dir_name(j) +' -name "*3.bcd_fp.fits" >  bcd_list.txt'
     command4 =  ' find '+ dir_name(j) +' -name "*4.bcd_fp.fits" >>  bcd_list.txt'

     command2 = 'cat bcd_list.txt | grep IRAC.1. > bcd_ch1.txt'
     command3 = 'cat bcd_list.txt | grep IRAC.2. | grep -v IRAC.1 > bcd_ch2.txt'
     
     commands = [command1, command4, command2, command3]
     for i = 0, n_elements(commands) -1 do spawn, commands(i)
     
     readcol,dir_name(j) + '/bcd_ch1.txt', bcdch1list_a, format="A", /silent
;     readcol,dir_name(j) + '/bcd_ch2.txt', bcdch2list_a, format="A", /silent

     bcdch1list = [bcdch1list , bcdch1list_a]
;     bcdch2list = [bcdch2list , bcdch2list_a]

  endfor
  
  print, 'servs',  bcdch1list

  xcen = [38,16]
  ycen = [111,116]
  flux1arr_0p6 =fltarr(n_elements(bcdch1list))
  flux2arr_0p6 =fltarr(n_elements(bcdch1list))
  time_0p6 = fltarr(n_elements(bcdch1list))
  count_0p6 = 0
  flux1arr_12 =fltarr(n_elements(bcdch1list))
  flux2arr_12 =fltarr(n_elements(bcdch1list))
  time_12 = fltarr(n_elements(bcdch1list))
  count_12 = 0
  flux1arr_100 =fltarr(n_elements(bcdch1list))
  flux2arr_100 =fltarr(n_elements(bcdch1list))
  time_100 = fltarr(n_elements(bcdch1list))
  count_100 = 0
  
  for i = 43, n_elements(bcdch1list) - 1  do begin
     fits_read, bcdch1list(i), data, imheader
     a = where(finite(data) lt 1, count)
     if count gt 0 then data[a] = 0
                                ;do aperture photometry on star position 
     object = data[xcen(0)-2:xcen(0)+2,ycen(0)-2:ycen(0)+2]
     nobject = n_elements(object)
     sky = total(data[xcen(0)-23:xcen(0)+23, ycen(0)-23:ycen(0)+23])- total(data[xcen(0)-15:xcen(0)+15,ycen(0)-15:ycen(0)+15])  
     nannulus = n_elements(data[xcen(0)-23:xcen(0)+23, ycen(0)-23:ycen(0)+23]) - n_elements(data[xcen(0)-15:xcen(0)+15,ycen(0)-15:ycen(0)+15])    
     sky = sky / nannulus       ;mean sky flux
     flux1 =  total( object)  - sky*nobject
     
     object = data[xcen(1)-2:xcen(1)+2,ycen(1)-2:ycen(1)+2]
     nobject = n_elements(object)
     sky = total(data[xcen(1)-13:xcen(1)+33, ycen(1)-23:ycen(1)+23])- total(data[xcen(1)-13:xcen(1)+15,ycen(1)-15:ycen(1)+15])  
     nannulus = n_elements(data[xcen(1)-13:xcen(1)+33, ycen(1)-23:ycen(1)+23]) - n_elements(data[xcen(1)-13:xcen(1)+15,ycen(1)-15:ycen(1)+15])    
     sky = sky / nannulus       ;mean sky flux
     flux2 =  total( object)  - sky*nobject

                                ;keep track of time of observation
     time= fxpar(imheader, 'SCLK_OBS')
     framtime = fxpar(imheader, 'FRAMTIME')
     
    
;     if framtime eq 200 then  begin
;        flux1arr_12(count_12) = flux1
;        flux2arr_12(count_12) = flux2
;        time_12(count_12) = time
;        count_12 = count_12+ 1
        
;     endif
;     if framtime eq 100 then  begin
        flux1arr_100(count_100) = flux1
        flux2arr_100(count_100) = flux2
        time_100(count_100) = time
        count_100 = count_100+ 1
;     endif
     
     
  endfor
  
;   flux1arr_12 = flux1arr_12[0:count_12 - 1]
;  flux2arr_12 = flux2arr_12[0:count_12 - 1]
;  time_12 = time_12[0:count_12 - 1]
  flux1arr_100 = flux1arr_100[0:count_100 - 1]
  flux2arr_100 = flux2arr_100[0:count_100 - 1]
  time_100 = time_100[0:count_100 - 1]

;print, 'x', (time_100 - time_100(0))/3600.
;print, 'y', flux1arr_100/ flux1arr_100(0)
  oplot, (time_100 - time_100(0))/3600., flux1arr_100/ flux1arr_100(0), psym = 5
;  oplot, (time_12 - time_12(0))/3600., flux1arr_12/ flux1arr_12(0), psym = 1, color = redcolor
  
  oplot, (time_100 - time_100(0))/3600., flux1arr_100/ flux1arr_100(0)
;  oplot, (time_12 - time_12(0))/3600., flux1arr_12/ flux1arr_12(0),color = redcolor
  
;  oplot, (time_100 - time_100(0))/3600., flux2arr_100/ flux2arr_100(0), psym = 5
;  oplot, (time_12 - time_12(0))/3600., flux2arr_12/ flux2arr_12(0), psym = 5, color = redcolor
  
;  oplot, (time_100 - time_100(0))/3600., flux2arr_100/ flux2arr_100(0)
;  oplot, (time_12 - time_12(0))/3600., flux2arr_12/ flux2arr_12(0),  color = redcolor
 
legend, ['5mag, 12sframtime, 48s exposed,  July 2','1mag, 2sframtime, 102s exposed(1 dither), July 2', '3.5mag, 12sframtime,  112.6s exposed, July 15', '4.5mag, 12sframtime, 112.6s exposed, July 15', '4.7mag, 100sframtime, 600s exposed(small dithers) July 22'], psym = [4,1, 2,6,5], color = [redcolor,  greencolor, bluecolor, bluecolor, blackcolor]
;-------------------------------------------------------------
;looking at the 5th mag star as part of the final setpoint test
  dir_name =['/Users/jkrick/iwic/iwicDAA/IRAC020700/bcd/0036840192', '/Users/jkrick/iwic/iwicDAA/IRAC020700/bcd/0036862976' ,'/Users/jkrick/iwic/iwicDAA/IRAC020700/bcd/0036839424']
  bcdch1list = ' '
  bcdch2list = ' '

  for j = 0, n_elements(dir_name) - 1 do begin
     CD,dir_name(j)
     
     command1 =  ' find '+ dir_name(j) +' -name "*bcd_fp.fits" >  bcd_list.txt'
     command2 = 'cat bcd_list.txt | grep IRAC.1. > bcd_ch1.txt'
     command3 = 'cat bcd_list.txt | grep IRAC.2. | grep -v IRAC.1 > bcd_ch2.txt'
     
     commands = [command1, command2, command3]
     for i = 0, n_elements(commands) -1 do spawn, commands(i)
     
     readcol,dir_name(j) + '/bcd_ch1.txt', bcdch1list_a, format="A", /silent
     readcol,dir_name(j) + '/bcd_ch2.txt', bcdch2list_a, format="A", /silent

     bcdch1list = [bcdch1list , bcdch1list_a]
     bcdch2list = [bcdch2list , bcdch2list_a]

  endfor
  
  ;print, bcdch1list

  xcen = [85,190]
  ycen = [105,150]
  flux1arr_0p6 =fltarr(n_elements(bcdch1list))
  flux2arr_0p6 =fltarr(n_elements(bcdch1list))
  time_0p6 = fltarr(n_elements(bcdch1list))
  count_0p6 = 0
  flux1arr_12 =fltarr(n_elements(bcdch1list))
  flux2arr_12 =fltarr(n_elements(bcdch1list))
  time_12 = fltarr(n_elements(bcdch1list))
  count_12 = 0
  flux1arr_100 =fltarr(n_elements(bcdch1list))
  flux2arr_100 =fltarr(n_elements(bcdch1list))
  time_100 = fltarr(n_elements(bcdch1list))
  count_100 = 0
  
  for i = 4, n_elements(bcdch1list) - 1  do begin
     fits_read, bcdch1list(i), data, imheader

     time= fxpar(imheader, 'SCLK_OBS')
     framtime = fxpar(imheader, 'FRAMTIME')

     if framtime eq 12 then  begin
 
        a = where(finite(data) lt 1, count)
        if count gt 0 then data[a] = 0
        
                                ;do aperture photometry on star position 
        object = data[xcen(0)-2:xcen(0)+2,ycen(0)-2:ycen(0)+2]
        nobject = n_elements(object)
        sky = total(data[xcen(0)-23:xcen(0)+23, ycen(0)-23:ycen(0)+23])- total(data[xcen(0)-15:xcen(0)+15,ycen(0)-15:ycen(0)+15])  
        nannulus = n_elements(data[xcen(0)-23:xcen(0)+23, ycen(0)-23:ycen(0)+23]) - n_elements(data[xcen(0)-15:xcen(0)+15,ycen(0)-15:ycen(0)+15])    
        sky = sky / nannulus    ;mean sky flux
        flux1 =  total( object)  - sky*nobject
        
        object = data[xcen(1)-2:xcen(1)+2,ycen(1)-2:ycen(1)+2]
        nobject = n_elements(object)
        sky = total(data[xcen(1)-13:xcen(1)+33, ycen(1)-23:ycen(1)+23])- total(data[xcen(1)-13:xcen(1)+15,ycen(1)-15:ycen(1)+15])  
        nannulus = n_elements(data[xcen(1)-13:xcen(1)+33, ycen(1)-23:ycen(1)+23]) - n_elements(data[xcen(1)-13:xcen(1)+15,ycen(1)-15:ycen(1)+15])    
        sky = sky / nannulus    ;mean sky flux
        flux2 =  total( object)  - sky*nobject
        
                                ;keep track of time of observation
     
    
;     if framtime eq 200 then  begin
;        flux1arr_12(count_12) = flux1
;        flux2arr_12(count_12) = flux2
;        time_12(count_12) = time
;        count_12 = count_12+ 1
        
;     endif
        flux1arr_100(count_100) = flux1
        flux2arr_100(count_100) = flux2
        time_100(count_100) = time
        count_100 = count_100+ 1
     endif
     
     
  endfor
  
;   flux1arr_12 = flux1arr_12[0:count_12 - 1]
;  flux2arr_12 = flux2arr_12[0:count_12 - 1]
;  time_12 = time_12[0:count_12 - 1]
  flux1arr_100 = flux1arr_100[0:count_100 - 1]
  flux2arr_100 = flux2arr_100[0:count_100 - 1]
  time_100 = time_100[0:count_100 - 1]

;print, 'x', (time_100 - time_100(0))/3600.
;print, 'y', flux1arr_100/ flux1arr_100(0)
  oplot, (time_100 - time_100(0))/3600., flux1arr_100/ flux1arr_100(0), psym = 4, color = redcolor
;  oplot, (time_12 - time_12(0))/3600., flux1arr_12/ flux1arr_12(0), psym = 1, color = redcolor
  
  oplot, (time_100 - time_100(0))/3600., flux1arr_100/ flux1arr_100(0), color = redcolor
;  oplot, (time_12 - time_12(0))/3600., flux1arr_12/ flux1arr_12(0),color = redcolor
  
;  oplot, (time_100 - time_100(0))/3600., flux2arr_100/ flux2arr_100(0), psym = 1, color = redcolor
;  oplot, (time_12 - time_12(0))/3600., flux2arr_12/ flux2arr_12(0), psym = 5, color = redcolor
  
;  oplot, (time_100 - time_100(0))/3600., flux2arr_100/ flux2arr_100(0), color = redcolor
;  oplot, (time_12 - time_12(0))/3600., flux2arr_12/ flux2arr_12(0),  color = redcolor
;-------------------------------------------------------------
;looking at the 1st mag star as part of the bright star test in DAA
  dir_name =['/Users/jkrick/iwic/iwicDAA/IRAC020700/bcd/0036894208', '/Users/jkrick/iwic/iwicDAA/IRAC020700/bcd/0036888320' ]
  bcdch1list = ' '
  bcdch2list = ' '

  for j = 0, n_elements(dir_name) - 1 do begin
     CD,dir_name(j)
     
     command1 =  ' find '+ dir_name(j) +' -name "*2.bcd_fp.fits" >  bcd_list.txt'
     command2 = 'cat bcd_list.txt | grep IRAC.1. > bcd_ch1.txt'
     command3 = 'cat bcd_list.txt | grep IRAC.2. | grep -v IRAC.1 > bcd_ch2.txt'
     
     commands = [command1, command2, command3]
     for i = 0, n_elements(commands) -1 do spawn, commands(i)
     
     readcol,dir_name(j) + '/bcd_ch1.txt', bcdch1list_a, format="A", /silent
     readcol,dir_name(j) + '/bcd_ch2.txt', bcdch2list_a, format="A", /silent

     bcdch1list = [bcdch1list , bcdch1list_a]
     bcdch2list = [bcdch2list , bcdch2list_a]

  endfor
  
  print, '1mag',bcdch1list

  xcen = [11,190]
  ycen = [38,150]
  flux1arr_0p6 =fltarr(n_elements(bcdch1list))
  flux2arr_0p6 =fltarr(n_elements(bcdch1list))
  time_0p6 = fltarr(n_elements(bcdch1list))
  count_0p6 = 0
  flux1arr_12 =fltarr(n_elements(bcdch1list))
  flux2arr_12 =fltarr(n_elements(bcdch1list))
  time_12 = fltarr(n_elements(bcdch1list))
  count_12 = 0
  flux1arr_100 =fltarr(n_elements(bcdch1list))
  flux2arr_100 =fltarr(n_elements(bcdch1list))
  time_100 = fltarr(n_elements(bcdch1list))
  count_100 = 0
  
  for i = 4, n_elements(bcdch1list) - 1  do begin
     fits_read, bcdch1list(i), data, imheader

     time= fxpar(imheader, 'SCLK_OBS')
     framtime = fxpar(imheader, 'FRAMTIME')

     if framtime eq 2 then  begin
 
        a = where(finite(data) lt 1, count)
        if count gt 0 then data[a] = 0
        
                                ;do aperture photometry on star position 
        object = data[xcen(0)-2:xcen(0)+2,ycen(0)-2:ycen(0)+2]
        nobject = n_elements(object)
        sky = total(data[xcen(0)-10:xcen(0)+33, ycen(0)-23:ycen(0)+23])- total(data[xcen(0)-10:xcen(0)+15,ycen(0)-15:ycen(0)+15])  
        nannulus = n_elements(data[xcen(0)-10:xcen(0)+33, ycen(0)-23:ycen(0)+23]) - n_elements(data[xcen(0)-10:xcen(0)+15,ycen(0)-15:ycen(0)+15])    
        sky = sky / nannulus    ;mean sky flux
        flux1 =  total( object)  - sky*nobject
        
        object = data[xcen(1)-2:xcen(1)+2,ycen(1)-2:ycen(1)+2]
        nobject = n_elements(object)
        sky = total(data[xcen(1)-13:xcen(1)+33, ycen(1)-23:ycen(1)+23])- total(data[xcen(1)-13:xcen(1)+15,ycen(1)-15:ycen(1)+15])  
        nannulus = n_elements(data[xcen(1)-13:xcen(1)+33, ycen(1)-23:ycen(1)+23]) - n_elements(data[xcen(1)-13:xcen(1)+15,ycen(1)-15:ycen(1)+15])    
        sky = sky / nannulus    ;mean sky flux
        flux2 =  total( object)  - sky*nobject
        
                                ;keep track of time of observation
     
    
;     if framtime eq 200 then  begin
;        flux1arr_12(count_12) = flux1
;        flux2arr_12(count_12) = flux2
;        time_12(count_12) = time
;        count_12 = count_12+ 1
        
;     endif
        flux1arr_100(count_100) = flux1
        flux2arr_100(count_100) = flux2
        time_100(count_100) = time
        count_100 = count_100+ 1
     endif
     
     
  endfor
  
;   flux1arr_12 = flux1arr_12[0:count_12 - 1]
;  flux2arr_12 = flux2arr_12[0:count_12 - 1]
;  time_12 = time_12[0:count_12 - 1]
  flux1arr_100 = flux1arr_100[0:count_100 - 1]
  flux2arr_100 = flux2arr_100[0:count_100 - 1]
  time_100 = time_100[0:count_100 - 1]

  oplot, (time_100 - time_100(0))/3600., flux1arr_100/ flux1arr_100(0), psym = 1, color = greencolor
;  oplot, (time_12 - time_12(0))/3600., flux1arr_12/ flux1arr_12(0), psym = 1, color = redcolor
  
  oplot, (time_100 - time_100(0))/3600., flux1arr_100/ flux1arr_100(0), color = greencolor
;  oplot, (time_12 - time_12(0))/3600., flux1arr_12/ flux1arr_12(0),color = redcolor
  
;  oplot, (time_100 - time_100(0))/3600., flux2arr_100/ flux2arr_100(0), psym = 1, color = redcolor
;  oplot, (time_12 - time_12(0))/3600., flux2arr_12/ flux2arr_12(0), psym = 5, color = redcolor
  
;  oplot, (time_100 - time_100(0))/3600., flux2arr_100/ flux2arr_100(0), color = redcolor
;  oplot, (time_12 - time_12(0))/3600., flux2arr_12/ flux2arr_12(0),  color = redcolor

 ; ps_close, /noprint,/noid
  ps_end, /png


end
