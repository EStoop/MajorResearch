pro latent
  ;look at the latent tests in pc2
;there are three of the same sets of observations
;want to see if the latents are 100% repeatable.
  ps_start,filename='/Users/jkrick/iwic/latent_repeatability.ps'
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


;read in the flats for later
  fits_read, '//Users/jkrick/irac_pc/pc2/IRAC022300/irac_b1_fa_superskyflat_070709.fits', ch1flat, ch1flathead
  fits_read, '//Users/jkrick/irac_pc/pc2/IRAC022300/irac_b2_fa_superskyflat_070709.fits', ch2flat, ch2flathead


  for test = 0, 2 do begin  ;for each instance of the latent test
   
     if test eq 0 then begin
     ;dir_name =[ '/Users/jkrick/irac_pc/pc2/IRAC022300/bcd/'] ;,
        dir_name = ['/Users/jkrick/irac_pc/pc2/IRAC022300/bcd/0037837312']
        colorname = blackcolor
     endif

     if test eq 1 then begin
        dir_name =[ '/Users/jkrick/irac_pc/pc2/IRAC022300/bcd/0037837568']
        colorname = redcolor
     endif

     if test eq 2 then begin
        dir_name =[ '/Users/jkrick/irac_pc/pc2/IRAC022300/bcd/0037837824']
        colorname = bluecolor
     endif

     bcdch1list = ' '
     bcdch2list = ' '

     for j = 0, n_elements(dir_name) - 1 do begin
        CD,dir_name(j)
     
        command1 =  ' find '+ dir_name(j) +' -name "*bcd_fp.fits" >  bcd_list.txt'
        command2 = 'cat bcd_list.txt | grep IRAC.1. > bcd_ch1.txt'
        command3 = 'cat bcd_list.txt | grep IRAC.2. | grep -v IRAC.1 > bcd_ch2.txt'
     
        commands = [command1, command2, command3]
        for i = 0, n_elements(commands) -1 do spawn, commands(i)
        
        readcol,dir_name(j) + '/bcd_ch1.txt', bcdch1list_a, format="A", /silent, skipline = 2
        readcol,dir_name(j) + '/bcd_ch2.txt', bcdch2list_a, format="A", /silent, skipline = 2 ;skip the .6s and 12s frames

        bcdch1list = [bcdch1list , bcdch1list_a]
        bcdch2list = [bcdch2list , bcdch2list_a]

     endfor
  
 ; print, 'bcd1', bcdch1list

     xcenarr = [213,210]
     ycenarr = [46,46]
  
     for ch = 0, 1 do begin   ;for each channel
        if ch eq 0 then begin
           bcdlist = bcdch1list
           xcen = xcenarr[0]
           ycen = ycenarr[0]
           start = 4
        endif

        if ch eq 1 then begin
           bcdlist = bcdch2list
           xcen = xcenarr[1]
           ycen = ycenarr[1]
           start = 1
        endif
     
        flux1arr_100 =fltarr(n_elements(bcdlist))
        flux2arr_100 =fltarr(n_elements(bcdlist))
        time_100 = dblarr(n_elements(bcdlist))
        count_100 = 0

  
        for i = start, n_elements(bcdlist) - 1  do begin
           fits_read, bcdlist(i), data, imheader
     ;print, 'working on image', bcdlist(i)
  

     ;make the background positive

           if ch eq 0 then begin
              data = data * ch1flat
              data = data + 0.2 + 0.04  ;from zody_est in header + additional to bring positive
              data = data / ch1flat
           endif
          if ch eq 1 then begin
              data = data * ch2flat
              data = data + 0.17 ;from zody_est in header
              data = data / ch2flat
           endif

     ;get rid of NaN's
           data[where(finite(data) lt 1)] = 0


      ;do aperture photometry on star position 
           object = data[xcen-2:xcen+2,ycen-2:ycen+2]
           nobject = n_elements(object)
           sky = total(data[xcen-23:xcen+23, ycen-23:ycen+23])- total(data[xcen-15:xcen+15,ycen-15:ycen+15])  
           nannulus = n_elements(data[xcen-23:xcen+23, ycen-23:ycen+23]) - n_elements(data[xcen-15:xcen+15,ycen-15:ycen+15])    
           sky = sky / nannulus ;mean sky flux
           mmm, data[170:255, 1:90], skymode, sigma, skw
           flux1 =  total( object)  - skymode*nobject

     
     ;keep track of time of observation
           time= fxpar(imheader, 'SCLK_OBS')
           framtime = fxpar(imheader, 'FRAMTIME')
     
   
           if framtime eq 100 then  begin
              flux1arr_100(count_100) = flux1
              time_100(count_100) = time
              count_100 = count_100+ 1
           endif
     
     
        endfor

        flux1arr_100 = flux1arr_100[0:count_100 - 1]
        time_100 = time_100[0:count_100 - 1]

        print, 'flux1arr_100', flux1arr_100
        print, 'time_100', (time_100 - time_100(0))/3600.

        ;plot
        if ch eq 0  and test eq 0 then plot, (time_100 - time_100(0))/3600., flux1arr_100, psym = 2,$ 
                                             xtitle = 'hours', ytitle = 'Flux (counts in 100s)', title = 'hd166780 Kmag = 4.0',   $
                                             color = colorname, xrange=[.01,2.0], yrange= [0.1,1000],/ylog,/xlog
        if ch eq 0 and test ne 0 then oplot, (time_100 - time_100(0))/3600., flux1arr_100, psym = 2, color = colorname ;/ max(flux1arr_100)
        if ch eq 1 then oplot, (time_100 - time_100(0))/3600., flux1arr_100, psym = 1, color = colorname ;/ max(flux1arr_100)



        ;keep track of everything so I can use it later.
        ;this code is too ugly to deal with every time

        if ch eq 0 and test eq 0 then begin
           ch1_test0_x = (time_100 - time_100(0))/3600.
           ch1_test0_y = flux1arr_100
        endif
        if ch eq 0 and test eq 1 then begin
           ch1_test1_x = (time_100 - time_100(0))/3600.
           ch1_test1_y = flux1arr_100
        endif
        if ch eq 0 and test eq 2 then begin
           ch1_test2_x = (time_100 - time_100(0))/3600.
           ch1_test2_y = flux1arr_100
        endif
;--
        if ch eq 1 and test eq 0 then begin
           ch2_test0_x = (time_100 - time_100(0))/3600.
           ch2_test0_y = flux1arr_100
        endif
        if ch eq 1 and test eq 1 then begin
           ch2_test1_x = (time_100 - time_100(0))/3600.
           ch2_test1_y = flux1arr_100
        endif
        if ch eq 1 and test eq 2 then begin
           ch2_test2_x = (time_100 - time_100(0))/3600.
           ch2_test2_y = flux1arr_100
        endif


 
     endfor

;got this by hand from median combine below
     x = 1.73
     xerr = 0.39
     y1 = [0.25911969,0.29547089,0.29771935]
     y2 = [-0.33238834,-0.32438448,  -0.36428360] ;don't know how to deal with negative
     xyouts, x, y1(0), '*'
     xyouts, x, y1(1),'*', color = redcolor
     xyouts, x, y1(2), '*', color = bluecolor

     legend, ['test1','test2','test3'], color = [blackcolor, redcolor, bluecolor],linestyle = [0,0,0], /top, /right
     legend, ['ch1', 'ch2'] , psym = [2, 1],/bottom, /left
  endfor



;average together the 3 tests
ch1_y = fltarr(n_elements(ch1_test0_y))
ch2_y = fltarr(n_elements(ch2_test0_y))

 for i = 0, n_elements(ch1_test0_y) - 1 do ch1_y(i) = (ch1_test0_y(i) + ch1_test1_y(i) + ch1_test2_y(i)) / 3.
 for i = 0, n_elements(ch2_test0_y) - 1 do ch2_y(i) = (ch2_test0_y(i) + ch2_test1_y(i) + ch2_test2_y(i)) / 3.

 ch1_x = ch1_test2_x
 ch2_x = ch2_test2_x

;add in the median combine value
ch1_x = [ch1_x, 1.73]
ch2_x = [ch2_x, 1.73]
ch1_y = [ch1_y, mean(y1)]
ch2_y = [ch2_y, mean(y2)]

save, filename = '/Users/jkrick/IWIC/latent.sav', ch1_x, ch2_x, ch1_y, ch2_y

;--------------------------------------------------------------------------
;now measure flux on median combine of the skydarks after the latent
;test
med_names = ['/Users/jkrick/irac_pc/pc2/IRAC022300/bcd/0037838080','/Users/jkrick/irac_pc/pc2/IRAC022300/bcd/0037838336','/Users/jkrick/irac_pc/pc2/IRAC022300/bcd/0037838592']

for i = 0, n_elements(med_names) - 1 do begin
   cd, med_names(i)
   fits_read, 'meddark_ch1.fits', data1, header1
   fits_read, 'meddark_ch2.fits', data2, header2

     ;get rid of NaN's
     data1[where(finite(data1) lt 1)] = 0
     data1 = data1 + .25
     object = data1[xcen-2:xcen+2,ycen-2:ycen+2]
     nobject = n_elements(object)
     mmm, data1[170:255, 1:90], skymode, sigma, skw
     flux1 =  total( object)  - skymode*nobject
     print, 'flux1', flux1
     ;keep track of time of observation
     time= fxpar(imheader, 'SCLK_OBS')

     print, 'time 1', (time - time_100(0))/3600.

     ;get rid of NaN's
     data2[where(finite(data2) lt 1)] = 0
     data2 = data2 + .25
     object = data2[xcen-2:xcen+2,ycen-2:ycen+2]
     nobject = n_elements(object)
     mmm, data2[170:255, 1:90], skymode, sigma, skw
     flux2 =  total( object)  - skymode*nobject
     print, 'flux2', flux2
     ;keep track of time of observation
     time= fxpar(imheader, 'SCLK_OBS')
     print, 'time 2', (time - time_100(0))/3600.
  endfor

   

;-------------------------------------------------------------

 ; ps_close, /noprint,/noid
  ps_end,/png


end


;  ch1_test0_x = fltarr(100)
;  ch1_test0_y = fltarr(100)
;  ch2_test0_x = fltarr(100)
;  ch2_test0_y = fltarr(100)
;  ch1_test1_x = fltarr(100)
;  ch1_test1_y = fltarr(100)
;  ch2_test1_x = fltarr(100)
;  ch2_test1_y = fltarr(100)
;  ch1_test2_x = fltarr(100)
;  ch1_test2_y = fltarr(100)
;  ch2_test2_x = fltarr(100)
;  ch2_test2_y = fltarr(100)
