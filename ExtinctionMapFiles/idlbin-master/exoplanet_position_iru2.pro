pro exoplanet_position_iru2

  colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY'    ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'TAN', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET' , 'WHEAT'    , 'WHITE_SMOKE'   , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU', ' PINK',  'PLUM',  ' POWDER_BLUE' ]

  ra_ref =     286.97079
  dec_ref=    46.868328


  dir =[ '/Users/jkrick/iracdata/flight/IWIC/IRAC028400/bcd/0042620160'  ,'/Users/jkrick/irac_warm/gyros/r39525376/ch2/bcd']
  for d = 0, 1 do begin  ;look at both aors
     
     CD, dir[d]                    ; change directories to the correct AOR directory

     if d eq 0 then command = 'ls *bcd_fp.fits > /Users/jkrick/irac_warm/gyros/cbcdlist_42620160.txt'
     if d eq 1 then command = 'ls *cbcd.fits > /Users/jkrick/irac_warm/gyros/cbcdlist_39525376.txt'
     spawn, command
     
     if d eq 0 then readcol,'/Users/jkrick/irac_warm/gyros/cbcdlist_42620160.txt',fitsname, format = 'A', /silent
     if d eq 1 then readcol,'/Users/jkrick/irac_warm/gyros/cbcdlist_39525376.txt',fitsname, format = 'A', /silent
     yarr =   fltarr(n_elements(fitsname))
     xarr =  fltarr(n_elements(fitsname))
     timearr =  fltarr(n_elements(fitsname))
     c = 0
     
     for i =0, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
                                ;print, 'working on ', fitsname(i)         
                                ;now were did it really point to?
        
        get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5],ra = ra_ref, dec = dec_ref, /silent
        
        header = headfits(fitsname(i)) ;
        
        xarr[i] = x_center
        yarr[i] = y_center
        timearr[i] = sxpar(header, 'SCLK_OBS')
        
     endfor                     ; for each fits file in the AOR
     
     timearr = timearr - timearr[0] 
                                ;convert to hrs
     timearr = timearr / 60. / 60.
     
                                ;print, 'timearr', timearr
                                ;print, 'xarr', xarr
                                ;print, 'yarr', yarr
     if d eq 0 then begin
        s = plot(timearr, yarr,'6r1+' ,xtitle = 'Time(hrs)',ytitle = 'Y pix', yrange = [128.0, 128.3], xrange = [0,8], title = 'Ch2 6s Koi13_1 Aug_2011')
        t = plot(timearr, xarr,'6r1+' ,xtitle = 'Time(hrs)',ytitle = 'X pix', yrange = [126.5, 126.7], xrange = [0,8], title = 'Ch2 6s Koi13_1 Aug_2011')
     endif
   if d eq 1 then begin
        s = plot(timearr, yarr,'6r1+' ,xtitle = 'Time(hrs)',ytitle = 'Y pix' , yrange = [127.7, 128.0], xrange = [0,8],title = 'Ch2 30s Koi13_1 Aug_2010')
        t = plot(timearr, xarr,'6r1+' ,xtitle = 'Time(hrs)',ytitle = 'X pix', yrange = [126.75, 126.95], xrange = [0,8], title = 'Ch2 30s Koi13_1 Aug_2010')
     endif

  endfor ; for each AOR


 end
  
