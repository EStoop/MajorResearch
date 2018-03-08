pro calc_test_bkgd, planetname

  ;;have run different background regions on WASP14b
  ;;read in the save files and calculate
  ;;differences in the photometry between
  ;;each different version

  ;;for difference = maybe the mean stddev amongst the different
  ;;snapshots?

  planetinfo = create_planetinfo()
  aorname= planetinfo[planetname, 'aorname_ch2']

  if planetname eq 'WASP-14b' then begin
     dirname = '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/'
     bs = 5E-6
     ;;-----------
     restore, dirname + 'WASP-14b_phot_ch2_2.25000_selfbkgd.sav'
     selfbkgdarr = fltarr(n_elements(aorname))
     for a =5, n_elements(aorname) - 1 do begin ;should be skipping the stares 
        selfbkgdarr(a) = stddev(planethash[aorname(a),'corrflux'],/nan)
     endfor
     plothist, selfbkgdarr, xhist, yhist, /noplot, bin = bs
     p1 = plot(xhist, yhist, xtitle = 'standard deviation', ytitle = 'number', color = 'blue', thick =2, xrange = [0.00018, 0.0003])
     print, 'mean stddev selfbkgd', mean(selfbkgdarr)
     
     
     ;;-----------
     restore, dirname + 'WASP-14b_phot_ch2_2.25000_nomb.sav'
     nombarr = fltarr(n_elements(aorname))
     for a =5, n_elements(aorname) - 1 do begin ;should be skipping the stares 
        nombarr(a) = stddev(planethash[aorname(a),'corrflux'],/nan)
     endfor
     plothist, nombarr, xhist, yhist, /noplot, bin = bs
     p1 = plot(xhist, yhist, color = 'red', thick =2, overplot = p1)
     print, 'mean stddev nomb', mean(nombarr)
     
     
     
     ;;-----------
     ;;using the pmap dataset which also has 3_15 background
     restore, dirname + 'WASP-14b_phot_ch2_2.25000_160126.sav'
     threefifteen = fltarr(n_elements(aorname))
     for a =5, n_elements(aorname) - 1 do begin ;should be skipping the stares 
        threefifteen(a) = stddev(planethash[aorname(a),'corrflux'],/nan)
     endfor
     plothist, threefifteen, xhist, yhist, /noplot, bin =bs
     p1 = plot(xhist, yhist, color = 'green', thick =2, overplot = p1)
     print, 'mean stddev 3_15', mean(threefifteen)
     
     
     ;;-----------
     restore, dirname + 'WASP-14b_phot_ch2_2.25000_715.sav'
     sevenfifteen = fltarr(n_elements(aorname))
     for a =5, n_elements(aorname) - 1 do begin ;should be skipping the stares 
        sevenfifteen(a) = stddev(planethash[aorname(a),'corrflux'],/nan)
     endfor
     plothist, sevenfifteen, xhist, yhist, /noplot, bin = bs
     p1 = plot(xhist, yhist, color = 'purple', thick =2, overplot = p1)
     print, 'mean stddev 7_15', mean(sevenfifteen)
     
     
     ;;-----------
     restore, dirname + 'WASP-14b_phot_ch2_2.25000_150723_newtime.sav'
     threeseven = fltarr(n_elements(aorname))
     for a =5, n_elements(aorname) - 1 do begin ;should be skipping the stares 
        threeseven(a) = stddev(planethash[aorname(a),'corrflux'],/nan)
     endfor
     plothist, threeseven, xhist, yhist, /noplot, bin = bs
     p1 = plot(xhist, yhist, color = 'black', thick =2, overplot = p1)
     print, 'mean stddev 3-7', mean(threeseven)
     
     print, 'delta', (mean(threeseven) - mean(threefifteen))/ mean(threefifteen)
  endif


    if planetname eq 'HD158460' then begin
     dirname = '/Users/jkrick/irac_warm/pcrs_planets/HD158460/'
     bs = 1E-4
     ;;-----------
     restore, dirname + 'HD158460_phot_ch2_2.25000_150723.sav'
     threeseven = fltarr(n_elements(aorname))
     for a =5, n_elements(aorname) - 1 do begin ;should be skipping the stares 
        threeseven(a) = stddev(planethash[aorname(a),'corrflux'],/nan)
     endfor
          print, 'threeseven', threeseven

     plothist, threeseven, xhist, yhist, /noplot, bin = bs
     p1 = plot(xhist, yhist, xtitle = 'standard deviation', ytitle = 'number', color = 'blue', thick =2)
     print, 'mean stddev threeseven', mean(threeseven)
     
     ;;-----------
     restore, dirname + 'HD158460_phot_ch2_2.25000_160126.sav'
     threefifteen = fltarr(n_elements(aorname))
     for a =5, n_elements(aorname) - 1 do begin ;should be skipping the stares 
        threefifteen(a) = stddev(planethash[aorname(a),'corrflux'],/nan)
     endfor
     print, 'threefifteen', threefifteen
     plothist, threefifteen, xhist, yhist, /noplot, bin = bs
     p1 = plot(xhist, yhist, color = 'red', thick =2 , overplot = p1)
     print, 'mean stddev threefifteen', mean(threefifteen)

     print, 'delta', (mean(threeseven) - mean(threefifteen))/ mean(threefifteen)
   
  endif
    
  end
  
