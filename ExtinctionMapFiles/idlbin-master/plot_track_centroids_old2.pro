pro plot_track_centroids, run_data = run_data

   if keyword_set(run_data) then begin

     savenames = ['/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/track_centroids_pixval.sav', '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/track_centroids_pixval_2.sav']

     for s = 1, 1 do begin
        restore, savenames(s)
        print, 'planethhash at beginning'
        help, planethash
        aorlist = planethash.keys()
        sigmax = fltarr(n_elements(aorlist))*alog10(-1)
        sigmay = sigmax
        xjd = sigmax
        xdrift = sigmax         ; *alog10(-1)
        ydrift = sigmax         ;*alog10(-1)
        dpa = sigmax
        pa = sigmax
        exptimearr = sigmax
        short_drift = sigmax
        slope_drift = sigmax
        startyear = sigmax
        pkperiod = list(length=2) ;n_elements(aorlist))
        pkstrength =  list(length=2) 
        pktime = list(length=2)
        
        for n = 0,  n_elements(aorlist) - 1 do begin
           print, '--------------'
           print, 'working on ',n, ' ',aorlist(n), n_elements(planethash[aorlist(n)].xcen)
           timearr = planethash[aorlist(n)].timearr
           bmjdarr = planethash[aorlist(n)].bmjdarr
           short_drift[n] =  planethash[aorlist(n)].short_drift
           ;;slope_drift[n] =  planethash[aorlist(n)].slope_drift
           if planethash[aorlist(n)].haskey('slope_drift') gt 0 then slope_drift[n] =  planethash[aorlist(n)].slope_drift else slope_drift[n] = alog10(-1)
           print, 'slope drift', slope_drift[n]
           xjd[n] = bmjdarr(0) + 2400000.5
           CALDAT, xjd[n], Month, Day, year
           startyear[n] = year
           
           time0 = timearr(0)
           timearr = (timearr - time0)/60./60. ; now in hours instead of sclk
           exptimearr[n] = planethash[aorlist(n)].exptime
           
           ;;-------------------------------------
           ;;sigmax & sigmay &sigmaxy vs. time
           ;;not sure what sigmaxy is?
           if max(timearr) gt 1.2 then begin
              sigmax[n] = stddev(planethash[aorlist(n)].xcen,/nan)
              sigmay[n] = stddev(planethash[aorlist(n)].ycen,/nan)
              print, 'sigma x, y ', sigmax[n], sigmay[n]
              
              ;;-------------------------------------
              ;;long term xdrift vs.& y drift
              start = [.05,15.0]
              ;;don't have errors in position, instead, fake it.
              noise = fltarr(n_elements(planethash[aorlist(n)].xcen))
              noise = noise + 1.
              xcenfit= MPFITFUN('linear',timearr, planethash[aorlist(n)].xcen, noise, start,/Quiet)
              xdrift[n] = xcenfit(0)
              ycenfit= MPFITFUN('linear',timearr, planethash[aorlist(n)].ycen, noise, start,/Quiet)
              ydrift[n] = ycenfit(0)
              ;;do some quick paring down of the data
              xnum = findgen(n_elements(timearr))
              i = where(xnum mod 10 lt 1) ;pick out the odd numbers only
              tx = timearr(i)
              ycen = planethash[aorlist(n)].ycen
              ty = ycen(i)
              ;;pl = plot(timearr, ycen, title = aorlist(n), xtitle = 'time(hrs)', yrange = [mean(planethash[aorlist(n)].ycen,/nan) -0.5, mean(planethash[aorlist(n)].ycen,/nan) +0.5])
              ;;pl = plot(timearr, ycenfit(0)*timearr + ycenfit(1), color = 'red', overplot = pl)
              ;;XX don't want to keep this value if dithered.
           endif
           
           ;;-------------------------------------
           ;;short term drift
           ;;want duration and slope.
           
           
           
           ;;-------------------------------------
           ;;delta pitch angle
           prepitchangle = planethash[aorlist(n)].prepitchangle
           dpa[n] = prepitchangle(n_elements(prepitchangle) - 2) - planethash[aorlist(n)].pitchangle
           pa[n] = planethash[aorlist(n)].pitchangle
           
           ;;-------------------------------------
           ;;width of the peak in the power spectrim at 30min
           ;;don't need to do this for pre-AORs
           ;;periodogram
           if max(timearr) gt 1.2 then begin
              xday = timearr*60. ;in minutes             
              ycen = planethash[aorlist(n)].ycen
              bad = where(finite(ycen) lt 1,nbad)
              print, 'number of nans ', nbad
              if nbad eq n_elements(ycen) then CONTINUE ; no good data points
              if nbad gt 0  then remove, bad, xday, ycen
              result = LNP_TEST(xday, ycen,/double, WK1 = wk1, WK2 = wk2, JMAX = jmax)
              
              ;;find the peaks above N* the random level
              ;;but N appears to be a function of the number of data points.
              d0 = wk2 - shift(wk2, 1)
              d1 = wk2 - shift(wk2, -1)
              pk = where(d0 gt 0 and d1 gt 0, npk)
              peakheight = wk2[pk]
              peakfreq = wk1[pk]
              peakperiod = 1/peakfreq
              ;;print, 'npk', npk
              
              ;;define significance relative to other peaks between 5 -  60 minutes
              short = where(peakperiod gt 5 and peakperiod lt 60)
              peakheight = peakheight(short)
              peakfreq = peakfreq(short)
              peakperiod = peakperiod(short)
              mn = robust_mean(peakheight, 4)
              sig = robust_sigma(peakheight)
              nsig = 5
              realpk = where(peakheight gt mn + nsig*sig and peakperiod gt 30. and peakperiod le 50., nrealpk)
              
              ;;concatenating arrays- messy, but not sure how else to do it
              if n eq 0 then begin
                 if nrealpk gt 0 then begin
                    pkperiod = peakperiod(realpk)
                    pkstrength =(peakheight(realpk) -mn) /sig ;in units of some sort of significance
                    pt = fltarr(nrealpk)
                    pt(*) = xjd(n)
                    pktime = pt
                 endif else begin
                    pkperiod = 0
                    pkstrength = 0
                    pktime = xjd(n)
                 endelse
              endif else begin
                 if nrealpk gt 0 then begin
                    pkperiod = [pkperiod, peakperiod(realpk)]
                    pkstrength =[pkstrength, (peakheight(realpk) -mn) /sig ] ;in units of some sort of significance
                    pt = fltarr(nrealpk)
                    pt(*) = xjd(n)
                    pktime = [pktime, pt]
                 endif else begin
                    pkperiod =[pkperiod, 0]
                    pkstrength = [pkstrength,0]
                    pktime = [pktime, xjd(n)]
                 endelse
                 
              endelse
              
           endif
           
        endfor
        save, sigmax, sigmay, xjd, xdrift, ydrift, dpa, pa, exptimearr, short_drift, slope_drift, startyear, filename = '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/plot_track_centroids.sav'
        help, planethash
        print, planethash.remove(/all)
        help, planethash
     endfor                     ; for each save file
     
  endif else begin ;;keyword_set run_data
     print, 'restoring data'
     restore, '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/plot_track_centroids.sav'
     
  endelse
  xjd = xjd[0:n-1]
  sigmax = sigmax[0:n-1]
  sigmay = sigmay[0:n-1]
  xdrift = xdrift[0:n-1]
  ydrift = ydrift[0:n-1]
  dpa = dpa[0:n-1]
  pa = pa[0:n-1]
  exptimearr = exptimearr[0:n-1]
  short_drift = short_drift[0:n-1]
  slope_drift = slope_drift[0:n-1]
  startyear = startyear[0:n-1]

  ;;set up color coding by exposure time
  colorarr = intarr(3, n_elements(exptimearr))
  zerop02 = where(exptimearr lt 0.05 and exptimearr gt 0., n0p02)
  zerop1 = where(exptimearr lt 0.2 and exptimearr gt 0.05,n0p1)
  zerop4 = where(exptimearr lt 0.4 and exptimearr gt 0.3,n0p4)
  twop0 = where(exptimearr gt 1.0 and exptimearr lt 2.0, n2p0)
  sixp0 = where(exptimearr gt 4.0 and exptimearr lt 7.0, n6p0)
  twelve = where(exptimearr gt 9.0 and exptimearr lt 15,n12)
  thirty = where(exptimearr gt 20 and exptimearr lt 90,n30)
  hundred = where(exptimearr gt 50 and exptimearr lt 200, n100)
  for c = 0, n0p02 - 1 do colorarr[*,zerop02(c)] = [128,0,0];'maroon'
  for c = 0, n0p1 - 1 do colorarr[*,zerop1(c)] = [255,0,0];'red'
  for c = 0, n0p4 - 1 do colorarr[*,zerop4(c)] = [255,69,0];'orange_red'
  for c = 0, n2p0 - 1 do colorarr[*,twop0(c)] = [238,118,33];'dark_orange'
  for c = 0, n6p0 - 1 do colorarr[*,sixp0(c)] = [127,255,0];'lime_green'
  for c = 0, n12 - 1 do colorarr[*,twelve(c)] =[64,224,208]; 'aqua'
  for c = 0, n30 - 1 do colorarr[*,thirty(c)] = [0,0,255];'blue'
  for c = 0, n100 - 1 do colorarr[*,hundred(c)] = [155,48,255];'Purple'

  print, 'set up colorarr'
  ;;set up color coding by year of observation
  ;;coloryear = bytarr(3, n_elements(exptimearr))
  coloryear = strarr( n_elements(exptimearr))
  twenty10 = where(startyear eq 2010, n2010)
  twenty11 = where(startyear eq 2011, n2011)
  twenty12 = where(startyear eq 2012, n2012)
  twenty13 = where(startyear eq 2013, n2013)
  twenty14 = where(startyear eq 2014, n2014)
  twenty15 = where(startyear eq 2015, n2015)
  twenty16 = where(startyear eq 2016, n2016)
  twenty17 = where(startyear eq 2017, n2017)
  twenty18 = where(startyear eq 2018, n2018)
  
  ;;for c = 0, n2010 - 1 do coloryear[*,twenty10(c)] = [128,0,0]   ;'maroon'
  ;;for c = 0, n2011 - 1 do coloryear[*,twenty11(c)] = [255,0,0];'red'
  ;;for c = 0, n2012 - 1 do coloryear[*,twenty12(c)] = [255,69,0];'orange_red'
  ;;for c = 0, n2013 - 1 do coloryear[*,twenty13(c)] = [238,118,33];'dark_orange'
  ;;for c = 0, n2014 - 1 do coloryear[*,twenty14(c)] = [127,255,0];'lime_green'
  ;;for c = 0, n2015 - 1 do coloryear[*,twenty15(c)] =[64,224,208]; 'aqua'
  ;;for c = 0, n2016 - 1 do coloryear[*,twenty16(c)] = [0,0,255];'blue'
  ;;for c = 0, n2017 - 1 do coloryear[*,twenty17(c)] = [155,48,255];'Purple'
  ;;for c = 0, n2018 - 1 do coloryear[*,twenty18(c)] = [51,0,102];'Purple'
  coloryear(twenty10) = 'black'
  coloryear(twenty11) = 'blue'
  for ci = 0, n_elements(startyear) - 1 do print, 'year, color ', startyear(ci), xjd(ci), coloryear(ci)
  ;;------------------------------------------------
  ;;sigmax & sigmay &sigmaxy vs. time
  ;;------------------------------------------------

  ;;psx = plot(xjd, sigmax,'1s', sym_size = 0.5,   /sym_filled , ytitle = 'X Size of cloud (stddev in position)', $
  ;;           XTICKFORMAT='(C(CMoA,1x,CYI))', xtickunits = ['Time'], xminor =11, yrange = [0, 0.2],$`
  ;;           vert_colors = colorarr)
  ;;psy = plot(xjd, sigmay,'1s', sym_size = 0.5,   sym_filled = 1, ytitle = 'Y Size of cloud (stddev in position)', $
  ;;           XTICKFORMAT='(C(CDI,1x,CMoA,1x,CYI))', xtickunits = ['Time'], xminor = 11, yrange = [0, 0.2],$
  ;;          vert_colors = colorarr)
  
  ;;------------------------------------------------
  ;;Long term xdrift vs. ydrift
  ;;------------------------------------------------
  
  maxdrift = 0.05
  help, xdrift
  help, ydrift
  pxydrift = plot(xdrift, ydrift,  '1s', sym_size = 0.5,   sym_filled = 1,  ytitle = 'Long Term Y drift (px/hr)', $
                  xtitle = 'Long Term X drift (px/hr)', yrange = [-0.04, 0.04], xrange =[-0.04,0.04],$;[-1*maxdrift,maxdrift],$
                  vert_colors = colorarr)
  print, 'done plot'
  t1 =text(0.02,0.034, '0.02s', color = [128,0,0],/data)
  t1 =text(0.02,0.029, '0.1s',  color= [255,0,0],/data)    ;'red'
  t1 =text(0.02,0.024, '0.4s',  color= [255,69,0] ,/data)  ;'orange_red'
  t1 =text(0.02,0.019, '2.0s',  color= [238,118,33],/data) ;'dark_orange'
  t1 =text(0.02,0.014, '6.0s',  color= [127,255,0],/data)  ;'lime_green'
  t1 =text(0.02,0.009, '12s',  color=[64,224,208] ,/data) ; 'aqua'
  t1 =text(0.02,0.004, '30s',  color= [0,0,255] ,/data)   ;'blue'
  t1 =text(0.02,-0.001, '100s',  color= [155,48,255],/data) ;'Purple'
  print, 'done text'
  pl1 = polyline([0.0,0.0], [-0.04, 0.04], /data, target = pxydrift, color = 'black')
  pl1 = polyline( [-0.04, 0.04], [0.0,0.0], /data, target = pxydrift, color = 'black')
  print, 'done long term xdrift vs. ydrift'
  ;;------------------------------------------------
  ;;Long term xdrift & y drift vs. time
  ;;------------------------------------------------

  ;;first get rid of outliers
  ;print, 'xdrift', xdrift
  meanclip,xdrift,  mean, meansig;3 Sigma clipping
  gx = where(xdrift gt (mean - 3*meansig) and xdrift lt (mean + 3*meansig), ngood)
  print, 'num good out of total',ngood, n_elements(xdrift)
  xjdgx = xjd(gx)
  xdriftgx = xdrift(gx)
  meanclip,ydrift,  mean, meansig ;3 Sigma clipping
  gy = where(ydrift gt (mean - 3*meansig) and ydrift lt (mean + 3*meansig), ngood)
  print, 'num, ngood',ngood, n_elements(ydrift)
  xjdgy = xjd(gy)
  ydriftgy = ydrift(gy)

  
  pdrift = plot(xjdgx, xdriftgx, '1s', sym_size = 0.7,   sym_filled = 1,  ytitle = 'Long Term Drift (px/hr)', $
                yrange = [-0.04, 0.04],$;[-1*maxdrift, maxdrift+0.01],$
                color = 'blue', XTICKFORMAT='(C(CMoA,1x,CYI))', xtickunits = ['Time'], xminor =11, name='xdrift')
  pdrift2 = plot(xjdgy, ydriftgy, '1s', sym_size = 0.7,   sym_filled = 1,  overplot = pdrift, color = 'red', name = 'ydrift')

  ;;fit a linear function to these data.
  start = [.0005,0.00001]
        ;;don't have errors in drift, instead, fake it.
  xnoise = fltarr(n_elements(xdriftgx))
  xnoise = xnoise + 1.
  xdriftfit= MPFITFUN('linear',xjdgx, xdriftgx, xnoise, start);,/Quiet)
  pdrift3 = plot(xjdgx, xdriftfit(0)*xjdgx + xdriftfit(1), color = 'blue', overplot = pdrift)
  ynoise = fltarr(n_elements(ydriftgy))
  ynoise = ynoise + 1.
  ydriftfit= MPFITFUN('linear',xjdgy, ydriftgy, ynoise, start);,/Quiet)
  pdrift4 = plot(xjdgy, ydriftfit(0)*xjdgy + ydriftfit(1), color = 'red', overplot = pdrift)

  lpd = legend(target = [pdrift, pdrift2], /auto_text_color, position = [xjd[20], maxdrift-0.01],/data)
  ;;------------------------------------------------
  ;;delta pitch angle
  ;;------------------------------------------------

  ;;pdpa = plot( dpa, xdrift, '1s', sym_size = 0.5,   sym_filled = 1, xtitle = 'Change in pitch angle from previous observation', $
  ;;            ytitle = 'Long Term Drift (px/hr)', color = 'blue',yrange = [-1*maxdrift,maxdrift], name = 'xdrift')
  ;;pdpa2 = plot(dpa, ydrift,'1s', sym_size = 0.5,   sym_filled = 1, overplot = pdpa, color = 'red', name = 'ydrift')
  ;;ldpa = legend(target = [pdpa, pdpa2], /auto_text_color, position = [dpa[20], maxdrift-0.01],/data)
  
  ;;periodogram fun
  ;;how do I collapse a list into a single array?
 ;; pperiod = bubbleplot(pktime, pkperiod, /shaded, magnitude = pkstrength , exponent = 0.5, $
 ;;                      ytitle = 'Period of the power spectrum peaks (min)',XTICKFORMAT='(C(CMoA,1x,CYI))', $
 ;;                      xtickunits = ['Time'], xminor =11 , yrange = [20,60])
     

  ;;------------------------------------------------
  ;;short term drift plotting
  ;;------------------------------------------------
  ;;worry about zero vs. nan vs. negative.
  drift_dist = slope_drift * short_drift  ;;oops this is what Carl plots

  ;;testplot = plot(xjd, startyear, '1D', xtitle = 'xjd', ytitle = 'startyear')
  timeshortdrift = plot(xjd(twenty10), short_drift(twenty10), '1D', sym_size = 1.0,   /sym_filled , ytitle = 'Duration (hours)', $
                        XTICKFORMAT='(C(CMoA,1x,CYI))', xtickunits = ['Time'], xminor =11,ytickinterval = 0.4,$
                        xshowtext =0, position = [0.2,0.65,0.9,0.9], title = 'Ycen Short Term Drift',color = 'sky blue' )
  timeshortdrift = plot(xjd(twenty11), short_drift(twenty11), '1D', sym_size = 1.0,   /sym_filled , color = 'deep sky blue', $
                        overplot = timeshortdrift )
 timeshortdrift = plot(xjd(twenty12), short_drift(twenty12), '1D', sym_size = 1.0,   /sym_filled , color = 'royal blue', $
                        overplot = timeshortdrift )
  
  timedriftdist = plot(xjd(twenty10), drift_dist(twenty10), '1D', sym_size = 1.0,   /sym_filled , ytitle = 'Length (pixels)', $
                       XTICKFORMAT='(C(CMoA,1x,CYI))', xtickunits = ['Time'], xminor =11,$
                       xshowtext=0,/current, position = [0.2, 0.38, 0.9, 0.63],ytickinterval = 0.3,color = 'sky blue')
  timedriftdist = plot(xjd(twenty11), drift_dist(twenty11), '1D', sym_size = 1.0,   /sym_filled ,overplot = timedriftdist, $
                       color = 'deep sky blue')
  timedriftdist = plot(xjd(twenty12), drift_dist(twenty12), '1D', sym_size = 1.0,   /sym_filled ,overplot = timedriftdist, $
                       color = 'royal blue')
  
  timeslopedrift = plot(xjd(twenty10), slope_drift(twenty10), '1D', sym_size = 1.0,   /sym_filled , ytitle = 'Slope (pix/hr)', $
                        XTICKFORMAT='(C(CMoA,1x,CYI))', xtickunits = ['Time'], xminor =11,$
                        /current, position = [0.2,0.11, 0.9,0.36])
  timeslopedrift = plot(xjd(twenty11), slope_drift(twenty11), '1D', sym_size = 1.0,   /sym_filled ,overplot = timeslopedrift, $
                        color = 'deep sky blue')
  timeslopedrift = plot(xjd(twenty12), slope_drift(twenty12), '1D', sym_size = 1.0,   /sym_filled ,overplot = timeslopedrift, $
                        color = 'royal blue')
 
  ;;-------------
  timeshortdrift = plot(pa(twenty10), short_drift(twenty10), '1D', sym_size = 1.0,   /sym_filled , ytitle = 'Duration (hours)', $
                        xshowtext =0, position = [0.2,0.65,0.9,0.9], ytickinterval = 0.4,title = 'Ycen Short Term Drift', $
                        color = 'black')
  timeshortdrift = plot(pa(twenty11), short_drift(twenty11), '1D', sym_size = 1.0,   /sym_filled, overplot = timeshortdrift, $
                        color = 'deep sky blue' )
  timeshortdrift = plot(pa(twenty12), short_drift(twenty12), '1D', sym_size = 1.0,   /sym_filled, overplot = timeshortdrift, $
                        color = 'royal blue' )
 
  timedriftdist = plot(pa(twenty10), drift_dist(twenty10), '1D', sym_size = 1.0,   /sym_filled , ytitle = 'Length (pixels)', $
                       xshowtext=0,/current, position = [0.2, 0.38, 0.9, 0.63],ytickinterval = 0.3, color = 'black')
  timedriftdist = plot(pa(twenty11), drift_dist(twenty11), '1D', sym_size = 1.0,   /sym_filled , overplot = timedriftdist, $
                       color = 'deep sky blue')
   timedriftdist = plot(pa(twenty12), drift_dist(twenty12), '1D', sym_size = 1.0,   /sym_filled , overplot = timedriftdist, $
                       color = 'royal blue')
 
  timeslopedrift = plot(pa(twenty10), slope_drift(twenty10), '1D', sym_size = 1.0,   /sym_filled , ytitle = 'Slope (pix/hr)', $
                        /current, position = [0.2,0.11, 0.9,0.36], xtitle = 'Pitch Angle', color = 'black')
  timeslopedrift = plot(pa(twenty11), slope_drift(twenty11), '1D', sym_size = 1.0,   /sym_filled ,overplot = timeslopedrift, $
                        color = 'deep sky blue')
  timeslopedrift = plot(pa(twenty12), slope_drift(twenty12), '1D', sym_size = 1.0,   /sym_filled ,overplot = timeslopedrift, $
                        color = 'royal blue')
 
  ;;--------------
  
  timeshortdrift = plot(dpa(twenty10), short_drift(twenty10), '1D', sym_size = 1.0,   /sym_filled , ytitle = 'Duration (hours)', $
                        xshowtext =0, position = [0.2,0.65,0.9,0.9], ytickinterval = 0.4, title = 'Ycen Short Term Drift', $
                        color = 'black')
  timeshortdrift = plot(dpa(twenty11), short_drift(twenty11), '1D', sym_size = 1.0,   /sym_filled, overplot = timeshortdrift, $
                        color = 'deep sky blue')
  timeshortdrift = plot(dpa(twenty12), short_drift(twenty12), '1D', sym_size = 1.0,   /sym_filled, overplot = timeshortdrift, $
                        color = 'royal blue')
 
  timedriftdist = plot(dpa(twenty10), drift_dist(twenty10), '1D', sym_size = 1.0,   /sym_filled , ytitle = 'Length (pixels)', $
                       xshowtext=0,/current, position = [0.2, 0.38, 0.9, 0.63],ytickinterval = 0.3, color = 'black')
  timedriftdist = plot(dpa(twenty11), drift_dist(twenty11), '1D', sym_size = 1.0,   /sym_filled , overplot = timedriftdist,$
                       color = 'deep sky blue')
  timedriftdist = plot(dpa(twenty12), drift_dist(twenty12), '1D', sym_size = 1.0,   /sym_filled , overplot = timedriftdist,$
                       color = 'royal blue')
  
  timeslopedrift = plot(dpa(twenty10), slope_drift(twenty10), '1D', sym_size = 1.0,   /sym_filled , ytitle = 'Slope (pix/hr)', $
                        /current, position = [0.2,0.11, 0.9,0.36], xtitle = 'Delta Pitch Angle', color = 'black')
  timeslopedrift = plot(dpa(twenty11), slope_drift(twenty11), '1D', sym_size = 1.0,   /sym_filled ,overplot = timeslopedrift, $
                        color = 'deep sky blue')
   timeslopedrift = plot(dpa(twenty12), slope_drift(twenty12), '1D', sym_size = 1.0,   /sym_filled ,overplot = timeslopedrift, $
                        color = 'royal blue')
  
  
  ;;sdrift = plot(dpa, short_drift, '1s', sym_size = 0.5,   /sym_filled , ytitle = 'duration of ycen short term drift', $
  ;;              xtitle = 'Change in pitch angle from previous observation', vert_colors = colorarr)
  ;;sdrift = plot(pa, short_drift, '1s', sym_size = 0.5,   /sym_filled , ytitle = 'duration of ycen short term drift', $
  ;;              xtitle = 'Pitch angle', vert_colors = colorarr)
  ;;sdrift = plot(dpa, slope_drift, '1s', sym_size = 0.5,   /sym_filled , ytitle = 'slope of ycen short term drift', $
  ;;              xtitle = 'Change in pitch angle from previous observation', vert_colors = colorarr)
  ;;sdrift = plot(pa, slope_drift, '1s', sym_size = 0.5,   /sym_filled , ytitle = 'slope of ycen short term drift', $
  ;;              xtitle = 'Pitch angle', vert_colors = colorarr)

  ;;sdrift = plot(dpa, drift_dist, '1s', sym_size = 0.5,   /sym_filled , ytitle = 'total drift (pixels)', $
  ;;              xtitle = 'Change in pitch angle from previous observation', vert_colors = colorarr)
  ;;sdrift = plot(pa, drift_dist, '1s', sym_size = 0.5,   /sym_filled , ytitle = 'total drift (pixels)', $
  ;;              xtitle = 'Pitch angle', vert_colors = colorarr)

  
  
  ;;if keyword_set(run_data) then save, /variables, filename = '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/plot_track_centroids.sav'
  
end



    ;;plothist, planethash[aorlist(n)].xcen, xhist, yhist, bin = 0.05,/noplot
     ;;pg = barplot(xhist, yhist, title = aorlist(n))
     ;;start = [0.,10., 2.]
     ;;noise = fltarr(n_elements(yhist))
     ;;noise[*] = 1                                              ;equally weight the values
     ;;result= MPFITFUN('mygauss',xhist,yhist, noise, start) ;/quiet   ; fit a gaussian to the histogram sorted data


       ;; xnum = findgen(n_elements(xday))
       ;; i = where(xnum mod 10 lt 1) ;pick out the odd numbers only
;;        xday = xday(i)
;;        ycen = ycen(i)
     ;;pl = plot(tx, ty, title = aorlist(n), xtitle = 'time(hrs)', yrange = [mean(planethash[aorlist(n)].ycen,/nan) -0.5, mean(planethash[aorlist(n)
       ;;b = plot(1/wk1, wk2, xtitle = 'Frequency(minutes)', ytitle = 'Power', xrange =[0,100],$
        ;;         thick = 2, color = 'red',name = 'Y centroids', title = 'Y centroids', yrange = [0, max(wk2[10:100])])
        ;;print, 'LNP result', result
        ;;test = where(1/wk1 gt 40 and 1/wk1 lt 50,ntest)
        ;;if ntest gt 0 then print, 'values around 45', wk2[test]
        ;;want to know if there is a real peak between 30 - 50 minutes
        ;;try monte carloing the noise periodogram
        ;;randomize the order of the photometry while keeping the time
        ;;information the same.
        ;;nwk2 = n_elements(wk2)  ;which is bin_flux cut down
        ;;ny2= n_elements(ycen)
        ;;make a randomly ordered array with nel elements
        ;;nmc = 3
        ;;final_max = fltarr(nmc)

        ;;for rcount = 0, nmc - 1 do begin
        ;;   rand = randomu(seed, ny2)
           ;;now use that as the pointer/orderer for the photometry
        ;;   randflux = ycen[sort(rand)]
        ;;   rr = LNP_TEST(xday, randflux, /double,  WK1 = rrwk1, WK2 = rrwk2, JMAX = jmax)
           ;;print, 'random result', rr
        ;;   final_max(rcount) = rr(0)
           ;;randp = plot(wk1, wk2,overplot=b, color = 'gray')
        ;;endfor
;;print, 'max finalmax', max(final_max)  ;; this is the same as
;;max(rr(0))
 ;realpk = where(peakheight gt
                                ;5.*max(final_max)  and peakperiod gt
                                ;30. and peakperiod le 50.,nrealpk)
 ;;plothist, peakheight, xhist, yhist, bin = 0.1, /noplot
        ;;pb = barplot(xhist, yhist, xtitle = 'peakheight', xrange = [0, 10])
        ;;print, 'mean ,sig', mn , sig
        ;;print, 'nrealpk', nrealpk
        ;;print, 'periods'
        ;;print,peakperiod(realpk)
        ;;print, 'strengths'
        ;;print,(peakheight(realpk) -mn) /sig
        ;;want to keep the period and strength of these peaks
