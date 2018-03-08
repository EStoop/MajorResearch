pro plot_rocky_planets, figure1 = figure1, figure2 = figure2
  ;readcol, '/Users/jkrick/external/irac_warm/senior_review_2016/test.txt',  NAME,MSINI,A,PER,ECC,year,DEPTH,KS,MASS,SEP,R,TEFF, format = '(A, F10, F10, F10, F10, I10, F10, F10, F10, F10, F10, I10 )'

  readcol, '/Users/jkrick/external/irac_warm/senior_review_2016/exoplanetsprg_10_22_15.txt', NAME,PER,FIRSTREF,DEPTH,DEPTHUPPER,DEPTHLOWER,UDEPTH,KS,MASS,SEP,R,DIST,TEFF, Rstar, format = '(A, F10, I10, F10, F10, F10, F10, F10, F10, F10, F10,F10, I10 , F10)'

  ;make sure they are transiting
  good = where(depth gt 0, goodcount)
  print, 'transiting', goodcount
  sep = sep(good)
  R = R(good)
  firstref = firstref(good)
  depth = depth(good)
  name = name(good)
  TEFF = TEFF(good)
  Tstar = Rstar(good)
  PER = PER(good)
  Rstar = Rstar(good)
  Ks = Ks(good)
  
  ;;make sure they have stellar temperatures
  good = where(TEFF gt 0, goodcount)
  print, 'TEFF', goodcount
  sep = sep(good)
  R = R(good)
  firstref = firstref(good)
  depth = depth(good)
  name = name(good)
  TEFF = TEFF(good)
  Tstar = Rstar(good)
  PER=PER(good)
  Rstar = Rstar(good)
  Ks = Ks(good)
  
  ;;make a cut on magnitude
  ;;right now choosing Ks = 11 as a best case scenario for JWST
  ;;NIRSPEC
  print, 'mag', Ks
  maggood = where(KS le 11., goodcount)
  print, 'KS', goodcount
  sep = sep(maggood)
  R = R(maggood)
  firstref = firstref(maggood)
  depth = depth(maggood)
  name = name(maggood)
  TEFF = TEFF(maggood)
  Tstar = Rstar(maggood)
  PER=PER(maggood)
  Rstar = Rstar(maggood)
  Ks = Ks(maggood)
  
  ;;calculate equilibrium temperature
  ;; is it scary that I am using wikipedia for this?
  ;;https://en.wikipedia.org/wiki/Planetary_equilibrium_temperature#Calculation_for_extrasolar_planets
  a = 0                         ; assume zero albedo
  Rstar = Rstar * 695500.       ; convert to km
  D = sep * 1.49460E8            ; convert to Km
  
  Teq = TEFF*(sqrt(Rstar/(2*D)))  

  ;;not a good idea but...
  ;;we don't have the Teff for GJ1132, but we do have Teq
  Teq[3] = 400

  ;;print, 'Teq', Teq

  ;;keep only the non-nan's - why are there nan's?
  good = where(finite(Teq) gt 0, goodcount) 
  print, 'Teq', goodcount
  sep = sep(good)
  R = R(good)
  firstref = firstref(good)
  depth = depth(good)
  name = name(good)
  TEFF = TEFF(good)
  Tstar = Rstar(good)
  Teq = Teq(good)
  PER = PER(good)
  Rstar = Rstar(good)
  Ks = Ks(good)
  
  ;;sort by year of discovery
  k2 = where (strmid(name, 0, 2) eq 'K2', complement = old)
  new = where(firstref gt 2014 , complement = old)
  GJ = where(strmid(name, 0, 6) eq 'GJ1132', complement = dm)
  if keyword_set(figure1) then begin

     ;;scale the bubble sizes by transit depth
     p1 = bubbleplot(Teq(old), R(old), magnitude = depth(old), /shaded, exponent = 0.5, xtitle = 'Equilibrium Temperature (K)', $
                     ytitle ='Planet Radius ($R_{Jup}$)', /xlog,/ylog, yrange = [.03,5], xrange = [100, 4000], axis_style = 1, /nodata,$
                    xthick = 2, ythick = 2, xtickfont_style = 1, ytickfont_style = 1, margin = 0.12)
     
     ;;add a region for liquid water
     ;;making this up a bit
     xf = [270,370,370,270]
     yf = [0.03, 0.03, 5, 5]
     poly = polygon( xf, yf,/data, /fill_background, fill_color = 'light blue',transparency = 2, overplot = p1)
     
     p1 = bubbleplot(Teq(old), R(old), magnitude = depth(old), /shaded, exponent = 0.5,  /xlog,/ylog, $
                     color = 'light salmon', overplot = p1)
     p2 = bubbleplot(Teq(new), R(new),magnitude = depth(new), /shaded, exponent = p1.exponent, color = 'light salmon', /xlog, /ylog,  overplot = p1, max_value= p1.max_value) ; labels = name(new),
     p2 = bubbleplot(Teq(k2), R(k2),magnitude = depth(k2), /shaded, exponent = p1.exponent, color = 'blue', /xlog, /ylog,  overplot = p1, label_position = "top", max_value= p1.max_value) ; labels = name(k2),

     p2 = bubbleplot(Teq(GJ), R(GJ),magnitude = depth(GJ), /shaded, exponent = p1.exponent, color = 'red', /xlog, /ylog,  overplot = p1, label_position = "top", max_value= p1.max_value) ; labels = name(k2),
     
     yaxis = axis('Y', location = 'right', Title = 'Planet Radius ($R_{Earth}$)',coord_transform=[0,10.97], thick = 2, tickfont_style = 1)
     xaxis = axis('X', location = 'top', showtext = 0, thick = 2)
     
     ;;maybe an image of earth at the "right" location
     ;;not working
;;  read_jpeg, '/Users/jkrick/external/irac_warm/senior_review_2016/earth_smaller.jpg',earth
;;  help, earth
     ;;i1 = image(earth,/xlog, /ylog, overplot = 1)
     s = text(270,0.1, '$\oplus$', /data, overplot = p1, font_style = 1)
     s2 = text(205, 0.11, 'Earth', /data, overplot = p1, font_style = 1)
     s3 = text(310, 0.08, 'GJ1132b', /data, overplot = p1, font_style = 1)
     ;;s = symbol(300, 0.1, 's', /data, sym_color = 'yellow', overplot =
     ;;p1)
     
  endif
  
;----------------------------------------------------------------------------------------
;;spitzer
  S_period=[  0.2016,0.2659,1.8723,3.9890,9.6799,0.8932,0.5222,1.3414, 41.]
  S_rad =[  0.3445, 0.3705,0.6045,0.7308,0.730,0.5009,0.4377, 0.5566,0.73]
  
  S_rad = S_rad(sort(S_period))
  S_period = S_period(sort(S_period))
  
;;Mearth
  M_period =[0.2033,0.3078,0.3803,0.5265,0.6611,1.0173,1.3634,1.6438,2.1499,3.2026,4.2920,6.0889,8.2939,9.4466, 15.0]
  M_rad = [0.6825,0.7494,0.7828,0.8423,0.8868,0.9723,1.0391,1.0837,1.1509,1.2549,1.3400,1.4443,1.5480,1.5889, 1.8]
  
;;Harps
  H_period=[0.2016,0.3285,0.4775,0.7653,1.3414,2.2030,3.3627,5.3896,9.7590, 40.]
  H_rad = [0.9760,1.0206,1.0577,1.1023,1.1617,1.2137,1.2620,1.3177,1.3883, 1.6]
  
  if keyword_set(figure2) then begin
     p2 = plot(S_period, S_rad, color = 'deep sky blue', /xlog, thick = 2, xrange = [0.2, 10], yrange = [0.1,1.8],$
               xtitle = 'Period(Days)', ytitle = 'Planet Radius (Rearth)')

     ;;make polygon for habitable zone
     xf = [2.0,5.3,5.3,2.0]
     yf = [0.1, 0.1, 1.8, 1.8]
     poly = polygon( xf, yf,/data, /fill_background, fill_color = 'light blue',transparency = 20, overplot = p2)
     p2 = plot(S_period, S_rad, color = 'deep sky blue', /xlog, thick = 2,overplot = p2)
     
     t2 = text(4.3, 0.74,'Spitzer',/data, target = p2, color = 'deep sky blue', font_style = 1)
     p2 = plot(M_period, M_rad, color = 'brown', /xlog, thick = 2,overplot = p2)
     t2 = text(0.3,0.65,'MEarth',/data, target = p2, color = 'brown', font_style = 1)
     p2 = plot(H_period, H_rad, color = 'chocolate', /xlog, thick = 2,overplot = p2)
     t2 = text(0.3, 1.08,'HARPS',/data, target = p2, color = 'chocolate', font_style = 1)
     
;;add known planets
     x = [1.6];,3.1, 9.3, 0.85]
     y = [1.2];,1.6, 1.6, 1.68]
     short_depth = [0.0028,0.000359,0.0007 ]
     locations = ['top', 'top','top']
     names = ['GJ1132b', 'HD219134b','K2-21b']
     p2 = bubbleplot(x,y,magnitude = short_depth,  exponent = 0.5, color = 'blue', /shaded, max_value = max(depth), $
                      label_position = locations, labels = names, /xlog, label_vertical_alignment = 0., overplot = p2) ;max_value = max(depth(new)), label_position = locations,
;     p2 = bubbleplot(PER,R*10.97,magnitude = depth,  color = 'blue', /shaded, $
;                      labels = name, /xlog, overplot = p2) ;max_value = max(depth(new)), label_position = locations,
     
     
;levels for mars & moon?
     l2 = polyline([0.2, 10.], [0.53,0.53], /data, target = p2, color = 'gray', linestyle =2, thick = 2)
     t2 = text(5.0, 0.43, 'Mars', /data, target = p2, color = 'gray', font_style = 1)
     l2 = polyline([0.2, 10.], [0.27,0.27], /data, target = p2, color = 'gray', linestyle =2, thick = 2)
     t2 = text(5.0, 0.19, 'Moon', /data, target = p2,  color = 'gray', font_style = 1)
     l2 = polyline([0.2, 10.], [0.413,0.413], /data, target = p2, color = 'gray', linestyle =2, thick = 2)
     t2 = text(5.0, 0.31, 'Ganymede', /data, target = p2,  color = 'gray', font_style = 1)

    
  endif
  
end

