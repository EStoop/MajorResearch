pro rms_binsize_v3, planetname
if planetname eq 'WASP-14b' then begin
   colorarr = ['black', 'gray','red','orange','blue','burlywood','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia']
endif else begin
 colorarr = ['black', 'blue', 'red', 'green','burlywood','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia']
endelse

  planetinfo = create_planetinfo()
  aorname= planetinfo[planetname, 'aorname_ch2'] 
  basedir = planetinfo[planetname, 'basedir']
  startaor = planetinfo[planetname, 'stareaor']
  ;;restore, strcompress(basedir + planetname + '/' + planetname + '_phot_ch2_2.25000_150723_newtime.sav',/remove_all)
  restore, strcompress(basedir + planetname + '/' + planetname + '_phot_ch2_2.25000_160126.sav',/remove_all)

  ;;need a 2D array to keep track of all rms arrays, then be able
  ;;toaverage them.  binning scales are always the same, so can just
  ;;keep the last one.
  bin_scale = 200              ; largest bin_scale = 1/3 of the number of data points

  allrms =fltarr(n_elements(aorname), bin_scale )
  
;  readcol, '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/wong/ch2.dat', model_phase, model_flux
  print, 'startaor', startaor
  ;startaor = 36
  for a =startaor , n_elements(aorname) -1 do begin
     corrflux =planethash[aorname(a),'corrflux_d'];'corrflux_d'
     corrfluxerr =planethash[aorname(a),'corrfluxerr'];'corrfluxerr'
     phasearr = planethash[aorname(a),'phase']

     ;;get rid of NANs
     b = where(finite(corrflux) gt 0, good) ;and phasearr lt -0.002 and phasearr gt -0.012

     ;;only consider those AORs with more than 20% of their data in
     ;;the sweet spot.
     if good gt 0.1*n_elements(corrflux) and a ne 36 then begin ;8
        print, 'a', a, good, 0.2*n_elements(corrflux), ' ', colorarr(a), mean(corrflux, /nan), phasearr[10]

        corrflux = corrflux(b)
        corrfluxerr = corrfluxerr(b)
        phasearr = phasearr(b)

        rmsarr = fltarr(bin_scale + 1 )
        numberarr = findgen(n_elements(corrflux))
     
        for n = 1, bin_scale - 1  do begin
           
           ;;bin up the snapshot data
           h = histogram(numberarr, OMIN=om, binsize = n, reverse_indices = ri)
           ;;setup some arrays
           bin_corrflux = dblarr(n_elements(h))
           if n eq 1 then begin  ; no binning
              rmsarr(n) = stddev(corrflux)                
           endif else begin  ; for the binning scales
              c = 0
              for j = 0L, n_elements(h) - 1 do begin
                 
                 if (ri[j+1] gt ri[j] + 1)  then begin ;require 2 elements in the bin
                    ;;what is the mean flux in each bin
                    meanclip, corrflux[ri[ri[j]:ri[j+1]-1]], meanflux, sigmaflux
                    bin_corrflux[c] = meanflux
                    c = c + 1
                 endif
              endfor
              bin_corrflux = bin_corrflux[0:c-1]

;              if a eq 1 and n lt 65 and n mod 3 eq 0 then begin
;                 plothist, bin_corrflux, xhist, yhist, bin = 1E-4, /noprint, /noplot
;                 ph = barplot(xhist, yhist, fill_color = 'sky blue')
;              endif

              rmsarr(n) =  stddev(bin_corrflux)
;              print,'binning numbers', n, n_elements(bin_corrflux), rms(bin_corrflux), stddev(bin_corrflux)
;              print, 'bin_corrflux', bin_corrflux

           endelse
           allrms(a,n) = rmsarr(n)

        endfor
 
        p2 = plot(findgen(n_elements(rmsarr)), rmsarr/ rmsarr(1), /xlog, /ylog, overplot = p2, $   ;
                  color = colorarr(a),  xrange =[1, bin_scale], xtitle = 'Number of frames',yrange = [1E-2,1], $
                   ytitle = 'Normalized Residual RMS') ; axis_style = 1

     endif  ; on the sweet spot

  endfor ; for all AORs

  ;;what is the average or median of the allrms array?
  medianrms = median(allrms, dimension = 1)

  ;;try normalizing
  medianrms = medianrms/medianrms(1)
  p2 = plot(findgen(n_elements(rmsarr)), medianrms, thick =3,color = 'black', overplot = p2)
;------------------------
;now add straight root N
     
  source_mjy = 68.29            ;mJy  from Star-pet
  exptime = 2.
     
;ch2
  gain = 3.71
  pixel_scale = 1.22
  flux_conv = .1469
  xmax = max(p2.xrange)
  bin_scale = findgen(xmax) + 1
  root_n = sqrt(bin_scale)
  source_electrons = mjy_to_electron( source_mjy, pixel_scale, gain, exptime, flux_conv)
  sigma_poisson = sqrt(source_electrons)
  y =  (sigma_poisson / root_n) / source_electrons
  
  p2 = plot( bin_scale, y/y(0), thick = 3,  overplot = p2)  ;
;now add root2 * root N
  y2 =  1.1 * y/y(0)
;  p2 = plot( bin_scale, y2, thick = 3, linestyle = 2, color = 'cyan', overplot = p2)
  


;------------------------
;add WOng et al. 4.5 micron plot
  xw = [1.018, 1.414, 8.061, 23.332, 56.334, 69.524, 84.906, 106.094, 131.287, 160.696, 198.176, 383.034, 473.986, 587.036, 737.974, 923.027, 1.161e3, 1.455e3, 1.846e3]

  yw = [3.866e-3, 3.280e-3, 1.378e-3, 8.197e-4, 5.484e-4, 4.971e-4, 4.452e-4, 4.047e-4, 3.668e-4, 3.300e-4, 2.996e-4, 2.278e-4, 2.091e-4, 1.920e-4, 1.799e-4, 1.684e-4, 1.541e-4, 1.433e-4, 1.341e-4]

  p2 = plot(xw, yw/ yw(0), thick = 3, linestyle = 2, overplot = p2, color = 'silver') ;


  ;;-----------------------------------
  xaxis = axis('x', location = [0,max(p2.yrange)], coord_transform = [0, .0337],target = p2, $
               textpos = 1, title = 'Binning Scale (Min.)')
  xaxis = axis('y', location = [max(p2.xrange),0], target = p2, tickdir = 0, textpos = 0, showtext = 0)
  
end


       
;             mp = where(bin_phasearr(im) gt (model_phase - 5E-5) and bin_phasearr(im) lt (model_phase + 5E-5),num_mp)
;              if num_mp eq 0 then begin
;                 print, 'hmmm no phase matches'
;                 print, 'n corr', n_elements(corrflux)
;                 print, bin_phasearr(im)
;                 ;;ok but I should be able to find the nearest in
;                 ;;phase and use that, the gaps aren't that big;;
;
;                 ;;use function closest - do this for all of these in
;                 ;;                       fact....
;              endif
;              
;              if num_mp eq 1 then bin_modelarr[im] = model_flux(mp)
;              if num_mp gt 1 then begin
;                 bin_modelarr[im] = mean(model_flux(mp))
;              endif
              
             ;;only fit the function for n = 1, t`hen for the other
              ;;binnings, re-scale
;              start = [1E-7,mean(corrflux, /nan)]
;              result= MPFITFUN('linear',bin_timearr,bin_corrflux, bin_corrfluxerr, start,/Quiet) 
