pro rms_binsize_v2
 colorarr = ['gray', 'gray','gray','gray','gray','burlywood','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia']
 
  planetinfo = create_planetinfo()
  aorname= planetinfo['WASP-14b', 'aorname_ch2'] 

  restore, '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/WASP-14b_phot_ch2_2.25000_150226_bcdsdcorr.sav'
  readcol, '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/wong/ch2.dat', model_phase, model_flux

  for a =0, 5 do begin; n_elements(aorname) -1 do begin
     corrflux =planethash[aorname(a),'corrflux_d']
     corrfluxerr =planethash[aorname(a),'corrfluxerr']
     phasearr = planethash[aorname(a),'phase']

     ;;get rid of NANs
     b = where(finite(corrflux) gt 0, good)

     ;;only consider those AORs with more than 20% of their data in
     ;;the sweet spot.
     if good gt 0.2*n_elements(corrflux) then begin
        print, 'a', a, good, 0.2*n_elements(corrflux), ' ', colorarr(a), mean(corrflux, /nan)

        corrflux = corrflux(b)
        corrfluxerr = corrfluxerr(b)
        phasearr = phasearr(b)

        ;; make a model light curve with the same *unbinned* phase as the
        ;; snapshot data
        modelarr = fltarr(n_elements(corrflux))
        for im = 0, n_elements(corrflux) -1 do begin
           index = closest(model_phase, phasearr(im))
           modelarr[im] = model_flux(index)
        endfor

        ;;subtract model from the snapshot corrfluxes
        delta = corrflux - modelarr

        bin_scale = 189         ; largest bin_scale = 1/3 of the number of data points
        rmsarr = fltarr(bin_scale + 1 )
        numberarr = findgen(n_elements(corrflux))
     
        for n = 1, bin_scale  do begin

           ;;bin up the snapshot data
           h = histogram(numberarr, OMIN=om, binsize = n, reverse_indices = ri)
           ;;setup some arrays
           bin_delta = dblarr(n_elements(h))
           
           if n eq 1 then begin  ; no binning
              rmsarr(n) = stddev(delta)                
           endif else begin  ; for the binning scales
              c = 0
              for j = 0L, n_elements(h) - 1 do begin
                 
                 if (ri[j+1] gt ri[j] + 1)  then begin ;require 2 elements in the bin
                    ;;what is the mean flux in each bin
                    meanclip, delta[ri[ri[j]:ri[j+1]-1]], meanflux, sigmaflux
                    bin_delta[c] = meanflux
                    c = c + 1
                 endif
              endfor
              bin_delta = bin_delta[0:c-1]
              rmsarr(n) =  stddev(bin_delta)

           endelse

 
        endfor
 

        p2 = plot(findgen(n_elements(rmsarr)), rmsarr, /xlog, /ylog, overplot = p2, $   ;/ rmsarr(1)
                  color = colorarr(a),  xrange =[1, 200], xtitle = 'Number of frames',$  ;yrange = [1E-2,1],
                  axis_style = 1, ytitle = 'Normalized Residual RMS')
        xaxis = axis('x', location = [0,max(p2.yrange)], coord_transform = [0, .0337],target = p2, $
        textpos = 1, title = 'Binning Scale (Min.)')
        xaxis = axis('y', location = [max(p2.xrange),0], target = p2, tickdir = 0, textpos = 0, showtext = 0)

     endif  ; on the sweet spot

  endfor ; for all AORs
     
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
  
  p2 = plot( bin_scale, y, thick = 3,  overplot = p2)  ;/y(0)
;now add root2 * root N
  y2 =  1.5 * y/y(0)
;  p2 = plot( bin_scale, y2, thick = 3, linestyle = 2, color = 'grey', overplot = p2)
  


;------------------------
;add WOng et al. 4.5 micron plot
  xw = [1.018, 1.414, 8.061, 23.332, 56.334, 69.524, 84.906, 106.094, 131.287, 160.696, 198.176, 383.034, 473.986, 587.036, 737.974, 923.027, 1.161e3, 1.455e3, 1.846e3]

  yw = [3.866e-3, 3.280e-3, 1.378e-3, 8.197e-4, 5.484e-4, 4.971e-4, 4.452e-4, 4.047e-4, 3.668e-4, 3.300e-4, 2.996e-4, 2.278e-4, 2.091e-4, 1.920e-4, 1.799e-4, 1.684e-4, 1.541e-4, 1.433e-4, 1.341e-4]

  p2 = plot(xw, yw, thick = 3, linestyle = 2, overplot = p2) ;/ yw(0)

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
