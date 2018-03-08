pro stare_snap

;plot the snaps vs. stare photometry to see if the difference is
;additive or multiplicative


;need to associate stares to snaps
;use closest in phase

  exosystem = 'WASP-14 b'
  planetname = 'WASP-14b'
;run code to read in all the input planet parameters
  planetinfo = create_planetinfo()
  aorname= planetinfo[planetname, 'aorname_ch2'] 
  basedir = planetinfo[planetname, 'basedir']
  period =  planetinfo[planetname, 'period']
  stareaor = planetinfo[planetname, 'stareaor']

  dirname = strcompress(basedir + planetname +'/')                                                            ;+'/hybrid_pmap_nn/')
  savefilename = strcompress(dirname + planetname +'_phot_ch2_2.25000.sav',/remove_all) ;
  restore, savefilename

  startaor =  0;0                 ;  n_elements(aorname) -29
  stopaor =     n_elements(aorname) - 1
  time_0 = (planethash[aorname(startaor),'timearr']) 
  time_0 = time_0(0)
  stareaor = 5

  ;;make one array of all the stare photometry
  for a = 0, stareaor - 1 do begin
     a_phase =  planethash[aorname(a),'phase']
     a_flux = planethash[aorname(a),'corrflux_d']
     g = where(finite(a_flux) gt 0)
     a_phase = a_phase(g)
     a_flux = a_flux(g)
     if a eq 0 then begin
        stare_phase = a_phase
        stare_flux = a_flux
     endif else begin
        stare_phase = [stare_phase, a_phase]
        stare_flux = [stare_flux, a_flux]
     endelse

  endfor

  ;;make an array of all the snap photometry
  for a = stareaor, stopaor do begin
     a_flux = planethash[aorname(a),'corrflux_d']
     a_phase = planethash[aorname(a),'phase']
     g = where(finite(a_flux) gt 0)
     a_flux = a_flux(g)
     a_phase = a_phase(g)

     meanclip, a_flux, meanflux, sigmaflux
     meanclip, a_phase, meanphase, sigmaphase

     if a eq stareaor then begin
        snap_phase = meanphase; a_phase
        snap_flux =meanflux; a_flux
        
     endif else begin
        snap_phase = [snap_phase,  a_phase]
        snap_flux = [snap_flux, a_flux]
     endelse

  endfor

  ;;for each snapshot data point
  ;;find the closest in phase point in the stare dataset
  stare_close = fltarr(n_elements(snap_phase))
  for i = 0, n_elements(snap_phase) -1 do begin
     stare_close[i] = closest(stare_phase, snap_phase(i))
  endfor

  stare_close_flux = stare_flux(stare_close)
  min =.0565;0.97
  max = .0590;1.03
  norm = 1.0 ; median(stare_flux)
  factor =  .00012; 1.002
  p = plot(stare_close_flux/norm, (snap_flux/norm) + factor, '1.',sym_filled =1, sym_size  = 0.1,  $
           xtitle = 'Stare Corrflux', ytitle = 'Snap Corrflux', title = 'additive' ,xrange = [min, max], yrange = [min, max])
  pl = plot([min, max], [min, max], overplot = p, color = 'red')
end
