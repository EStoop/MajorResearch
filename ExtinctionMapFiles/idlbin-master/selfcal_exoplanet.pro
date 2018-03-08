pro selfcal_exoplanet, planetname, bin_level, apradius, chname, binning=binning, sine=sine
  
;example call: selfcal_exoplanet, 'wasp15', 10*63L, /binning

;run code to read in all the input planet parameters
planetinfo = create_planetinfo()
if chname eq '2' then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
basedir = planetinfo[planetname, 'basedir']
intended_phase = planetinfo[planetname, 'intended_phase']
utmjd_center = planetinfo[planetname, 'utmjd_center']
period = planetinfo[planetname, 'period']
dirname = strcompress(basedir + planetname +'/')
savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'.sav',/remove_all)
restore, savefilename

;==========================================
exoplanet_data_file = '/Users/jkrick/idlbin/exoplanets.csv'
exosystem = strmid(planetname, 0, 7) + ' b' ;'HD 209458 b' ;
if planetname eq 'WASP-52b' then teq_p = 1315
if chname eq '2' then lambdaname  = '4.5'
if chname eq '1' then lambdaname  = '3.6'
get_exoplanet_data,EXOSYSTEM=exosystem,MSINI=msini,MSTAR=mstar,RSTAR = rstar, TRANSIT_DEPTH=transit_depth,RP_RSTAR=rp_rstar,$
                   AR_SEMIMAJ=ar_semimaj,$
                   TEQ_P=1315,TEFF_STAR=teff_star,SECONDARY_DEPTH=secondary_depth,SECONDARY_LAMBDA=lambdaname,$
                   INCLINATION=inclination,MJD_TRANSIT=mjd_transit,P_ORBIT=p_orbit,EXODATA=exodata,RA=ra,DEC=dec,VMAG=vmag,$
                   DISTANCE=distance,ECC=ecc,T14=t14,F36=f36,F45=f45,FP_FSTAR0=fp_fstar0,VERBOSE=verbose
ra_ref = ra*15.                 ; comes in hours!
dec_ref = dec
;utmjd_center = mjd_transit
;period = p_orbit
semimaj = ar_semimaj
dstar = 2.*rstar
m_star = mstar
t_dur =  t14  ; in days
;==========================================


for a = 0,  n_elements(aorname) - 1 do begin
   print, 'working on aorname', aorname(a)

 ;for chopping off some initial part of the light curve
   startnum = intarr(n_elements(aorname))   ;empty now ; not using
   
   xarr = (planethash[aorname(a),'xcen'])[startnum(a):*]
   yarr = (planethash[aorname(a),'ycen'])[startnum(a):*]
   timearr = (planethash[aorname(a),'timearr'])[startnum(a):*]
   fluxerr = (planethash[aorname(a),'fluxerr'])[startnum(a):*]
   flux = (planethash[aorname(a),'flux'])[startnum(a):*]
   phase = (planethash[aorname(a),'phase'])[startnum(a):*]
   
  
;==========================================
;get rid of position outliers
   print, ' fluxarr', n_elements(flux), n_elements(planethash[aorname(a),'flux'])

   print, 'positions', mean(xarr) + 3.0*stddev(xarr),  mean(xarr) -3.0*stddev(xarr),  mean(yarr) +3.0*stddev(yarr), mean(yarr) - 3.0*stddev(yarr)
   
   good = where(xarr lt mean(xarr) + 3.0*stddev(xarr) and xarr gt mean(xarr) -3.0*stddev(xarr) and yarr lt mean(yarr) +3.0*stddev(yarr) and yarr gt mean(yarr) - 3.0*stddev(yarr),ngood_pmap, complement=bad) 
   print, 'bad position',n_elements(bad), n_elements(good)
   
   xarr = xarr[good]
   yarr = yarr[good]
   timearr = timearr[good]
   fluxerr = fluxerr[good]
   flux = flux[good]
   phase = phase[good]
   print, 'middle flux', n_elements(flux)
   
  ;try getting rid of flux outliers.
  ;do some running mean with clipping
  start = 0
  print, 'nflux', n_elements(flux)
  for ni = 100, n_elements(flux) -1,100 do begin
     meanclip,flux[start:ni], m, s, subs = subs;,/verbose
    ; print, 'good', subs+start
     ;now keep a list of the good ones
     if ni eq 100 then good_ni = subs+start else good_ni = [good_ni, subs+start]
     start = ni + 1
  endfor

  ;see if it worked
  xarr = xarr[good_ni]
  yarr = yarr[good_ni]
  timearr = timearr[good_ni]
  fluxerr = fluxerr[good_ni]
  flux = flux[good_ni]
  phase = phase[good_ni]
  print, 'nflux after flux clip', n_elements(flux)

;==========================================
 ;try just cutting out the transit and eclipse
  ; later can make this somehow fit back in.
  transit_phase = t_dur / period ; what fraction of the phase is this in transit (or eclipse)
  print, 'test transit phase', transit_phase
  print, 'phase going into flatline', phase(0), n_elements(phase)
  flatline_1 = where(phase gt -0.5 + transit_phase /2. and phase lt 0.D - transit_phase/2., n_flatline_1)
  flatline_2 = where(phase gt 0.D + transit_phase /2. and phase lt 0.5 - transit_phase /2., n_flatline_2)
  ;put all the out of transit and eclipse phases together
 
  if n_flatline_1 gt 0 and n_flatline_2 gt 0 then flatline = [flatline_1, flatline_2]
  if n_flatline_1 gt 0 and n_flatline_2 lt 1 then flatline = flatline_1
  if n_flatline_1 lt 1 and n_flatline_2 gt 0 then flatline = flatline_2

  xarr = xarr[flatline]
  yarr = yarr[flatline]
  timearr = timearr[flatline]
  fluxerr = fluxerr[flatline]
  flux = flux[flatline]
  phase = phase[flatline]
  print, 'nflux after flatline', n_elements(flux)


;put time in hours
  timearr = (timearr - timearr(0))/60./60.
  print, 'timearr', min(timearr), max(timearr)

;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
;now work on the data fitting with a polynomial for comparison:
;can I get rid of pixel phase effect?

    x = timearr
    y =flux
    yerr = fluxerr
    xcenarr = xarr
    ycenarr = yarr
    amx = median(xcenarr)         
    amy = median(ycenarr)        
    adx = xcenarr - amx           
    ady = ycenarr - amy    
    phasearr = phase
    
;    utcsarr=utcs
    nseg = 1

; Initial guesses	
  if keyword_set(sine) then begin  ;use a fitting function with a sine curve or not
     known_period = period*24.; hours 
     pa0 = [median(y), 1.,1.,1.,1.,1.,0.001,known_period*(2*!Pi),0.9]
     func = 'fpa1_xfunc3'
  endif else begin
     pa0 = [median(y), 1.,1.,1.,1.,1.]
     func = 'fpa1_xfunc2'

  endelse

  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(pa0))
   ;limit the range of the sin cuve phase
  if keyword_set(sine) then begin  ;use a fitting function with a sine curve
     parinfo[7].fixed = 1          
     parinfo[6].limited[0] = 1
     parinfo[6].limits[0] = 0.0
  endif


  afargs = {FLUX:y, DX:adx, DY:ady, T:x, ERR:yerr}
  pa = mpfit('fpa1_xfunc2', pa0, FUNCTARGS=afargs, PERROR=spa, BESTNORM=achi, DOF=adof, COVAR = COV, status = status, $
             errmsg = errmsg, parinfo = parinfo,/quiet)
  print, 'status', status
  print, errmsg
  achi = achi / adof
  print, 'reduced chi squared', achi
  if keyword_set(sine) then begin  ;use a fitting function with a sine curve or not
     model_fit = pa[0] * (1+ pa[1] * adx + pa[2] * ady + pa[3] * adx * ady + pa[4] * adx * adx + pa[5] * ady * ady)+ pa[6]*sin(x/pa[7] + pa[8]) 
  endif else begin
     model_fit = pa[0] * (1+ pa[1] * adx + pa[2] * ady + pa[3] * adx * ady + pa[4] * adx * adx + pa[5] * ady * ady)
  endelse

  ;is there some way to overplot the y vs. time plot with the model?
  ;XXX not working
  model_y = pa[0] * (1+  pa[2] * ady + pa[3] * adx * ady  + pa[5] * ady * ady)

  
  ;plot residuals
  sub = y - model_fit + pa[0]  ; add back in the overall level in Jy

 ; f1 = plot(x,   model  , '6b1.',/overplot)
  ;f2 = plot(x, sub+add, '1sr',/overplot , sym_size = 0.4, sym_filled = 1.)

  print, 'nsub', n_elements(sub)

;then put it all back together again
  if nseg eq 1 then begin
     full_xarr = xcenarr
     full_yarr = ycenarr
     full_flux = y
     full_fluxerr = yerr
     full_phase = phasearr
     full_model = model_fit
     full_model_y = model_y
     full_sub = sub
     full_time = x
  endif else begin
     full_xarr = [full_xarr, xcenarr]
     full_yarr = [full_yarr, ycenarr]
     full_flux = [full_flux,y]
     full_fluxerr = [full_fluxerr, yerr]
      full_phase = [full_phase, phasearr]
     full_model = [full_model,model_fit]
     full_model_y = [full_model_y, model_y]
     full_sub = [full_sub,sub]
     full_time = [full_time, x]

  endelse
  print, 'nfull', n_elements(full_xarr)
;endfor  ;for nseg


;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------

 ;binning after fitting
  if keyword_set(binning) then begin
     
     nfits = n_elements(full_xarr) / bin_level
     print, 'nfits',  nfits

     numberarr = findgen(n_elements(full_xarr))

     h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
     print, 'omin', om, 'nh', n_elements(h)

     bin_flux = dblarr(n_elements(h))
     bin_fluxerr = dblarr(n_elements(h))
     bin_timearr = dblarr(n_elements(h))
     bin_phasearr = dblarr(n_elements(h))
     bin_sub = dblarr(n_elements(h))
     bin_model = dblarr(n_elements(h))
     bin_model_y = dblarr(n_elements(h))

     c = 0
     for j = 0L, n_elements(h) - 1 do begin

;get rid of the bins with no values and low numbers, meaning low overlap
        if (ri[j+1] gt ri[j])  then begin
        
           meanclip, full_flux[ri[ri[j]:ri[j+1]-1]], meanflux, sigmaflux
           bin_flux[c] = meanflux ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

           meanclip, full_time[ri[ri[j]:ri[j+1]-1]], meantimearr, sigmatimearr
           bin_timearr[c]=meantimearr

           meanclip, full_phase[ri[ri[j]:ri[j+1]-1]], meanphasearr, sigmaphasearr
           bin_phasearr[c]= meanphasearr

           idataerr = full_fluxerr[ri[ri[j]:ri[j+1]-1]]
           bin_fluxerr[c] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))

           meanclip, full_sub[ri[ri[j]:ri[j+1]-1]], meansub, sigmasub
           bin_sub[c] = meansub

           meanclip, full_model[ri[ri[j]:ri[j+1]-1]], meanmodel, sigmamodel
           bin_model[c] = meanmodel

           meanclip, full_model_y[ri[ri[j]:ri[j+1]-1]], meanmodely, sigmamodely
           bin_model_y[c] = meanmodely

           c = c + 1
        endif
     endfor
  bin_flux = bin_flux[0:c-1]
  bin_fluxerr = bin_fluxerr[0:c-1]
  bin_timearr = bin_timearr[0:c-1]
  bin_phasearr = bin_phasearr[0:c-1]
  bin_sub = bin_sub[0:c-1]
  bin_model = bin_model[0:c-1]
  bin_model_y = bin_model_y[0:c-1]

  print, 'nbin', n_elements(bin_timearr)
;-------------------------------------------------------------------------------------

     x = bin_phasearr

     y = (bin_sub )/mean(bin_sub)
     yerr = bin_fluxerr/ mean(bin_sub)
     print, 'n x, y, yerr', n_elements(x), n_elements(y), n_elements(yerr)
     
  endif                         ;keyword_set binning
  selfcal_timearr = full_time
  selfcal_flux = full_flux
  
  save, bin_phasearr, y, bin_timearr, yerr, filename=strcompress(dirname + 'selfcal'+ aorname(a) +string(apradius)+'.sav',/remove_all)

  print, 'testing phase in selfcal_exoplanet', bin_phasearr
endfor

end


function fpa1_xfunc2, p, FLUX=flux, DX=dx, DY=dy, T=t, ERR=err
; p is parameters, f0, a1..a5
; t is time
    n = n_elements(dx)
    model = dblarr(n) + p[0]

    scale = p[1] * dx + p[2] * dy + p[3] * dx * dy + p[4] * dx * dx + $
            p[5] * dy * dy 
    dt = t - min(t)
    if (n_elements(p) ge 7) then scale2 = 1.0 + p[6] * dt else scale2 = 1.0
    if (n_elements(p) eq 8) then scale2 = scale2 + p[7] * dt * dt
    
    model = model * (1. + scale) * scale2
    model = (flux - model) / err

return, model
end

function fpa1_xfunc3, p, FLUX=flux, DX=dx, DY=dy, T=t, ERR=err
; p is parameters, f0, a1..a5
; t is time
    n = n_elements(dx)
    model = dblarr(n) + p[0]

    scale = p[1] * dx + p[2] * dy + p[3] * dx * dy + p[4] * dx * dx + $
            p[5] * dy * dy 
    dt = t - min(t)
    
    model = model * (1. + scale) 
    
    ;now add in the sin curve 
    model = model + p[6]*sin(t/p[7] + p[8]) 
 

   model = (flux - model) / err

return, model
end
