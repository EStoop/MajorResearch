pro position_outliers, planetname, apradius, chname

;run code to read in all the input planet parameters
  planetinfo = create_planetinfo()
  if chname eq '2' then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
  basedir = planetinfo[planetname, 'basedir']
  
  dirname = strcompress(basedir + planetname +'/')
  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'.sav',/remove_all)
  print, 'restoring ', savefilename
  restore, savefilename
  print, 'aorname', aorname(0)

  exoplanet_data_file = '/Users/jkrick/idlbin/exoplanets.csv'
  exosystem = strmid(planetname, 0, 7) + ' b' ;'HD 209458 b' ;
;  if planetname eq 'WASP-52b' then 
  teq_p = 1315
  if chname eq '2' then lambdaname  = '4.5'
  if chname eq '1' then lambdaname  = '3.6'
  get_exoplanet_data,EXOSYSTEM=exosystem,MSINI=msini,MSTAR=mstar,TRANSIT_DEPTH=transit_depth,RP_RSTAR=rp_rstar,AR_SEMIMAJ=ar_semimaj,$
                     TEQ_P=teq_p,TEFF_STAR=teff_star,SECONDARY_DEPTH=secondary_depth,SECONDARY_LAMBDA=lambdaname,$
                     INCLINATION=inclination,MJD_TRANSIT=mjd_transit,P_ORBIT=p_orbit,EXODATA=exodata,RA=ra,DEC=dec,VMAG=vmag,$
                     DISTANCE=distance,ECC=ecc,T14=t14,F36=f36,F45=f45,FP_FSTAR0=fp_fstar0,VERBOSE=verbose
  ra_ref = ra*15.               ; comes in hours!
  dec_ref = dec

  
  for a = 0, n_elements(aorname) - 1 do begin
     xoutlier = where(planethash[aorname(a),'xcen'] gt 15.25, xcount)
     print, 'xcount', xcount
     for c = 0, xcount - 1 do begin
        
        print,xoutlier(c) / 64.,  ',', (planethash[aorname(a),'xcen'])(xoutlier(c)), ',',  (planethash[aorname(a),'ycen'])(xoutlier(c)), ',',$
               (planethash[aorname(a),'timearr'])(xoutlier(c)), format = '(F10.2,A, F0, A, F0, A, F0)'
     endfor

     dir = dirname+ string(aorname(a) ) 
     CD, dir                    ; change directories to the correct AOR directory
     command  = strcompress( 'find ch'+chname+"/bcd -name '*bcd.fits' > "+dirname+'bcdlist.txt')
     spawn, command
     command2 =  strcompress('find ch'+chname+"/bcd -name '*bunc.fits' > "+dirname + 'bunclist.txt')
     spawn, command2
     readcol,strcompress(dirname +'bcdlist.txt'),fitsname, format = 'A', /silent
     readcol,strcompress(dirname+'bunclist.txt'),buncname, format = 'A', /silent
     
     for j = 0, 70 do begin
        fits_read, fitsname(j), im, h
        fits_read, buncname(j), unc, hunc
        get_centroids_for_calstar_jk,im, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                     x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                     x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                     xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                     xfwhm, yfwhm, /WARM
        for i = 0, 63 do begin
           print, 'frame',j,  i, x3(i), y3(i)
        endfor   ; for each subframe

     endfor                     ; for each fits image
     
  endfor  ; for each AOR

  
end
