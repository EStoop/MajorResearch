pro phot_exoplanet_ML, planetname, apradius,chname, columntrack = columntrack, breatheap = breatheap, hybrid = hybrid, rawfile = rawfile, simulated = simulated
;do photometry on any IRAC staring mode exoplanet data
;now with hashes of hashes!

  t1 = systime(1)
  journal = '/Users/jkrick/external/irac_warm/HD75289/journal.txt'
;convert aperture radius in pixels into what get_centroids_for_calstar_jk uses 
case apradius of
;[ 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25]
   1.5: apval = 0
   2.0: begin
      apval = 0
      if chname eq '2' then pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_rmulti_s3_7_0p1s_x4_150226.sav'
;      if chname eq '2' then pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_r2p00_s3_7_0p1s_x4_140314.sav'
      if chname eq '1' then pmapfile =  '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch1_r2p00_s3_7_0p4s_140314.sav'
   end
   2.25: begin
      apval = 1
      if chname eq '2' then begin
         pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_rmulti_s3_15_0p1s_x4_160126.sav'
         ;;pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_rmulti_s3_7_0p1s_x4_150723.sav'
         
         ;restore, '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_indices_150309.sav'
      endif

      if chname eq '1' then pmapfile =  '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch1_rmulti_s3_7_sd_0p4s_sdark_150722.sav'
;pmap_data_ch1_r2p25_s3_7_0p4s_140716.sav'
   end
   2.5: begin
      apval = 2
      if chname eq '2' then pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_rmulti_s3_7_0p1s_x4_sdark_150223.sav'
;      if chname eq '2' then pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_r2p50_s3_7_sd_0p1s_x4_sdark_150127.sav'
      if chname eq '1' then pmapfile =  '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch1_r2p50_s3_7_0p4s_140314.sav'
   end
   2.75: apval = 2
   3.25: apval = 2
   3: apval = 3

   Else: apval = 2              ; if no decision, then choose an apradius = 2.5 pixels
endcase

if keyword_set(simulated) then pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_sim_data_ch2_rmulti_s3_7_0p1sx4_151030.sav'
print, 'using pmap file', pmapfile, apval

;run code to read in all the planet parameters
planetinfo = create_planetinfo()
;ra_ref = planetinfo[planetname, 'ra']
;dec_ref = planetinfo[planetname, 'dec']
if chname eq '2'  then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
basedir = planetinfo[planetname, 'basedir']
;utmjd_center = planetinfo[planetname, 'utmjd_center']
;period = planetinfo[planetname, 'period']
intended_phase = planetinfo[planetname, 'intended_phase']
;ra_ref = planetinfo[planetname, 'ra']
;dec_ref = planetinfo[planetname, 'dec']
exosystem = strmid(planetname, 0, 8 )+ ' b' ;'HD 209458 b' ;

exoplanet_data_file = '/Users/jkrick/idlbin/exoplanets.csv'
;exosystem = planetname
if planetname eq 'WASP-13b' then exosystem = 'WASP-13 b'
if planetname eq 'WASP-14b' then exosystem = 'WASP-14 b'
if planetname eq 'WASP-15b' then exosystem = 'WASP-15 b'
if planetname eq 'WASP-16b' then exosystem = 'WASP-16 b'
if planetname eq 'WASP-38b' then exosystem = 'WASP-38 b'
if planetname eq 'WASP-62b' then exosystem = 'WASP-62 b'
if planetname eq 'WASP-52b' then exosystem = 'WASP-52 b'
if planetname eq 'WASP-76b' then exosystem = 'WASP-76 b'
if planetname eq 'HAT-P-22' then exosystem = 'HAT-P-22 b'
if planetname eq 'GJ1214' then exosystem = 'GJ 1214 b'
if planetname eq '55Cnc' then exosystem = '55 Cnc e'
if planetname eq 'HD209458' then exosystem = 'HD 209458 b'
if planetname eq 'Kepler-5' then exosystem = 'Kepler-5 b'
if planetname eq 'Kepler-17' then exosystem = 'Kepler-17 b'
if planetname eq 'upsandb' then exosystem = 'upsilon And b'
if planetname eq 'XO3' then exosystem = 'XO-3 b' 
if planetname eq 'simul_XO3' then exosystem = 'XO-3 b' 
if planetname eq 'tres2' then exosystem = 'TrES-2 b'
if planetname eq 'WASP-63b' then exosystem = 'WASP-63 b'
if planetname eq 'KOI069' then exosystem = 'Kepler-93 b'

;print, exosystem, 'exosystem'
if planetname eq 'WASP-52b' then teq_p = 1315
if planetname eq 'HD 7924 b' then begin
   inclination = 85.
   rp_rstar = 0.001
endif
if planetname eq 'upsandb' then begin
   inclination = 90.
   rp_rstar = 0.001  ;made up
   mjd_transit = 0.
endif

if chname eq '2' then lambdaname  = '4.5'
if chname eq '1' then lambdaname  = '3.6'
get_exoplanet_data,EXOSYSTEM=exosystem,MSINI=msini,MSTAR=mstar,TRANSIT_DEPTH=transit_depth,RP_RSTAR=rp_rstar,AR_SEMIMAJ=ar_semimaj,$
                       TEQ_P=1315,TEFF_STAR=teff_star,SECONDARY_DEPTH=secondary_depth,SECONDARY_LAMBDA=lambdaname,$
                       INCLINATION=inclination,MJD_TRANSIT=mjd_transit,P_ORBIT=p_orbit,EXODATA=exodata,RA=ra,DEC=dec,VMAG=vmag,$
                       DISTANCE=distance,ECC=ecc,T14=t14,F36=f36,F45=f45,FP_FSTAR0=fp_fstar0,/verbose
;ra = 1000.

  print, 'parameters from create_planetinfo'
   ra_ref = planetinfo[planetname, 'ra']
   dec_ref = planetinfo[planetname, 'dec']
   utmjd_center = planetinfo[planetname, 'utmjd_center']
   period = planetinfo[planetname, 'period']
   ra = ra_ref
print, 'ra', ra_ref, 'dec', dec_ref
if ra lt 400. then begin  ; that means get_exoplanet_data actually found the target
  ; ra_ref = double(ra)*15.D       ; comes in hours!;
  ; help, ra_ref
  ; help, double(ra)
 ;  print, ra, ra_ref, double(ra) * 15.D, double(ra) * 15.D/ 15.D
 ;  ra_ref = ra
 ;  dec_ref = dec
;   utmjd_center = mjd_transit
;   period = p_orbit
endif else begin
   print, 'parameters from create_planetinfo'
   ra_ref = planetinfo[planetname, 'ra']
   dec_ref = planetinfo[planetname, 'dec']
   utmjd_center = planetinfo[planetname, 'utmjd_center']
   period = planetinfo[planetname, 'period']
endelse

if planetname eq 'upsandb' then begin
   ;proper motions must be the problem, this is the easy fix
    ra_ref = 24.198512
   dec_ref = 41.404099
endif

if planetname eq '55Cnc' then begin
   ra_ref = 133.14771
   dec_ref = 28.33032
endif

;---------------

;print, 'ut_mjd',utmjd_center
dirname = strcompress(basedir + planetname +'/')
planethash = hash()


if chname eq '2' then occ_filename =  '/Users/jkrick/irac_warm/pmap/pmap_fits/pmap_ch2_0p1s_x4_500x500_0043_120827_occthresh.fits'$
else occ_filename = '/Users/jkrick/irac_warm/pmap/pmap_fits/pmap_ch1_500x500_0043_120828_occthresh.fits'
fits_read,occ_filename, occdata, occheader
startaor =1
stopaor =   n_elements(aorname) - 1
for a =startaor,  stopaor do begin
   print, '----------------------------'
   print, 'working on ',a, ' ', aorname(a)
   dir = dirname+ string(aorname(a) ) 
   CD, dir                      ; change directories to the correct AOR directory
   pwd
   if( n_elements(aorname) eq 3) and (a eq 2) then chname = '2'  ;;xxx only for 2mass pre_dither observations
   print, 'chname in use', chname
;   command  = strcompress( 'find ch'+chname+"/bcd -name 'sdcorrSPITZER*_bcd.fits' > "+dirname+'bcdlist.txt')
   command  = strcompress( 'find ch'+chname+"/bcd -name 'SPITZER*_bcd.fits' > "+dirname+'bcdlist.txt')
;   print, 'command', command
   spawn, command
   command2 =  strcompress('find ch'+chname+"/bcd -name '*bunc.fits' > "+dirname + 'bunclist.txt')
   spawn, command2
   command3 =  strcompress('find ch'+chname+"/raw -name '*dce.fits' > "+dirname + 'dcelist.txt')
   spawn, command3
 
   readcol,strcompress(dirname +'bcdlist.txt'),fitsname, format = 'A', /silent
   readcol,strcompress(dirname+'bunclist.txt'),buncname, format = 'A', /silent
   readcol,strcompress(dirname+'dcelist.txt'),rawname, format = 'A', /silent

   print,'n_elements(fitsname)', n_elements(fitsname)
;     aparr = dblarr(n_elements(fitsname))  ;keep the aperture sizes used
   
   startfits = 0L


   for i =startfits,  n_elements(fitsname) - 1  do begin ;read each cbcd file, find centroid, keep track
;       print, 'working on ', fitsname(i)         
      header = headfits(fitsname(i)) ;
      sclk_obs= sxpar(header, 'SCLK_OBS')
      frametime = sxpar(header, 'FRAMTIME')
      bmjd_obs = sxpar(header, 'BMJD_OBS')
       ;utcs_obs = sxpar(header, 'UTCS_OBS')
      ch = sxpar(header, 'CHNLNUM')
      ronoise = sxpar(header, 'RONOISE') ; these are zeros
      gain = sxpar(header, 'GAIN')
      fluxconv = sxpar(header, 'FLUXCONV')
      exptime = sxpar(header, 'EXPTIME')
      aintbeg = sxpar(header, 'AINTBEG')
      atimeend = sxpar(header, 'ATIMEEND')
      ;;print, 'aintbeg, atimeend', aintbeg, atimeend
      naxis = sxpar(header, 'NAXIS')
      framedly = sxpar(header, 'FRAMEDLY')
      if ra_ref gt 359. then begin
        print, 'inside ra_ref gt 359'
        ra_ref = sxpar(header, 'RA_REF')
        dec_ref = sxpar(header, 'DEC_REF')
     endif

      if ch eq '2' and frametime eq 2 then ronoise = 12.1
      if i eq startfits then begin
         print, 'ronoise', ronoise, gain, fluxconv, exptime, ra_ref, dec_ref, naxis
         sclk_0 = sclk_obs
      endif



      if i eq startfits and naxis eq 3 then begin
         xarr = fltarr(63*(n_elements(fitsname)))
         yarr = xarr
         xerr = xarr
         yerr = xarr
         fluxarr = xarr
         fluxerrarr = xarr
         corrfluxarr = xarr
         corrfluxerrarr = xarr
         timearr =  dblarr(63*(n_elements(fitsname)))
         bmjd = dblarr(63*(n_elements(fitsname)))
         backarr = xarr
         backerrarr = xarr
         nparr = xarr
         npcentroidsarr = xarr
         xfwhmarr = xarr
         yfwhmarr = xarr
         peakpixDNarr = xarr
         piarr = fltarr(7,7,63L*n_elements(fitsname))
      endif
      if i eq startfits and naxis ne 3 then begin
         xarr = fltarr(n_elements(fitsname))
         yarr = xarr
         xerr = xarr
         yerr = xarr
         fluxarr = xarr
         fluxerrarr = xarr
         corrfluxarr = xarr
         corrfluxerrarr = xarr
         timearr =dblarr(n_elements(fitsname))
         bmjd = dblarr(n_elements(fitsname))
         backarr = xarr
         backerrarr = xarr
         nparr = xarr
         npcentroidsarr = xarr
         xfwhmarr = xarr
         yfwhmarr = xarr
         peakpixDNarr = xarr
         piarr =fltarr(7,7,n_elements(fitsname))
      endif
      fdarr = fltarr(n_elements(fitsname))
      fdarr[i] = framedly

      if naxis eq 3 then begin
         deltatime = (atimeend - aintbeg) / 64.D ; real time for each of the 64 frames
         nt = dindgen(64)
         sclkarr = sclk_obs  + (deltatime*nt);/60./60./24.D ; 0.5*frametime + frametime*nt
         bmjdarr= bmjd_obs + (deltatime*nt)/60./60./24.D   ; 0.5*(frametime/60./60./24.) + (frametime/60./60./24.)*nt
      endif else begin          ; full array, so don't need to expand out the times
         sclkarr = sclk_obs
         bmjdarr = bmjd_obs
      endelse
        ;; print, 'sclk split up', sclkarr(0), ' ' ,sclkarr(1), ' ' ,sclkarr(60),  format = '(A, D0, A, D0, A, D0)'

      ;read in the files
      fits_read, fitsname(i), im, h
      fits_read, buncname(i), unc, hunc
      test = sxpar(h, 'RAWFILE')
      
    ; fits_read, covname, covdata, covheader

      ;apply mask file if necessary
      if planetname eq 'hd189733' then  im[13:16, 4:7, *] = !Values.F_NAN ;mask region with nan set for bad regions
      if planetname eq 'WASP-14b' then  begin
         ra_off = 218.27655; 218.27718
         dec_off = 21.897652; 21.89808
         adxy, h, ra_off, dec_off, x_off, y_off
         x_off = round(x_off)
         y_off = round(y_off)
         im[(x_off-3):(x_off + 3), (y_off - 3):(y_off+3), *] = !Values.F_NAN ;mask region with nan set for bad regions
;         one = strmid(fitsname(i), 0,8)
;         three = strmid(fitsname(i), 8)
;         newname = strcompress( one + 'test' + three, /remove_all)
;         fits_write, newname, im, h
      endif

      if planetname eq 'HAT-P-22' and chname eq '2' then  im[19:25, 9:15, *] = !Values.F_NAN  ;mask region with nan set for bad regions
      if planetname eq 'HAT-P-22' and chname eq '1' then  im[5:11, 11:17, *] = !Values.F_NAN  ;mask region with nan set for bad regions
      if planetname eq 'HD93385' then im[17:21, 21:27, *] = !Values.F_NAN ;mask region with nan set for bad regions                      
      
      ;;calculate a background manually
 ;     bkgd = calc_bkgd(im, h, ra_ref, dec_ref)

      ;run the centroiding and photometry
      get_centroids_for_calstar_jk,im, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                   x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                   x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                   xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                   xfwhm, yfwhm,  /WARM
      nanfound = where(FINITE(np) lt 1, nancount)
      if nancount gt 0 then print, 'NAN: ', fitsname(i), nanfound
      
      x_center = temporary(x3)
      y_center = temporary(y3)
      xcenerr = temporary(x3s)
      ycenerr = temporary(y3s)
     ;choose the requested pixel aperture
      abcdflux = f[*,apval]      
      fs = fs[*,apval]
     ; 3-7 pixel background
      back = b[*,0]
      backerr = bs[*,0]
      npcentroids = np   ;keep track of np from get_centroids
;      print, 'testing noise pixel', npcentroids
      ;calculate noise pixel
;      print, mean(x_center), mean(y_center), naxis
;      np = noisepix(im, x_center, y_center, ronoise, gain, exptime, fluxconv, naxis, 3.5, 6.5, 12.5)
      np = findgen(64)
      ; if changing the apertures then use this to calculate photometry
;;      if keyword_set(breatheap) then begin
;;         abcdflux = betap(im, x_center, y_center, ronoise, gain, exptime, fluxconv,np, chname)
         ;XXXfake these for now
;;         fs = abcdflux
;;         back = abcdflux
;;         backerr = abcdflux
;;      endif 
      
      if keyword_set(rawfile) then begin
;read in the raw data file and get the DN level of the peakpixel      
         fits_read, rawname(i), rawdata, rawheader
        
         barrel = fxpar(rawheader, 'A0741D00')
         fowlnum = fxpar(rawheader, 'A0614D00')
         pedsig = fxpar(rawheader, 'A0742D00')
         ichan = fxpar (rawheader, 'CHNLNUM')
              
         ;;use Bill's code to conver to DN
         dewrap2, rawdata, ichan, barrel, fowlnum, pedsig, 0, rawdata
     
         ;;or use Jim's code to convert to DN
         ;;rawdata = irac_raw2dn(rawdata,ichan,barrel,fowlnum)
         rawdata = reform(rawdata, 32, 32, 64)
         peakpixDN = abcdflux   ; set up the array to be the same size as abcdflux
         if naxis eq 3 then begin
            for pp = 0, 63 do begin
               peakpixDN[pp] = max(rawdata[13:16,13:16,pp])
            endfor
         endif
      endif

;track the values of the 7x7 pixel box around the centroid
;        pi = track_box(im, x_center, y_center)   ; tb now a 25 x 64 element array
      if naxis eq 3 then pi = im[12:18, 12:18,*] ; now a 7x7x64 element array
      if naxis eq 2 then pi = im[(fix(x_center) - 3):(fix(x_center+3)), (fix(y_center) - 3):(fix(y_center+3))] ; now a 7x7x64 element array
;;      print, 'im[15,15]', im[15,15,*]
;track the value of a column
      if keyword_set(columntrack) then begin 
         centerpixval1 = findgen(64)
         centerpixval2 = findgen(64)
         centerpixval3 = findgen(64)
         centerpixval4 = findgen(64)
         centerpixval5 = findgen(64)
         centerpixval6 = findgen(64)
         sigmapixval1 = findgen(64)
         sigmapixval2 = findgen(64)
         sigmapixval3 = findgen(64)
         sigmapixval4 = findgen(64)
         sigmapixval5 = findgen(64)
         sigmapixval6 = findgen(64)
         for nframes = 0, 64 - 1 do begin
            meanclip, im[13, 12:18,nframes], meancol1, sigmacol
            centerpixval1[nframes] = meancol1
            sigmapixval1[nframes] = sigmacol
            meanclip, im[13, 18:24,nframes], meancol2, sigmacol
            centerpixval2[nframes] = meancol2   
            sigmapixval2[nframes] = sigmacol     
            meanclip, im[14, 12:18,nframes], meancol3, sigmacol
            centerpixval3[nframes] = meancol3           
            sigmapixval3[nframes] = sigmacol
            meanclip, im[14, 18:24,nframes], meancol4, sigmacol
            centerpixval4[nframes] = meancol4
            sigmapixval4[nframes] = sigmacol
            meanclip, im[13, 10:22,nframes], meancol5, sigmacol
            centerpixval5[nframes] = meancol5
            sigmapixval5[nframes] = sigmacol
            meanclip, im[14, 10:22,nframes], meancol6, sigmacol
            centerpixval6[nframes] = meancol6
            sigmapixval6[nframes] = sigmacol
         endfor
      endif
      
      
      
;---------------------------------
;use pmap data to find nearest neighbors in pmap dataset and find a
;correction based on those neighbors.  
      if keyword_set(hybrid) then begin
                                ;use hybrid technique with pmap dataset and nn techniques
;         print, 'starting corrflux', x_center[0], y_center[0],
;         abcdflux[0], ch, xfwhm[0], yfwhm[0], npcentroids[0]
         ;;read in a list of indices for rejecting some of the pmap
         ;;data points
         if naxis eq 3 then begin
            corrflux = pmap_correct(x_center,y_center,abcdflux,ch,xfwhm,yfwhm,NP = npcentroids,$
                                    FUNC=fs,CORR_UNC=corrfluxerr, DATAFILE=pmapfile,NNEAREST=nn, $
                                    R_USE = apradius, USE_PMAP = IMAIN) ;,/VERBOSE
         endif else begin
            corrflux = pmap_correct(x_center,y_center,abcdflux,ch,xfwhm,yfwhm,NP = npcentroids,$
                                    FUNC=fs,CORR_UNC=corrfluxerr, DATAFILE=pmapfile,NNEAREST=nn, $
                                    R_USE = apradius, USE_PMAP = IMAIN,/full) ;,/VERBOSE
         endelse
         

            
;         print, 'corrflux', corrflux[0:10]
                                ;corrflux = pmap_correct(x_center,y_center,abcdflux,ch,npcentroids,occdata, corr_unc = corrfluxerr, func = fs,$
         ;                       datafile =pmapfile,/threshold_occ,/use_np) 
      endif else begin
                                ;make up some junk for corrflux, so I can keep it as a variable in the rest of this
         corrflux = fltarr(n_elements(x_center))
         corrfluxerr = corrflux
               ;correct for pixel phase effect based on pmaps from Jim
                                ;file_suffix =
                                ;['500x500_0043_120828.fits','0p1s_x4_500x500_0043_121120.fits']
;         corrflux = iracpc_pmap_corr(abcdflux,x_center,y_center,ch,/threshold_occ, threshold_val = 20) ;, file_suffix = file_suffix)
;         corrflux = abcdflux
;         corrfluxerr = fs       ;leave out the pmap err for now
      endelse


;---------------------------------

      if naxis eq 3 then begin  ; and i eq 0 then begin
         xarr[i*63] = x_center[1:*]
         yarr[i*63] = y_center[1:*]
         xerr[i*63] = xcenerr[1:*]
         yerr[i*63] = ycenerr[1:*]
         fluxarr[i*63] = abcdflux[1:*]
         fluxerrarr[i*63] = fs[1:*]
         corrfluxarr[i*63] = corrflux[1:*]
         corrfluxerrarr[i*63] = corrfluxerr[1:*]
         timearr[i*63] = sclkarr[1:*]        
         bmjd[i*63] = bmjdarr[1:*]
         backarr[i*63] = back[1:*]
         backerrarr[i*63] = backerr[1:*]
         nparr[i*63] = np[1:*]
         npcentroidsarr[i*63] = npcentroids[1:*]
         xfwhmarr[i*63] = xfwhm[1:*]
         yfwhmarr[i*63] = yfwhm[1:*]
         if keyword_set(rawfile) then peakpixDNarr[i*63] = peakpixDN[1:*]
         ;;print, 'testing pi', pi
         ;;need to loop through each image individually
         for pj = 0, 62 do begin
            ;;print, 'adding to piarr', pi[*,*,pj]
            piarr(*,*,i*63 + pj) = pi[*,*,pj]
         endfor
         
         ;;piarr(i*63L) = pi ;;this line is my problem
         ;;pi
         ;;print, 'piarr', piarr[*,*,i*63L]
        ;; print, 'setting up max piarr, min piarr', max(piarr), min(piarr)

         if keyword_set(columntrack) then begin 
                                ; I think I deleted more parts of this than I may have intended, so if it is not working, that may be why
            centerpixarr1 = centerpixval1[1:*]
            centerpixarr2 = centerpixval2[1:*]
            centerpixarr3 = centerpixval3[1:*]
            centerpixarr4 = centerpixval4[1:*]
            centerpixarr5 = centerpixval5[1:*]
            centerpixarr6 = centerpixval6[1:*]
            sigmapixarr1 = sigmapixval1[1:*]
            sigmapixarr2 = sigmapixval2[1:*]
            sigmapixarr3 = sigmapixval3[1:*]
            sigmapixarr4 = sigmapixval4[1:*]
            sigmapixarr5 = sigmapixval5[1:*]
            sigmapixarr6 = sigmapixval6[1:*]
         endif
      endif 
      if naxis eq 2 then begin; and i eq 0 then begin
         xarr[i] = x_center
         yarr[i]  =  y_center
         yerr[i] = ycenerr
         xerr[i] = xcenerr
         fluxarr[i]  =  abcdflux
         fluxerrarr[i]  =  fs
         corrfluxarr[i]  = corrflux
         corrfluxerrarr[i]  =  corrfluxerr
         timearr[i]  = sclkarr
         bmjd[i]  = bmjdarr
         backarr[i]  =  back
         backerrarr[i]  = backerr
         nparr[i]  = npcentroids
         npcentroidsarr[i] = npcentroids
         xfwhmarr[i] = xfwhm
         yfwhmarr[i] = yfwhm
         if keyword_set(rawfile) then peakpixDNarr[i] = peakpixDN
         piarr[*,*,i] = pi
      endif

      if a eq startaor and i eq startfits then begin
         time_0 = bmjdarr(0)
      endif

   endfor; for each fits file in the AOR


   
   ;---------------------------------

;   print, 'testing bmjd, utmjd', bmjd(0), ' ' ,bmjd(20), ' ' ,bmjd(60), ' ' ,bmjd(100) , format = '(A, D0, A, D0, A, D0, A, D0)'
;get phasing out of the way here
   bmjd_dist = bmjd - utmjd_center ; how many UTC away from the transit center
   phase =( bmjd_dist / period )- fix(bmjd_dist/period)
   
;   low = where(phase lt-0.5 and phase ge -1.0)
;   phase(low) = phase(low) + 1.
   phase = phase + (phase lt -0.5 and phase ge -1.0)
   
;   high = where(phase gt 0.5 and phase le 1.0)
;   phase(high) = phase(high) - 1.0
   phase = phase- (phase gt 0.5 and phase le 1.0)
   
;instead of wrapping the phase at 0.5, this makes the secondary eclipse be in the middle
   if intended_phase gt 0.4 and intended_phase lt 0.6 then begin
      n = where(phase lt 0, ncount)
      if ncount gt 0 then phase(n) = phase(n) + 1
   endif
  

;   if intended_phase gt 0.4 and intended_phase lt 0.6 then begin ;secondary eclipse
;      print, 'secondary eclipse intended'
;      phase = temporary(phase)+0.5
;   endif 
   
;print, 'testing phase', phase[0:200]
;   print, 'end phase', phase[n_elements(phase) - 1]

  ;---------------------------------

   ;; make a correction for flux degradation as a function of time
   ;;since the beginning of this targets AORs.


   if chname eq '2' then factor = 0.0005
   if chname eq '1' then factor = 0.0010
   ;;0.05% per year or 0.0005*(time in
   ;;years since start of observation)
   deltatime = bmjd - time_0 ; in days
   deltatime = deltatime /365. ; in years
   degrade = factor*deltatime
   corrflux_d =  corrfluxarr +(corrfluxarr* degrade)
   
;--------------------------------
;fill in that hash of hases
;--------------------------------


   if keyword_set(columntrack) then begin 
      
      keys =['ra', 'dec', 'xcen', 'ycen', 'flux','fluxerr', 'corrflux', 'corrfluxerr', 'sclktime_0', 'timearr', 'aor', 'bmjdarr',  'bkgd', 'bkgderr','np', 'centerpixarr1','centerpixarr2','centerpixarr3','centerpixarr4','centerpixarr5','centerpixarr6','sigmapixarr1','sigmapixarr2','sigmapixarr3','sigmapixarr4','sigmapixarr5','sigmapixarr6','phase']
      values=list(ra_ref,  dec_ref, xarr, yarr, fluxarr, fluxerrarr, corrfluxarr, corrfluxerrarr, sclk_0, timearr, aorname(a), bmjd,  backarr, backerrarr, nparr, centerpixarr1, centerpixarr2, centerpixarr3, centerpixarr4, centerpixarr5, centerpixarr6, sigmapixarr1, sigmapixarr2, sigmapixarr3, sigmapixarr4, sigmapixarr5, sigmapixarr6, phase)
      planethash[aorname(a)] = HASH(keys, values)
   endif else begin
      keys =['ra', 'dec', 'xcen', 'ycen', 'xcenerr', 'ycenerr', 'flux','fluxerr', 'corrflux', 'corrfluxerr', 'sclktime_0', 'timearr', 'aor', 'bmjdarr', 'bkgd', 'bkgderr','np','phase', 'npcentroids','exptime','xfwhm', 'yfwhm','framedly','pixvals','corrflux_d']
      values=list(ra_ref,  dec_ref, xarr, yarr, xerr, yerr, fluxarr, fluxerrarr, corrfluxarr, corrfluxerrarr, sclk_0, timearr, aorname(a), bmjd,  backarr, backerrarr, nparr, phase, npcentroidsarr, exptime, xfwhmarr, yfwhmarr, fdarr,piarr, corrflux_d)
      planethash[aorname(a)] = HASH(keys, values)
;      keys =['ra', 'dec', 'xcen', 'ycen', 'flux','fluxerr', 'corrflux', 'corrfluxerr', 'sclktime_0', 'timearr', 'aor', 'bmjdarr', 'bkgd', 'bkgderr','np','phase', 'npcentroids','exptime','xfwhm', 'yfwhm','peakpixDN', 'framedly','pixvals','corrflux_d']
;      values=list(ra_ref,  dec_ref, xarr, yarr, fluxarr, fluxerrarr, corrfluxarr, corrfluxerrarr, sclk_0, timearr, aorname(a), bmjd,  backarr, backerrarr, nparr, phase, npcentroidsarr, exptime, xfwhmarr, yfwhmarr, peakpixDNarr, fdarr,piarr, corrflux_d)
      planethash[aorname(a)] = HASH(keys, values)
   endelse


;;make a 2D array to hold variables instead of hash to output to csv file
  pix2_2 = reform(piarr[2, 2, *])
  pix2_3 = reform(piarr[2, 3, *])
  pix2_4 = reform(piarr[2, 4, *])
  pix3_2 = reform(piarr[3, 2, *])
  pix3_3 = reform(piarr[3, 3, *])
  pix3_4 = reform(piarr[3, 4, *])
  pix4_2 = reform(piarr[4, 2, *])
  pix4_3 = reform(piarr[4, 3, *])
  pix4_4 = reform(piarr[4, 4, *])

  help, pix2_2
  help, xarr
  
  all_arr = [[xarr],[xerr], [yarr],[yerr], [fluxarr], [fluxerrarr], [npcentroidsarr], [xfwhmarr], [yfwhmarr],  [bmjd],  [backarr], [backerrarr], [pix2_2] , [pix2_3], [pix2_4], [pix3_2] , [pix3_3], [pix3_4], [pix4_2] , [pix4_3], [pix4_4]]
  help, all_arr

  trans_arr = transpose(all_arr)
  help, trans_arr
  header = ['xpos','xerr','ypos','yerr','flux', 'fluxerr','np', 'xfwhm','yfwhm', 'bmjd',  'bg_flux', 'sigma_bg_flux','pix1','pix2', 'pix3','pix4','pix5','pix6','pix7','pix8','pix9']

  ;;use write_csv
  csvname = strcompress('/Users/jkrick/external/irac_warm/XO3/XO3_'+ aorname(a)+'.csv')
  print, csvname
  write_csv, csvname, trans_arr, header = header



   
endfor                          ;for each AOR


if keyword_set(breatheap) then begin
   savename = strcompress(dirname + planetname +'_phot_ch'+chname+'_varap.sav')
endif else begin
   pmapname = strmid(pmapfile, 9, 6,/reverse_offset)
   savename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'_' + pmapname + '.sav',/remove_all)
endelse

save, planethash, filename=savename
print, 'saving planethash', savename
print, 'time check', systime(1) - t1

  print, 'max piarr, min piarr', max(piarr), min(piarr)
;;  b = barplot(xhist, yhist, title = 'my data')
  a = where(piarr gt 22000, na, ncomplement = nc)
  print, 'number gt 22000',na, nc


 



  
  journal

end


;function to calcluate noise pixel
function noisepix, im, xcen, ycen, ronoise, gain, exptime, fluxconv,naxis
  apradnp = 3;4
  skyradin = 3;10
  skyradout = 6;12
  convfac = gain*exptime/fluxconv

  if naxis gt 2 then begin  ; for subarray
     np = fltarr(64,/NOZERO)
     for npj = 0, 63 do begin
        indim = im[*,*,npj]
        indim = indim*convfac
        ; aper requires a 2d array
        aper, indim, xcen[npj], ycen[npj], topflux, topfluxerr, xb, xbs, 1.0, apradnp,[skyradin,skyradout],/flux,/exact, /silent, /nan, readnoise = ronoise, setskyval = 0
        aper, indim^2, xcen[npj], ycen[npj], bottomflux, bottomfluxerr, xb, xbs, 1.0,apradnp,[skyradin,skyradout],/flux,/exact, /silent, /nan, readnoise = ronoise, setskyval = 0
        
        np[npj] = topflux^2 / bottomflux
     endfor
  endif

  if naxis lt 3 then begin ; for full array
     indim = im*convfac
        
        aper, indim, xcen, ycen, topflux, topfluxerr, xb, xbs, 1.0, apradnp,[skyradin,skyradout],/flux,/exact, /silent, /nan, readnoise = ronoise, setskyval = 0
        aper, indim^2, xcen, ycen, bottomflux, bottomfluxerr, xb, xbs, 1.0,apradnp,[skyradin,skyradout],/flux,/exact, /silent, /nan, readnoise = ronoise, setskyval = 0
        
        np = topflux^2 / bottomflux
  endif

     return, np               
end


function betap, im, xcen, ycen, ronoise, gain, exptime, fluxconv, np, chname
 ;this function does aperture photometry based on an aperture size that is allowed to vary
  ; as a function of noise pixel
  ; ref: Knutson et al. 2012

  backap = [3.,7.] 
  convfac = gain*exptime/fluxconv
  eim = im * convfac

  varap = sqrt(np)  + 0.2
  ;print, 'testing vartap', varap

  ;XXX add some way of keeping track of varap
  ;don't know how to return that

  abcdflux = fltarr(64,/NOZERO)
  badpix = [-9., 9.] * 1.D8
  pxscal1 = [-1.22334117768332D, -1.21641835430637D, -1.22673962032422D, -1.2244968325831D]
  pxscal1 = abs(pxscal1)
  pxscal2 = [1.22328355209902D, 1.21585676679388D, 1.22298117494211D, 1.21904995758086D]
  pscale2 = pxscal1[long(chname) - 1] * pxscal2[long(chname) - 1]
  scale = pscale2 * !DPI * !DPI / (3600.D * 3600.D * 180.D * 180.D) * 1.0D+06

  for s= 0, 63 do begin
     eslice = eim[*,*,s]
     aper, eslice, xcen[s], ycen[s], xf, xfs, xb, xbs, 1.0, varap[s], backap, $
			      badpix, /FLUX, /EXACT, /SILENT, /NAN, $;
			      READNOISE=ronoise
     f = xf/ convfac
     f = f * scale
     abcdflux[s] =f 
;     print, 'varap, abcdflux', varap[s], abcdflux[s]
  endfor


  return, abcdflux
end

         ;slice up image
;           	for s = 0, 63 do begin
;                   eslice = eim[*, *, s];;;

