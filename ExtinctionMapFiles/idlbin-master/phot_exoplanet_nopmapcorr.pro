pro phot_exoplanet_nopmapcorr, planetname, apradius,chname, columntrack = columntrack, breatheap = breatheap
;do photometry on any IRAC staring mode exoplanet data
;now with hashes of hashes!
;print, 'waiting'
;wait, 3600  ; waiting for another code to finish
;print, 'done waiting'
 t1 = systime(1)

;convert aperture radius in pixels into what get_centroids_for_calstar_jk uses 
case apradius of
;[ 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25]
   1.5: apval = 0
   1.75: begin
      apval = 0
      if chname eq '2' then pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_r1p75_s3_7_0p1s_x4_140314.sav'
      if chname eq '1' then pmapfile =  '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch1_r1p75_s3_7_0p4s_140314.sav'
   end
   2.0: begin
      apval = 1
      if chname eq '2' then pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_r2p00_s3_7_0p1s_x4_140314.sav'
      if chname eq '1' then pmapfile =  '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch1_r2p00_s3_7_0p4s_140314.sav'
   end
   2.25: begin
      apval = 2
      if chname eq '2' then pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_r2p25_s3_7_0p1s_x4_140716.sav'
      if chname eq '1' then pmapfile =  '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch1_r2p25_s3_7_0p4s_140314.sav'
   end
   2.5: begin
      apval = 3
      if chname eq '2' then pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_r2p50_s3_7_0p1s_x4_140314.sav'
      if chname eq '1' then pmapfile =  '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch1_r2p50_s3_7_0p4s_140314.sav'
   end
   2.75: apval = 5
   3.25: apval = 7
   3: apval = 6

   Else: apval = 2              ; if no decision, then choose an apradius = 2.5 pixels
endcase


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

exoplanet_data_file = '/Users/jkrick/idlbin/exoplanets.csv'
exosystem = strmid(planetname, 0, 8 )+ ' b' ;'HD 209458 b' ;
;exosystem = planetname
if planetname eq 'WASP-13b' then exosystem = 'WASP-13 b'
if planetname eq 'WASP-15b' then exosystem = 'WASP-15 b'
if planetname eq 'WASP-14b' then exosystem = 'WASP-14 b'
if planetname eq 'WASP-16b' then exosystem = 'WASP-16 b'
if planetname eq 'WASP-38b' then exosystem = 'WASP-38 b'
if planetname eq 'WASP-62b' then exosystem = 'WASP-62 b'
if planetname eq 'WASP-52b' then exosystem = 'WASP-52 b'
if planetname eq 'HAT-P-22' then exosystem = 'HAT-P-22 b'
if planetname eq 'HAT-P-8' then exosystem = 'HAT-P-8 b'
if planetname eq 'HD209458' then exosystem = 'HD209458 b'


print, exosystem, 'exosystem'
if planetname eq 'WASP-52b' then teq_p = 1315
if planetname eq 'HD 7924 b' then begin
   inclination = 85.
   rp_rstar = 0.001
endif

if chname eq '2' then lambdaname  = '4.5'
if chname eq '1' then lambdaname  = '3.6'
get_exoplanet_data,EXOSYSTEM=exosystem,MSINI=msini,MSTAR=mstar,TRANSIT_DEPTH=transit_depth,RP_RSTAR=rp_rstar,AR_SEMIMAJ=ar_semimaj,$
                       TEQ_P=1315,TEFF_STAR=teff_star,SECONDARY_DEPTH=secondary_depth,SECONDARY_LAMBDA=lambdaname,$
                       INCLINATION=inclination,MJD_TRANSIT=mjd_transit,P_ORBIT=p_orbit,EXODATA=exodata,RA=ra,DEC=dec,VMAG=vmag,$
                       DISTANCE=distance,ECC=ecc,T14=t14,F36=f36,F45=f45,FP_FSTAR0=fp_fstar0,/verbose
if ra lt 400 then begin  ; that means get_exoplanet_data actually found the target
   ra_ref = ra*15.       ; comes in hours!;
   dec_ref = dec
   utmjd_center = mjd_transit
   period = p_orbit
endif else begin
   ;warning these could be junk as well
   ra_ref = planetinfo[planetname, 'ra']
   dec_ref = planetinfo[planetname, 'dec']
   utmjd_center = planetinfo[planetname, 'utmjd_center']
   period = planetinfo[planetname, 'period']
endelse

;---------------

print, 'ra, dec', ra_ref, dec_ref
dirname = strcompress(basedir + planetname +'/')
planethash = hash()


if chname eq '2' then occ_filename =  '/Users/jkrick/irac_warm/pmap/pmap_fits/pmap_ch2_0p1s_x4_500x500_0043_120827_occthresh.fits'$
                                      else occ_filename = '/Users/jkrick/irac_warm/pmap/pmap_fits/pmap_ch1_500x500_0043_120828_occthresh.fits'
fits_read,occ_filename, occdata, occheader
startaor =0
stopaor =  n_elements(aorname) - 1
for a =startaor, stopaor do begin
   print, 'working on ',aorname(a)
   dir = dirname+ string(aorname(a) ) 
   CD, dir                      ; change directories to the correct AOR directory
   command  = strcompress( 'find ch'+chname+"/bcd -name 'SPITZER*bcd.fits' > "+dirname+'bcdlist.txt')
   print, 'command', command
   spawn, command
   command2 =  strcompress('find ch'+chname+"/bcd -name '*bunc.fits' > "+dirname + 'bunclist.txt')
   spawn, command2
   command3 =  strcompress('find ch'+chname+"/raw -name '*dce.fits' > "+dirname + 'dcelist.txt')
   spawn, command3
 
   readcol,strcompress(dirname +'bcdlist.txt'),fitsname, format = 'A', /silent
   readcol,strcompress(dirname+'bunclist.txt'),buncname, format = 'A', /silent
   readcol,strcompress(dirname+'dcelist.txt'),rawname, format = 'A', /silent

   print,'n_elements(fitsname)', n_elements(fitsname), n_elements(buncname), n_elements(rawname)
;     aparr = dblarr(n_elements(fitsname))  ;keep the aperture sizes used
   
   startfits = 0.D


   for i =startfits,  n_elements(fitsname) - 1  do begin ;read each cbcd file, find centroid, keep track
 ;      print, 'working on ', fitsname(i)         
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
      naxis = sxpar(header, 'NAXIS')
      ;XXXX remove this later, for now I only want subarray
      if naxis eq 2 then goto, SKIPAOR

      if ch eq '2' and frametime eq 2 then ronoise = 12.1
      if i eq startfits then print, 'ronoise', ronoise, gain, fluxconv, exptime
      if i eq startfits then sclk_0 = sclk_obs

      if i eq startfits and naxis eq 3 then begin
         xarr = fltarr(63*(n_elements(fitsname)))
         yarr = xarr
         fluxarr = xarr
         fluxerrarr = xarr
         timearr =  dblarr(63*(n_elements(fitsname)))
         bmjd = dblarr(63*(n_elements(fitsname)))
         backarr = xarr
         backerrarr = xarr
         nparr = xarr
         npcentroidsarr = xarr
         xfwhmarr = xarr
         yfwhmarr = xarr
         peakpixDNarr = xarr
      endif
      if i eq startfits and naxis ne 3 then begin
         xarr = fltarr(n_elements(fitsname))
         yarr = xarr
         fluxarr = xarr
         fluxerrarr = xarr
          timearr =dblarr(n_elements(fitsname))
         bmjd = dblarr(n_elements(fitsname))
         backarr = xarr
         backerrarr = xarr
         nparr = xarr
         npcentroidsarr = xarr
         xfwhmarr = xarr
         yfwhmarr = xarr
         peakpixDNarr = xarr
      endif


      if naxis eq 3 then begin
         deltatime = (atimeend - aintbeg) / 64.D ; real time for each of the 64 frames
         nt = dindgen(64)
         sclkarr = sclk_obs  + (deltatime*nt)/60./60./24.D ; 0.5*frametime + frametime*nt
         bmjdarr= bmjd_obs + (deltatime*nt)/60./60./24.D   ; 0.5*(frametime/60./60./24.) + (frametime/60./60./24.)*nt
      endif else begin          ; full array, so don't need to expand out the times
         sclkarr = sclk_obs
         bmjdarr = bmjd_obs
      endelse
;         print, 'bmjdarr split up', bmjdarr(0), ' ' ,bmjdarr(20), ' ' ,bmjdarr(60),  format = '(A, D0, A, D0, A, D0)'

      ;read in the files
      fits_read, fitsname(i), im, h
      fits_read, buncname(i), unc, hunc
      test = sxpar(h, 'RAWFILE')
      
    ; fits_read, covname, covdata, covheader

      ;apply mask file if necessary
      if planetname eq 'hd189733' then  im[13:16, 4:7, *] = !Values.F_NAN ;mask region with nan set for bad regions
      if planetname eq 'WASP-14b' then  begin
         im[4:7, 13:16, *] = !Values.F_NAN ;mask region with nan set for bad regions      
      endif

      if planetname eq 'HAT-P-22' and chname eq '2' then  im[19:25, 9:15, *] = !Values.F_NAN                              ;mask region with nan set for bad regions
      if planetname eq 'HAT-P-22' and chname eq '1' then  im[5:11, 11:17, *] = !Values.F_NAN                              ;mask region with nan set for bad regions

      if planetname eq 'HD93385' then im[17:21, 21:27, *] = !Values.F_NAN ;mask region with nan set for bad regions                         
      ;run the centroiding and photometry
      get_centroids_for_calstar_jk,im, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                   x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                   x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                   xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                   xfwhm, yfwhm, /WARM
      nanfound = where(FINITE(np) lt 1, nancount)
      if nancount gt 0 then print, 'NAN: ', fitsname(i), nanfound

      x_center = temporary(x3)
      y_center = temporary(y3)
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
      if keyword_set(breatheap) then begin
         abcdflux = betap(im, x_center, y_center, ronoise, gain, exptime, fluxconv,npcentroids, chname)
         ;XXXfake these for now
         fs = abcdflux
         back = abcdflux
         backerr = abcdflux
      endif 
      
;read in the raw data file and get the DN level of the peakpixel 
      ;only do this if the star is on the image
      if x_center(0) gt 0 and x_center(0) lt 32 and y_center(0) gt 0 and y_center(0) lt 32 then begin
        fits_read, rawname(i), rawdata, rawheader
        
        barrel = fxpar(rawheader, 'A0741D00')
        fowlnum = fxpar(rawheader, 'A0614D00')
        pedsig = fxpar(rawheader, 'A0742D00')
        ichan = fxpar (rawheader, 'CHNLNUM')
              
                                ;use Bill's code to conver to DN
;        dewrap2, rawdata, ichan, barrel, fowlnum, pedsig, 0, rawdata
     
                                ;or use Jim's code to convert to DN
        rawdata = irac_raw2dn(rawdata,ichan,barrel,fowlnum)
        rawdata = reform(rawdata, 32, 32, 64)
        peakpixDN = abcdflux  ; set up the array to be the same size as abcdflux
        if naxis eq 3 then begin
           for pp = 0, 63 do begin
           ;   print, 'x, y, f', x_center(pp), y_center(pp), rawdata[15,16,pp]
           ;   if x_center(pp) gt 0 and x_center(pp) lt 32 and y_center(pp) gt 0 and y_center(pp) lt 32 then begin
           ;      peakpixDN[pp] = rawdata[fix(x_center(pp)),fix(y_center(pp)),pp]
              peakpixDN[pp] = max(rawdata[13:16,13:16,pp])
           ;   endif else begin
          ;       print, 'dont understand', x_center(pp), y_center(pp)
          ;    endelse
           endfor
        endif
     endif

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

      if naxis eq 3 then begin  ; and i eq 0 then begin
         xarr[i*63] = x_center[1:*]
         yarr[i*63] = y_center[1:*]
         fluxarr[i*63] = abcdflux[1:*]
         fluxerrarr[i*63] = fs[1:*]
         timearr[i*63] = sclkarr[1:*]        
         bmjd[i*63] = bmjdarr[1:*]
         backarr[i*63] = back[1:*]
         backerrarr[i*63] = backerr[1:*]
         nparr[i*63] = np[1:*]
         npcentroidsarr[i*63] = npcentroids[1:*]
         xfwhmarr[i*63] = xfwhm[1:*]
         yfwhmarr[i*63] = yfwhm[1:*]
         peakpixDNarr[i*63] = peakpixDN[1:*]
 ;        help, bmjd
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
         fluxarr[i]  =  abcdflux
         fluxerrarr[i]  =  fs
          timearr[i]  = sclkarr
         bmjd[i]  = bmjdarr
         backarr[i]  =  back
         backerrarr[i]  = backerr
         nparr[i]  = np
         npcentroidsarr[i] = npcentroids
         xfwhmarr[i] = xfwhm
         yfwhmarr[i] = yfwhm
         peakpixDNarr[i] = peakpixDN
      endif

   endfor; for each fits file in the AOR


;    p1 = plot(xarr, yarr, '1.')
;    p2 = plot(timearr - timearr(0), yarr, '1.')
;    p3 = plot(timearr- timearr(0), corrfluxarr, '1.')
; print, 'fluxarr', fluxarr[0:10]
   
   
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
      n = where(phase lt 0)
      phase(n) = phase(n) + 1
   endif
  

;   if intended_phase gt 0.4 and intended_phase lt 0.6 then begin ;secondary eclipse
;      print, 'secondary eclipse intended'
;      phase = temporary(phase)+0.5
;   endif 
   
;print, 'testing phase', phase[0:200]
;print, 'end phase', phase[n_elements(phase) - 1]
;--------------------------------
;fill in that hash of hases
;--------------------------------


   if keyword_set(columntrack) then begin 
      
      keys =['ra', 'dec', 'xcen', 'ycen', 'flux','fluxerr', 'corrflux', 'corrfluxerr', 'sclktime_0', 'timearr', 'aor', 'bmjdarr',  'bkgd', 'bkgderr','np', 'centerpixarr1','centerpixarr2','centerpixarr3','centerpixarr4','centerpixarr5','centerpixarr6','sigmapixarr1','sigmapixarr2','sigmapixarr3','sigmapixarr4','sigmapixarr5','sigmapixarr6','phase']
      values=list(ra_ref,  dec_ref, xarr, yarr, fluxarr, fluxerrarr, corrfluxarr, corrfluxerrarr, sclk_0, timearr, aorname(a), bmjd,  backarr, backerrarr, nparr, centerpixarr1, centerpixarr2, centerpixarr3, centerpixarr4, centerpixarr5, centerpixarr6, sigmapixarr1, sigmapixarr2, sigmapixarr3, sigmapixarr4, sigmapixarr5, sigmapixarr6, phase)
      planethash[aorname(a)] = HASH(keys, values)
   endif else begin
      keys =['ra', 'dec', 'xcen', 'ycen', 'flux','fluxerr', 'sclktime_0', 'timearr', 'aor', 'bmjdarr', 'bkgd', 'bkgderr','np','phase', 'npcentroids','exptime','xfwhm', 'yfwhm','peakpixDN']
      values=list(ra_ref,  dec_ref, xarr, yarr, fluxarr, fluxerrarr, sclk_0, timearr, aorname(a), bmjd,  backarr, backerrarr, nparr, phase, npcentroidsarr, exptime, xfwhmarr, yfwhmarr, peakpixDNarr)
      planethash[aorname(a)] = HASH(keys, values)
   endelse
SKIPAOR: print, 'skipped this AOR'
endfor                          ;for each AOR


if keyword_set(breatheap) then begin
   savename = strcompress(dirname + planetname +'_phot_ch'+chname+'_varap.sav')
endif else begin
   savename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'.sav',/remove_all)
endelse

save, planethash, filename=savename
print, 'saving planethash', savename
print, 'time check', systime(1) - t1


;testplot = plot(timearr, xarr, '1s', sym_size = 0.1, sym_filled = 1, xtitle = 'time', ytitle = 'xcen')
plothist, fluxarr, xhist, yhist, /noplot
testplot = plot(xhist, yhist,  xtitle = 'Flux', ytitle = 'Number', thick =3, color = 'blue')

                                ; print, planethash.keys()
 ; print, planethash[aorname(0)].keys()
 ; print, 'testing (planethash[aorname(0),phase])[0:10]',(planethash[aorname(0),'phase'])[0:10]
;  print, 'n_elements in hash', n_elements(planethash[aorname(1),'xcen'])


;histogram the aperture sizes if breathing the aperture
;if keyword_set(breathap) then begin
;    plot_hist, aparr, xhist, yhist, bin = 0.05, /noplot
;   b = plot(xhist, yhist, title = 'aperture sizes')
; endif

;testplot the first AOR of raw raw flux before I start analyzing.

;tp = plot(planethash[aorname(0),'phase'], planethash[aorname(0),'flux'], '1rs')


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

  varap = sqrt(np)  +1.0  ; sqrt(np)  + 0.2

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

;  print, 'testing varap', varap
;  print, 'testing abcdflux', abcdflux
  return, abcdflux
end

         ;slice up image
;           	for s = 0, 63 do begin
;                   eslice = eim[*, *, s];;;

