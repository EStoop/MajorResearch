pro phot_exoplanet_print, planetname, chname
;do photometry on any IRAC staring mode exoplanet data in order to
;print output.  DOn't correct the photometry for pixel phase effect.
;now with hashes of hashes!
;print, 'waiting'
;wait, 3600  ; waiting for another code to finish
;print, 'done waiting'
 t1 = systime(1)
; aps1 = [ 2.0, 2.25, 2.5, 2.75, 3.0] - from get_centroids...

;case apradius of
;   2.0: apval = 0
;   2.25: apval = 1
;   2.5:  apval = 23
;   2.75: apval = 3
;   3.0: apval = 4;;

;   Else: apval = 2              ; if no decision, then choose an apradius = 2.5 pixels
;endcase

;run code to read in all the planet parameters
planetinfo = create_planetinfo()
;ra_ref = planetinfo[planetname, 'ra']
;dec_ref = planetinfo[planetname, 'dec']
if chname eq '2'  then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
basedir = planetinfo[planetname, 'basedir']
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
if planetname eq 'HAT-P-22' then exosystem = 'HAT-P-22 b'
if planetname eq 'GJ1214' then exosystem = 'GJ 1214 b'
if planetname eq '55Cnc' then exosystem = '55 Cnc e'
if planetname eq 'HD209458' then exosystem = 'HD 209458 b'
if planetname eq 'Kepler-5' then exosystem = 'Kepler-5 b'
if planetname eq 'Kepler-17' then exosystem = 'Kepler-17 b'
if planetname eq 'upsandb' then exosystem = 'upsilon And b'
if planetname eq 'XO3' then exosystem = 'XO-3 b' 
if planetname eq 'simul_XO3' then exosystem = 'XO-3 b' 

print, exosystem, 'exosystem'
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

print, 'parameters from create_planetinfo'
ra_ref = planetinfo[planetname, 'ra']
dec_ref = planetinfo[planetname, 'dec']
utmjd_center = planetinfo[planetname, 'utmjd_center']
period = planetinfo[planetname, 'period']
ra = ra_ref
print, 'ra', ra_ref, 'dec', dec_ref
if ra lt 400. then begin        ; that means get_exoplanet_data actually found the target
;   ra_ref = double(ra)*15.D       ; comes in hours!;
;   help, ra_ref
;   help, double(ra)
;   print, ra, ra_ref, double(ra) * 15.D, double(ra) * 15.D/ 15.D
;   dec_ref = dec
;   utmjd_center = mjd_transit
;   period = p_orbit
endif else begin
   print, 'parameters from create_planetinfo'
   ra_ref = planetinfo[planetname, 'ra']
   dec_ref = planetinfo[planetname, 'dec']
   utmjd_center = planetinfo[planetname, 'utmjd_center']
   period = planetinfo[planetname, 'period']
endelse
print, 'ra, dec', ra_ref, dec_ref
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

startaor = 0
stopaor =    n_elements(aorname) - 1
for a =startaor,  stopaor do begin
   print, 'working on ',aorname(a)
   dir = dirname+ string(aorname(a) ) 
   CD, dir                      ; change directories to the correct AOR directory
   command  = strcompress( 'find ch'+chname+"/bcd -name 'SPITZER*_bcd.fits' > "+dirname+'bcdlist.txt')
   spawn, command
   command2 =  strcompress('find ch'+chname+"/bcd -name '*bunc.fits' > "+dirname + 'bunclist.txt')
   spawn, command2
   command3 =  strcompress('find ch'+chname+"/raw -name '*dce.fits' > "+dirname + 'dcelist.txt')
   spawn, command3
 
   readcol,strcompress(dirname +'bcdlist.txt'),fitsname, format = 'A', /silent
   readcol,strcompress(dirname+'bunclist.txt'),buncname, format = 'A', /silent
   readcol,strcompress(dirname+'dcelist.txt'),rawname, format = 'A', /silent

   print,'n_elements(fitsname)', n_elements(fitsname)
   
   startfits = 0L


   for i =startfits, n_elements(fitsname) - 1  do begin ;read each cbcd file, find centroid, keep track
;;      print, 'working on ', fitsname(i)         
      header = headfits(fitsname(i)) ;
      sclk_obs= sxpar(header, 'SCLK_OBS')
      frametime = sxpar(header, 'FRAMTIME')
      bmjd_obs = sxpar(header, 'BMJD_OBS')
      ch = sxpar(header, 'CHNLNUM')
      ronoise = sxpar(header, 'RONOISE') ; these are zeros
      gain = sxpar(header, 'GAIN')
      fluxconv = sxpar(header, 'FLUXCONV')
      exptime = sxpar(header, 'EXPTIME')
      aintbeg = sxpar(header, 'AINTBEG')
      atimeend = sxpar(header, 'ATIMEEND')
      naxis = sxpar(header, 'NAXIS')
      framedly = sxpar(header, 'FRAMEDLY')
      if ra_ref gt 359. then begin
        print, 'inside ra_ref gt 359'
        ra_ref = sxpar(header, 'RA_REF')
        dec_ref = sxpar(header, 'DEC_REF')
     endif

      if ch eq '2' and frametime eq 2 then ronoise = 12.1
      if i eq startfits then begin
         print, 'ronoise', ronoise, gain, fluxconv, exptime, ra_ref, dec_ref
         sclk_0 = sclk_obs
      endif



      if i eq startfits and naxis eq 3 then begin
         xarr = fltarr(64*(n_elements(fitsname)))
         yarr = xarr
         fluxarr_2 = xarr
         fluxerrarr_2 = xarr
         fluxarr_2p25 = xarr
         fluxerrarr_2p25 = xarr
         fluxarr_2p5 = xarr
         fluxerrarr_2p5 = xarr
         fluxarr_2p75 = xarr
         fluxerrarr_2p75 = xarr
         fluxarr_3 = xarr
         fluxerrarr_3 = xarr
         bcdname = strarr(64*(n_elements(fitsname)))
;         framenum = intarr(64*(n_elements(fitsname)))
         timearr =  dblarr(64*(n_elements(fitsname)))
         bmjd = dblarr(64*(n_elements(fitsname)))
         backarr = xarr
         backerrarr = xarr
         nparr = xarr
         npcentroidsarr = xarr
         xfwhmarr = xarr
         yfwhmarr = xarr
         peakpixDNarr = xarr
         piarr = fltarr(64*n_elements(fitsname),9*9)

 ;        help, framenum

      endif
      if i eq startfits and naxis ne 3 then begin
         xarr = fltarr(n_elements(fitsname))
         yarr = xarr
         fluxarr_2 = xarr
         fluxerrarr_2 = xarr
         fluxarr_2p25 = xarr
         fluxerrarr_2p25 = xarr
         fluxarr_2p5 = xarr
         fluxerrarr_2p5 = xarr
         fluxarr_2p75 = xarr
         fluxerrarr_2p75 = xarr
         fluxarr_3 = xarr
         fluxerrarr_3 = xarr
         timearr =dblarr(n_elements(fitsname))
         bmjd = dblarr(n_elements(fitsname))
         backarr = xarr
         backerrarr = xarr
         nparr = xarr
         npcentroidsarr = xarr
         xfwhmarr = xarr
         yfwhmarr = xarr
         peakpixDNarr = xarr
         piarr = fltarr(n_elements(fitsname),9*9)
         bcdname = strarr((n_elements(fitsname)))

      endif
      fdarr = fltarr(n_elements(fitsname))
      fdarr[i] = framedly


      if naxis eq 3 then begin
         deltatime = (atimeend - aintbeg) / 64.D ; real time for each of the 64 frames
         nt = dindgen(64)
         sclkarr = sclk_obs  + (deltatime*nt)/60./60./24.D ; 0.5*frametime + frametime*nt
         bmjdarr= bmjd_obs + (deltatime*nt)/60./60./24.D   ; 0.5*(frametime/60./60./24.) + (frametime/60./60./24.)*nt
      endif else begin          ; full array, s`o don't need to expand out the times
         sclkarr = sclk_obs
         bmjdarr = bmjd_obs
      endelse

      ;read in the files
      fits_read, fitsname(i), im, h
      fits_read, buncname(i), unc, hunc
      test = sxpar(h, 'RAWFILE')
      

      ;apply mask file if necessary
      if planetname eq 'hd189733' then  im[13:16, 4:7, *] = !Values.F_NAN ;mask region with nan set for bad regions
      if planetname eq 'WASP-14b' then  begin
         ra_off = 218.27655; 218.27718
         dec_off = 21.897652; 21.89808
         adxy, h, ra_off, dec_off, x_off, y_off
         x_off = round(x_off)
         y_off = round(y_off)
         im[(x_off-3):(x_off + 3), (y_off - 3):(y_off+3), *] = !Values.F_NAN ;mask region with nan set for bad regions
      endif

      if planetname eq 'HAT-P-22' and chname eq '2' then  im[19:25, 9:15, *] = !Values.F_NAN   ;mask region with nan set for bad regions
      if planetname eq 'HAT-P-22' and chname eq '1' then  im[5:11, 11:17, *] = !Values.F_NAN   ;mask region with nan set for bad regions
      if planetname eq 'HD93385' then im[17:21, 21:27, *] = !Values.F_NAN ;mask region with nan set for bad regions                      
      

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
     ;choose the requested pixel aperture
;      abcdflux2 = f[*,0]      
;      fs2 = fs[*,0]
     ; 3-7 pixel background
      back = b[*,0]
      backerr = bs[*,0]
      npcentroids = np   ;keep track of np from get_centroids
      np = findgen(64)
      
      ;;read in the raw data file and get the DN level of the peakpixel      
      fits_read, rawname(i), rawdata, rawheader
        
      barrel = fxpar(rawheader, 'A0741D00')
      fowlnum = fxpar(rawheader, 'A0614D00')
      pedsig = fxpar(rawheader, 'A0742D00')
      ichan = fxpar (rawheader, 'CHNLNUM')
              
      ;;use Bill's code to convery to DN
      dewrap2, rawdata, ichan, barrel, fowlnum, pedsig, 0, rawdata
     
      ;;or use Jim's code to convert to DN
      ;;rawdata = irac_raw2dn(rawdata,ichan,barrel,fowlnum)
      rawdata = reform(rawdata, 32, 32, 64)
      peakpixDN = back      ; set up the array to be the same size as abcdflux
      if naxis eq 3 then begin
         for pp = 0, 63 do begin
            peakpixDN[pp] = max(rawdata[13:16,13:16,pp])
         endfor
      endif

;track the values of the 9x9 pixel box around the centroid
;        pi = track_box(im, x_center, y_center)   ; tb now a 25 x 64 element array
      if naxis eq 3 then begin
         if i eq startfits then nframe = 0
         for subframe = 0, 63 do begin
            earr = fltarr(81)
            c = 0
            
            for xcoord = 11, 19 do begin
               for ycoord = 11, 19 do begin
                  ;;if xcoord eq 11 and ycoord eq 11 and subframe eq 0 then print, '[11, 11, 0]', im[xcoord, ycoord, subframe]
                  earr[c] =  im[xcoord, ycoord, subframe]
                  c++
               endfor
            endfor
;            if subframe eq 0 then print, 'earr', earr
            piarr[nframe,*] =earr
            nframe++

         endfor

      endif
      
;      help, piarr
;      print, piarr[0,*]


      if naxis eq 2 then pi = im[(fix(x_center) - 4):(fix(x_center+4)), (fix(y_center) - 4):(fix(y_center+4))] ; now a 9x9x64 element array
      
;---------------------------------

      if naxis eq 3 then begin  ; and i eq 0 then begin
         xarr[i*64] = x_center
         yarr[i*64] = y_center
         fluxarr_2[i*64] = f[*,0]   
         fluxerrarr_2[i*64] = fs[*,0]
         fluxarr_2p25[i*64] = f[*,1]   
         fluxerrarr_2p25[i*64] = fs[*,1]
         fluxarr_2p5[i*64] = f[*,2]   
         fluxerrarr_2p5[i*64] = fs[*,2]
         fluxarr_2p75[i*64] = f[*,3]   
         fluxerrarr_2p75[i*64] = fs[*,3]
         fluxarr_3[i*64] = f[*,4]   
         fluxerrarr_3[i*64] = fs[*,4]
         timearr[i*64] = sclkarr      
         bmjd[i*64] = bmjdarr
         backarr[i*64] = back
         backerrarr[i*64] = backerr
         nparr[i*64] = np
         npcentroidsarr[i*64] = npcentroids
         xfwhmarr[i*64] = xfwhm
         yfwhmarr[i*64] = yfwhm
         if keyword_set(rawfile) then peakpixDNarr[i*64] = peakpixDN
;         piarr[i*64] = pi[*,*,*]
         bcdname[i*64:(i*64) + 63] = fitsname[i]
 ;        framenum[i*64] = indgen(64)
         peakpixdnarr[i*64] = peakpixDN
      endif 
      if naxis eq 2 then begin; and i eq 0 then begin
         xarr[i] = x_center
         yarr[i]  =  y_center
         fluxarr_2[i] = f[0]   
         fluxerrarr_2[i] = fs[0]
         fluxarr_2p25[i] = f[1]   
         fluxerrarr_2p25[i] = fs[1]
         fluxarr_2p5[i] = f[2]   
         fluxerrarr_2p5[i] = fs[2]
         fluxarr_2p75[i] = f[3]   
         fluxerrarr_2p75[i] = fs[3]
         fluxarr_3[i] = f[4]   
         fluxerrarr_3[i] = fs[4]
         timearr[i]  = sclkarr
         bmjd[i]  = bmjdarr
         backarr[i]  =  back
         backerrarr[i]  = backerr
         nparr[i]  = npcentroids
         npcentroidsarr[i] = npcentroids
         bcdname[i] = fitsname[i]

         xfwhmarr[i] = xfwhm
         yfwhmarr[i] = yfwhm
         if keyword_set(rawfile) then peakpixDNarr[i] = peakpixDN
         piarr[i,*] = pi
      endif

      if a eq startaor and i eq startfits then begin
         time_0 = bmjdarr(0)
      endif

   endfor; for each fits file in the AOR

   
;--------------------------------
;fill in that hash of hases
;--------------------------------
;make a framenum array
   if naxis eq 3 then begin
      framenum = intarr(64*n_elements(fitsname))
      for i = 0, n_elements(fitsname) - 1 do begin
         framenum[i*64] = indgen(64)
      endfor
   endif else begin
       framenum = indgen(n_elements(fitsname))
    endelse

   keys =[ 'bcdname', 'framenum','xcen', 'ycen', 'flux_2','fluxerr_2','flux_2p25','fluxerr_2p25','flux_2p5','fluxerr_2p5','flux_2p75','fluxerr_2p75','flux_3','fluxerr_3', 'BMJD', 'bkgd', 'bkgderr','np','xfwhm', 'yfwhm','peakpixDN'];,'pixvals']
   values=list(bcdname, framenum, xarr, yarr, fluxarr_2, fluxerrarr_2,fluxarr_2p25, fluxerrarr_2p25,fluxarr_2p5, fluxerrarr_2p5,fluxarr_2p75, fluxerrarr_2p75,fluxarr_3, fluxerrarr_3,  bmjd,  backarr, backerrarr, npcentroidsarr, xfwhmarr, yfwhmarr, peakpixDNarr);, piarr)

   help, bcdname, framenum, xarr, yarr, fluxarr_2, bmjd, backarr, bacjerrarr, npcentroidsarr, xfwhmarr, yfwhmarr, peakpixDNarr

   ;;create a structure to house all the variables except pixvals
   allstruct = { bcdname:bcdname, framenum:framenum, xcen:xarr, ycen:yarr, flux_2:fluxarr_2, fluxerr_2:fluxerrarr_2, flux_2p25:fluxarr_2p25, fluxerr_2p25:fluxerrarr_2p25, flux_2p5:fluxarr_2p5, fluxerr_2p5:fluxerrarr_2p5, flux_2p75:fluxarr_2p75, fluxerr_2p75:fluxerrarr_2p75, flux_3:fluxarr_3, fluxerr_3:fluxerrarr_3, bmjd:bmjd, bkgd:backarr, bkgderr:backerrarr,np:npcentroidsarr,xfwhm:xfwhmarr, yfwhm:yfwhmarr,peakpix:peakpixDNarr }

   outname = strcompress(dirname + planetname +'_'+aorname(a) +'_phot_ch'+chname+'.txt',/remove_all)
   write_csv, outname,allstruct, header = keys

   ;;create a structure to house all the variables including pixvals
   keys_long =[ 'bcdname', 'framenum','xcen', 'ycen', 'flux_2','fluxerr_2','flux_2p25','fluxerr_2p25','flux_2p5','fluxerr_2p5','flux_2p75','fluxerr_2p75','flux_3','fluxerr_3', 'BMJD', 'bkgd', 'bkgderr','np','xfwhm', 'yfwhm','peakpixDN','pixval0-0', 'pixval0-1', 'pixval0-2', 'pixval0-3', 'pixval0-4', 'pixval0-5', 'pixval0-6', 'pixval0-7', 'pixval0-8', 'pixval1-0', 'pixval1-1', 'pixval1-2', 'pixval1-3', 'pixval1-4', 'pixval1-5', 'pixval1-6', 'pixval1-7', 'pixval1-8', 'pixval2-0', 'pixval2-1', 'pixval2-2', 'pixval2-3', 'pixval2-4', 'pixval2-5', 'pixval2-6', 'pixval2-7', 'pixval2-8', 'pixval3-0', 'pixval3-1', 'pixval3-2', 'pixval3-3', 'pixval3-4', 'pixval3-5', 'pixval3-6', 'pixval3-7', 'pixval3-8', 'pixval4-0', 'pixval4-1', 'pixval4-2', 'pixval4-3', 'pixval4-4', 'pixval4-5', 'pixval4-6', 'pixval4-7', 'pixval4-8', 'pixval5-0', 'pixval5-1', 'pixval5-2', 'pixval5-3', 'pixval5-4', 'pixval5-5', 'pixval5-6', 'pixval5-7', 'pixval5-8', 'pixval6-0', 'pixval6-1', 'pixval6-2', 'pixval6-3', 'pixval6-4', 'pixval6-5', 'pixval6-6', 'pixval6-7', 'pixval6-8', 'pixval7-0', 'pixval7-1', 'pixval7-2', 'pixval7-3', 'pixval7-4', 'pixval7-5', 'pixval7-6', 'pixval7-7', 'pixval7-8', 'pixval8-0', 'pixval8-1', 'pixval8-2', 'pixval8-3', 'pixval8-4', 'pixval8-5', 'pixval8-6', 'pixval8-7', 'pixval8-8']

   allstruct_long = { bcdname:bcdname, framenum:framenum, xcen:xarr, ycen:yarr, flux_2:fluxarr_2, fluxerr_2:fluxerrarr_2, flux_2p25:fluxarr_2p25, fluxerr_2p25:fluxerrarr_2p25, flux_2p5:fluxarr_2p5, fluxerr_2p5:fluxerrarr_2p5, flux_2p75:fluxarr_2p75, fluxerr_2p75:fluxerrarr_2p75, flux_3:fluxarr_3, fluxerr_3:fluxerrarr_3, bmjd:bmjd, bkgd:backarr, bkgderr:backerrarr,np:npcentroidsarr,xfwhm:xfwhmarr, yfwhm:yfwhmarr,peakpix:peakpixDNarr, piarr0:piarr[*,0], piarr1:piarr[*,1], piarr2:piarr[*,2], piarr3:piarr[*,3], piarr4:piarr[*,4], piarr5:piarr[*,5], piarr6:piarr[*,6], piarr7:piarr[*,7], piarr8:piarr[*,8], piarr9:piarr[*,9], piarr10:piarr[*,10], piarr11:piarr[*,11], piarr12:piarr[*,12], piarr13:piarr[*,13], piarr14:piarr[*,14], piarr15:piarr[*,15], piarr16:piarr[*,16], piarr17:piarr[*,17], piarr18:piarr[*,18], piarr19:piarr[*,19], piarr20:piarr[*,20], piarr21:piarr[*,21], piarr22:piarr[*,22], piarr23:piarr[*,23], piarr24:piarr[*,24], piarr25:piarr[*,25], piarr26:piarr[*,26], piarr27:piarr[*,27], piarr28:piarr[*,28], piarr29:piarr[*,29], piarr30:piarr[*,30], piarr31:piarr[*,31], piarr32:piarr[*,32], piarr33:piarr[*,33], piarr34:piarr[*,34], piarr35:piarr[*,35], piarr36:piarr[*,36], piarr37:piarr[*,37], piarr38:piarr[*,38], piarr39:piarr[*,39], piarr40:piarr[*,40], piarr41:piarr[*,41], piarr42:piarr[*,42], piarr43:piarr[*,43], piarr44:piarr[*,44], piarr45:piarr[*,45], piarr46:piarr[*,46], piarr47:piarr[*,47], piarr48:piarr[*,48], piarr49:piarr[*,49], piarr50:piarr[*,50], piarr51:piarr[*,51], piarr52:piarr[*,52], piarr53:piarr[*,53], piarr54:piarr[*,54], piarr55:piarr[*,55], piarr56:piarr[*,56], piarr57:piarr[*,57], piarr58:piarr[*,58], piarr59:piarr[*,59], piarr60:piarr[*,60], piarr61:piarr[*,61], piarr62:piarr[*,62], piarr63:piarr[*,63], piarr64:piarr[*,64], piarr65:piarr[*,65], piarr66:piarr[*,66], piarr67:piarr[*,67], piarr68:piarr[*,68], piarr69:piarr[*,69], piarr70:piarr[*,70], piarr71:piarr[*,71], piarr72:piarr[*,72], piarr73:piarr[*,73], piarr74:piarr[*,74], piarr75:piarr[*,75], piarr76:piarr[*,76], piarr77:piarr[*,77], piarr78:piarr[*,78], piarr79:piarr[*,79], piarr80:piarr[*,80]}

   outname = strcompress(dirname + planetname +'_'+aorname(a) +'_phot_ch'+chname+'_pixvals.txt',/remove_all)
   write_csv, outname,allstruct_long, header = keys_long

;--------------------------------
;make an output text file
;--------------------------------

;   outname = strcompress(dirname + planetname +'_'+aorname(a) +'_phot_ch'+chname+'.txt',/remove_all)
;   openw, outlun, outname, /GET_LUN
;   printf, outlun, 'bcdname, framenum, xarr, yarr, bmjd, flux_2, fluxerr_2, flux_2p25, fluxerr_2p25, fluxarr_2p5, fluxerr_2p5, flux_2p75, fluxerr_2p75, flux_3, fluxerr_3, back, backerr, np, xfwhm, yfwhm, peakpixDN', format = '(A)'

;   outname_long = strcompress(dirname + planetname +'_'+aorname(a) +'_phot_ch'+chname+'_long.txt',/remove_all)
;   openw, outlun_long, outname_long, /GET_LUN

;   for j = 0,64*n_elements(fitsname) - 1 do begin

;      printf, outlun, bcdname(j), framenum(j), xarr(j), yarr(j), bmjd(j), fluxarr_2(j), fluxerrarr_2(j), fluxarr_2p25(j), fluxerrarr_2p25(j), fluxarr_2p5(j), fluxerrarr_2p5(j), fluxarr_2p75(j), fluxerrarr_2p75(j), fluxarr_3(j), fluxerrarr_3(j), backarr(j), backerrarr(j), npcentroidsarr(j), xfwhmarr(j), yfwhmarr(j), '  ', peakpixDNarr(j), format = '(A, I10,F10.5, F10.5, D, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.3, F10.3, F10.3, A, F10.5)'

;  printf, outlun_long, bcdname(j), framenum(j), xarr(j), yarr(j), bmjd(j), fluxarr_2(j), fluxerrarr_2(j), fluxarr_2p25(j), fluxerrarr_2p25(j), fluxarr_2p5(j), fluxerrarr_2p5(j), fluxarr_2p75(j), fluxerrarr_2p75(j), fluxarr_3(j), fluxerrarr_3(j), backarr(j), backerrarr(j), npcentroidsarr(j), xfwhmarr(j), yfwhmarr(j), '  ', peakpixDNarr(j), piarr[j,0],piarr[j,1],piarr[j,2],piarr[j,3],piarr[j,4],piarr[j,5],piarr[j,6],piarr[j,7],piarr[j,8],piarr[j,9],piarr[j,10],piarr[j,11],piarr[j,12],piarr[j,13],piarr[j,14],piarr[j,15],piarr[j,16],piarr[j,17],piarr[j,18],piarr[j,19],piarr[j,20],piarr[j,21],piarr[j,22],piarr[j,23],piarr[j,24],piarr[j,25],piarr[j,26],piarr[j,27],piarr[j,28],piarr[j,29],piarr[j,30],piarr[j,31],piarr[j,32],piarr[j,33],piarr[j,34],piarr[j,35],piarr[j,36],piarr[j,37],piarr[j,38],piarr[j,39],piarr[j,40],piarr[j,41],piarr[j,42],piarr[j,43],piarr[j,44],piarr[j,45],piarr[j,46],piarr[j,47],piarr[j,48],piarr[j,49],piarr[j,50],piarr[j,51],piarr[j,52],piarr[j,53],piarr[j,54],piarr[j,55],piarr[j,56],piarr[j,57],piarr[j,58],piarr[j,59],piarr[j,60],piarr[j,61],piarr[j,62],piarr[j,63],piarr[j,64],piarr[j,65],piarr[j,66],piarr[j,67],piarr[j,68],piarr[j,69],piarr[j,70],piarr[j,71],piarr[j,72],piarr[j,73],piarr[j,74],piarr[j,75],piarr[j,76],piarr[j,77],piarr[j,78],piarr[j,79],piarr[j,80], format = '(A, I10,F10.5, F10.5, D, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.3, F10.3, F10.3, A, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5 )'


;   endfor
;   free_lun, outlun
;   free_lun, outlun_long
endfor                          ;for each AOR


savename = strcompress(dirname + planetname +'_phot_ch'+chname+'.sav',/remove_all)

save, planethash, filename=savename
print, 'saving planethash', savename
print, 'time check', systime(1) - t1




 




end

