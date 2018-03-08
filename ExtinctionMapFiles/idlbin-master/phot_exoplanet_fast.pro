pro phot_exoplanet, planetname, columntrack = columntrack, breatheap = breatheap
;do photometry on any IRAC staring mode exoplanet data
;now with hashes of hashes!
 t1 = systime(1)
;run code to read in all the planet parameters
planetinfo = create_planetinfo()
chname =1
ra_ref =245.15145
dec_ref = 41.0479777777779 
aorname = '0046475776'
basedir = '/Users/jkrick/irac_warm/pcrs_test/'

;---------------
dirname = basedir; strcompress(basedir + planetname +'/')
planethash = hash()


  for a = 0,   n_elements(aorname) - 1 do begin
     print, 'working on ',aorname(a)
     dir = dirname+ string(aorname(a) ) 
     CD, dir                    ; change directories to the correct AOR directory
     command  = strcompress( 'ls *bcd_fp.fits > '+dirname+'bcdlist.txt')
     spawn, command
     command2 =  strcompress('ls *sig_dntoflux.fits > '+dirname + 'bunclist.txt')
     spawn, command2

     readcol,strcompress(dirname +'bcdlist.txt'),fitsname, format = 'A', /silent
     readcol,strcompress(dirname+'bunclist.txt'),buncname, format = 'A', /silent

;     aparr = dblarr(n_elements(fitsname))  ;keep the aperture sizes used

     for i =0.D,  n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
        ;print, 'working on ', fitsname(i)         
        header = headfits(fitsname(i)) ;
        sclk_obs= sxpar(header, 'SCLK_OBS')
        frametime = sxpar(header, 'FRAMTIME')
        bmjd_obs = sxpar(header, 'BMJD_OBS')
        ;utcs_obs = sxpar(header, 'UTCS_OBS')
        ch = sxpar(header, 'CHNLNUM')
        ronoise = sxpar(header, 'RONOISE')
        gain = sxpar(header, 'GAIN')
        fluxconv = sxpar(header, 'FLUXCONV')
        exptime = sxpar(header, 'EXPTIME')
        aintbeg = sxpar(header, 'AINTBEG')
        atimeend = sxpar(header, 'ATIMEEND')

        if i eq 0 then sclk_0 = sxpar(header, 'SCLK_OBS')

        sclkarr = dblarr(64)
        bmjdarr = dblarr(64)
        ;utcsarr = dblarr(64)
        deltatime = (atimeend - aintbeg) / 64.D  ; real time for each of the 64 frames
        for nt = 0., 64 - 1 do begin
           sclkarr[nt] = sclk_obs  + (deltatime*nt)/60./60./24.D ; 0.5*frametime + frametime*nt
          ; utcsarr[nt]= utcs_obs + (deltatime*nt)/60./60./24.D;+ 0.5*(frametime/60./60./24.)  + (frametime/60./60./24.) *nt

 ;          print, 'deltattime in seconds', (deltatime*nt)/60./60./24.D, format = '(A, F0)'
;          print, 'bmjdobs', bmjd_obs, format =  '(A, F0)'
;           print, 'addition',  bmjd_obs + (deltatime*nt)/60./60./24.D, format =  '(A, F0)'
           bmjdarr[nt]= bmjd_obs + (deltatime*nt)/60./60./24.D; 0.5*(frametime/60./60./24.) + (frametime/60./60./24.)*nt
;           print, 'time test', bmjdarr[nt], format = '(A, F0)'
        endfor
         
        ;read in the files
        fits_read, fitsname(i), im, h
        fits_read, buncname(i), unc, hunc

      
                                ; fits_read, covname, covdata, covheader
        get_centroids_for_calstar_jk,im, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                     x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                     x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                     xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                     xfwhm, yfwhm, /WARM
        
        x_center = x3
        y_center = y3

        ;calculate noise pixel
        np = noisepix(im, x_center, y_center, ronoise, gain, exptime, fluxconv)
        
        ; if changing the apertures then use this to calculate photometry
        if keyword_set(breatheap) then begin
           abcdflux = betap(im, x_center, y_center, ronoise, gain, exptime, fluxconv,np, chname)
           ;XXXfake these for now
           fs = abcdflux
           back = abcdflux
           backerr = abcdflux
        endif else begin

      ;choose 3 pixel aperture, 3-7 pixel background
           abcdflux = f[*,9]      
           fs = fs[*,9]
        ;choose [10, 12-20]
;        abcdflux = f[*,5]
;        fs = fs[*,5]
           
           back = b[*,2]
           backerr = bs[*,2]
        endelse


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

           
;--------------------------------
      ;correct for pixel phase effect based on pmaps from Jim
        file_suffix = ['500x500_0043_120828.fits','0p1s_x4_500x500_0043_121120.fits'];121120 ['500x500_0043_120409.fits','0p1s_x4_500x500_0043_120124.fits']
        corrflux = iracpc_pmap_corr(abcdflux,x_center,y_center,ch,/threshold_occ, threshold_val = 20, file_suffix = file_suffix)
        corrfluxerr = fs        ;leave out the pmap err for now
     
 ;---------------------------------
        if i eq 0 then begin
           xarr = x_center[1:*]
           yarr = y_center[1:*]
           fluxarr = abcdflux[1:*]
           fluxerrarr = fs[1:*]
           corrfluxarr = corrflux[1:*]
           corrfluxerrarr = corrfluxerr[1:*]
           timearr = sclkarr[1:*]        
           bmjd = bmjdarr[1:*]
;           utcs = utcsarr[1:*]
           backarr = back[1:*]
           backerrarr = backerr[1:*]
           nparr = np[1:*]

           if keyword_set(columntrack) then begin 
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


        endif else begin
           xarr = [xarr, x_center[1:*]]
           yarr = [yarr, y_center[1:*]]
           fluxarr = [fluxarr, abcdflux[1:*]]
           fluxerrarr = [fluxerrarr, fs[1:*]]
           corrfluxarr = [corrfluxarr, corrflux[1:*]]
           corrfluxerrarr = [corrfluxerrarr, corrfluxerr[1:*]]
           timearr = [timearr, sclkarr[1:*]]
           bmjd = [bmjd, bmjdarr[1:*]]
;           utcs = [utcs, utcsarr[1:*]]
           backarr = [backarr, back[1:*]]
           backerrarr = [backerrarr, backerr[1:*]]
           nparr = [nparr, np[1:*]]
          if keyword_set(columntrack) then begin 
              centerpixarr1 = [centerpixarr1, centerpixval1[1:*]]
              centerpixarr2 = [centerpixarr2, centerpixval2[1:*]]
              centerpixarr3 = [centerpixarr3, centerpixval3[1:*]]
              centerpixarr4 = [centerpixarr4, centerpixval4[1:*]]
              centerpixarr5 = [centerpixarr5, centerpixval5[1:*]]
              centerpixarr6 = [centerpixarr6, centerpixval6[1:*]]
              sigmapixarr1 = [sigmapixarr1, sigmapixval1[1:*]]
              sigmapixarr2 = [sigmapixarr2, sigmapixval2[1:*]]
              sigmapixarr3 = [sigmapixarr3, sigmapixval3[1:*]]
              sigmapixarr4 = [sigmapixarr4, sigmapixval4[1:*]]
              sigmapixarr5 = [sigmapixarr5, sigmapixval5[1:*]]
              sigmapixarr6 = [sigmapixarr6, sigmapixval6[1:*]]
           endif


         endelse
        
     endfor                     ; for each fits file in the AOR
 ;    p1 = plot(xarr, yarr, '1.')
 ;    p2 = plot(timearr - timearr(0), yarr, '1.')
 ;    p3 = plot(timearr- timearr(0), corrfluxarr, '1.')
    ; print, 'fluxarr', fluxarr[0:10]
                                ;fill in that hash of hases

          if keyword_set(columntrack) then begin 

             keys =['ra', 'dec', 'xcen', 'ycen', 'flux','fluxerr', 'corrflux', 'corrfluxerr', 'sclktime_0', 'timearr', 'aor', 'bmjdarr',  'bkgd', 'bkgderr','np', 'centerpixarr1','centerpixarr2','centerpixarr3','centerpixarr4','centerpixarr5','centerpixarr6','sigmapixarr1','sigmapixarr2','sigmapixarr3','sigmapixarr4','sigmapixarr5','sigmapixarr6']
             values=list(ra_ref,  dec_ref, xarr, yarr, fluxarr, fluxerrarr, corrfluxarr, corrfluxerrarr, sclk_0, timearr, aorname(a), bmjd,  backarr, backerrarr, nparr, centerpixarr1, centerpixarr2, centerpixarr3, centerpixarr4, centerpixarr5, centerpixarr6, sigmapixarr1, sigmapixarr2, sigmapixarr3, sigmapixarr4, sigmapixarr5, sigmapixarr6)
             planethash[aorname(a)] = HASH(keys, values)
          endif else begin
             keys =['ra', 'dec', 'xcen', 'ycen', 'flux','fluxerr', 'corrflux', 'corrfluxerr', 'sclktime_0', 'timearr', 'aor', 'bmjdarr', 'bkgd', 'bkgderr','np']
             values=list(ra_ref,  dec_ref, xarr, yarr, fluxarr, fluxerrarr, corrfluxarr, corrfluxerrarr, sclk_0, timearr, aorname(a), bmjd,  backarr, backerrarr, nparr)
             planethash[aorname(a)] = HASH(keys, values)
          endelse

  endfor                        ;for each AOR
  
  ;test
  if keyword_set(breatheap) then begin
     savename = strcompress(dirname + planetname +'_phot_ch'+chname+'_varap.sav')
  endif else begin
     savename = strcompress(dirname + planetname +'_phot_ch'+chname+'.sav')
  endelse
  save, planethash, filename=savename
  print, 'saving planethash', savename
  print, 'time check', systime(1) - t1

 ; print, planethash.keys()
;  print, planethash[aorname(0)].keys()
 ; print, 'testing (planethash[aorname(0),flux])[0:10]',(planethash[aorname(0),'flux'])[0:10]
;  print, 'n_elements in hash', n_elements(planethash[aorname(1),'xcen'])
  

  ;histogram the aperture sizes if breathing the aperture
  ;if keyword_set(breathap) then begin
 ;    plot_hist, aparr, xhist, yhist, bin = 0.05, /noplot
  ;   b = plot(xhist, yhist, title = 'aperture sizes')
 ; endif

end




;function to calcluate noise pixel
function noisepix, im, xcen, ycen, ronoise, gain, exptime, fluxconv

  convfac = gain*exptime/fluxconv
   np = fltarr(64)
  for npj = 0, 63 do begin
     indim = im[*,*,npj]
     indim = indim*convfac
     
     aper, indim, xcen[npj], ycen[npj], topflux, topfluxerr, xb, xbs, 1.0, 4,[10,12],/flux,/exact, /silent, /nan, readnoise = ronoise, setskyval = 0
     aper, indim^2, xcen[npj], ycen[npj], bottomflux, bottomfluxerr, xb, xbs, 1.0,4,[10,12],/flux,/exact, /silent, /nan, readnoise = ronoise, setskyval = 0
     
     beta = topflux^2 / bottomflux
;     print, npj, beta
     np[npj] = beta
  endfor

     return, np                 ;this should be a 64 element array
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

  abcdflux = fltarr(64)
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

