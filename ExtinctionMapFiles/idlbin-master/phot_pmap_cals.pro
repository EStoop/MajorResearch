pro phot_pmap_cals
  dirname = '/Users/jkrick/iracdata/flight/IWIC/'
  aorname = ['IRAC029900/bcd/0045217536', 'IRAC029900/bcd/0045217792', 'IRAC030000/bcd/0045255424', 'IRAC030000/bcd/0045255680','IRAC030100/bcd/0045286656', 'IRAC030100/bcd/0045286912','IRAC030200/bcd/0045381376', 'IRAC030200/bcd/0045381632','IRAC030300/bcd/0045423872', 'IRAC030300/bcd/0045424128','IRAC030400/bcd/0045528576', 'IRAC030400/bcd/0045528832','IRAC030500/bcd/0045567232', 'IRAC030500/bcd/0045567488']
  colorarr = ['deep_pink', 'magenta', 'medium_purple', 'hot_pink', 'light_pink', 'rosy_brown', 'chocolate', 'saddle_brown', 'maroon', 'orange_red', 'dark_orange', 'gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green', 'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki','black', 'light_cyan', 'lavender', 'peach_puff', 'pale_goldenrod','red']

  nfits = 242L*63L    
  print, 'nfits',nfits
  pmap = replicate({pmapob, xcen:dblarr(nfits),ycen:dblarr(nfits),flux:dblarr(nfits),fluxerr:dblarr(nfits), corrflux:dblarr(nfits), corrfluxerr:dblarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' ', bmjdarr:dblarr(nfits), utcsarr:dblarr(nfits)},n_elements(aorname))
  
  ra_ref =269.7283 
  dec_ref = 67.793528

  for a = 0,  n_elements(aorname) - 1 do begin
     print, 'working on ',aorname(a)
     dir = dirname+ string(aorname(a) ) 
     CD, dir                    ; change directories to the correct AOR directory
     command  =  "ls *bcd_fp.fits > /Users/jkrick/irac_warm/pmap/bcdlist.txt"
     spawn, command
     command2 =  "ls *sig_dntoflux.fits > /Users/jkrick/irac_warm/pmap/bunclist.txt"
     spawn, command2

     readcol,'/Users/jkrick/irac_warm/pmap/bcdlist.txt',fitsname, format = 'A', /silent
     readcol,'/Users/jkrick/irac_warm/pmap/bunclist.txt',buncname, format = 'A', /silent

     for i =0.D,  n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
        ;print, 'working on ', fitsname(i)         
                                ;now were did it really point to?
        header = headfits(fitsname(i)) ;
        sclk_obs= sxpar(header, 'SCLK_OBS')
        frametime = sxpar(header, 'FRAMTIME')
        bmjd_obs = sxpar(header, 'BMJD_OBS')
        utcs_obs = sxpar(header, 'UTCS_OBS')

        if i eq 0 then sclk_0 = sxpar(header, 'SCLK_OBS')

        sclkarr = dblarr(64)
        bmjdarr = fltarr(64)
        utcsarr = fltarr(64)
        for nt = 0., 64 - 1 do begin
           sclkarr[nt] = (sclk_obs ) + 0.5*frametime + frametime*nt
           bmjdarr[nt]= bmjd_obs + 0.5*frametime + frametime*nt
           utcsarr[nt]= utcs_obs + 0.5*frametime + frametime*nt
        endfor
         
        ;do the photometry
   
        fits_read, fitsname(i), im, h
        fits_read, buncname(i), unc, hunc
                                ; fits_read, covname, covdata, covheader
        get_centroids_for_calstar_jk,im, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                     x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                     x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                     xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                     xfwhm, yfwhm, /WARM
        
      ;choose 3 pixel aperture, 3-7 pixel background
        abcdflux = f[*,9]       ;put it back in same nomenclature
        fs = fs[*,9]
        x_center = x3
        y_center = y3
  
;--------------------------------
      ;correct for pixel phase effect based on pmaps from Jim
        
        corrflux = iracpc_pmap_corr(abcdflux, x_center, y_center, 2, /threshold_occ, threshold_val = 20)
        corrfluxerr = fs        ;leave out the pmap err for now
     
 ;--------------------------------

        if  i eq 0 then begin
           xarr = x_center[1:*]
           yarr = y_center[1:*]
           fluxarr = abcdflux[1:*]
           fluxerrarr = fs[1:*]
           corrfluxarr = corrflux[1:*]
           corrfluxerrarr = corrfluxerr[1:*]
           timearr = sclkarr[1:*]        
           bmjd = bmjdarr[1:*]
           utcs = utcsarr[1:*]

        endif else begin
           xarr = [xarr, x_center[1:*]]
           yarr = [yarr, y_center[1:*]]
           fluxarr = [fluxarr, abcdflux[1:*]]
           fluxerrarr = [fluxerrarr, fs[1:*]]
           corrfluxarr = [corrfluxarr, corrflux[1:*]]
           corrfluxerrarr = [corrfluxerrarr, corrfluxerr[1:*]]
           timearr = [timearr, sclkarr[1:*]]
           bmjd = [bmjd, bmjdarr[1:*]]
           utcs = [utcs, utcsarr[1:*]]
         endelse
        
     endfor                     ; for each fits file in the AOR
     if a eq 0 then begin
        p1 = plot(xarr, yarr, '1.', xtitle = 'X pix', ytitle = 'Y pix', xrange = [14.5, 15.5], yrange = [14.5, 15.5],color = colorarr[a]) 
     endif else begin
        p1 = plot(xarr, yarr, '1.',/overplot,color = colorarr[a]) 
     endelse

    print, 'n xarr, fluxarr, nfits', n_elements(xarr), n_elements(fluxarr), nfits

    if n_elements(xarr) lt nfits then begin
       xarr2 = dblarr(nfits)
       yarr2 = dblarr(nfits)
       fluxarr2= dblarr(nfits)
       fluxerrarr2= dblarr(nfits)
       corrfluxarr2= dblarr(nfits)
       corrfluxerrarr2= dblarr(nfits)
       timearr2= dblarr(nfits)
       bmjd2 = dblarr(nfits)
       utcs2 = dblarr(nfits)

       xarr2[0:n_elements(xarr) - 1] = xarr
       yarr2[0:n_elements(xarr) - 1] = yarr
       fluxarr2[0:n_elements(xarr) - 1] = fluxarr
       fluxerrarr2[0:n_elements(xarr) - 1] = fluxerrarr
       corrfluxarr2[0:n_elements(xarr) - 1] = corrfluxarr
       corrfluxerrarr2[0:n_elements(xarr) - 1] = corrfluxerrarr
       timearr2[0:n_elements(xarr) - 1] = timearr
       bmjd2[0:n_elements(xarr) - 1] = bmjd
       utcs2[0:n_elements(xarr) - 1] = utcs

       pmap[a] ={pmapob,  xarr2,yarr2,fluxarr2,fluxerrarr2, corrfluxarr2, corrfluxerrarr2, sclk_0, timearr2,aorname(a), bmjd2,utcs2}
    endif else begin
       pmap[a] ={pmapob,  xarr,yarr,fluxarr,fluxerrarr, corrfluxarr, corrfluxerrarr, sclk_0, timearr,aorname(a),bmjd, utcs}
    endelse
  endfor                        ;for each AOR
  
 save, pmap, filename='/Users/jkrick/pmap/pmap_cals.sav'


;now plot them all on the image 
;could also put this in its own plot_* code
fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p1s_x4_500x500_0043_120124/pmap_ch2_0p1s_x4_500x500_0043_120124.fits', pmapdata, pmapheader

 c = contour(pmapdata, /fill, n_levels = 21, rgb_table = 0, xtitle = 'X (pixel)', ytitle = 'Y (pixel)', aspect_ratio = 1, xrange = [0,500])

pmap.xcen = 500.* (pmap.xcen - 14.5)
pmap.ycen = 500.* (pmap.ycen - 14.5)

  for a = 0, n_elements(aorname) - 1 do begin
     if a eq 0 then begin
;        an = plot(bin_snaps[a].xcen, bin_snaps[a].ycen, '6r1s', sym_size = 0.2,   sym_filled = 1, xrange = [14.5, 15.5], yrange = [14.5, 15.5], xtitle = 'X (pixel)', ytitle = 'Y (pixel)', color = colorarr[a], aspect_ratio = 1)
       an = plot(pmap[a].xcen, pmap[a].ycen, '6r1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a],/overplot)
     endif
;
     if a gt 0 then begin
        an = plot(pmap[a].xcen, pmap[a].ycen, '6r1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a],/overplot)
     endif
  endfor

  xsweet = 15.120
  ysweet = 15.085  
  box_x = [xsweet-0.1, xsweet-0.1, xsweet + 0.1, xsweet + 0.1, xsweet -0.1]
  box_y = [ysweet-0.1, ysweet +0.1, ysweet +0.1, ysweet - 0.1,ysweet -0.1]

  
end
