pro hd158460_snap_v2
  colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' , ' red',  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'green', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET' , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU',   ' POWDER_BLUE' ]

 ; aorname = ['r44497408','r44497152','r44497664','r44497920','r44498176','r44498432','r44498688','r44498944','r44499200','r44499456','r44499712', 'r44499968','r44500224','r44500480','r44500736', 'r44500992', 'r44501248', 'r44501504' ]

;remove those with no pmap coverage
  ;aorname = ['r44497152','r44497664','r44497920','r44498176','r44498432','r44498688','r44498944','r44499200','r44499456','r44500224','r44500480','r44500736', 'r44500992', 'r44501248', 'r44501504' ]

;44497152,44497664,44497920,44498176,44498432,44498688,44498944,44499200,44499456,44500224,44500480,44500736, 44500992, 44501248, 44501504 

;for the pmap data
;aorname=['r44464128','r44463872','r44463616','r44463360','r44463104','r44462848','r44462592','r44462336','r44462080','r44461824']

;44464128,44463872,44463616,44463360,44463104,44462848,44462592,44462336,44462080,44461824
;these are the staring mode AORs
 ; aorname = ['r42506496','r42051584']

;January 2012 snapshot take two:
aorname = ['0045184256','0045184512','0045184768','0045185024','0045185280','0045185536','0045185792','0045186048','0045186304','0045186560','0045186816' ,'0045187072','0045188352','0045187584','0045187840','0045187328','0045188096','0045188608']



;m = pp_multiplot(multi_layout=[n_elements(aorname),2]);,global_xtitle='Test X axis title',global_ytitle='Test Y axis title')
  ;nfits = 75*63  ; for snapshots
;nfits = 242.*63.
 ;nfits = 1330*63.
nfits = 210*63.
print, 'nfits', nfits
  snapshots = replicate({snapob, ra:0D,dec:0D,xcen:fltarr(nfits),ycen:fltarr(nfits),flux:fltarr(nfits),fluxerr:fltarr(nfits), corrflux:fltarr(nfits), corrfluxerr:fltarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' ', peakpix:fltarr(nfits)},n_elements(aorname))

  ;read in all the pixel phase correction data
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/pmap_ch2_0p4s_500x500_0043_111206.fits', pmapdata, pmapheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/psigma_ch2_0p4s_500x500_0043_111206.fits', psigmadata, psigmaheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/occu_ch2_0p4s_500x500_0043_111206.fits', occdata, occheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/xgrid_ch2_0p4s_500x500_0043_111206.fits', xgriddata, xgridheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/ygrid_ch2_0p4s_500x500_0043_111206.fits', ygriddata, ygridheader

;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p1s_x4_500x500_0043_120124/pmap_ch2_0p1s_x4_500x500_0043_120124.fits', pmapdata, pmapheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p1s_x4_500x500_0043_120124/psigma_ch2_0p1s_x4_500x500_0043_120124.fits', psigmadata, psigmaheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p1s_x4_500x500_0043_120124/occu_ch2_0p1s_x4_500x500_0043_120124.fits', occdata, occheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p1s_x4_500x500_0043_120124/xgrid_ch2_0p1s_x4_500x500_0043_120124.fits', xgriddata, xgridheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p1s_x4_500x500_0043_120124/ygrid_ch2_0p1s_x4_500x500_0043_120124.fits', ygriddata, ygridheader

  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p1s_x4_1000x1000_0020_120124/pmap_ch2_0p1s_x4_1000x1000_0020_120124.fits', pmapdata, pmapheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p1s_x4_1000x1000_0020_120124/psigma_ch2_0p1s_x4_1000x1000_0020_120124.fits', psigmadata, psigmaheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p1s_x4_1000x1000_0020_120124/occu_ch2_0p1s_x4_1000x1000_0020_120124.fits', occdata, occheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p1s_x4_1000x1000_0020_120124/xgrid_ch2_0p1s_x4_1000x1000_0020_120124.fits', xgriddata, xgridheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p1s_x4_1000x1000_0020_120124/ygrid_ch2_0p1s_x4_1000x1000_0020_120124.fits', ygriddata, ygridheader

  for a =0,  n_elements(aorname) - 1 do begin
     print, 'working on ',aorname(a)
     bcddir = '/Users/jkrick/iracdata/flight/IWIC/IRAC029700/bcd/' + string(aorname(a)) 
     CD, bcddir                    ; change directories to the correct AOR directory
     command  =  " ls *bcd_fp.fits > /Users/jkrick/irac_warm/pmap/cbcdlist.txt"
     spawn, command
     command2 =  "ls *sig_dntoflux.fits > /Users/jkrick/irac_warm/pmap/bunclist.txt" ;this really should be the bunc file
     spawn, command2
     
     rawdir = '/Users/jkrick/iracdata/flight/IWIC/IRAC029700/raw/' + string(aorname(a)) 
     CD, rawdir                                                                 ; change directories to the correct AOR directory
     command  =  " ls *mipl.fits > /Users/jkrick/irac_warm/pmap/rawlist.txt" ;this really should be the raw file
     spawn, command
 
     readcol,'/Users/jkrick/irac_warm/pmap/rawlist.txt',rawname, format = 'A', /silent
     readcol,'/Users/jkrick/irac_warm/pmap/cbcdlist.txt',fitsname, format = 'A', /silent
     readcol,'/Users/jkrick/irac_warm/pmap/bunclist.txt',buncname, format = 'A', /silent

     
     yarr = fltarr(n_elements(fitsname)*63)
     xarr = fltarr(n_elements(fitsname)*63)
     timearr = fltarr(n_elements(fitsname)*63)
     fluxarr = fltarr(n_elements(fitsname)*63)
     fluxerrarr = fltarr(n_elements(fitsname)*63)
     corrfluxarr = fltarr(n_elements(fitsname)*63)
     corrfluxerrarr = fltarr(n_elements(fitsname)*63)
     peakarr = fltarr(n_elements(fitsname)*63)
     print, 'xarr', n_elements(xarr)

     cd, bcddir
     for i =0.D, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
      ;  print, 'working on ', fitsname(i)         
                                ;now were did it really point to?
        header = headfits(fitsname(i)) ;
        sclk_obs= sxpar(header, 'SCLK_OBS')
        if i eq 0 then sclk_0 = sxpar(header, 'SCLK_OBS')

        for nt = 0, 63 - 1 do begin
           timearr[i*63+nt] = ((sclk_obs - sclk_0) + 0.4*nt)/60./60.
           timearr[i*63+nt] = ((sclk_obs - sclk_0) + 0.4*nt)/60./60.
        endfor

        if i eq 0 then begin
           ra_ref = sxpar(header, 'RA_REF')
           dec_ref = sxpar(header, 'DEC_REF')
           print, 'ra, dec', ra_ref, dec_ref
        endif
        
        ;do the photometry
        ;get_centroids,fitsname(i), testt, testdt, testx_center, testy_center, testabcdflux, testxs, testys, testfs, testb, /WARM, /APER, APRAD = [3],ra = ra_ref, dec = dec_ref, /silent    

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
  
        ;print,'abcdflux', abcdflux
      ;correct for pixel phase effect based on pmaps from Jim
        ;corrflux = pixel_phase_correct_gauss(abcdflux,x_center,y_center,2,  '3_12_20')
        corrresult = interp2d(pmapdata,xgriddata,ygriddata,x_center,y_center,/regular)
        pmaperr = interp2d(psigmadata,xgriddata,ygriddata,x_center,y_center,/regular)
        occresult = interp2d(occdata,xgriddata,ygriddata,x_center,y_center,/regular)
               
        ;print out some summaries
;        print, 'x, y, flux, correction', mean(x_center), mean(y_center), mean(abcdflux), mean(corrresult), format = '(A, F10.4, F10.4, F10.4, F10.4)'

        ;start by applying the correction to all 65 subarray images per frame
        corrflux = abcdflux / corrresult
        corrfluxerr = sqrt(fs^2 + pmaperr^2)
     
        ;then nan out those which do not have a high enough occupation
        for ncor = 0, n_elements(abcdflux) - 1 do begin
        
           if occresult(ncor) lt 12. then begin ; 1E6 then begin
              ;print, 'out of range', fitsname(i), ncor
              corrflux(ncor) = abcdflux(ncor) / alog10(-1)
              corrfluxerr(ncor)= abcdflux(ncor)/ alog10(-1)
           endif
           
        endfor

        ;----------------------------
        ;add info about the peak pixel flux
        ;print, 'working on raw', rawname(i)
        cd, rawdir
        fits_read, rawname(i), rawdata, rawheader
        mc = reform(rawdata, 32, 32, 64)
        img = mc
; convert to real
        img=img*1.
; fix the problem with unsigned ints
        fix=where((img LT 0),count)
        if (count GT 0) then img[fix]=img[fix]+65536.
; flip the InSb
        img=65535.-img

        peaks = img[15, 16,*]
       ;array of 64 peak pixel values is  img[15,16]
        ;----------------------------
        cd, bcddir
    ; print, 'finished raw'

       if i eq 0 then begin
           xarr = x_center[1:*]
           yarr = y_center[1:*]
           fluxarr = abcdflux[1:*]
           fluxerrarr = fs[1:*]
           corrfluxarr = corrflux[1:*]
           corrfluxerrarr = corrfluxerr[1:*]
           peakarr = peaks[1:*]
        endif else begin
           xarr = [xarr, x_center[1:*]]
           yarr = [yarr, y_center[1:*]]
           fluxarr = [fluxarr, abcdflux[1:*]]
           fluxerrarr = [fluxerrarr, fs[1:*]]
           corrfluxarr = [corrfluxarr, corrflux[1:*]]
           corrfluxerrarr = [corrfluxerrarr, corrfluxerr[1:*]]
           ;print, 'filling peakarr'
           peakarr = [peakarr, peaks[1:*]]
           ;print, 'done filling peakarr'
        endelse
      ;print, 'finished filling arrays'

     endfor                     ; for each fits file in the AOR
     ;fill in that structure
     print,'xarr', n_elements(xarr)
     snapshots[a] ={snapob,  ra_ref,dec_ref,xarr,yarr,fluxarr,fluxerrarr, corrfluxarr, corrfluxerrarr, sclk_0, timearr,aorname(a), peakarr}
     ;print, 'finished filling structure'
  endfor                        ;for each AOR
  
  ;test
  print, snapshots[0].ra, snapshots[0].xcen[26]
  
  save, snapshots, filename='/Users/jkrick/irac_warm/pmap/hd158460_pmap1000.sav'
  undefine, snapshots
  

end
