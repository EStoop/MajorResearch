pro virgo_match_catalogs
;measure IRAC colors of
;ICL plume in the virgo cluster
;ps_open, filename='/Users/jkrick/Virgo/IRAC/background.ps',/portrait,/square,/color

openw, outlun, '/Users/jkrick/virgo/irac/matched_ch1_ch2.txt', /get_lun


readcol, '/Users/jkrick/Virgo/irac/ch1/ch1_Combine-mosaic/test.cat', NUMBER_ch1                 , ra_ch1                , dec_ch1                , X_ch1                , Y_ch1                , MAG_AUTO_ch1               , FLUX_BEST_ch1              , FLUXERR_BEST_ch1           , MAG_BEST_ch1               , MAGERR_BEST_ch1            , BACKGROUND_ch1             , FLUX_MAX_ch1               , ISOAREA_IMAGE_ch1          , ALPHA_J2000_ch1            , DELTA_J2000_ch1            , A_IMAGE_ch1                , B_IMAGE_ch1                , THETA_IMAGE_ch1            , MU_THRESHOLD_ch1           , MU_MAX_ch1                 , FLAGS_ch1                  , FWHM_IMAGE_ch1             , CLASS_STAR_ch1             , ELLIPTICITY_ch1  , format = "I10, D10.11, D10.11, F10.5, F10.5,      F10.5, F10.5,    F10.5, F10.5,    F10.5, F10.5,      F10.5, F10.5,    F10.5, F10.5, F10.5, F10.5,      F10.5, F10.5,    F10.5, I10, F10.5, F10.5, F10.5"

readcol, '/Users/jkrick/Virgo/irac/ch2/ch2_Combine-mosaic/test.cat', NUMBER_ch2                 , ra_ch2                , dec_ch2                , X_ch2                , Y_ch2                , MAG_AUTO_ch2               , FLUX_BEST_ch2              , FLUXERR_BEST_ch2           , MAG_BEST_ch2               , MAGERR_BEST_ch2            , BACKGROUND_ch2             , FLUX_MAX_ch2               , ISOAREA_IMAGE_ch2          , ALPHA_J2000_ch2            , DELTA_J2000_ch2            , A_IMAGE_ch2                , B_IMAGE_ch2                , THETA_IMAGE_ch2            , MU_THRESHOLD_ch2           , MU_MAX_ch2                 , FLAGS_ch2                  , FWHM_IMAGE_ch2             , CLASS_STAR_ch2             , ELLIPTICITY_ch2       , format = "I10, D10.11, D10.11, F10.5, F10.5,      F10.5, F10.5,    F10.5, F10.5,    F10.5, F10.5,      F10.5, F10.5,    F10.5, F10.5, F10.5, F10.5,      F10.5, F10.5,    F10.5, I10, F10.5, F10.5, F10.5"      


  
; match ch1 to ch2
   i36=n_elements(x_ch1)
   i36match=fltarr(i36)
   i36match[*]=-999
   print,'Matching '
   print,"Starting at "+systime()
   dist=i36match
   dist[*]=0
   
   for q=long(0),i36-1 do begin
      
      dist=sphdist(ra_ch1[q],dec_ch1[q],ra_ch2,dec_ch2, /degrees)
      sep=min(dist,ind)
      if (sep LT 0.0008) then begin
         i36match[q]=ind
                                ;iracobject[ind].iracmatch = q
                                ;iracobject[ind].iracmatchdist = sep
         printf, outlun, NUMBER_ch1[q]                 , ra_ch1[q]                , dec_ch1[q]                , X_ch1[q]                , Y_ch1[q]                , MAG_AUTO_ch1               , FLUX_BEST_ch1[q]              , FLUXERR_BEST_ch1[q]           , MAG_BEST_ch1[q]               , MAGERR_BEST_ch1[q]            , BACKGROUND_ch1[q]             , FLUX_MAX_ch1[q]               , ISOAREA_IMAGE_ch1[q]          , ALPHA_J2000_ch1[q]            , DELTA_J2000_ch1[q]            , A_IMAGE_ch1[q]                , B_IMAGE_ch1[q]                , THETA_IMAGE_ch1[q]            , MU_THRESHOLD_ch1[q]           , MU_MAX_ch1[q]                 , FLAGS_ch1[q]                  , FWHM_IMAGE_ch1[q]             , CLASS_STAR_ch1[q]             , ELLIPTICITY_ch1[q] ,NUMBER_ch2[ind]                 , ra_ch2[ind]                , dec_ch2[ind]                , X_ch2[ind]                , Y_ch2[ind]                , MAG_AUTO_ch2[ind]               , FLUX_BEST_ch2[ind]              , FLUXERR_BEST_ch2[ind]           , MAG_BEST_ch2[ind]               , MAGERR_BEST_ch2[ind]            , BACKGROUND_ch2[ind]             , FLUX_MAX_ch2[ind]               , ISOAREA_IMAGE_ch2[ind]          , ALPHA_J2000_ch2[ind]            , DELTA_J2000_ch2[ind]            , A_IMAGE_ch2[ind]                , B_IMAGE_ch2[ind]                , THETA_IMAGE_ch2[ind]            , MU_THRESHOLD_ch2[ind]           , MU_MAX_ch2[ind]                 , FLAGS_ch2[ind]                  , FWHM_IMAGE_ch2[ind]             , CLASS_STAR_ch2[ind]             , ELLIPTICITY_ch2[ind]   

      endif

   endfor
   
   print,"Finished at "+systime()
   matched=where(i36match GE 0)
   nonmatched = where(i36match lt 0)
   print, n_elements(matched),"matched"
   print, n_elements(nonmatched),"nonmatched"
   
close, outlun
free_lun, outlun

end
