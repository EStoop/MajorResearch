PRO pmap_correct_interpolate_single,x,y,ch,func,f_interp,f_interp_unc,xfwhm,yfwhm,NNEAREST=nnearest,POSITION_ONLY=position_only,MAX_DIST=max_dist,$
                                    NP=np
   COMMON pmap_data,x_pmap,y_pmap,f_pmap,np_pmap,xfwhm_pmap,yfwhm_pmap,func_pmap,scale,sigscale
   ;;; Find the nearest neighbors in the pmap dataset
   IF N_ELEMENTS(nnearest) EQ 0 THEN nnearest = 50
   IF N_ELEMENTS(max_dist) EQ 0 THEN max_dist = 0.05
   IF N_ELEMENTS(np) NE 0 THEN BEGIN
      np1d = SQRT(np)
      np_pmap_1d = SQRT(np_pmap)
   ENDIF
   
   ;;; (Squared) Distance of pmap dataset points to current point
   b = FIX(ch) EQ 1 ? 0.8 : 1.0  ;;; Different scaling for y, according to Lewis et al.
   d2_spatial = (x_pmap-x)^2 + ((y_pmap-y)/b)^2
   IF KEYWORD_SET(POSITION_ONLY) THEN BEGIN
      d2 = d2_spatial
   ENDIF ELSE BEGIN
      IF N_ELEMENTS(NP) NE 0 THEN d2 = d2_spatial + (np_pmap_1d-np1d)^2  ELSE d2 = d2_spatial + (xfwhm_pmap-xfwhm)^2 + (yfwhm_pmap-yfwhm)^2
   ENDELSE 
   ;;; Sort by distance
   d_increasing = SORT(d2)
   ;;; Pick the NNEAREST entries
   use_index = d_increasing[0:nnearest-1]
   ;;; Figure out whether we have enough nearby data points 
   IF d2[use_index[nnearest-1]] GT MAX_DIST THEN BEGIN
      f_interp = !values.F_NAN
      f_interp_unc = !values.F_NAN
   ENDIF ELSE BEGIN
   ;;; Compute the spread in the pmap data for the nearest neighbor points 
      sig_x = STDDEV(x_pmap[use_index],/NAN,/DOUBLE)
      sig_y = STDDEV(y_pmap[use_index],/NAN,/DOUBLE)
      occ_vec_spatial = EXP( -(x-x_pmap[use_index])^2/(2*sig_x^2)  $
        -(y-y_pmap[use_index])^2/(2*sig_y^2) )

      IF ~KEYWORD_SET(POSITION_ONLY) THEN BEGIN
         IF N_ELEMENTS(NP) THEN BEGIN
             sig_np_1d = STDDEV(np_pmap_1d[use_index],/NAN,/DOUBLE) 
             occ_vec = occ_vec_spatial * $
               EXP(-(np1d-np_pmap_1d[use_index])^2/(2*sig_np_1d^2))
         ENDIF ELSE BEGIN
             sig_x_fwhm = STDDEV(xfwhm_pmap[use_index],/NAN,/DOUBLE)
             sig_y_fwhm = STDDEV(yfwhm_pmap[use_index],/NAN,/DOUBLE)
             occ_vec = occ_vec_spatial * $
               EXP( -(xfwhm-xfwhm_pmap[use_index])^2/(2*sig_x_fwhm^2) $
                    -(yfwhm-yfwhm_pmap[use_index])^2/(2*sig_y_fwhm^2) )
         ENDELSE
      ENDIF ELSE BEGIN
        occ_vec = occ_vec_spatial
      ENDELSE
   
   ;;; Compute the weighted occupation vector at the data point (for the nearest neighbor points)
      kernel = occ_vec / func_pmap[use_index]^2
   ;;; Compute the weighted average (interpolated) pmap flux at point (X,Y,NP)  
      f_interp = TOTAL(f_pmap[use_index] * kernel,/NAN,/DOUBLE) / TOTAL(kernel,/NAN,/DOUBLE)
      f_interp_unc = SQRT(TOTAL(occ_vec * kernel,/NAN,/DOUBLE) / TOTAL(kernel,/NAN,/DOUBLE)^2)
   ENDELSE
               
RETURN
END
FUNCTION pmap_correct,x,y,f,ch,xfwhm,yfwhm,POSITION_ONLY=position_only,NP=NP,FUNC=func,CORR_UNC=corr_unc,$
                      FULL=full,DATAFILE=datafile,NNEAREST=nnearest,MAX_DIST=max_dist,Verbose=verbose
;;  Correct photometric data using the pmap measurements, with independent variables (X,Y,XFWHM,YFWHM)
;;  Set NNEAREST= the number of nearest neighbors to use in the multidimensional space when computing the correction.  Default is 50.                        
;;  Use "NP=" keyword and input the values of noise pixels to use (X,Y,NP) to determine correction, instead of (X,Y,XFWHM,YFWHM)
;;  Set /POSITION_ONLY to use only X,Y positions instead of (X,Y,XFWHM,YFWHM) or (X,Y,NP)  
;;  Set MAX_DIST= the maximum radial (X,Y) distance of NNEAREST nearest neighbors (in pixels).  Default is 0.05.  If any of the
;;                neighbors is further away, then that value of the corrected flux is set to NaN.  
;;; Common block for pmap data set:
   COMMON pmap_data,x_pmap,y_pmap,f_pmap,np_pmap,xfwhm_pmap,yfwhm_pmap,func_pmap,scale,sigscale
   IF KEYWORD_SET(FULL) EQ 1 THEN BEGIN
      xfov = 8.
      yfov = 216.
   ENDIF ELSE BEGIN
      xfov = 0.
      yfov = 0.
   ENDELSE
   
   IF N_ELEMENTS(datafile) EQ 0 THEN BEGIN
      savfiles = FILE_SEARCH('pmap_data_ch'+STRTRIM(STRING(ch),2)+'*.sav',count=nfiles)
      use_file = savfiles[nfiles-1]   ;;; Choose the latest file
      IF keyword_set(vERBOSE) THEN print,'Using pmap data file '+use_file 
   ENDIF ELSE use_file = datafile
   RESTORE,use_file
   
   ndata = N_ELEMENTS(f)
   
   IF N_ELEMENTS(func) EQ 0 THEN func = REPLICATE(1.0,ndata)
   
   ;; Worry about this later
   ;IF N_ELEMENTS(NP) GT 0 THEN BEGIN
   ;  IF N_ELEMENTS(x) NE ndata OR N_ELEMENTS(y) NE ndata OR N_ELEMENTS(np) NE ndata OR N_ELEMENTS(func) NE ndata THEN BEGIN
   ;    result = DIALOG_MESSAGE('X,Y,F,FUNC,NP must all have the same number of elements.  Exiting.',/ERROR,/CENTER)
   ;    RETURN,-1
   ;  ENDIF
   ;ENDIF ELSE BEGIN
   ;  IF N_ELEMENTS(x) NE ndata OR N_ELEMENTS(y) NE ndata OR N_ELEMENTS(xfwhm) NE ndata OR N_ELEMENTS(yfwhm) NE ndata OR N_ELEMENTS(func) NE ndata THEN BEGIN
   ;    result = DIALOG_MESSAGE('X,Y,F,FUNC,XFWHM,YFWHM must all have the same number of elements.  Exiting.',/ERROR,/CENTER)
   ;    RETURN,-1
   ;  ENDIF
   ;ENDELSE
   
   IF N_ELEMENTS(ch) NE ndata THEN channel = REPLICATE(ch[0],ndata) ELSE channel = ch
      
   f_corr = DBLARR(ndata)
   corr_unc = DBLARR(ndata)
   f_pmap_interp = DBLARR(ndata)
   f_pmap_interp_unc = DBLARR(ndata)
   
;;; Loop through the input data
   IF keyword_set(verbose) THEN Print,'Correcting '+STRN(ndata)+' data points'
   IF KEYWORD_SET(POSITION_ONLY) THEN BEGIN
      FOR i = 0,ndata-1 DO BEGIN
         IF KEYWORD_SET(verbose) THEN Iterwait,i,10,100,Text='Completed '
         pmap_correct_interpolate_single,x[i]-xfov,y[i]-yfov,channel[i],func[i],fi,fiu,NNEAREST=nnearest,MAX_DIST=max_dist,/POSITION_ONLY
         f_pmap_interp[i] = fi
         f_pmap_interp_unc[i]  = fiu
      ENDFOR
   ENDIF ELSE BEGIN
      IF N_ELEMENTS(NP) EQ ndata THEN BEGIN
        FOR i = 0,ndata-1 DO BEGIN
          IF KEYWORD_SET(verbose) THEN Iterwait,i,10,100,Text='Completed '
          pmap_correct_interpolate_single,x[i]-xfov,y[i]-yfov,channel[i],func[i],fi,fiu,NP=np[i],NNEAREST=nnearest,MAX_DIST=max_dist
          f_pmap_interp[i] = fi
          f_pmap_interp_unc[i]  = fiu
        ENDFOR
      ENDIF ELSE BEGIN
         FOR i = 0,ndata-1 DO BEGIN
            IF KEYWORD_SET(verbose) THEN Iterwait,i,10,100,Text='Completed '
            pmap_correct_interpolate_single,x[i]-xfov,y[i]-yfov,channel[i],func[i],fi,fiu,xfwhm[i],yfwhm[i],NNEAREST=nnearest,MAX_DIST=max_dist
            f_pmap_interp[i] = fi
            f_pmap_interp_unc[i]  = fiu
         ENDFOR
      ENDELSE
   ENDELSE
   
;;; Correct the measured fluxes by dividing by the normalized interpolated pmap flux
   f_corr = f * scale /f_pmap_interp
   corr_unc = f_corr * SQRT( (f_pmap_interp_unc/f_pmap_interp)^2 + (func/f)^2 ) 
   
   RETURN,f_corr
END
