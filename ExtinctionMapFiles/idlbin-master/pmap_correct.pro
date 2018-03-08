PRO pmap_correct_interpolate_single,x,y,ch,func,f_interp,f_interp_unc,xfwhm,yfwhm,NNEAREST=nnearest,POSITION_ONLY=position_only,MAX_DIST=max_dist,$
                                    NP=np,OCC_MIN=occ_min,OCCUPATION=occupation
   COMMON pmap_data,x_pmap,y_pmap,f_pmap,np_pmap,xfwhm_pmap,yfwhm_pmap,func_pmap,scale,sigscale
   ;;; Find the nearest neighbors in the pmap dataset
   IF N_ELEMENTS(nnearest) EQ 0 THEN nnearest = 50
   IF N_ELEMENTS(max_dist) EQ 0 THEN max_dist = 0.0025d0
   IF N_ELEMENTS(occ_min) EQ 0 THEN occ_min = 20
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
      occupation = 0d0
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
      occupation = TOTAL(occ_vec_spatial)
      ;;; Figure out whether the constraint on occupation number is satisfied - this gives the effective number of data points 
      ;;; impinging on the value of flux at the point (X,Y)
      IF occ_min GT 0 AND occupation LT occ_min THEN BEGIN
        f_interp = !values.F_NAN
        f_interp_unc = !values.F_NAN
      ENDIF ELSE BEGIN
   ;;; Compute the weighted occupation vector at the data point (for the nearest neighbor points)
         kernel = occ_vec / func_pmap[use_index]^2
   ;;; Compute the weighted average (interpolated) pmap flux at point (X,Y,NP)  
         f_interp = TOTAL(f_pmap[use_index] * kernel,/NAN,/DOUBLE) / TOTAL(kernel,/NAN,/DOUBLE)
         f_interp_unc = SQRT(TOTAL(occ_vec * kernel,/NAN,/DOUBLE) / TOTAL(kernel,/NAN,/DOUBLE)^2)
      ENDELSE
   ENDELSE
               
RETURN
END
FUNCTION pmap_correct,x,y,f,ch,xfwhm,yfwhm,POSITION_ONLY=position_only,NP=NP,FUNC=func,CORR_UNC=corr_unc,$
                      FULL=full,DATAFILE=datafile,NNEAREST=nnearest,MAX_DIST=max_dist,Verbose=verbose,$
                      OCC_MIN=occ_min,OCCUPATION=occupation,USE_PMAP=use_pmap,R_USE=r_use
;;  Correct photometric data using the pmap measurements, with independent variables (X,Y,XFWHM,YFWHM)
;;  Set NNEAREST= the number of nearest neighbors to use in the multidimensional space when computing the correction.  Default is 50.                        
;;  Use "NP=" keyword and input the values of noise pixels to use (X,Y,NP) to determine correction, instead of (X,Y,XFWHM,YFWHM)
;;  Set /POSITION_ONLY to use only X,Y positions instead of (X,Y,XFWHM,YFWHM) or (X,Y,NP)  
;;  Set MAX_DIST= the maximum radial (X,Y) distance of NNEAREST nearest neighbors (in pixels).  Default is 0.0025.  If any of the
;;                neighbors is further away, then that value of the corrected flux is set to NaN.  Set to a large number to turn off this 
;;                constraint.  
;;  Set OCC_MIN = the minimum (spatial) occupation number to be obtained by the NNEAREST points.  This is a further constraint on the proximity
;;                of pmap data points. Default is 20.  Setting keyword OCCUPATION will return the array of occupation values.  Set to 0 to 
;;                turn off this constraint.
;;  Set R_USE = the radius of the photometric aperture used on the dataset.  This will pick out an aperture in the pmap dataset 
;;               (if multiple apertures exist).  
;;  Set USE_PMAP = the indices to the pmap dataset that you would like to use to correct the data.
;;
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
      IF keyword_set(VERBOSE) THEN print,'Using pmap data file '+use_file 
   ENDIF ELSE use_file = datafile
   RESTORE,use_file
   n_aper_pmap = N_ELEMENTS(r_aper)
   IF KEYWORD_SET(R_USE) THEN BEGIN
      iaper_pmap = WHERE(r_use[0] EQ r_aper,/NULL)
      IF ~N_ELEMENTS(iaper_pmap) THEN BEGIN
         PRINT,'Aperture '+STRING(r_use[0],format='(f4.2)')+' px not found in pmap dataset.  Exiting.'
         RETURN,!NULL
      ENDIF
      iaper_pmap = iaper_pmap[0]
   ENDIF ELSE BEGIN
      iaper_pmap = 0
   ENDELSE
   IF N_ELEMENTS(USE_PMAP) NE 0 THEN BEGIN
      x_pmap = x_pmap[use_pmap]
      y_pmap = y_pmap[use_pmap]
      f_pmap = f_pmap[use_pmap,iaper_pmap[0]]
      np_pmap = np_pmap[use_pmap]
      xfwhm_pmap = xfwhm_pmap[use_pmap]
      yfwhm_pmap = yfwhm_pmap[use_pmap]
      func_pmap = func_pmap[use_pmap,iaper_pmap[0]]
      scale = scale[iaper_pmap[0]]
   ENDIF ELSE BEGIN
      f_pmap = f_pmap[*,iaper_pmap]
      func_pmap = func_pmap[*,iaper_pmap]
      scale = scale[iaper_pmap]
   ENDELSE
   
   ndata = N_ELEMENTS(f)
   sz = SIZE(f)
   IF sz[0] EQ 0 THEN sz = [1,1,sz[1],sz[2]]
   IF N_ELEMENTS(func) EQ 0 THEN func = MAKE_ARRAY(SIZE=sz,VALUE=1.0)
   
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
      
   f_corr = MAKE_ARRAY(SIZE=sz,VALUE=!VALUES.F_NAN)
   occupation = MAKE_ARRAY(SIZE=sz,VALUE=!VALUES.F_NAN)
   corr_unc = MAKE_ARRAY(SIZE=sz,VALUE=!VALUES.F_NAN)
   f_pmap_interp = DBLARR(ndata)
   f_pmap_interp_unc = DBLARR(ndata)
   
;;; Loop through the input data
   igood = WHERE(FINITE(F) AND FINITE(x) AND FINITE(y),ngood)
   IF keyword_set(verbose) THEN Print,'Attempting to correct '+STRN(ngood)+' out of '+STRN(ndata)+' data points'
   CASE 1 OF 
      KEYWORD_SET(POSITION_ONLY): BEGIN
         FOR i = 0,ngood-1 DO BEGIN
            IF KEYWORD_SET(verbose) THEN Iterwait,i,100,1000,Text='Completed '
            pmap_correct_interpolate_single,x[igood[i]]-xfov,y[igood[i]]-yfov,channel[igood[i]],func[igood[i]],fi,fiu,NNEAREST=nnearest,$
                                            MAX_DIST=max_dist,/POSITION_ONLY,OCCUPATION=occi,OCC_MIN=occ_min
            f_pmap_interp[igood[i]] = fi
            f_pmap_interp_unc[igood[i]]  = fiu
            occupation[igood[i]] = occi
         ENDFOR
      END
      N_ELEMENTS(NP) EQ ndata: BEGIN
         FOR i = 0,ngood-1 DO BEGIN
            IF KEYWORD_SET(verbose) THEN Iterwait,i,100,1000,Text='Completed '
            pmap_correct_interpolate_single,x[igood[i]]-xfov,y[igood[i]]-yfov,channel[igood[i]],func[igood[i]],fi,fiu,NP=np[igood[i]],$
                                            NNEAREST=nnearest,MAX_DIST=max_dist,OCCUPATION=occi,OCC_MIN=occ_min
            f_pmap_interp[igood[i]] = fi
            f_pmap_interp_unc[igood[i]]  = fiu
            occupation[igood[i]] = occi
         ENDFOR
      END
      ELSE: BEGIN
         FOR i = 0,ngood-1 DO BEGIN
            IF KEYWORD_SET(verbose) THEN Iterwait,i,100,1000,Text='Completed '
            pmap_correct_interpolate_single,x[igood[i]]-xfov,y[igood[i]]-yfov,channel[igood[i]],func[igood[i]],fi,fiu,xfwhm[igood[i]],yfwhm[igood[i]],$
                                            NNEAREST=nnearest,MAX_DIST=max_dist,OCCUPATION=occi,OCC_MIN=occ_min
            f_pmap_interp[igood[i]] = fi
            f_pmap_interp_unc[igood[i]]  = fiu
            occupation[igood[i]] = occi
         ENDFOR
      END
   ENDCASE
   ;;; Report how many values were corrected, if less than the input.
   IF KEYWORD_SET(VERBOSE) THEN BEGIN
      igood_interp = WHERE(FINITE(f_pmap_interp[igood]),ngood_interp)
      IF ngood_interp LT ngood THEN PRINT,'Out of '+STRN(ngood)+' good input data points, '+STRN(ngood_interp)+' were successfully corrected.'
   ENDIF
   
;;; Correct the measured fluxes by dividing by the normalized interpolated pmap flux
   f_corr[igood] = f[igood] * scale /f_pmap_interp[igood]
   corr_unc[igood] = f_corr[igood] * SQRT( (f_pmap_interp_unc[igood]/f_pmap_interp[igood])^2 + (func[igood]/f[igood])^2 ) 
   corr_unc = REFORM(corr_unc)
   occupation = REFORM(occupation)
   RETURN,REFORM(f_corr)
END
