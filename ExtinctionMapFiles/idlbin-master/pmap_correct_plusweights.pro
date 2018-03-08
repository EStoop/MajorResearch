PRO pmap_correct_interpolate_single_plusweights,x,y,f,ch,np,func,f_interp,f_interp_unc,NNEAREST=nnearest, USE_NP = use_np
  COMMON pmap_data_plusweights,x_pmap,y_pmap,f_pmap,np_pmap,func_pmap,scale,sigscale
  ;;; Find the nearest neighbors in the pmap dataset
   IF N_ELEMENTS(nnearest) EQ 0 THEN nnearest = 50
   np1d = SQRT(np)
   np_pmap_1d = SQRT(np_pmap)
   
   ;;; (Squared) Distance of pmap dataset points to current point
   b = FIX(ch) EQ 1 ? 0.8 : 1.0
   d2 = (x_pmap-x)^2 + ((y_pmap-y)/b)^2 + keyword_set(use_np)*(np_pmap_1d-np1d)^2  ; only add the np to distance if called
   ;;; Sort by distance
   d_increasing = SORT(d2)
   ;;; Pick the NNEAREST entries
   use_index = d_increasing[0:nnearest-1]
   ;;;What is the distribution of those neighbors
;   plothist, d2(use_index), xhist, yhist, /noplot,/autobin
;   dist = plot(xhist, yhist, xtitle = 'X, Y, NP distance', /overplot )
   ;;; Compute the spread in the pmap data for the nearest neighbor points 
   sig_x = STDDEV(x_pmap[use_index],/NAN,/DOUBLE)
   sig_y = STDDEV(y_pmap[use_index],/NAN,/DOUBLE)
   sig_np_1d = STDDEV(np_pmap_1d[use_index],/NAN,/DOUBLE)
   
   ;;; Compute the weighted occupation vector at point (X,Y,NP) (for the nearest neighbor points)
   ;notice this uses addition in the weights, not multiplication
   occ_vec = ( EXP(-(x-x_pmap[use_index])^2/(2*sig_x^2)) + $
              EXP(-(y-y_pmap[use_index])^2/(2*sig_y^2)) + $
              keyword_set(use_np)*EXP(-(np1d-np_pmap_1d[use_index])^2/(2*sig_np_1d^2)) $
              ) 
   kernel = occ_vec / func_pmap[use_index]^2
   ;;; Compute the weighted average (interpolated) pmap flux at point (X,Y,NP)  
   f_interp = TOTAL(f_pmap[use_index] * kernel,/NAN,/DOUBLE) / TOTAL(kernel,/NAN,/DOUBLE)
   f_interp_unc = TOTAL(occ_vec * kernel,/NAN,/DOUBLE) / TOTAL(kernel,/NAN,/DOUBLE)^2
               
RETURN
END

FUNCTION pmap_correct_plusweights,x,y,f,ch,np,occdata,FUNC=func,CORR_UNC=corr_unc,FULL=full,DATAFILE=datafile,NNEAREST=nnearest,Verbose=verbose, Use_np = use_np, THRESHOLD_OCC = threshold_occ, THRESHOLD_VAL = threshold_val

;;; Common block for pmap data set:
   COMMON pmap_data_plusweights,x_pmap,y_pmap,f_pmap,np_pmap,func_pmap,scale,sigscale

   if n_elements(threshold_val) eq 0 then threshold_val = 20

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
   IF N_ELEMENTS(x) NE ndata OR N_ELEMENTS(y) NE ndata OR N_ELEMENTS(np) NE ndata OR N_ELEMENTS(func) NE ndata THEN BEGIN
      result = DIALOG_MESSAGE('X,Y,F,FUNC,NP must all have the same number of elements.  Exiting.',/ERROR,/CENTER)
      RETURN,-1
   ENDIF
   
   IF N_ELEMENTS(ch) NE ndata THEN channel = REPLICATE(ch[0],ndata) ELSE channel = ch
      
   f_corr = DBLARR(ndata)
   corr_unc = DBLARR(ndata)
   f_pmap_interp = DBLARR(ndata)
   f_pmap_interp_unc = DBLARR(ndata)
   
;;; Loop through the input data
   IF keyword_set(verbose) THEN Print,'Correcting '+STRN(ndata)+' data points'
   FOR i = 0,ndata-1 DO BEGIN
      IF KEYWORD_SET(verbose) THEN Iterwait,i,10,100,Text='Completed '
      pmap_correct_interpolate_single_plusweights,x[i]-xfov,y[i]-yfov,f[i],channel[i],np[i],func[i],fi,fiu,NNEAREST=nnearest, USE_NP = use_np
      f_pmap_interp[i] = fi
      f_pmap_interp_unc[i]  = fiu
      
;now check to see if we are in a region of the pmap with sufficient data
      IF KEYWORD_SET(THRESHOLD_OCC) THEN BEGIN
         IF (finite(x[i]) lt 1)|| (occdata[x[i]-xfov, y[i]-yfov] EQ !VALUES.F_NAN) THEN f_pmap_interp[i] = !VALUES.F_NAN
      endif

   ENDFOR
   
;;; Correct the measured fluxes by dividing by the normalized
;;; interpolated pmap flux
   f_corr = f * scale /f_pmap_interp
   corr_unc = f_corr * SQRT( (f_pmap_interp_unc/f_pmap_interp)^2 + (func/f)^2 ) 
   



   RETURN,f_corr
END
