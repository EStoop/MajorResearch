;+
; NAME: IRACPC_PMAP_CORR
;
;
;
; PURPOSE:     Correct IRAC observed aperture photometry-derived fluxes from the post-cryogenic, or warm 
;              mission (IRACPC) for the intra-pixel response, using the pixel maps of the 3.6 micron (channel 1) 
;              and 4.5 micron (channel 2) subarray "sweet spots".  
;              
;             The correction is derived by interpolating high resolution maps of flux as a function of centroid location over the photometric
;              pixel, which is the nominal subarray center pixel.  The quality of the maps varies with position, but in a region of approximately
;              0.5 x 0.5 pixel around the peak of the intrapixel response (the "sweet spot") the pixel gain map ("pmap") is estimated
;              to have formal errors in the gain of better than 1e-3 (where the average gain across the pixel is 1.0).  The effective
;              sweet spot is defined by an "occupation" image which is a weighted sum of the number of photometric data points that
;              contribute to a given grid point in the pmap.
;
; CATEGORY:    Photometry, corrections
;
;
; CALLING SEQUENCE:   corrected_flux = IRACPC_PMAP_CORR ( OBSERVED_FLUX, X, Y, CHANNEL, $
;                                               [,FILE_SUFFIX=file_suffix] [,FILE_DIR=file_dir] [,FLUX_SIGMA=flux_sigma] [,CORRECTED_SIGMA=corrected_sigma] $
;                                               [, /THRESHOLD_OCC ] [,THRESHOLD_VAL=threshold_val] [,/FULL] [,MISSING=missing] [,COMBINE_GAUSS=combine_gauss]
;
; INPUTS:   OBSERVED_FLUX: array of observed aperture fluxes
;                          corresponding to (X,Y)
;                       X: array of X centroids (float), defined such
;                          that the center of the bottom left pixel is
;                          (0.0,0.0).  The bottom left corner of the array is
;                          therefore (-0.5,-0.5) and the top right corner is
;                          (255.5,255.5).  Same number of elements as
;                          OBSERVED_FLUX.
;                       Y: array of Y centroids (float). Same number of
;                          elements as OBSERVED_FLUX.
;                 CHANNEL: IRAC channel number (1 or 2).  Same number of elements
;                          as OBSERVED_FLUX, or a single-element array or scalar
;                          (if only one element, the program assumes all elements of 
;                          OBSERVED_FLUX are for the same channel).
;
;
; OPTIONAL INPUTS: NONE
;
;
; KEYWORD PARAMETERS: 
;           /THRESHOLD_OCC: Set this keyword to use only data in the pixel gain map for which the "occupation" is greater than 
;                           THRESHOLD_VAL.  The provided gain map images were built by gaussian binning photometric measurements in the vicinity
;                           of each bin.  The "occupation" is a map giving the gaussian kernel-weighted number of data points that contribute to
;                           each grid location in the gain map.  You may want to set /THRESHOLD_OCC if the correction seems to add 
;                           noise to your data. We find that formal errors in the gain map are better than 10^-3 when
;                           occupation > 20 (this value of occupation partially defines the "sweet spot").  
;           THRESHOLD_VAL: The minimum occupation to accept.  The default value is THRESHOLD_VAL=20.
;           MISSING:       The value to assign to a photometric data point that has occupancy less than THRESHOLD_VAL (if /THRESHOLD_OCC
;                          is set).  Default is !Values.F_NAN
;           /FULL :        Set this keyword if your observations were taken in Full Array mode.
;           FILE_SUFFIX: Suffix of pixel gain map ("pmap") data files to be used.  For most cases, the default values should not be changed.
;           FILE_DIR: Directory where pixel gain map ("pmap") data files are located.  The default location is the subdirectory "pmap_fits"
;                     of the same directory that holds the iracpc_pmap_corr.pro code.
;           FLUX_SIGMA:      Array of errors on the aperture photometry (input). Should have the same number of elements as OBSERVED_FLUX.
;                            [CURRENTLY INACTIVE]
;           CORRECTED_SIGMA: Array of errors on the corrected fluxes (output). [CURRENTLY INACTIVE]  
;           /COMBINE_GAUSS:  Combine the gain map in the sweet spot with a double-gaussian model for data outside the sweet spot.  
;                            This technique is currently in an experimental stage and we do not recommend using it. [CURRENTLY INACTIVE for CH1]
;
;
; OUTPUTS:  This function returns the corrected flux, which will be either a scalar or 
;           array, depending on the number of elements of OBSERVED_FLUX.
;
;
;
; OPTIONAL OUTPUTS: NONE
;
;
; COMMON BLOCKS: NONE
;
;
;
; SIDE EFFECTS:
;
;
; RESTRICTIONS: This code requires that the pixel map fits files be located on the machine running the program.  By default this is the 
;               subdirectory "pmap_fits/" of the same directory that holds the iracpc_pmap_corr.pro code.    
; 
; The pixel gain map was derived from aperture photometry performed with a 3 pixel
;                aperture (3-7 pixel background annulus), and so is best used to correct data from the 
;                same aperture. (A larger aperture will give a less pronounced variation in pixel phase 
;                response, in which case this program may overcorrect your data.)
;                
;                This correction is only relevant for warm (post-cryogenic) data.
;                
;                As the highest quality portion of the gain map is in regions where occupancy > 20, we do not 
;                vouch for the effectiveness of this correction outside of the "sweet spot".
;
; INSTALLATION: 1) Download the tar.gz file containing the program file the iracpc_pmap_corr.pro and "pmap_fits/" fits file subdirectory.
;               2) Unpack the file in a directory of your choosing
;               3) Ensure that the directory in step (2) is in your IDL Path.
;
; PROCEDURE: 1) Measure X and Y pixel centroid position(s) for stellar image(s) (recommended method: first moment).
;               X and Y are decimal values defined such that the center of the bottom left pixel of the image has (X,Y) = (0.0,0.0). 
;            2) Perform aperture photometry on the observation(s) to derive flux(es) (recommended routine: aper.pro from IDL Astronomy User's Library).  
;               ***NOTE: This correction is optimized for the 3 pixel radius aperture, and 3-7 pixel sky background annulus.
;            3) Correct measured flux(es) using this function. 
;
;
; EXAMPLE:   IDL> corrected_flux = IRACPC_PMAP_CORR( OBSERVED_FLUX, X, Y, CHANNEL, /THRESHOLD_OCC )
;
;
;
; MODIFICATION HISTORY: 2011 12 13 - Created.                                          J. Ingalls (SSC)
;                       2012 02 14 - Added /COMBINE_GAUSS keyword for experimentation  J. Ingalls (SSC)
;                       2012 04 11 - Added CH1 fits file default prefix.               J. Ingalls (SSC)
;                       2012 05 02 - Fixed bug which caused inability to find CH2 Fits file
;                                    Added /FULL keyword for Full Array Observing mode J. Ingalls (SSC)
;                                                
;-
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;Here is interp2d.pro and associated subroutines.  iracpc_pmap_corr will follow afterwards.
;+
; NAME:
;   interp2d
; PURPOSE:
;   Perform bilinear 2d interpolation using the IDL intrinsic 
;   interpolate procedure
; CALLING SEQUENCE:
;   result = interp2d(A,x0,y0,x1,y1)
;   result = interp2d(A,x0,y0,x1,y1,/grid)
;   result = interp2d(A,x0,y0,x1,y1,/regular,/cubic)
;   result = interp2d(A,x0,y0,x1,y1,missing=missing)
; INPUTS:
;   A = 2d array to interpolate
;   x0  = Values that correspond to A(0,0), A(1,0), ...
;   y0  = Values that correspond to A(0,0), A(0,1), ...
;   x1  = New X values at which A should be interpolated
;   y1  = New Y values at which A should be interpolated
; OPTIONAL INPUTS:
;   nxny = [nx,ny] Vector of length 2 which specifies the size of
;         the regular linearized grid produced with trigrid.  The
;         default is nxny = [51,51].  If the size of A is much larger
;   than 51 by 51, greater accuracy may be obtained by having
;         nxny = [n_elements(A(*,0),n_elements(A(0,*))]
; OPTIONAL INPUT KEYWORDS:
;   grid= If set, return an n_elements(X1) by n_elements(y1) grid
;   missing = Value to points which have X1 gt max(X0) or X1 lt min(X0)
;   and the same for Y1.
;   quintic = If set, use smooth interpolation in call to trigrid
;   regular = If set, do not call trigrid -- x0 and y0 must be linear.
;   cubic   = If set, use cubic convolution
;   extrapolate = If set, then extrapolate beyond boundary points
;   bin = set to bin data prior to interpolation.
;         (e.g. bin=2 interpolate every second pixel)
; Returned:
;   result = a vector N_elements(X1) long 
;      or, if /grid is set
;   result = an array that is N_elements(X1) by N_elements(Y1)
;
; PROCEDURE:
;   First call the IDL intrinsic routines TRIANGULATE & TRIGRID to make
; sure that X0 and Y0 are linear (if /regular is not set).
;   Then call the IDL intrinsic INTERPOLATE to do bilinear interpolation.
; RESTRICTIONS:
;   X0 and Y0 must be linear functions.
;   A must be a 2-d array
; HISTORY:
;    9-mar-94, J. R. Lemen LPARL, Written.
;   20-Jan-95, JRL, Added the REGULAR & CUBIC keywords
;   6-Sept-97, Zarro, GSFC, allowed for 2-d (X-Y) coordinate inputs
;  22-Apri-99, Zarro, SM&A/GSFC - added /triangulate and made /regular
;              the default (much faster).
;----------------------------------------------------------------------------
function exist,var

return,n_elements(var) ne 0

end
;
;
; Project     : SOHO - CDS
;
; Name        : WHERE_VECTOR
;
; Purpose     : WHERE function for vectors
;
; Category    : Utility
;
; Explanation :
;
; Syntax      : IDL> ok=where_vector(vector,array,count)
;
; Inputs      : VECTOR = vector with with search elements
;               ARRAY = array to search for each element
;
; Opt. Inputs : None
;
; Outputs     : OK = subscripts of elements in ARRAY that match elements in vector

; Opt. Outputs: COUNT = total # of matches found
;
; Keywords    : TRIM = trim inputs if string inputs
;               CASE = make case sensitive if string inputs 
;               REST = indicies in ARRAY that don't match VECTOR
;               RCOUNT = # of non-matching elements
;               NOSORT = skip sorting input search vector
;
; Common      : None
;
; Restrictions: None
;
; Side effects: None
;
; History     : Version 1,  25-Dec-1995,  D.M. Zarro.  Written
; History     :             03-Feb-2011,  J.G. Ingalls.  changed call to DATATYPE to TYPENAME

function where_vector,vector,array,count,nosort=nosort,$
       trim_string=trim_string,case_sens=case_sens,rest=rest,rcount=rcount

if not exist(vector) or not exist(array) then return,-1
count=0
np=n_elements(array)
rcount=np
rest=lindgen(np)

;-- protect inputs and modify

trim_string=keyword_set(trim_string)
case_sens=keyword_set(case_sens)

svec=vector & sarr=array

if not keyword_set(nosort) then begin
 rs=uniq([svec],sort([svec]))
 svec=svec(rs) 
endif

if typename(vector) eq 'STRING' then begin
 if trim_string then svec=strtrim(svec,2)
 if not case_sens then svec=strupcase(svec)
endif
if typename(array) eq 'STRING' then begin
 if trim_string then sarr=strtrim(sarr,2)
 if not case_sens then sarr=strupcase(sarr)
endif

state=''
nvecs=n_elements(svec)
pieces=strarr(nvecs)
v=svec & s=sarr
for i=0,nvecs-1 do begin
 index=strtrim(string(i),2)
 pieces(i)='(v('+index+') eq s)'
 if i eq 0 then pieces(i)='clook=where('+pieces(i)
 if i eq (nvecs-1) then pieces(i)=pieces(i)+',count)'
 if (nvecs eq 1) or (i eq 0) then conn='' else conn=' or '
 state=state+conn+pieces(i)
endfor

status=execute(strcompress(strtrim(state,2)))

if count gt 0 then begin
 rest(clook)=-1
 rlook=where(rest gt -1,rcount)
 if rcount gt 0 then rest=rest(rlook) else rest=-1
endif

return,clook & end


 function interp2d, A, x0, y0, x1, y1, nxny, missing=missing,   $
  grid=grid, quintic=quintic, regular=regular, cubic=cubic,$
        extrapolate=extrapolate,bin=bin,trigrid=trigrid

;-- Check for self-consistent input dimensions

sz = size(a)

if sz(0) ne 2 then begin
  message,'input data array must be 2-D',/cont
  return,-1
endif

sx=size(x0) & sy=size(y0)
chk=where_vector(sx,sy,count)
if count eq 0 then begin
 message,'input X-Y coordinate arrays do not match in size',/cont
 return,-1
endif

twod=0
case 1 of
 sx(0) eq 1: begin
  nx=sx(1) & ny=sy(1)
 end
 sx(0) eq 2: begin
  twod=1
  nx=sx(1) & ny=sx(2)
 end
 else: begin
  message,'input X-Y coordinate arrays must 1- or 2-D',/cont
  return,-1
 end
endcase
 
if (sz(1) ne nx) or (sz(2) ne ny) then begin
 message,'Dimensions of Data, X0, Y0 are not consistent',/cont
 return,-1
endif

;-- Call triangulate and trigrid to get a regularly spaced grid

if exist(regular) then regular = (0 > regular < 1) else $
 if keyword_set(trigrid) then regular=0 else regular=1

if not keyword_set(regular) then begin
  message,'using trigrid option...',/cont
  if n_elements(nxny) eq 0 then nxny = [nx,ny]

  dprint,'% nx,ny: ',nxny(0),nxny(1)

  gs = [(max(X0)-min(X0))/(nxny(0)-1), (max(Y0)-min(Y0))/(nxny(1)-1)]

  if not twod then begin
   x0 = reform(temporary(x0),nx,ny)
   y0 = transpose(reform(temporary(y0),nx,ny))
  endif 

  triangulate, x0, y0, tr,bound

  if n_elements(quintic) eq 0 then quintic = 0  ; Make sure quintic is defined

  if keyword_set(extrapolate) then begin
   zz = trigrid(x0,y0, A, tr, gs, quintic=quintic,$
                extrapolate=bound)
  endif else begin
   zz = trigrid(x0,y0, A, tr, gs, quintic=quintic)
  endelse

  if not twod then begin
   x0=reform(temporary(x0),nx*ny)
   y0=reform(transpose(temporary(y0)),nx*ny)
  endif

  zz = temporary(zz(0:nxny(0)-1,0:nxny(1)-1)) ; Make sure the dimensions are matched
endif else zz = A       ; /regular was set -- x0 and y0 are linear

sz = size(zz)
xslope = (max(X0)-min(X0)) / (sz(1)-1)
yslope = (max(Y0)-min(Y0)) / (sz(2)-1)


; Map the coordinates

x2 = (x1 - min(x0)) / xslope
y2 = (y1 - min(y0)) / yslope

; Now interpolate

if n_elements(grid)    eq 0 then grid = 0
if n_elements(cubic)   eq 0 then cubic= 0
if n_elements(missing) eq 0 then $
    return,interpolate(zz,x2,y2,grid=grid,cubic=cubic) else $
    return,interpolate(zz,x2,y2,grid=grid,missing=missing,cubic=cubic)
end

;;;;;;;;;;;;;;;;;;;;;;;;HERE IS IRACPC_PMAP_CORR PROPER;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION iracpc_pmap_corr,f0,x0,y0,channel,FILE_SUFFIX=file_suffix,FILE_DIR=file_dir,FLUX_SIGMA=flux_sigma,CORRECTED_SIGMA=corrected_sigma,$
                          THRESHOLD_OCC=threshold_occ, THRESHOLD_VAL=threshold_val, COMBINE_GAUSS=combine_gauss, MISSING=missing,FULL=full

IF KEYWORD_SET(FULL) EQ 1 THEN BEGIN
   xfov = 8.
   yfov = 216.
ENDIF ELSE BEGIN
   xfov = 0.
   yfov = 0.
ENDELSE
IF N_ELEMENTS(FILE_DIR) EQ 0 THEN BEGIN
;; By default the fits files are located in pmap_dir/ under the location of IRACPC_PMAP_CORR (this program)
   this_struc = ROUTINE_INFO('iracpc_pmap_corr',/SOURCE,/FUNCTIONS)   ;;; Get the path to this program
   this_dir = FILE_DIRNAME(this_struc.path[0],/MARK_DIRECTORY)  ;;; Extract directory name
   file_dir = this_dir + 'pmap_fits'
ENDIF
IF N_ELEMENTS(FILE_SUFFIX) EQ 0 THEN file_suffix = ['500x500_0043_120409.fits','0p1s_x4_500x500_0043_120124.fits']
IF KEYWORD_SET(THRESHOLD_OCC) THEN BEGIN
   IF N_ELEMENTS(THRESHOLD_VAL) EQ 0 THEN threshold_val = 20
ENDIF

IF N_ELEMENTS(MISSING) EQ 0 THEN MISSING = !VALUES.F_NAN

IF N_ELEMENTS(f0) GT 1 THEN BEGIN
   corrected_flux= MAKE_ARRAY(SIZE=SIZE(f0),/DOUBLE) 
   occu_interp= MAKE_ARRAY(SIZE=SIZE(f0),/FLOAT)
ENDIF ELSE BEGIN
   corrected_flux = 0.d0
   occu_interp = 0.0
ENDELSE

IF N_ELEMENTS(channel) EQ 1 THEN BEGIN
   chan_uniq = channel[0]
   IF N_ELEMENTS(f0) GT 1 THEN chn = FIX(MAKE_ARRAY(SIZE=size(f0),value=channel[0]) ) ELSE chn = channel[0]
   nch = 1
ENDIF ELSE BEGIN
   iuniq = UNIQ(channel,SORT(CHANNEL))
   chan_uniq = channel[iuniq]
   nch = N_ELEMENTS(chan_uniq)
   chn = channel
ENDELSE

FOR i = 0,nch-1 DO BEGIN
   Index = WHERE(chn EQ chan_uniq[i],n)
   IF N NE 0 THEN BEGIN
      suffix = STRING(chan_uniq[i],format='(i1)')+'_'+file_suffix[chan_uniq[i]-1]
      xgrid = READFITS(file_dir + '/xgrid_ch'+suffix,/SILENT)
      xrange = [MIN(xgrid,max = mx),mx]
      ygrid = READFITS(file_dir + '/ygrid_ch'+suffix,/SILENT)
      yrange = [MIN(ygrid,max = my),my]
      IF keyword_set(COMBINE_GAUSS) THEN pmap = READFITS(file_dir + '/pmap_combined_ch'+suffix,/SILENT) ELSE pmap = READFITS(file_dir + '/pmap_ch'+suffix,/SILENT)
      psigma = READFITS(file_dir + '/psigma_ch'+suffix,/SILENT)
      occu = READFITS(file_dir + '/occu_ch'+suffix,/SILENT)
;;; Need to work on how to incorporate sigmas into the computation.  For now, just interpolate the pmap and divide.
      corrected_flux[index] = f0[Index]/interp2d(pmap,xgrid,ygrid,x0[Index]-xfov,y0[Index]-yfov,/regular)

      occu_interp[index] = interp2d(float(occu),xgrid,ygrid,x0[Index]-xfov,y0[Index]-yfov,/regular)
      IF KEYWORD_SET(THRESHOLD_OCC) THEN BEGIN
         iout = WHERE(occu_interp[Index] LT threshold_val,nout)
         IF nout NE 0 THEN BEGIN
            IF ~keyword_set(COMBINE_GAUSS) THEN corrected_flux[index[iout]] = missing
          ENDIF
      ENDIF
      ;iin = WHERE( x0[index] LE xrange[1] AND x0[index] GE xrange[0] AND y0[index] LE yrange[1] AND y0[index] GE yrange[0], nin, complement=iout, ncomplement=nout)
      ;IF nout NE 0 THEN corrected_flux[index[iout]] = missing
   ENDIF
ENDFOR

RETURN,corrected_flux

END
