; PURPOSE:
; This code will calculate the offsets to be entered into Spot fixed
; cluster target offsets when transforming coordinates from an
; observed position using guide-star peak-up to
; the sweet spot position.  This is only useful if you already have a
; guide-star peak-up observation of your target using the SSC
; recommended fixed cluster target offsets which did not land on the sweet spot, and
; you want to transform the observed coordinates to the sweet spot for
; designing another epoch of observations.
;
; EXAMPLE: irac_transformtosweet, 1, 'sub', 15.8, 15.8;
;
;
; INPUTS:
;  CH:  integer; Channel numbner, either 1 or 2
;  FOV: string; Are you observing with 'Full' or 'Sub' array?
;  X1: float; Observed X position of the target in pixels
;  Y1: float; Observed Y position of the target in pixels
;
; OUTPUTS:
;   The row and column offsets in arcseconds to input into Spot as a
;     Fixed Cluster target.
;   Further recommendations for designing sweet spot AORs are at:
;     http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/pcrs_obs.shtml
;
; NOTES: 
;    IDL coordinates start at [0.0, 0.0] for the lower left pixel.
;    Sweet spots are defined as:
;        Ch1 = [15.198,15.020], Ch2 = [15.120,15.085].
;    There is scatter in how well the telescope points to the sweet
;    spot, so the best results here will come from an average [X,Y]
;    centroid positions from your entire AOR of observations.
;
; INFORMATION: For stars which are either too faint or too bright to 
;  use PCRS Peak-Up on the target itself, it is necessary to use a 
;  catalog Peak-Up star.  This is accurate if the proper motion of 
;  the target star and catalog star are both well known.  Proper 
;  motions errors are provided for catalog stars when running 
;  pick_pcrs_catalog (contributed software provided by the SSC), and 
;  usually have low errors.  The only concern with a guide-star is if 
;  there is no acceptable catalog star near your target (this should 
;  be very rare).  If the proper motion of your target star is not 
;  well known, and your science requires the target to be on the 
;  sweet spot, then we recommend that users observe a first epoch of 
;  their catalog + target star pair using PCRS Peak-Up and then 
;  transform the coordinates from that first epoch to move the target 
;  onto the sweet spot for further science epochs.  This code will 
;  perform the coordinate transformation between the observed 
;  coordinates and the desired sweet spot coordinates.  

; MODIFICATION HISTORY:
;  V1: April 2015 JK
;-
pro irac_transformtosweet, CH, FOV, X1, Y1

  ;;do some error checking on the inputs
  if (N_params() lt 3) then begin
     print,'Wrong number of inputs - ' + $
           'irac_transformtosweet, CH, X1, Y1'
     return
  endif 
  if size(CH,/TYPE) gt 5 then begin
     print, 'CH must be an integer'
     return
  endif
  if (CH gt 2) or (CH lt 1) then begin
     print, 'CH must be the integer 1 or 2'
     return
  endif
;  if size(MODE,/TYPE) ne 7 then begin
;     print, 'MODE must be a string'
;     return
;  endif
  if size(FOV,/TYPE) ne 7 then begin
     print, 'FOV must be a string'
     return
  endif


  ;;set up all the input parameters
  case CH of
     1: begin
        linparams = [16.2,-1.09,-26.43,1.76]
        drow0 = [-0.344,130.932,-0.352] ; [full, full], [full, sub], [sub, sub]
        dcol0 = [0.171,127.429,0.064]
        X0 = 15.198             ;sweet spot center
        Y0 = 15.020
     end
     2: begin
        linparams = [23.01,-1.55,-24.60,1.63]
        drow0 = [-0.113, 126.649, -0.511]
        dcol0 = [0.398,124.529, 0.039]
        X0 = 15.120
        Y0 = 15.085
     end
  endcase
  
  ;;figure out which MODE and FOV to use
  ;; we are not using the full array sweet spot at the moment
  if  (FOV eq 'sub') then nmode = 2 
  if  (FOV eq 'full') then nmode = 1 
;  if (MODE eq 'full') and (FOV eq 'full') then nmode = 0
  
  ;;linear transformation
  drow = linparams(1)*(X0-X1) + drow0(nmode)
  dcol = linparams(3)*(Y0-Y1) + dcol0(nmode)
  
  ;;output
  print, 'Cluster offsets should be ', drow, ' ', dcol, format = '(A, F10.3, A, F10.3)'
  
end
