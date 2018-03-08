;**************************************************************************
;+
;*NAME: 
;    BINS     (General IDL Library 01) 6-JAN-83
;
;*CLASS:
;       measurement
;*CATEGORY
;
;*PURPOSE:  
;    To bin flux data on a specified wavelength grid with or without weights. 
;
;*CALLING SEQUENCE:
;    BINS,WAVE,FLUX,WEIGHT,WCENTER,WIDTHS,WMEAN,WSIGMA,WGT
;
;*PARAMETERS:
;    WAVE    (REQ) (I) (1) (I L F D)
;            Required input vector giving the wavelength scale for the
;            flux data which are to be binned.
;    FLUX    (REQ) (I) (1) (I L F D)
;            Required input vector giving the flux data to be binned.
;
;    WEIGHT  (OPT) (I) (1) (F D)
;            Optional input vector giving point by point weigths for the
;            flux vector. 
;
;    WCENTER (REQ) (I/O) (0 1) (I L F D)
;            Required input scalar or vector giving the centers for
;            the wavelength bin(s). The units MUST be the same as for
;            the WAVE vector. This vector is recomputed for output to the
;            actual value used.
;
;    WIDTHS  (REQ) (I/O) (0 1) (I L F D)
;            Required input scalar or vector giving the bin width in the
;            same units as the WAVE array for each bin.
;            This is recomputed to the actual value used.
;
;    WMEAN   (REQ) (O) (0 1) (F D)
;            Required scalar or vector output variable containing the
;            (weighted) mean value of the FLUX array in each bin.
;
;    WSIGMA  (REQ) (O) (0 1) (F D)
;            Required scalar or vector output variable containing the 
;            (weighted) rms standard deviation from the WMEAN in each bin.   
;
;    WGT     (REQ) (O) (0 1) (F D)
;            Required scalar or vector output variable containing the
;            new weight values in each bin. If weight is not specified,
;            WGT equals the number of points in each bin.
; 
;*EXAMPLES:
;    with WEIGHT assigned by procedure WEIGHT, 
;      BINS,WAVE,FLUX,WEIGHT,WCENTER,WIDTHS,WMEAN,WSIGMA,WGT
;
;    For uniform weighting:
;      BINS,WAVE,FLUX,WCENTER,WIDTHS,WMEAN,WSIGMA,WGT
;
;*SYSTEM VARIABLES USED:
;    !NOPRINT if=0, the results are printed in tabular form
;               >0, the results are not printed.
;
;*INTERACTIVE INPUT:
;    None.
;
;*SUBROUTINES CALLED:
;    TABINV
;    PARSHIFT
;
;*FILES USED:
;    None.
;
;*SIDE EFFECTS:
;    None.
;
;*NOTES:
;        Actually uses (WEIGHT>0) for weights.
;        All inputs may be vector or scalar.
;          If WEIGHT is a scalar then all points have this weight.
;          If WIDTHS is a scalar then all bins have this width.
;        The results are printed out if !NOPRINT = 0.
;        WCENTER and WIDTHS are recomputed to be the actual values used.
;        Typing BINS without parameters will display the procedure call
;          statement.
;
;*PROCEDURE:
;        Treats original data rigorously as binned data. 
;
;*MODIFICATION HISTORY:
;   Jan  6 1983 RJP GSFC Modified BINS to allow weights
;   Apr 17 1987 RWT GSFC Added PARCHECK, made WEIGHT optional
;   May  8 1987 RWT GSFC Fix PARSHIFT
;   Jul 23 1987 RWT GSFC Use !NOPRINT to suppress tabular printout 
;   Aug 25 1987 RWT GSFC correct calculation of N for scalar WCEN,
;                        and add listing of procedure call statement
;   Mar  7 1988 CAG GSFC Add VAX RDAF-style prolog
;   5-23-88 RWT fix scaling factor error in printout format
;-
;*****************************************************************************
PRO BINS,WAVE,FLUX,WEIGHT,WCENTER,WIDTHS,WMEAN,WSIGMA,WGT,flxmin,flxmax,npoints
;
; check parameters & redefine if optional parameters are missing
;
NPARM = N_PARAMS(0)
IF NPARM EQ 0 THEN BEGIN
   
PRINT,'BINS,WAVE,FLUX,WEIGHT,WCENTER,WIDTHS,WMEAN,WSIGMA,WGT,flxmin,flxmax,npoints
   RETALL & END
PARCHECK,NPARM,[10,11],'BINS'
WEI = WEIGHT
WCEN = WCENTER
IF NPARM EQ 11 THEN WID = WIDTHS
IF NPARM EQ 10 THEN BEGIN
   stop
   PARSHIFT,1,WEI,WCEN,WID   ; shift input parameters
   stop
   WEI = 1.0
   END
;
    S =SIZE(WCEN) 
    IF S(0) LT 1 THEN WCEN= FLTARR(1) + WCEN  ;convert scalar to array
    N = FIX(N_ELEMENTS(WCEN))
    SWID=SIZE(WID)
    IF SWID(0) LT 1 THEN WID=WID+FLTARR(N)
;
; RESCALE WAVELENGTH AND FLUX TO AVOID TROUBLE WITH NUMERICAL PRECISION
;
    WSCALE = TOTAL(WCEN) / TOTAL(WCEN*0+1)           ; MEAN WAVELENGTH
    WCENT  = WCEN-WSCALE
    W      = WAVE-WSCALE
    SCALE  = TOTAL(ABS(FLUX)) / TOTAL(FLUX*0+1)            ;NORMALIZATION VALUE
    SCALE  = SCALE+(SCALE EQ 0)
;
; WEIGHT = WEIGHT>0
;
    SWEI   = SIZE(WEI)
    IF (SWEI(0) LT 1) THEN WEI = 0.* WAVE + WEI
    WEI    = WEI>0
;
; SETUP QUANTITIES FOR INTEGRATION
;
    IND    = INDGEN(n_elements(W))
    WDIFF  = (w(IND+1) - w(IND-1)) / 2.
    AREA   = FLUX / SCALE * WDIFF * WEI
    SQRS   = AREA * FLUX / SCALE
;
; CREATE OUTPUT ARRAYS
;
    WMEAN  = FLTARR(N)
    WSIGMA = WMEAN
    WGT    = WMEAN
;
;  COMPUTE INTEGRATION LIMITS
;
;    LOW SIDE
;
    TABINV,W,WCENT-WID/2.,LIN
    IND    = FIX(LIN)  &  R= LIN-IND
    LW     = w(IND)*(1-R) + w(IND+1)*R
    LIND   = FIX(LIN+0.5)
    RL     = LW - (w(LIND) + w(LIND-1)) /2.
    tmp2=flux/scale
    FL     = tmp2(LIND)
    WGTL   = wei(LIND)
;
;   HIGH SIDE
;
    TABINV,W,WCENT+WID/2.,UIN
    IND    = FIX(UIN)  &  R= UIN-IND
    UW     = w(IND)*(1-R) + w(IND+1)*R
    UIND   = FIX(UIN+0.5)
    RU     = (w(UIND) + w(UIND+1)) /2. - UW
    tmp3=flux/scale
    FU     = tmp3(UIND)
    WGTU   = wei(UIND)
;
;  RECOMPUTE ACTUAL AREAS TO BE AVERAGED
;
    WCENT  = (UW+LW) /2.  
    WID = UW - LW
    NPOINTS= UIND - LIND+1
;
;  INTEGRATE BY SUMMATION 
;
    flxmin = wmean*0.
    flxmax = flxmin
    FOR I=0,N-1 DO BEGIN
        WMEAN(I)  = TOTAL(AREA(LIND(I):NPOINTS(I)+lind(i)-1))
        WSIGMA(I) = TOTAL(SQRS(LIND(I):NPOINTS(I)+lind(i)-1))
        TMP = WEI*WDIFF
        WGT(I)    = TOTAL(TMP(LIND(I):NPOINTS(I)+lind(i)-1))
	flxmin(i)     = min(flux(lind(i):npoints(i)+lind(i)-1))
        flxmax(i)     = max(flux(lind(i):npoints(i)+lind(i)-1))
        END  ; INTEGRATION BY SUMMATION
;
;  CORRECT FOR END EFFECTS
;
    ENDS   = WGTL * RL + WGTU * RU                         ; NON-NEGATIVE
    WGT    = WGT - ENDS
    WMEAN  = WMEAN - RL * FL * WGTL - RU * FU * WGTU
    WSIGMA = WSIGMA - FL * FL * RL * WGTL - FU * FU * RU * WGTU
;
;  COMPUTE STATISTICS
;
    WSIGMA = WSIGMA * WGT - WMEAN * WMEAN
    PWGT   = WGT + (WGT LE 0)                    ; INSURE POSITIVE FOR DIVISION
    WSIGMA = SQRT(WSIGMA +(WSIGMA LE 0)) * (WSIGMA GT 0)/ PWGT*SCALE*(WGT GT 0)
;
;          NOTE THAT WE AVOID ERROR SQRT(NEG OR ZERO)
;          WHICH CAN OCCUR WHEN NPOINTS IS SMALL
;          DUE TO THE NUMERICAL APPROXIMATIONS
;
    WMEAN  = WMEAN / PWGT * SCALE * (WGT GT 0)
    WGT    = WGT/(WID+(WID LE 0))*(WID GT 0)*NPOINTS $ ; AVG WT * NPTS
      * WGT/ (WGT + ENDS + (WGT LE 0))       ; ACCOUNTS ROUGHLY FOR END EFFECTS
;
; AVG WT * NPOINTS
;
; PRINTOUT OPTION
;
;IF !NOPRINT EQ 0 THEN BEGIN ; PRINTOUT
;  PRINT,'                            WEIGHTED VALUES'
;  IF NPARM EQ 11 THEN BEGIN
;     PRINT,'   WCENTER  WIDTHS    MEAN        SIGMA      WEIGHT'   
;     FOR I=0,N-1 DO $
;      PRINT,"$(1X,F10.2,F7.2,1P2E12.3,0PF10.2)",WCENT(I)+WSCALE,WID(I), $
;        WMEAN(I),WSIGMA(I),WGT(I)
;  END ELSE BEGIN
;     PRINT,'   WCENTER  WIDTHS    MEAN        SIGMA     NPNTS' 
;     FOR I=0,N-1 DO $
;      PRINT,"$(1X,F10.2,F7.2,1P2E12.3,I6)",WCENT(I)+WSCALE,WID(I), $
;        WMEAN(I),WSIGMA(I),NPOINTS(I)
;     END
;  END ; PRINTOUT
;
;  CONVERT OUTPUT TO SCALAR IF NECESSARY
;
    IF S(0) LT 1 THEN BEGIN ; CONVERT OUTPUT TO SCALAR
        WCENT= WCENT(0) & WID= WID(0)
        WMEAN= WMEAN(0) & WSIGMA= WSIGMA(0) & WGT=WGT(0)
        END
;
; RESCALE WAVELENGH SCALE
;
    WCEN=WCENT+WSCALE
;
; reset input parameters if WEIGHT was not specified
;
    IF NPARM EQ 10 THEN BEGIN
      
PARSHIFT,-1,WEI,WCEN,WID,WMEAN,WSIGMA,WGT
      WGT = NPOINTS
      END
WCENTER = WCEN
WIDTHS = WID
;
RETURN ; BINS
END
;**********************************************************************
;+
;*NAME:
;    TABINV     (General IDL Library 01) JUNE 15, 1981
; 
;*CLASS:
;    interpolation
;
;*CATEGORY:
;
;*PURPOSE:  
;    To find the effective index of a function value in the
;    domain of a vector array representing the function.
; 
;*CALLING SEQUENCE: 
;    TABINV,XARR,X,IEFF
; 
;*PARAMETERS:
;    XARR   (REQ) (I) (1) (I L F D)
;           Required input vector containing the data to be searched
;
;    X      (REQ) (I) (0 1) (I L F D)
;           Required input scalar or vector specifying the value(s)
;           whose effective index is sought.
;
;    IEFF   (REQ) (O) (0 1) (F)
;           Required output scalar or vector giving the effective
;           index or indices of X in XARR.
; 
; 
;*EXAMPLE:
;        Set all flux values of a spectrum (WAVE vs FLUX) to zero
;        for wavelengths less than 1150 Angstroms.
; 
;           TABINV,WAVE,1150.0,I
;           FOR J=0,FIX(I) DO FLUX(I)=0.
;
;*SYSTEM VARIABLES USED:
;      !ERR
;
;*INTERACTIVE INPUT:
;      None
;
;*SUBROUTINES CALLED:
;      PARCHECK
;
;*FILES USED:
;      None.
;
;*SIDE EFFECTS:
;      The procedure will abort with an error message if XARR is non-monotonic.
;      See the NOTES for more information.
;
;*RESTRICTIONS:
;      None
;
;*NOTES:
;      -  The binary search technique used in TABINV requires that
;        the input array XARR is monotonic (increasing or decreasing).
;        This means that input vectors with padded zeroes will
;        cause the routine to abort.
;      -  Note that IUE wavelength arrays which include the 
;        vacuum-to-air correction may not be monotonic.
;      -  If X is outside the range of XARR, IEFF is set to the 
;        indice of the closest endpoint (i.e. IEFF = IEFF >0.0 <NPT).
; 
;*PROCEDURE: 
;           The input array XARR is first tested for monotonicity
;        by evaluating the array XARR - SHIFT(XARR,1). If all
;        but one value are positive the array is considered to be
;        monotonicly increasing; if all but one are negative
;        the array is monotonicly decreasing. (Because the shifts are
;        circular, the difference between the first point and the 
;        last point is ignored.) Any other result will cause
;        TABINV to abort.
;           A binary search is then used to find the values XARR(I)
;        and XARR(I+1) where XARR(I) LE X LT XARR(I+1).
;        IEFF is then computed using linear interpolation 
;        between I and I+1.
; 
;                IEFF = I + (X-XARR(I)) / (XARR(I+1)-XARR(I))
; 
;        Let N = number of elements in XARR
; 
;                if x < XARR(0) then IEFF is set to 0
;                if x > XARR(N-1) then IEFF is set to N-1
; 
;*MODIFICATION HISTORY:
;      Feb. 17 1980   D.J. Lindler   initial program
;      Apr. 21 1982   FHS3    GSFC   CR#039 - monotonic functions.
;      Jan. 17 1984   RWT     GSFC   copied into IUER_GL01
;      Mar.  7 1984   RWT     GSFC   test for monotonicity added 
;                                    (see SMR #34)
;      Mar.  9 1984   RWT     GSFC   monotonicity test improved (now
;                                    aborts when 2 or more points have the
;                                    same value in XARR).
;      Mar. 15 1984   RWT     GSFC   more efficient coding of monotonicity test.
;      Jan   5 1988   RWT     GSFC   modify for large arrays (i.e. use longword
;                                    integers), use IDL SHIFT command to test
;                                    monotonicity, and add procedure call
;                                    listing.
;      Mar. 21 1988   CAG     GSFC   add VAX RDAF-style prolog, and PARCHECK.
;-
;***********************************************************************
PRO TABINV,XARR,X,IEFF
IF N_PARAMS(0) EQ 0 THEN BEGIN
   PRINT,' TABINV,XARR,X,IEFF'
   RETALL & END
PARCHECK,N_PARAMS(0),3,'TABINV'
;
; GET SIZE OF THE ARRAYS
;
  S=SIZE(XARR)
  SX=SIZE(X)
  NPOINTS= S(1) & NPT= NPOINTS-1
;
; MAKE X INTO A VECTOR IF IT IS A SCALAR
;
  IF SX(0) EQ 0 THEN BEGIN
     XSAMP=FLTARR(1)+X
     NX=1   
     END ELSE BEGIN
     XSAMP=X
     NX= SX(1)
  END; IF
;
; INITIALIZE BINARY SEARCH AREA AND COMPUTE NUMBER OF DIVISIONS NEEDED
;
  ILEFT=LONARR(NX) & IRIGHT=ILEFT
  NDIVISIONS = FIX( ALOG10(NPOINTS)/ALOG10(2.0)+1.0 )
;
; TEST FOR MONOTONICITY (CLEVER)
;
  I = XARR - SHIFT(XARR,1)
  A = WHERE(I GT 0)
  IF (!ERR EQ NPT) THEN BEGIN  ; increasing array
     IRIGHT = IRIGHT + NPT
  END ELSE BEGIN
     A = WHERE(I LT 0)  ; test for decreasing array
     IF (!ERR EQ NPT) THEN ILEFT=ILEFT+NPT ELSE BEGIN
        PRINT,'ABORTING TABINV: XARR VECTOR NOT MONOTONIC' & RETALL & END
     END          ; monotonicity 
;
; PERFORM BINARY SEARCH BY DIVIDING SEARCH INTERVAL IN
; HALF NDIVISION TIMES
;
  FOR I=1,NDIVISIONS DO BEGIN
    IDIV=(ILEFT+IRIGHT)/2   ;SPLIT INTERVAL IN HALF
    XVAL=xarr(IDIV)  ;FIND FUNCTION VALUES AT CENTER
    GREATER= LONG(XSAMP GT XVAL)  ;DETERMINE WHICH SIDE XSAMP IS ON
    LESS   = LONG(XSAMP LE XVAL)  
    ILEFT=ILEFT*LESS + IDIV*GREATER ;COMPUTE NEW SEARCH AREA
    IRIGHT=IRIGHT*GREATER + IDIV*LESS
    END; FOR I LOOP
;
; INTERPOLATE BETWEEN INTERVAL OF WIDTH = 1
;
  XLEFT =xarr(ILEFT)        ;VALUE ON LEFT SIDE
  XRIGHT=xarr(IRIGHT)        ;VALUE ON RIGHT SIDE
  IEFF=(XRIGHT-XSAMP)*ILEFT+(XSAMP-XLEFT)*IRIGHT + ILEFT*(XRIGHT EQ XLEFT)
  IEFF=IEFF/(XRIGHT-XLEFT+(XRIGHT EQ XLEFT)) ;INTERPOLATE
  IEFF=IEFF>0.0 < NPT               ;DO NOT ALLOW EXTRAPOLTION BEYOND ENDS
  IF SX(0) EQ 0 THEN IEFF=IEFF(0)   ;IF NEEDED CHANGE TO SCALAR
  RETURN
  END

;*******************************************************************
;+
;*NAME:
;    PARSHIFT     23 MAY 1987
;
;*CLASS:
;    i/o
;
;*CATEGORY:
; 
;*PURPOSE:  
;    TO PROPERLY ASSIGN PARAMETERS SPECIFIED IN A
;    procedure call statement containing unspecified optional
;    parameters.
; 
;*CALLING SEQUENCE: 
;    PARSHIFT,N,P1,P2,...P10
; 
;*PARAMETERS:
;       N   (REQ) (I) (0 1) (I)
;           A scalar or vector describing which of the parameters
;           P1 through P10 were unspecified in the procedure in which
;           PARSHIFT is executed. If the first element is > 0, then
;           each parameter is assigned the value of the closest lower
;           parameter not flagged as unspecified, and the unspecified
;           parameter(s) are set to 0. If the first element is < 0, then
;           each parameter is assigned the value of the closest higher
;           parameter not flagged as unspecified, and the right most
;           parameter(s) are set to 0.
; 
;       P1,...P10  (REQ) (I/O) (0 1 2) (B I L F D S)
;           Parameters specified in the procedure call statement
;           in the procedure in which PARSHIFT is executed which are to be
;           reassigned because of unspecified optional parameters. Since only
;           parameters to the right of the first unspecified optional
;           parameter need to be redefined, P1 is usually used as the first
;           unspecified optional parameter. The parameters must be specified
;           in the same order as the inital procedure call statement. Any
;           number of parameters from 2 to 10 can be specified.
; 
;*EXAMPLES:
;   MYPRO.PRO  allows 6 parameters (4 input and 2 output) of
;   which the 2nd & 3rd are optional (e.g. PRO MYFILE,A,B,C,D,E,F).
;   To allow the mode MYPRO,A,D,E,F the following lines can be used
;   at the beginning of the procedure :
;        IF (N_PARAMS(0) EQ 4) THEN BEGIN
;             C = 0          ;the 3rd & 4th parameters are now output
;             D = 0          ;parameters and therefore need to be defined
;             PARSHIFT,[1,2],B,C,D
;   Since A does not change it doesn't need to be specified in the call
;   to PARSHIFT. Also, since the output parameters  E & F are not defined,
;   they can not be included. The above line will shift parameter D
;   to the correct parameter position for the procedure. At the end of
;   the procedure, add the command:
;        IF (N_PARAMS(0) EQ 4) THEN PARSHIFT,[-1,2],B,C,D,E,F
;   to shift the parameters back to the positions they were originally
;   specified in.
;
;       Procedure RV has 4 parameters H,W,F,E of which H and E are
;       optional. All parameters are input parameters. To allow the
;       mode RV,W,F  PARSHIFT is executed at the beginning and end of
;       RV as follows:
;       beginning of procedure -
;       PAR = N_PARAMS(0)
;       IF PAR EQ 2 THEN BEGIN
;          PARSHIFT,1,H,W,F   ; E is undefined, don't include
;          E = 0              ; define E here (H already set to 0)
;          END
;       end of procedure - 
;       IF PAR EQ 2 THEN PARSHIFT,-1,H,W,F  
;       Without reassigning the parameters, RV (in the above example)
;       would interpret W as H (since H is the first parameter in the
;       procedure call statement) and F as W. The parameters must
;       be switched back to their original positions however before
;       the procedure ends.
; 
;       Consider MYPRO,IN1,IN2,OPIN1,OPIN2,IN3,OUT1,OUT2 which has
;       5 input parameters, 2 output parameters, and the third and
;       fourth input parameters are optional. Note that since the
;       total number of parameters is used to determine the mode
;       it is not generally possible to distinguish running
;                 MYPRO,IN1,IN2,OPIN1,IN3,OUT1,OUT2  from
;                 MYPRO,IN1,IN2,OPIN2,IN3,OUT1,OUT2  
;       and therefore can not be allowed.  
;       To allow the mode MYPRO,IN1,IN2,IN3,OUT1,OUT2 do the
;       following: 
;             at start of MYPRO -
;       IF (N_PARAMS(0) EQ 5) THEN BEGIN
;             OPIN2=0     ;not yet defined since the 4th param. is now output
;             IN3=0       ;also not yet defined but not really necessary here 
;             PARSHIFT,[1,2],OPIN1,OPIN2,IN3
;             END
;             at end of MYPRO -
;       IF N_PARAMS(0) EQ 5 THEN PARSHIFT,[-1,2],OPIN1,OPIN2,IN3,OUT1,OUT2
;
;*SYSTEM VARIABLES USED:
;    None
;
;*INTERACTIVE INPUT:
;    None
;
;*SUBROUTINES CALLED:
;    PARCHECK
;
;*FILES USED:
;    None
;
;*SIDE EFFECTS:
;    None
;
;*RESTRICTIONS:
;    None
; 
;*NOTES:
;    - PARSHIFT is intended to be used at the beginning of a
;    procedure with N > 0, and only the input parameters specified.
;    An error will occur if any undefined parameter is used. This
;    includes all output parameters AND any optional input parameters 
;    that are unspecified or reassigned to be output parameters.
;    If the optional parameters are all output parameters, PARSHIFT
;    doesn't need to be run till the end of the procedure as 
;    described below.
;    - PARSHIFT should be used at the end of a procedure with
;    the first element of N < 0, and with all input and output 
;    parameters (starting with the first unspecified optional
;    parameter) included in the parameter list.
;    - PARSHIFT is not necessary for procedures in which the
;    optional parameters occur at the end of the parameter list
;    since no required parameters need to be redefined.
;    - If N is a scalar, it is converted to a 1-element array.
;    - If the procedure aborts between the 1st and last execution of PARSHIFT,
;    the parameters will not be reassigned as they were originally specified.
;    Its a good idea to assign the input parameters to temporary variables
;    and use these as input to PARSHIFT.
; 
;*PROCEDURE: 
;       After converting scalar values of N to 1-elemet arrays, PARSHIFT
;       checks the sign of the first element. If > 0, the parameters are
;       shifted to the right (e.g., P10 = P9, P9 = P8, ...) for as many
;       elements as contained in N. The parameters specified by N are then
;       set to 0. If the first element is < 0, then the parameters are
;       shifted to the left (e.g., P1 = P2, P2 = P3, ...). The right most
;       parameters are then set to zero for each element of N (e.g., if
;       P1 to P6 are specified and N has 2 elements, P6 and P5 =0).
;
;
;*MODIFICATION HISTYORY:
;       Apr 17 1987  RWT  GSFC   initial program 
;       Mar 15 1988  CAG  GSFC   add VAX RDAF-style prolog
;-
;****************************************************************************
PRO PARSHIFT,N,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10
;
NPAR = N_PARAMS(0)
IF NPAR EQ 0 THEN BEGIN
    PRINT,' PARSHIFT,N,P1,P2,p3,p4,p5,p6,p7,p8,p9,p10'
    RETALL & END
PARCHECK,NPAR,[3,4,5,6,7,8,9,10,11],'PARSHIFT'
S = SIZE(N)
IF S(0) EQ 0 THEN BEGIN
   N = INTARR(1) + N                 ; N was scalar
   NOPT = 0
   END ELSE NOPT = (S(1) > 1) - 1    ; NOPT = # of optional parameters -1
;
; If N(0) > 0, shift parameters according to procedure call statement
;
IF N(0) GT 0 THEN BEGIN
    ; 
    ; shift NOPT parameters & set unspecified parameters to 0
    ;
    FOR I=0,NOPT DO BEGIN            ;redefine parameters
      IF (N(I)+1 GT NPAR) THEN BEGIN
         PRINT,'Error specifying parameters for PARSHIFT'
         PRINT,'NPAR = ',NPAR,'  N = ',N(I)
         RETALL & END
      IF (NPAR EQ 11) AND (11 GT N(I)) THEN P10 = P9
      IF (NPAR GE 10) AND (10 GT N(I)) THEN P9 = P8
      IF (NPAR GE 9) AND (9 GT N(I))   THEN P8 = P7
      IF (NPAR GE 8) AND (8 GT N(I))   THEN P7 = P6
      IF (NPAR GE 7) AND (7 GT N(I))   THEN P6 = P5
      IF (NPAR GE 6) AND (6 GT N(I))   THEN P5 = P4
      IF (NPAR GE 5) AND (5 GT N(I))   THEN P4 = P3
      IF (NPAR GE 4) AND (4 GT N(I))   THEN P3 = P2
      IF (NPAR GE 3) AND (3 GT N(I))   THEN P2 = P1
      END   ;redefining parameters               
    FOR I=0,NOPT DO BEGIN          ;set to 0
      IF N(I) EQ 10 THEN P10 = 0
      IF N(I) EQ 9  THEN P9  = 0
      IF N(I) EQ 8  THEN P8  = 0
      IF N(I) EQ 7  THEN P7  = 0
      IF N(I) EQ 6  THEN P6  = 0
      IF N(I) EQ 5  THEN P5  = 0
      IF N(I) EQ 4  THEN P4  = 0
      IF N(I) EQ 3  THEN P3  = 0
      IF N(I) EQ 2  THEN P2  = 0
      IF N(I) EQ 1  THEN P1  = 0
      END       ;setting parameters to 0
    END  $    ; N(0) > 0
 ELSE BEGIN     ; N(0) < 0
;
; If N(0) < 0, shift parameters back based on what user specified
;
    N = ABS(N)
    FOR I=0,NOPT DO BEGIN            ;redefine parameters
      IF (N(I)+1 GT NPAR) THEN BEGIN
           PRINT,'Error specifying parameters for PARSHIFT'
           PRINT,'NPAR = ',NPAR,'  N = ',N(I)
           RETALL & END
      IF (NPAR GE 3) AND (3 GT N(I))   THEN P1 = P2
      IF (NPAR GE 4) AND (4 GT N(I))   THEN P2 = P3
      IF (NPAR GE 5) AND (5 GT N(I))   THEN P3 = P4
      IF (NPAR GE 6) AND (6 GT N(I))   THEN P4 = P5
      IF (NPAR GE 7) AND (7 GT N(I))   THEN P5 = P6
      IF (NPAR GE 8) AND (8 GT N(I))   THEN P6 = P7
      IF (NPAR GE 9) AND (9 GT N(I))   THEN P7 = P8
      IF (NPAR GE 10) AND (10 GT N(I)) THEN P8 = P9
      IF (NPAR EQ 11) AND (11 GT N(I)) THEN P9 = P10
      END   ;redefining parameters               
    FOR I=0,NOPT DO BEGIN          ;set to 0
      IF (N(I) GE 1) AND (NPAR EQ 11)  THEN P10 = 0
      IF (N(I) LE 2) AND (NPAR LE 10)  THEN P9  = 0
      IF (N(I) LE 3) AND (NPAR LE 9)   THEN P8  = 0
      IF (N(I) LE 4) AND (NPAR LE 8)   THEN P7  = 0
      IF (N(I) LE 5) AND (NPAR LE 7)   THEN P6  = 0
      IF (N(I) LE 6) AND (NPAR LE 6)   THEN P5  = 0
      IF (N(I) LE 7) AND (NPAR LE 5)   THEN P4  = 0
      IF (N(I) LE 8) AND (NPAR LE 4)   THEN P3  = 0
      IF (N(I) LE 9) AND (NPAR LE 3)   THEN P2  = 0
      IF (N(I) LE 10) AND (NPAR LE 2)  THEN P1  = 0
      END       ;setting parameters to 0
    END     ; N(0) < 0
RETURN
END

;*********************************************************************
;+
;
;*NAME:
;    PARCHECK     (General IDL Library 01)  30-MAR-1987
;
;*CLASS:
;    Error Checking
;
;*CATEGORY:
; 
;*PURPOSE:
;    To check that a procedure has been called with the minimum of allowed
;    number of parameters. 
; 
;*CALLING SEQUENCE:
;    PARCHECK,NPARM,MINPARMS,CALLINGPRO
; 
;*PARAMETERS:
;    NPARM       (REQ) (I) (0) (I)
;                Required input scalar giving the number of parameters
;                in the procedure call (i.e. N_PARAMS(0)).
;    MINPARMS    (REQ) (I) (0 1) (I)
;                If scalar, the minimum number of parameters needed for the
;                procedure to execute properly.
;                If an array, it represents the allowed numbers of
;                parameters (e.g. if 3,4 or 6 parameters are allowed,
;                then set MINPARMS([0,1,2]) = [3,4,6] ).
;    CALLINGPRO  (REQ) (I) (0) (S)
;                Required string giving the name of the calling procedure.
;
;*EXAMPLES:
;     Determine if procedure PRO, which contains a
;     call to PARCHECK has the minimum number of parameters
;     (i.e. 4):
;             PARCHECK,N_PARAMS(0),4,'PRO'
;     If the same procedure can have 4,5,7, or 8 parameters
;     then use:
;             PARCHECK,N_PARAMS(0),[4,5,7,8],'PRO'
;
;*SYSTEM VARIABLES USED:
;     !ERR
;
;*INTERACTIVE INPUT:
;     None
;
;*SUBROUTINES CALLED:
;     PCHECK
;
;*FILES USED:
;     None 
;
;*SIDE EFFECTS:
;     None
;
;*RESTRICTIONS:
;     None
;
;*NOTES:
;     None
;
;*PROCEDURE:
;     The input parameters to PARCHECK are first checked themselves
;     using PCHECK. If MINPARMS is a scalar it is compared to NPARM.
;     If NPARM < MINPARMS, then an error message is printed and the
;     procedure returns to the main level. If MINPARMS is a vector,
;     then NPARM is subtracted from each value of MINPARMS and the
;     resulting vector is checked for zeroes. If no values are zero,
;     then error messages are printed and the program returns to the
;     main level.
; 
;*MODIFICATION HISTORY :
;     Mar 30 1987    CAG  GSFC  initial program
;     Apr    1987    RWT  GSFC  add vector input for parameters
;     Mar 15 1988    CAG  GSFC  add VAX RDAF-style prolog
;-
;***********************************************************************
PRO PARCHECK,NPARM,MINPARMS,CALLINGPRO
; check input parameters
;
PCHECK,NPARM,1,100,0111
PCHECK,MINPARMS,2,110,0111
PCHECK,CALLINGPRO,3,100,1000
CPRO = STRUPCASE(CALLINGPRO)
;
; check if number of parameters is > minimum (if MINPARMS is scalar)
;
S = SIZE(MINPARMS)
IF S(0) EQ 0 THEN BEGIN
  IF NPARM LT MINPARMS THEN BEGIN
     PRINT,"$(1X,3A,I2,A)",'Procedure ',CPRO,  $
       ' needs at least ',FIX(MINPARMS),' parameters to execute.'
     PRINT,"$(1X,A,I2,A)",'Only ',FIX(NPARM),' parameters were specified'
     PRINT,'ACTION: RETURNING TO THE MAIN PROGRAM LEVEL'
     RETALL
     END
;
; check if NPARM is an allowed number of parameters (if MINPARMS is a vector)
;
  END ELSE BEGIN
     IND = WHERE(MINPARMS-NPARM)   ; !ERR will be # of non-zero values
     IF S(1) EQ !ERR THEN BEGIN
        PRINT,'Invalid number of input parameters for procedure ',CPRO
        PRINT,'ACTION: RETURNING TO THE MAIN PROGRAM LEVEL'
        RETALL
        END
     END
RETURN
END

;*************************************************************************
;+
;
;*NAME:   
;    PCHECK    (General IDL Library 01)   DEC. 1980
;
;*CLASS:
;    Error Checking
;
;*CATEGORY:
;
;*PURPOSE: 
;     TO CHECK THE PARAMETERS OF A PROCEDURE FOR VALID TYPES.
;
;*CALLING SEQUENCE:
;     PCHECK,VARIABLE,POSITION,DIMENSION,TYPE
;
;*PARAMETERS: 
;        VARIABLE  (REQ) (I) (0 1 2) (B I L F D S)
;                  Required input parameter containing the variable 
;                  to be checked.
;        POSITION  (REQ) (I) (0)  (I)
;                  Required input scalar giving the variable position
;                  in the calling sequence of the procedure
;        DIMENSION (REQ) (I) (0)  (I)
;                  Required input scalar 3 digit integer giving the
;                  valid dimensions for VARIABLE. Each digit is either 0 or
;                  1, with 1 indicating that the dimension is valid for
;                  VARIABLE.
;                  1st digit - scalar
;                  2nd digit - one dimensional array (vector)
;                  3rd digit - two dimensional array
;        TYPE      (REQ) (I) (0) (I)
;                  Required input integer scalar giving the valid data
;                  types for VARIABLE. Each digit is either 0 or 1, with
;                  1 indicating that the type is valid.
;                        1st digit - string
;                        2nd       - byte
;                        3rd       - integer or longword integer
;                        4th       - floating point, double precision or
;                                    complex data type
;
;*EXAMPLES:
;        For a procedure DOIT requiring parameters;
;          WAVE- one dimensional floating point array
;          FLUX- one dimensional floating point or integer array
;          K - scalar (floating point, integer, or byte)
;        use:
;            PRO DOIT,WAVE,FLUX,K,FLUXOUT
;            ;
;            ;    initial comments
;            ;
;            PCHECK,WAVE,FLUX,1,010,0001
;            PCHECK,FLUX,2,010,0110
;            PCHECK,K,3,100,0111
;               ...
;
;*SYSTEM VARIABLES USED:
;      None
;
;*INTERACTIVE INPUT:
;      None
;
;*SUBROUTINES CALLED:
;      None
;
;*FILES USED:
;      None
;
;*SIDE EFFECTS:
;      None
;
;*RESTRICTIONS:
;      None
;
;*NOTES:
;         
;        To exit from the error condition it is necessary to enter either
;        a .CONT or a RETALL after the error message has been printed.
;        
;        This routine should be used at the begining of all user 
;        procedures
;
;*PROCEDURE:
;        If VARIABLE is an improper type or undefined, an appropriate
;        message is printed and the procedure stops.  The stop gives a 
;        indication of the line number and procedure in which the
;        error occurred.
;        
;*MODIFICATION HISTORY:
;        Dec    1980  D. Lindler       initial program
;        Nov 17 1981  FHS3       GSFC  indicate the calling procedure
;                                      name on an error condition (CR 016).
;        Jan 17 1984  RWT        GSFC  copied into IUER_GL01
;        Oct 24 1985  RWT        GSFC  modify to accept new data types
;                                      longword integer, double precision,
;                                      and complex. Modified to be compatible
;                                      with existing version (e.g. no new
;                                      input parameters).
;        Mar 15 1988  CAG        GSFC  add VAX RDAF-style prolog.
;-
;**************************************************************************
PRO PCHECK,VARIABLE,POSITION,DIMENSION,TYPE
;
; DECODE VALID DIMENSIONS
;
  VDIM=BYTARR(3)
  IDIM=DIMENSION
  FOR I=0,2 DO BEGIN
    ITEMP=IDIM/10*10
    VDIM(2-I)=IDIM-ITEMP
    IDIM=IDIM/10
  END
;
; DECODE VALID TYPES
;
  ITYPE=FIX(TYPE)
  VTYPE=BYTARR(4)
  FOR I=0,3 DO BEGIN
    ITEMP=ITYPE/10*10
    VTYPE(I)=ITYPE-ITEMP
    ITYPE=ITYPE/10
  END
;
; DETERMINE TYPE AND DIMENSION OF THE INPUT VARIABLE
;
  S=SIZE(VARIABLE)
  NDIM=S(0)
  NTYPE=S(NDIM+1)
;
; CHECK IF IT IS DEFINED
;
  IF NTYPE EQ 0 THEN BEGIN
  PRINT,' INPUT VARIABLE NUMBER',POSITION,' IS NOT DEFINED'
  GOTO,EXPLAIN
  END
;
; CHECK FOR VALID TYPE
;
  IF (NTYPE EQ 7 ) AND ( VTYPE(3) NE 0 ) THEN GOTO,CHECKD ;string
  IF (NTYPE EQ 1 ) AND ( VTYPE(2) NE 0 ) THEN GOTO,CHECKD ;byte
  IF (NTYPE EQ 2 ) AND ( VTYPE(1) NE 0 ) THEN GOTO,CHECKD ;integer
  IF (NTYPE EQ 3 ) AND ( VTYPE(1) NE 0 ) THEN GOTO,CHECKD ;longword integer
  IF (NTYPE EQ 4 ) AND ( VTYPE(0) NE 0 ) THEN GOTO,CHECKD ;floating pt.
  IF (NTYPE EQ 5 ) AND ( VTYPE(0) NE 0 ) THEN GOTO,CHECKD ;double precision
  IF (NTYPE EQ 6 ) AND ( VTYPE(0) NE 0 ) THEN GOTO,CHECKD ;complex
;
  PRINT,'THE PARAMETER IN POSITION',POSITION,' IS OF INVALID TYPE'
  GOTO,EXPLAIN
;
;
; CHECK FOR A VALID DIMENSION
;
  CHECKD:  IF VDIM(NDIM) NE 0 THEN RETURN
  PRINT,'INPUT VARIABLE IN POSITION',POSITION,' HAS WRONG DIMENSION'
;
; GIVE CORRECT TYPE(S) AND DIMENSIONS
;
  EXPLAIN:
  TMSG='     THE VALID TYPE(S) ARE: '
  IF VTYPE(3) NE 0 THEN TMSG=TMSG+' STRING '
  IF VTYPE(2) NE 0 THEN TMSG=TMSG+' BYTE '
  IF VTYPE(1) NE 0 THEN TMSG=TMSG+' INTEGER '
  IF VTYPE(0) NE 0 THEN TMSG=TMSG+' FLOATING POINT '
  PRINT,TMSG
;
  TMSG='  THE VALID DIMENSIONS ARE:'
  IF VDIM(2) NE 0 THEN TMSG=TMSG+' 2-D '
  IF VDIM(1) NE 0 THEN TMSG=TMSG+' 1-D '
  IF VDIM(0) NE 0 THEN TMSG=TMSG+' SCALAR'
  PRINT,TMSG
  PRINT,'TYPE .CON OR RETALL TO CONTINUE'
  STOP
  RETALL
END


pro qmset,num,pos,noeras

npar=n_params(0)
if npar lt 1 then num=1
if npar lt 2 then pos='T'
if npar lt 3 then noeras=0

!NOERAS=NOERAS
!COLOR=2

if num eq 1 then set_viewport,.24,.9,.20,.88

if num eq 2 then begin
CASE STRUPCASE(POS) OF
 'T':      set_viewport,.24,.90,.6,.9
 'B':      set_viewport,.24,.90,.15,.45
 'L':      set_viewport,.26,.52,.23,.875
 'R':      set_viewport,.64,.90,.23,.875
ELSE: BEGIN
       PRINT,'INVALID PAGE POSITION CODE FOR N=2: RETURNING'
       RETURN
      END
ENDCASE
END ;NUM=4

if num eq 3 then begin
CASE STRUPCASE(POS) OF
 'T': SET_VIEWPORT,.235,.905,.695,.875
 'M': SET_VIEWPORT,.235,.905,.445,.625
 'B': SET_VIEWPORT,.235,.905,.195,.375
ELSE: BEGIN
       PRINT,'INVALID PAGE POSITION CODE FOR N=3: RETURNING'
       RETURN
      END
ENDCASE
end ; num=3

if num eq 4 then begin
CASE STRUPCASE(POS) OF
 'UL':      set_viewport,.19,.53,.575,.875
 'UR':      set_viewport,.59,.93,.575,.875
 'LL':      set_viewport,.19,.53,.20,.50
 'LR':      set_viewport,.59,.93,.20,.50
ELSE: BEGIN
       PRINT,'INVALID PAGE POSITION CODE FOR N=4: RETURNING'
       RETURN
      END
ENDCASE
END ;NUM=4

RETURN
END

pro newfield,fdescp,field


vrybng: field=' '
fdescp=' '
print,' '
esc=string(27b)
print,esc,'[2J'
print,'$(6(/),a)',' '
print,'$(5x,a,t65,a)','Attitude Control Field','enter a'
print,'$(5x,a,t65,a)','Command And Data Handling Field','enter c'
print,'$(5x,a,t65,a)','Power Field','enter p'
print,'$(5x,a,t65,a)','Thermal Field','enter t'
print,'$(5x,a,t65,a)','Dewar Field','enter d'
print,'$(5x,a,t65,a)','Transmitter Field','enter x'
print,'$(5x,a,t65,a)','Exit from Routine','enter e'
print,'$(6(/),a)',' '
ans=get_kbrd(1)
ans=strupcase(ans)
case ans of
'A': goto,acs
'C': goto,cdh
'P': goto,powr
'T': goto,therm
'D': goto,dewar
'X': goto,trans
'E': goto,vrynd
else: begin
	print,'  Not a choice, try again'
	wait,2
	goto,vrybng
	end
endcase

acs: attcf,fdescp,field
if strlen(field) lt 10 then goto,vrybng else goto,nd

cdh: comdhf,fdescp,field
if strlen(field) lt 10 then goto,vrybng else goto,nd

powr: power,fdescp,field
if strlen(field) lt 10 then goto,vrybng else goto,nd

therm: thermal,fdescp,field
if strlen(field) lt 10 then goto,vrybng else goto,nd


dewar: print,esc,'[2J'
begin
print,'$(5x,a,t65,a)','WAPCOHXT  <4-7>SCU HRM # 9  APCO HEAT EXCHANGER 4 MSB','enter a'
print,'$(5x,a,t65,a)',"WIVCSHDT  <4-7>SCU HRM #24  INNER VAPOR COOLED 4 MSB'S",'enter b'
print,'$(5x,a,t65,a)',"WPPINTT   <4-7>SCU HRM #15  POROUS PLUG 4 MSB'S",'enter c' 
print,'$(5x,a,t65,a)',"Return to Main Menu",'enter d'
ans=get_kbrd(1)
ans=strupcase(ans) 
case ans of
'A': begin
	field='nsb_scdb_eus.wapcohxt'
	fdescp='WAPCOHXT  <4-7>SCU HRM # 9  APCO HEAT EXCHANGER 4 MSB'
	end
'B': begin
	field='nsb_scdb_eus.wivcshdt'
	fdescp="WIVCSHDT  <4-7>SCU HRM #24  INNER VAPOR COOLED 4 MSB'S"
	end
'C': begin
	field='nsb_scdb_eus.wppintt'
	fdescp="WPPINTT   <4-7>SCU HRM #15  POROUS PLUG 4 MSB'S"
	end
'D': goto,vrybng
else: begin
	print,' Not a choice, try again'
	wait,1
	goto, dewar
	end
endcase
end
goto, nd

trans: print,esc,'[2J'
begin
print,'$(5x,a,t65,a)','X1RFPOUT  XP1 RF POWER OUTPUT','enter a'
print,'$(5x,a,t65,a)','X2RFPOUT  XP2 RF POWER OUT','enter b'
print,'$(5x,a,t65,a)','XM1ST     < 3 >XMTR 1 ON (TRANSMITTING) 0=OFF 1=ON','enter c'
print,'$(5x,a,t65,a)','XM2ST     < 3 >XMIT 2 ON (TRANSMITTING) 0=OFF 1=ON','enter d'
print,'$(5x,a,t65,a)','Return to Main Menu','enter e'
ans=get_kbrd(2)
ans=strupcase(ans)
case ans of
'A': begin
	field='nsb_scdb_eus.x1rfpout'
	fdescp='X1RFPOUT  XP1 RF POWER OUTPUT'
	end
'B': begin
	field='nsb_scdb_eus.x2rfpout'
	fdescp='X2RFPOUT  XP2 RF POWER OUT'
	end
'C': begin
	field='nsb_scdb_eus.xm1st'
	fdescp='XM1ST     < 3 >XMTR 1 ON (TRANSMITTING) 0=OFF 1=ON'
	end
'D': begin
	field='nsb_scdb_eus.xm2st'
	fdescp='XM2ST     < 3 >XMIT 2 ON (TRANSMITTING) 0=OFF 1=ON'
	end
'E': goto,vrybng
else: begin
	print,' Not a choice, try again'
	wait,1
	goto, trans
	end
endcase
end
nd:
print,esc,'[2J'
print,'$(6(/),a)',' '
print,'$(5x,a)','You have selected the following field'
print,'$(5x,a,2x,a)',field,fdescp
print,'$(5x,a,t25,a)','Is this correct?','<cr> means yes, any other key means no'
print,'$(6(/),a)',' '
ans=get_kbrd(1)
ans=strupcase(ans)
ret=string("15b)
if ans eq ret then goto,vrynd else goto,vrybng
vrynd:
return
end

pro attcf,fdescp,field
field=' '
fdescp=' '
vrybng: print,' '
esc=string(27b)
print,esc,'[2J'
begin
print,'$(5(/),a)',' '
print,'$(5x,a,t65,a)','Earth Sensor/Gyro/Current','enter a'
print,'$(5x,a,t65,a)','Wheel/Sun Sensor/TAM','enter b'
print,'$(5(/),a)',' '
case get_kbrd(1) of
'a': goto, egc
'b': goto, wst
else: begin
	print, ' Not a choice, try again'
	wait,1
	goto, vrybng
	end
endcase
end

egc: print,esc,'[2J'
begin
print,'$(5x,a,t65,a)','AESARPER  MODULO 64 EARTH SENSOR A-AXIS RESOLVED','enter a'
print,'$(5x,a,t65,a)','AESBRPER  MODULO 64 EARTH SENSOR  B-AXIS RESOLVED','enter b'     
print,'$(5x,a,t65,a)','AESCRPER  MODULO 64 EARTH SENSOR  C-AXIS RESOLVED','enter c'     
print,'$(5x,a,t65,a)','AGYABPT   GYRO A (1) BASEPLATE TEMP','enter d'
print,'$(5x,a,t65,a)','AGYAORB   GYRO A  (A-AXIS) ORBIT RATE BIAS','enter e'
print,'$(5x,a,t65,a)','AGYARRER  MODULO 64 GYRO A  (A-AXIS) RESOLVED RATE','enter f' 
print,'$(5x,a,t65,a)','AGYBBPT   GYRO B (3) BASEPLATE TEMP','enter g'
print,'$(5x,a,t65,a)','AGYBORB   GYRO B  (B-AXIS) ORBIT RATE BIAS','enter h'
print,'$(5x,a,t65,a)','AGYBRRER  MODULO 64 GYRO B  (B-AXIS) RESOLVED RATE','enter i' 
print,'$(5x,a,t65,a)','AGYCBPT   GYRO C (5) BASEPLATE TEMP','enter j'
print,'$(5x,a,t65,a)','AGYCORB   GYRO C   (C-AXIS) ORBIT RATE BIAS','enter k' 
print,'$(5x,a,t65,a)','AGYCRRER  MODULO 64 GYRO C  (C-AXIS) RESOLVED RATE','enter l'
print,'$(5x,a,t65,a)','AMABI     MODULO 64 MMA A BAR CURRENT','enter m'
print,'$(5x,a,t65,a)','AMACE1ER  ACE CE 1 SPIN RATE COMMAND','enter n'     
print,'$(5x,a,t65,a)','AMACE2ER  ACE CE 2 SPIN RATE ERROR','enter o'      
print,'$(5x,a,t65,a)','AMAXBI    MODULO 64 MMA AX BAR CURRENT','enter p'        
print,'$(5x,a,t65,a)','AMBBI     MODULO 64 MMA B BAR CURRENT','enter q'    
print,'$(5x,a,t65,a)','AMBXBI    MODULO 64 MMA BX BAR CURRENT','enter r'     
print,'$(5x,a,t65,a)','AMCBI     MODULO 64 MMA C BAR CURRENT','enter s'    
print,'$(5x,a,t65,a)','AMCXBI    MODULO 64 MMA CX BAR CURRENT','enter t'     
print,'$(5x,a,t65,a)','Return to Main Menu','enter u'
ans=get_kbrd(1)
ans=strupcase(ans)     
case ans of
'A': begin
	field='nsb_scdb_eus.aesarper'
	fdescp='AESARPER  MODULO 64 EARTH SENSOR A-AXIS RESOLVED'
	end
'B': begin
	field='nsb_scdb_eus.aesbrper'
	fdescp='AESBRPER  MODULO 64 EARTH SENSOR B-AXIS RESOLVED'
	end
'C': begin
	field='nsb_scdb_eus.aescrper'
	fdescp='AESCRPER  MODULO 64 EARTH SENSOR C-AXIS RESOLVE'
	end
'D': begin
	field='nsb_scdb_eus.agyabpt'
	fdescp='AGYABPT  GYRO A (1) BASEPLATE TEMP'
	end
'E': begin
	field='nsb_scdb_eus.agyaorb'
	fdescp='AGYAORB   GYRO A  (A-AXIS) ORBIT RATE BIAS'
	end
'F': begin
	field='nsb_scdb_eus.agyarrer' 
	fdescp='AGYARRER  MODULO 64 GYRO A  (A-AXIS) RESOLVED RATE'
	end 
'G': begin
	field='nsb_scdb_eus.agybbpt'
	fdescp='AGYBBPT   GYRO B (3) BASEPLATE TEMP'
	end
'H': begin
	field='nsb_scdb_eus.agyborb' 
	fdescp='AGYBORB   GYRO B  (B-AXIS) ORBIT RATE BIAS'
	end
'I': begin
	field='nsb_scdb_eus.agybrrer' 
	fdescp='AGYBRRER  MODULO 64 GYRO B  (B-AXIS) RESOLVED RATE'
	end
'J': begin
	field='nsb_scdb_eus.agycbpt' 
	fdescp='AGYCBPT   GYRO C (5) BASEPLATE TEMP'
	end
'K': begin
	field='nsb_scdb_eus.agycorb'
	fdescp='AGYCORB   GYRO C   (C-AXIS) ORBIT RATE BIAS'
	end
'L': begin
	field='nsb_scdb_eus.agycrrer'
	fdescp='AGYCRRER  MODULO 64 GYRO C  (C-AXIS) RESOLVED RATE'
	end
'M': begin
	field='nsb_scdb_eus.amabi'
	fdescp='AMABI     MODULO 64 MMA A BAR CURRENT'
	end
'N': begin
	field='nsb_scdb_eus.amace1er'
	fdescp='AMACE1ER  ACE CE 1 SPIN RATE COMMAND'
	end
'O': begin
	field='nsb_scdb_eus.amace2er'
	fdescp='AMACE2ER  ACE CE 2 SPIN RATE ERROR'
	end
'P': begin
	field='nsb_scdb_eus.amaxbi'
	fdescp='AMAXBI    MODULO 64 MMA AX BAR CURRENT'
	end
'Q': begin
	field='nsb_scdb_eus.ambbi'
	fdescp='AMBBI     MODULO 64 MMA B BAR CURRENT'
	end
'R': begin
	field='nsb_scdb_eus.ambxbi'
	fdescp='AMBXBI    MODULO 64 MMA BX BAR CURRENT'
	end
'S': begin
	field='nsb_scdb_eus.amcbi'
	fdescp='AMCBI     MODULO 64 MMA C BAR CURRENT'
	end
'T': begin
	field='nsb_scdb_eus.amcxbi'
	fdescp='AMCXBI    MODULO 64 MMA CX BAR CURRENT'
	end
'U': goto, nd
else: begin
	print,' Not a choice, try again'
	wait,1
	goto, egc
	end
endcase
end
goto,nd

wst: print,esc,'[2J'
begin
print,'$(5x,a,t65,a)','AMW1INTG  MWEA 1 WHEEL INTEGRATOR','enter a'
print,'$(5x,a,t65,a)','AMW1SPER  MWEA 1 WHEEL SPEED ERROR INTEGRATORS','enter b'
print,'$(5x,a,t65,a)','AMW1SRV   MWEA 1 SERVO','enter c'
print,'$(5x,a,t65,a)','AMW1WSP   MWEA 1 WHEEL SPEED','enter d'         
print,'$(5x,a,t65,a)','AMW2INTG  MWEA 2 WHEEL INTEGRATOR ','enter e'         
print,'$(5x,a,t65,a)','AMW2SPER  MWEA 2 WHEEL SPEED ERROR INTEGRATORS','enter f'
print,'$(5x,a,t65,a)','AMW2SRV   MWEA 2 SERVO','enter g'
print,'$(5x,a,t65,a)','AMW2WSP   MWEA 2 WHEEL SPEED','enter h'           
print,'$(5x,a,t65,a)','ARWASP    MODULO 32 RWA A SPEED','enter i'           
print,'$(5x,a,t65,a)','ARWBSP    MODULO 32 RWA B SPEED','enter j'           
print,'$(5x,a,t65,a)','ARWCSP    MODULO 32 RWA C SPEED','enter k'           
print,'$(5x,a,t65,a)','ASSCEER   MODULO 64 SUN SENSOR CORRECTED ELEVATION','enter l'
print,'$(5x,a,t65,a)','ASSEOBS   SUN SENSOR ELEVATION OFFSET BIAS','enter m'
print,'$(5x,a,t65,a)','ATAMA     MODULO 32 TAM A  FIELD','enter n'   
print,'$(5x,a,t65,a)','ATAMAP    MODULO 32 TAM AP FIELD  (PERPENDICULAR)','enter o'
print,'$(5x,a,t65,a)','ATAMAX    MODULO 32 TAM AX FIELD','enter p'
print,'$(5x,a,t65,a)','ATAMB     MODULO 32 TAM B  FIELD','enter q'
print,'$(5x,a,t65,a)','ATAMBP    MODULO 32 TAM BP FIELD','enter r' 
print,'$(5x,a,t65,a)','ATAMBX    MODULO 32 TAM BX FIELD','enter s' 
print,'$(5x,a,t65,a)','ATAMC     MODULO 32 TAM C  FIELD','enter t'
print,'$(5x,a,t65,a)','ATAMCP    MODULO 32 TAM CP FIELD','enter u'     
print,'$(5x,a,t65,a)','ATAMCX    MODULO 32 TAM CX FIELD','enter v'
Print,'$(5x,a,t65,a)','Return to Main Menu','enter w'
ans=get_kbrd(1)
ans=strupcase(ans)
case ans of
'A': begin
	field='nsb_scdb_eus.amw1intg'
	fdescp='AMW1INTG  MWEA 1 WHEEL INTEGRATOR'
	end
'B': begin
	field='nsb_scdb_eus.amw1sper'
	fdescp='AMW1SPER  MWEA 1 WHEEL SPEED ERROR INTEGRATORS'
	end
'C': begin
	field='nsb_scdb_eus.amw1srv'
	fdescp='AMW1SRV   MWEA 1 SERVO'
	end
'D': begin
	field='nsb_scdb_eus.amw1wsp'
	fdescp='AMW1WSP   MWEA 1 WHEEL SPEED'
	end
'E': begin
	field='nsb_scdb_eus.amw2intg'
	fdescp='AMW2INTG  MWEA 2 WHEEL INTEGRATOR '
	end
'F': begin
	field='nsb_scdb_eus.amw2sper'
	fdescp='AMW2SPER  MWEA 2 WHEEL SPEED ERROR INTEGRATORS'
	end
'G': begin
	field='nsb_scdb_eus.amw2srv'
	fdescp='AMW2SRV   MWEA 2 SERVO'
	end
'H': begin
	field='nsb_scdb_eus.amw2wsp'
	fdescp='AMW2WSP   MWEA 2 WHEEL SPEED'
	end
'I': begin
	field='nsb_scdb_eus.arwasp'
	fdescp='ARWASP    MODULO 32 RWA A SPEED'
	end
'J': begin
	field='nsb_scdb_eus.arwbsp'
	fdescp='ARWBSP    MODULO 32 RWA B SPEED'
	end
'K': begin
	field='nsb_scdb_eus.arwcsp'
	fdescp='ARWCSP    MODULO 32 RWA C SPEED'
	end
'L': begin
	field='nsb_scdb_eus.assceer'
	fdescp='ASSCEER   MODULO 64 SUN SENSOR CORRECTED ELEVATION'
	end
'M': begin
	field='nsb_scdb_eus.asseobs'
	fdescp='ASSEOBS   SUN SENSOR ELEVATION OFFSET BIAS'
	end
'N': begin
	field='nsb_scdb_eus.atama'
	fdescp='ATAMA     MODULO 32 TAM A  FIELD'
	end
'O': begin
	field='nsb_scdb_eus.atamap'
	fdescp='ATAMAP    MODULO 32 TAM AP FIELD  (PERPENDICULAR)'
	end
'P': begin
	field='nsb_scdb_eus.atamax'
	fdescp='ATAMAX    MODULO 32 TAM AX FIELD'
	end
'Q': begin
	field='nsb_scdb_eus.atamb'
	fdescp='ATAMB     MODULO 32 TAM B  FIELD'
	end
'R': begin
	field='nsb_scdb_eus.atambp'
	fdescp='ATAMBP    MODULO 32 TAM BP FIELD'
	end
'S': begin
	field='nsb_scdb_eus.atambx'
	fdescp='ATAMBX    MODULO 32 TAM BX FIELD'
	end
'T': begin
	field='nsb_scdb_eus.atamc'
	fdescp='ATAMC     MODULO 32 TAM C  FIELD'
	end
'U': begin
	field='nsb_scdb_eus.atamcp'
	fdescp='ATAMCP    MODULO 32 TAM CP FIELD'
	end
'V': begin
	field='nsb_scdb_eus.atamcx'
	fdescp='ATAMCX    MODULO 32 TAM CX FIELD'
	end
'W': goto,nd
else: begin
	print,' Not a choice, try again'
	wait,1
	goto, wst
	end
endcase
end
nd:

return
end
pro power,fdescp,field
field=' '
fdescp=' '

esc=string(27b)

power: print,esc,'[2J'
begin
print,'$(5x,a,t65,a)','PEBI    ESSENTIAL BUS CURRENT','enter a'
print,'$(5x,a,t65,a)','PEBVC   MODULO 04 ESSENTIAL BUS VOLT 0-40 VOLTS','enter b'
print,'$(5x,a,t65,a)','PEBVF1  MODULO 04 ESSENTIAL BUS VOLT 27-29 VOLTS','enter c'
print,'$(5x,a,t65,a)','PEBVF2  MODULO 04 ESSENTIAL BUS VOLT 27-29 VOLTS','enter d'
print,'$(5x,a,t65,a)','PEBVF3  MODULO 04 ESSENTIAL BUS VOLT 27-29 VOLTS','enter e'
print,'$(5x,a,t65,a)','PNEBI   NON-ESSENTIAL BUS CURRENT','enter f'
print,'$(5x,a,t65,a)','PSAI    SOLAR ARRAY CURRENT','enter g'
print,'$(5x,a,t65,a)','PSDD1T  SHUNT DISSIPATOR PANEL D #1 TEMP.','enter h'
print,'$(5x,a,t65,a)','PSDD2T  SHUNT DISSIPATOR PANEL D #2 TEMP.','enter i'
print,'$(5x,a,t65,a)','PSDET   SHUNT DISSIPATOR PANEL E TEMP.','enter j'
print,'$(5x,a,t65,a)','PSDFT   SHUNT DISSIPATOR PANEL F TEMP.','enter k'
print,'$(5x,a,t65,a)','PSDGT   SHUNT DISSIPATOR PANEL G TEMP.','enter l'
print,'$(5x,a,t65,a)','PSDHT   SHUNT DISSIPATOR PANEL H TEMP.','enter m'
print,'$(5x,a,t65,a)','PSHI    SHUNT CURRENT','enter n'
print,'$(5x,a,t65,a)','Return to Main Menu','enter o'
ans=get_kbrd(1)
ans=strupcase(ans)
case ans of
'A': begin
	field='nsb_scdb_eus.pebi'
	fdescp='PEBI    ESSENTIAL BUS CURRENT'
	end
'B': begin
	field='nsb_scdb_eus.pebvc'
	fdescp='PEBVC   MODULO 04 ESSENTIAL BUS VOLT 0-40 VOLTS'
	end
'C': begin
	field='nsb_scdb_eus.pebvf1'
	fdescp='PEBVF1  MODULO 04 ESSENTIAL BUS VOLT 27-29 VOLTS'
	end
'D': begin
	field='nsb_scdb_eus.pebvf2'
	fdescp='PEBVF2  MODULO 04 ESSENTIAL BUS VOLT 27-29 VOLTS'
	end
'E': begin
	field='nsb_scdb_eus.pebvf3'
	fdescp='PEBVF3  MODULO 04 ESSENTIAL BUS VOLT 27-29 VOLTS'
	end
'F': begin
	field='nsb_scdb_eus.pnebi'
	fdescp='PNEBI   NON-ESSENTIAL BUS CURRENT'
	end
'G': begin
	field='nsb_scdb_eus.psai'
	fdescp='PSAI    SOLAR ARRAY CURRENT'
	end
'H': begin
	field='nsb_scdb_eus.psdd1t'
	fdescp='PSDD1T  SHUNT DISSIPATOR PANEL D #1 TEMP.'
	end
'I': begin
	field='nsb_scdb_eus.psdd2t'
	fdescp='PSDD2T  SHUNT DISSIPATOR PANEL D #2 TEMP.'
	end
'J': begin
	field='nsb_scdb_eus.psdet'
	fdescp='PSDET   SHUNT DISSIPATOR PANEL E TEMP.'
	end
'K': begin
	field='nsb_scdb_eus.psdft'
	fdescp='PSDFT   SHUNT DISSIPATOR PANEL F TEMP.'
	end
'L': begin
	field='nsb_scdb_eus.psdgt'
	fdescp='PSDGT   SHUNT DISSIPATOR PANEL G TEMP.'
	end
'M': begin
	field='nsb_scdb_eus.psdht'
	fdescp='PSDHT   SHUNT DISSIPATOR PANEL H TEMP.'
	end
'N': begin
	field='nsb_scdb_eus.pshi'
	fdescp='PSHI    SHUNT CURRENT'
	end
'O': goto,nd
else: begin
	print,' Not a choice, try again'
	wait,1
	goto, power
	end
endcase
end
nd:

return
end
pro thermal,fdescp,field
field=' '
fdescp=' '
vrybng:
print,' '
esc=string(27b)
print,esc,'[2J'

begin
print,'$(5x,a,t65,a)',"TARMZT    <4-7>SCU HRM # 5  ADAPTER RING (-Z) 4 MSB'S",'enter a'
print,'$(5x,a,t65,a)','TARPYT    <4-7>SCU HRM # 6  ADAPTER RING (+Y/30/-Z)4 M','enter b'
print,'$(5x,a,t65,a)','TBDB1T    BOTTOM DECK BAY 1 TEMP.','enter c'
print,'$(5x,a,t65,a)','TBDB2T    BOTTOM DECK BAY 2 TEMP.','enter d'
print,'$(5x,a,t65,a)','TBDB3T    BOTTOM DECK BAY 3 TEMP.','enter e'
print,'$(5x,a,t65,a)','TCLAMP53  SCU LRM #21  UPPER GIRTH RING 53 GHZ','enter f'
print,'$(5x,a,t65,a)','TCLAMPYG  SCU LRM #22  LOWER GIRTH RING 53 GHZ','enter g'
print,'$(5x,a,t65,a)','TDAIPDFT  DMR IPDU A FOOTPRINT TEMP. K3','enter h'
print,'$(5x,a,t65,a)','TDBIPDFT  DMR IPDU B FOOTPRINT TEMP. K3','enter i'
print,'$(5x,a,t65,a)','TDI2HTST  < 0 > DMR IPDU B HTR PWR K16','enter j'
print,'$(5x,a,t65,a)','TDM31SRT  SCU LRM #23  DMR 31 GHZ SUPPORT RING','enter k'
print,'$(5x,a,t65,a)','TDM53SRT  SCU LRM # 7  DMR 53 GHZ SUPPORT RING','enter l'
print,'$(5x,a,t65,a)','TDM90SRT  SCU LRM # 8  DMR 90 GHZ SUPPORT RING','enter m'
print,'$(5x,a,t65,a)','TFPDHTST  < 4 > FIRAS PREAMP AND DMR IPDU B HTR  K36','enter n'
print,'$(5x,a,t65,a)',"TIISHST   <4-7>SCU HRM # 7  IIS HEAT STRAP I/F4 MSB'S",'enter o'
print,'$(5x,a,t65,a)','TITSMMYT  <4-7>MODULO 64SCU HRM # 2 INNER THERMAL SHD.','enter p'
print,'$(5x,a,t65,a)','TITSMPYT  <4-7>MODULO 64SCU HRM#3 INNR THM SHD. MLI(+Y','enter q'
print,'$(5x,a,t65,a)','TSHCLT    SCU LRM # 3  THERMAL SHD. LOWER HONEYCOMB','enter r'
print,'$(5x,a,t65,a)','TSHCUT    SCU LRM # 2  THERMAL SHD. UPPER HONEYCOMB','enter s'
print,'$(5x,a,t65,a)','TSOMLIUT  SCU LRM #16  THERMAL SHIELD UPPER OUTER MLI','enter t'
print,'$(5x,a,t65,a)','TTDB1T    TOP DECK BAY 1 TEMP.','enter u'
print,'$(5x,a,t65,a)','TTDB2T    TOP DECK BAY 2 TEMP.','enter v'
print,'$(5x,a,t65,a)','TTDB3T    TOP DECK BAY 3 TEMP.','enter w'
print,'$(5x,a,t65,a)','Return to Main Menu','enter x'
ans=get_kbrd(1)
ans=strupcase(ans)
case ans of
'A': begin
	field='nsb_scdb_eus.tarmzt'
	fdescp="TARMZT    <4-7>SCU HRM # 5  ADAPTER RING (-Z) 4 MSB'S"
	end
'B': begin
	field='nsb_scdb_eus.tarpyt'
	fdescp='TARPYT    <4-7>SCU HRM # 6  ADAPTER RING (+Y/30/-Z)4 M'
	end
'C': begin
	field='nsb_scdb_eus.tbdb1t'
	fdescp='TBDB1T    BOTTOM DECK BAY 1 TEMP.'
	end
'D': begin
	field='nsb_scdb_eus.tbdb2t'
	fdescp='TBDB2T    BOTTOM DECK BAY 2 TEMP.'
	end
'E': begin
	field='nsb_scdb_eus.tbdb3t'
	fdescp='TBDB3T    BOTTOM DECK BAY 3 TEMP.'
	end
'F': begin
	field='nsb_scdb_eus.tclamp53'
	fdescp='TCLAMP53  SCU LRM #21  UPPER GIRTH RING 53 GHZ'
	end
'G': begin
	field='nsb_scdb_eus.tclampyg'
	fdescp='TCLAMPYG  SCU LRM #22  LOWER GIRTH RING 53 GHZ'
	end
'H': begin
	field='nsb_scdb_eus.tdaipdft'
	fdescp='TDAIPDFT  DMR IPDU A FOOTPRINT TEMP. K3'
	end
'I': begin
	field='nsb_scdb_eus.tdbipdft'
	fdescp='TDBIPDFT  DMR IPDU B FOOTPRINT TEMP. K3'
	end
'J': begin
	field='nsb_scdb_eus.tdi2htst'
	fdescp='TDI2HTST  < 0 > DMR IPDU B HTR PWR K16'
	end
'K': begin
	field='nsb_scdb_eus.tdm31srt'
	fdescp='TDM31SRT  SCU LRM #23  DMR 31 GHZ SUPPORT RING'
	end
'L': begin
	field='nsb_scdb_eus.tdm53srt'
	fdescp='TDM53SRT  SCU LRM # 7  DMR 53 GHZ SUPPORT RING'
	end
'M': begin
	field='nsb_scdb_eus.tdm90srt'
	fdescp='TDM90SRT  SCU LRM # 8  DMR 90 GHZ SUPPORT RING'
	end
'N': begin
	field='nsb_scdb_eus.tfpdhtst'
	fdescp='TFPDHTST  < 4 > FIRAS PREAMP AND DMR IPDU B HTR  K36'
	end
'O': begin
	field='nsb_scdb_eus.tiishst'
	fdescp="TIISHST   <4-7>SCU HRM # 7  IIS HEAT STRAP I/F4 MSB'S"
	end
'P': begin
	field='nsb_scdb_eus.titsmmyt'
	fdescp='TITSMMYT  <4-7>MODULO 64SCU HRM # 2 INNER THERMAL SHD.'
	end
'Q': begin
	field='nsb_scdb_eus.titsmpyt'
	fdescp='TITSMPYT  <4-7>MODULO 64SCU HRM#3 INNR THM SHD. MLI(+Y'
	end
'R': begin
	field='nsb_scdb_eus.tshclt'
	fdescp='TSHCLT    SCU LRM # 3  THERMAL SHD. LOWER HONEYCOMB'
	end
'S': begin
	field='nsb_scdb_eus.tshcut'
	fdescp='TSHCUT    SCU LRM # 2  THERMAL SHD. UPPER HONEYCOMB'
	end
'T': begin
	field='nsb_scdb_eus.tsomliut'
	fdescp='TSOMLIUT  SCU LRM #16  THERMAL SHIELD UPPER OUTER MLI'
	end
'U': begin
	field='nsb_scdb_eus.ttdb1t'
	fdescp='TTDB1T    TOP DECK BAY 1 TEMP.'
	end
'V': begin
	field='nsb_scdb_eus.ttdb2t'
	fdescp='TTDB2T    TOP DECK BAY 2 TEMP.'
	end
'W': begin
	field='nsb_scdb_eus.ttdb3t'
	fdescp='TTDB3T    TOP DECK BAY 3 TEMP.'
	end
'X': goto,nd
else: begin
	print,' Not a choice, try again'
	wait,1
	goto, vrybng
	end
endcase
end
nd:
return
end
pro comdhf,fdescp,field
field=' '
fdescp=' '


esc=string(27b)
print,esc,'[2J'

cdh: print,esc,'[2J'
begin
print,'$(5x,a,t65,a)','CICU1PWR  < 3 > INS CMD UNIT 1 PWR K3','enter a'
print,'$(5x,a,t65,a)','CICU2PWR  < 0 > INST CMD UNIT 2 PWR K8','enter b'
print,'$(5x,a,t65,a)','CITU1PWR  < 1 > INS TLM UNIT 1 PWR K1','enter c'
print,'$(5x,a,t65,a)','CITU2PWR  < 1 > INST TLM UNIT 2 PWR K17','enter d'
print,'$(5x,a,t65,a)','CSTU1PWR  < 0 > S/C TLM UNIT 1 PWR K0','enter e'
print,'$(5x,a,t65,a)','CSTU2PWR  < 0 > S/C  TLM UNIT 2 PWR K16','enter f'
print,'$(5x,a,t65,a)','CTU1PWR   < 6 > CENTRAL TLM UNIT 1 PWR K6','enter g'
print,'$(5x,a,t65,a)','CTU2PWR   < 5 > CENTRAL TLM UNIT 2 PWR K13','enter h'
print,'$(5x,a,t65,a)','Return to Main Menu','enter i'
ans=get_kbrd(1)
ans=strupcase(ans)
case ans of
'A': begin
	field='nsb_scdb_eus.cicu1pwr'
	fdescp='CICU1PWR  < 3 > INS CMD UNIT 1 PWR K3'
	end
'B': begin
	field='nsb_scdb_eus.cicu2pwr'
	fdescp='CICU2PWR  < 0 > INST CMD UNIT 2 PWR K8'
	end
'C': begin
	field='nsb_scdb_eus.citu1pwr'
	fdescp='CITU1PWR  < 1 > INS TLM UNIT 1 PWR K1'
	end
'D': begin
	field='nsb_scdb_eus.citu2pwr'
	fdescp='CITU2PWR  < 1 > INST TLM UNIT 2 PWR K17'
	end
'E': begin
	field='nsb_scdb_eus.cstu1pwr'
	fdescp='CSTU1PWR  < 0 > S/C TLM UNIT 1 PWR K0'
	end
'F': begin
	field='nsb_scdb_eus.cstu2pwr'
	fdescp='CSTU2PWR  < 0 > S/C  TLM UNIT 2 PWR K16'
	end
'G': begin
	field='nsb_scdb_eus.ctu1pwr'
	fdescp='CTU1PWR   < 6 > CENTRAL TLM UNIT 1 PWR K6'
	end
'H': begin
	field='nsb_scdb_eus.ctu2pwr'
	fdecsp='CTU2PWR   < 5 > CENTRAL TLM UNIT 2 PWR K13'
	end
'I': goto,nd
else: begin
	print,' Not a choice, try again'
	wait,1
	goto, cdh
	end
endcase
end

nd:

return
end
pro timeconv,asciitime,ftime,to

if N_PARAMS() eq 0 THEN BEGIN
	PRINT,'scdtt,asciitime,ftime,to'
	PRINT,"scdtt,timearray,ftimearray,to"
	RETALL
	END

to=asciitime(0)
yr=FLTARR(N_ELEMENTS(asciitime))
day=yr
hour=yr
minute=yr
msec=yr
FOR W=0,N_ELEMENTS(asciitime)-1 DO BEGIN
	yr(W)=(FLOAT((STRMID(asciitime(W),0,2))))
 	day(W)=(FLOAT((STRMID(asciitime(W),2,3))))
	hour(W)=(FLOAT((STRMID(asciitime(W),5,2))))
	minute(W)=(FLOAT((STRMID(asciitime(W),7,2))))
	msec(W)=(FLOAT((STRMID(asciitime(W),9,5))))
	ENDFOR

yr=yr-yr(0)
day=day-day(0)
ftime=FLTARR(N_ELEMENTS(asciitime))
ftime=((yr*365.25)+day)*24+hour+minute/60.+msec*1.E-3/3600
;ftime=ftime-ftime(0)

return
end

PRO  SCDTT
;+
;******************************************************************
; NAME:
;	SCDTT
;
; PURPOSE:
;	To make short term statistical analysis of Spacecraft data
;	simple and easy.
;
; DESCRIPTION:
;	SCDTT is an IDL version 2 routine.
;
;
;
; INVOCATION:
;	Once inside the IDL enviroment just type
;	SCDTT<cr>
;	The menus are easy follow and the user should be able to execute
;	the entire routine.
;
; OUTPUTS:
;	The user has several choices for output.  The user can create a 
;	hardfile to be plotted later of any plot displayed while running 
;	the routine.  The plot files have a PostScript format.
;	In addition the user is able to save most of the statistical 
;	values calculated in an IDL SAVE file. 
;
; COMMON BLOCKS:
;	None.
;
; RESTRICTIONS:
;	None.
; NOTES:
;	The routine assumes that the user has issued tha SET_PLOT 
;	command and the terminal is setup for plots.
;***********************************************************
;#
;*NAME:
;	SCDTT
;
;*CLASS:
;	Analysis and Measurement
;
;*CATEGORY:
;
;*PURPOSE:
;	To make doing statistical analysis of Spacecraft data simple 
;	and easy.
;
;*CALLING SEQUENCE:
;	SCDTT
;
;*PARAMETERS:
;
;*EXAMPLES:
;
;*SYSTEM KEYWORDS USED:
;	!D.NAME, used to reset terminal after a hardcopy is made.
;
;*SUBROUTINES CALLED:
;	BINS
;	NEWFIELD
;	ATTCF
;	POWER
;	THERMAL
;	TIMECONV
;	QMSET
;	COMDHF
;
;**********************************************************************
;-

print,' '
print,' '
start=' '
stop=' '

; PROMPT USER FOR START AND STOP TIMES TO BE ANALYZED

read,'      Enter start time  ',start
read,'      Enter stop time   ',stop


SETLOG,'CSDR$SPACE_ARCHIVE','CSDR$SPACE_EDIT'
trmtyp=!d.name

; THIS ROUTINE LETS THE USER CHOOSE ANY FIELD IN THE SPACECRAFT DATA 
; RECORD
newfield,fdescp,field

start=STRTRIM(start,2)
stop=STRTRIM(stop,2)
field=STRTRIM(field,2)

vrybng: print,' '
field=strupcase(field)
if strlen(field) lt 10 then goto,nd

print,'  Retrieving data field '+field+' from Archive.....'
vec=READ_ARCV(field,start,stop)                 
					;Retrieving data record
print,' '

timeaxis=READ_ARCV('nsb_scdb_eus.ascii_time',start,stop)
					;Retrieving time array
to=timeaxis(0)

if trmtyp ne 'X' then qmset,1
print,' '

; THIS ROUTINE CONVERTS THE ASC II TIME ARRAY INTO FLOATING POINT           
; NUMBERS FOR PLOTTING
timeconv,timeaxis,ntime,to
to=strtrim(to,2)
to=strmid(to,0,9)

print,' '
print,' '
num=N_ELEMENTS(vec)
sd=size(vec)
if sd(0) gt 1 then begin	; This section of code deals with fields
	data=vec(*,0)		; that are sampled more than once per
	p2=sd(2)		; major frame
	nntime=(findgen(num)+1.)/num
	nntime=nntime*(ntime(p2-1)-ntime(0))+ntime(0)
	ntime=nntime
	for k=1,p2-1 do begin
		data=[data,vec(*,k)]
		end
	end else data=vec

snum=strtrim(num,2)
;**************************************************************************
; This section of code prompts the user for the sample interval and sets 
; up the other arrays for calculating the statistics.
;
cint:
print,' There are '+snum+' data points in your timerange'
print,'$(3(/),a)',' '
READ,'  Enter the interval you would like to trend over ',int
int=Float(int)
sint=fix(int)
wcenter=(findgen((num/int)+0.5)+1.0)/(num/int)
wcenter=wcenter*(ntime(num-1)-ntime(0))
wcenter=wcenter+ntime(0)
width=(ntime(num-1)-ntime(0))/(num/int)
weight=(data*0.0)+1.0
;************************************************************************
sf=' '  			; Selected Field Flag
oras=0  			; Override Auto-Scaling Flag
XT='TIME (hours) starting at '+to
TT=fdescp
esc=string(27b)
;
; THIS ROUTINE CALCULATES ALL OF THE STATISTICS
bins,ntime,data,weight,wcenter,width,meandata,wsig,wgt,smplmin,smplmax,smplpoints
;
;***************************************************************************
; This section puts a menu on screen so the user may select which data 
; array he/she would like to view
begnng:
hc=0				; Hardcopy Flag
print,'$(6(/),a)',' '
print,esc,'[2J'
print,'$(5x,a,t65,a)','Plot of Field '+field,'Enter A'
print,'$(5x,a,t65,a)','Plot of Field with Mean overplotted','Enter B'
print,'$(5x,a,t65,a)','Plot of RMS Deviation from Mean','Enter C'
print,'$(5x,a,t65,a)','Plot of Min, Max per sample interval','Enter D'
print,'$(5x,a,t65,a)','Plot of Mean, with Min & Max overplotted','Enter G'
print,'$(5x,a,t65,a)','Plot of FFT of '+field,'Enter F'
print,'$(5x,a,t65,a)','Change the trend interval, current interval '+strtrim(sint,2),'Enter H'
print,'$(5x,a,t65,a)','Override auto-scale of plot','Enter S'
print,'$(5x,a,t65,a)','Go back to auto-scaling','Enter T'
print,'$(5x,a,t65,a)','Create an IDL save file','Enter I
print,'$(5x,a,t65,a)','Move to the next Menu','Enter M'
print,'$(5x,a,t65,a)','End','Enter E'
print,'$(5(/),a)',' '

;*************************************************************************
; THIS SECTION EXECUTES THE COMMAND MADE BY THE USER ABOVE
ans=get_kbrd(1)
ans=strupcase(ans)
pltlp: if hc eq 1 then ans=sf
case ans of 
'A': begin
	sf='A'
	if oras eq 0 then begin		; This section plots the raw data
		dxmn=min(ntime)		; only
		dxmx=max(ntime)
		!x.range(0)=dxmn-(dxmn*.025)
		!x.range(1)=dxmx+(dxmx*.025)
		dymn=total(smplmin)/n_elements(smplmin)
		dymx=total(smplmax)/n_elements(smplmax)
		!y.range(0)=dymn-(.025*dymn)
		!y.range(1)=dymx+(.025*dymx)
		end
	PLOT,ntime,data,xtitle=xt,title=tt,psym=3
	end
'B': begin 
	sf='B'
	if oras eq 0 then begin		; This section plots the raw 
		dxmn=min(ntime)		; data and overplots the mean
		dxmx=max(ntime)
		!x.range(0)=dxmn-(dxmn*.025)
		!x.range(1)=dxmx+(dxmx*.025)
		dymn=total(smplmin)/n_elements(smplmin)
		dymx=total(smplmax)/n_elements(smplmax)
		!y.range(0)=dymn-(.025*dymn)
		!y.range(1)=dymx+(.025*dymx)
		end
        yt=' Field with Mean overplotted, binsize '+strtrim(sint,2)
	PLOT,ntime,data,xtitle=xt,title=tt,ytitle=yt,psym=3
	OPLOT,wcenter,meandata,psym=-6
	end
'C': begin
	sf='C'
	if oras eq 0 then begin		; This section plots the RMS 
		dxmn=min(wcenter)	; deviation from the mean
		dxmx=max(wcenter)
		!x.range(0)=dxmn-(dxmn*.025)
		!x.range(1)=dxmx+(dxmx*.025)
		mwsig=total(wsig)/n_elements(wsig)
		!y.range(0)=-.25
		!y.range(1)=3*mwsig
		end
	yt=' RMS Deviation from the mean'
	plot,wcenter,wsig,ytitle=yt,xtitle=xt,title=tt,psym=-1
	end 
'D': begin
	sf='D'
	if oras eq 0 then begin		; This section plots the raw data
		dxmn=min(wcenter)	; with the minimum and maximum 
		dxmx=max(wcenter)	; of the sample interval overplotted
		!x.range(0)=dxmn-(dxmn*.025)
		!x.range(1)=dxmx+(dxmx*.025)
		dymn=total(smplmin)/n_elements(smplmin)
		dymx=total(smplmax)/n_elements(smplmax)
		!y.range(0)=dymn-(.025*dymn)
		!y.range(1)=dymx+(.025*dymx)
		end
	yt=' Min=X and Max=+ of Sample Interval 
	plot,wcenter,smplmin,psym=-7,ytitle=yt,xtitle=xt,title=tt
	oplot,wcenter,smplmax,psym=-1
	end
'G': begin
	sf='G'
	if oras eq 0 then begin		; This section plots the raw data
		dxmn=min(wcenter)	; and overplots the mean, minimum
		dxmx=max(wcenter)	; and maximum of the sample interval
		!x.range(0)=dxmn-(dxmn*.025)
		!x.range(1)=dxmx+(dxmx*.025)
		dymn=TOTAL(smplmin)/n_elements(smplmin)
		dymx=TOTAL(smplmax)/n_elements(smplmax)
		!y.range(0)=dymn-(.025*dymn)
		!y.range(1)=dymx+(.025*dymx)
		end
	yt=' Mean(squares) with Min(X) and Max(+) Overplotted'
	PLOT,ntime,data,psym=3,ytitle=yt,xtitle=xt,title=tt
	OPLOT,wcenter,meandata,psym=-6
	OPLOT,wcenter,smplmin,psym=-7
	OPLOT,wcenter,smplmax,psym=-1
	end		
'S': begin
	ans2='B'	; This section allows the user to set the scale
	oras=1		; of the plots
	read,' Which axis(ies) do you want to override (x/y/[b]) ',ans2
	ans2=strupcase(ans2)
	if ans2 eq 'B' then begin
		read,' Enter xmin,xmax,ymin,ymax ',xmn,xmx,ymn,ymx
		!x.range(0)=xmn 
		!x.range(1)=xmx 
		!y.range(0)=ymn 
		!y.range(1)=ymx
		end
	if ans2 eq 'X' then begin
		read,' Enter xmin,xmax ',xmn,xmx
		!x.range(0)=xmn 
		!x.range(1)=xmx 
		end
	if ans2 eq 'Y' then begin
		read,' Enter ymin,ymax ',ymn,ymx
		!y.range(0)=ymn 
		!y.range(1)=ymx
		end
	goto,begnng
	end
'F': begin
	sf='F'		; This section calculates the FFT of the raw data
	mnd=total(data)/n_elements(data)
	fd=data-mnd
	fftd=fft(fd,-1)
	if oras eq 0 then begin
		xray=n_elements(fftd)
		xry=indgen(xray)+1
		dxmn=min(xry)
		dxmx=max(xry)
		!x.range(0)=dxmn-(dxmn*.025)
		!x.range(1)=dxmx+(dxmx*.025)
		dymn=min(abs(fftd))
		dymx=max(abs(fftd))
		!y.range(0)=dymn-(dymx*.025)
		!y.range(1)=dymx+(dymx*.025)
		end
	TT='Field '+field
	PLOT,xry,abs(fftd),ytitle='Absolute Value of FFT',title=tt
	end
'T': begin
	oras=0		; This section returns the plots back to auto-scaling
	print,' Back to auto-scaling'
	wait,1
	goto,begnng
	end
'E': goto, nd
'H': goto,cint
'M': goto, mmenu
'I': begin
	ssf=strarr(7)		; This section allows the user to create
	rbh=0 & dbh=0 & xbh=0	; IDL save files for the data arrays
	mbh=0 & nsf=0
	chnt=0 & chmt=0
	rfile='' & dfile='' & mfile='' & xmnfile='' & xmxfile=''
	sel:print,esc,'[2J'
	print,'$(5(/),a)',' '
	print,'$(5x,a)','Select the files you want to save'
	print,' '
	print,'$(5x,a,t35,a)','Raw data','Enter R'
	print,'$(5x,a,t35,a)','Mean','Enter M'
	print,'$(5x,a,t35,a)','Min and Max','Enter X'
	print,'$(5x,a,t35,a)','RMS Deviation','Enter D'
	print,'$(5x,a,t35,a)','Quit','Enter Q'
	print,'$(6(/),a)',' '
	ans=get_kbrd(3)
	ans=strupcase(ans)
	case ans of
	'R': begin
		if rbh eq 1 then begin
			print,' Already selected'
			goto,sel
			end
		rbh=1
		print,' Please keep name 10 characters or less'
		read,' Enter name for raw data array ',rfile
		rfile=strupcase(rfile)
		ssf(nsf)=rfile
		nsf=nsf+1
		chr=execute(rfile+'=data')
		xrfile=strupcase(rfile+'_xarray')
		ssf(nsf)=xrfile
		chnt=execute(xrfile+'=ntime')
		nsf=nsf+1
		if nsf ge 7 then goto,al
		end
	'M': begin
		if mbh eq 1 then begin
			print,' Already selected'
			goto,sel
			end
		print,' Please keep name 10 characters or less'
		read,' Enter name for Mean data array ',mfile
		mfile=strupcase(mfile)
		ssf(nsf)=mfile
		nsf=nsf+1
		mbh=1
		chm=execute(mfile+'=meandata')
		if chmt eq 0 then begin
			ssf(nsf)='STAT_XARRAY'
			stat_xarray=wcenter
			chmt=1
			nsf=nsf+1
			end
		if nsf ge 7 then goto,al
		end
	'X': begin
		if xbh eq 1 then begin
			print,' Already selected'
			goto,sel
			end
		print,' Please keep name 10 characters or less'
		read,' Enter name for Min data array ',xmnfile
		xmnfile=strupcase(xmnfile)
		read,' Enter name for Max data array ',xmxfile
		xmxfile=strupcase(xmxfile)
		ssf(nsf)=xmnfile
		nsf=nsf+1
		ssf(nsf)=xmxfile
		nsf=nsf+1
		xbh=1
		chxmn=execute(xmnfile+'=smplmin')
		chxmx=execute(xmxfile+'=smplmax')
		if chmt eq 0 then begin
			ssf(nsf)='STAT_XARRAY'
			stat_xarray=wcenter
			chmt=1
			nsf=nsf+1
			end
		if nsf ge 7 then goto,al
		end
	'D': begin
		if dbh eq 1 then begin
			print,' Already selected'
			goto,sel
			end
		print,' Please keep name 10 characters or less'
		read,' Enter name for RMS Deviation array ',dfile
		dfile=strupcase(dfile)
		ssf(nsf)=dfile
		nsf=nsf+1
		dbh=1
		chd=execute(dfile+'=wsig')
		if chmt eq 0 then begin
			ssf(nsf)='STAT_XARRAY'
			stat_xarray=wcenter
			chmt=1
			nsf=nsf+1
			end
		if nsf ge 7 then goto,al
		end
	'Q': begin
		al:if nsf ge 1 then begin
			ssf=strtrim(ssf,2)
			sfile=''
			ss='SAVE'
			read,' Enter the name of the file where the arrays will be stored ',sfile
			sfile=strupcase(sfile)
			ss=ss+',filename='+'"'+sfile+'"'
			for g=0,nsf-1 do begin
				ss=ss+','+ssf(g)
				end
			print,' The following Save command is being used'
			print,ss
			wait,3
			t=execute(ss)
			end
		goto,begnng
		end
	else: begin
		print,' Not a selection, please try again'
		wait,1
		goto,sel
		end
	endcase
	goto,sel
	end
else: begin
	print,' Not a choice, try again'
	wait, 3
	goto, begnng
	end
endcase
;***************************************************************************
if hc eq 1 then begin
	device,/close
	set_plot,trmtyp
	goto,begnng
	end
			;This section gives the user the ability to make
			; a file for plotting later
print,'$(5x,a)','Would you like a hardcopy of this plot, enter y if yes'
ans3=get_kbrd(3)
ans3=strupcase(ans3)
case ans3 of 
'Y': begin
	print,esc,'2'
	goto, plt
	end
string("15b): begin
	print,esc,'2'
	goto,begnng
	end
else: begin
	print,esc,'2'
	goto,begnng
	end
endcase
;
;*******************************************************************
;
mmenu:print,esc,'[2J'
print,'$(6(/),a)',' '
print,'$(5x,a,t65,a)','Obtain another data field','enter a'
print,'$(5x,a,t65,a)','Quit the routine','enter q'
print,'$(7(/),a)',' '
ans4=get_kbrd(3)
ans4=strupcase(ans4)
case ans4 of 
'A': begin
	print,esc,'2'	; This section allows the user to choose another
	goto,newf	; data field or quit the routine
	end
'Q': goto,nd
else: begin
	print,esc,'2'
	print,' Not a choice, try again'
	wait, 2
	goto,mmenu
	end
endcase
;
;*********************************************************************
;
;
newf: begin
	newfield,fdescp,field
	goto, vrybng
	end
			; This section actually makes the file to be
			; plotted
plt: begin
	print,'$(5x,a)','Making a hardcopy'
	set_plot,'ps'
	device,filename=field,/landscape
	hc=1
	goto, pltlp
	end
nd:

return
end
;DISCLAIMER:
;
;This software was written at the Cosmology Data Analysis Center in
;support of the Cosmic Background Explorer (COBE) Project under NASA
;contract number NAS5-30750.
;
;This software may be used, copied, modified or redistributed so long
;as it is not sold and this disclaimer is distributed along with the
;software.  If you modify the software please indicate your
;modifications in a prominent place in the source code.  
;
;All routines are provided "as is" without any express or implied
;warranties whatsoever.  All routines are distributed without guarantee
;of support.  If errors are found in this code it is requested that you
;contact us by sending email to the address below to report the errors
;but we make no claims regarding timely fixes.  This software has been 
;used for analysis of COBE data but has not been validated and has not 
;been used to create validated data sets of any type.
;
;Please send bug reports to CGIS@ZWICKY.GSFC.NASA.GOV.


