pro read_covar,dsname,covar=covar,real_real_covar=real_real_covar, $
  imag_imag_covar=imag_imag_covar,imag_real_covar=imag_real_covar,$
  real_imag_covar=real_imag_covar,temps_covar=temps_covar,$
  real_temps_covar=real_temps_covar,imag_temps_covar=imag_temps_covar,$
  gigahz_scale=gigahz_scale,model_label=model_label,cov_label=cov_label,$
  chanscan=chanscan,nu_zero=nu_zero,delta_nu=delta_nu,num_freq=num_freq,$
  num_ifgs=num_ifgs,num_coadds=num_coadds,deg_freedom=deg_freedom,$
  bol_avg=bol_avg,cbias_avg=cbias_avg,volt_avg=volt_avg,s0=s0,tau=tau,$
  tbol=tbol,avg_rcalspec=avg_rcalspec,avg_icalspec=avg_icalspec,$
  bin_total=bin_total,wtd_bin_total=wtd_bin_total,eff_wt=eff_wt,$
  real_dispersion=real_dispersion,imag_dispersion=imag_dispersion,temp=temp

;+NAME/ONE-LINE DESCRIPTION:
;      READ_COVAR: used to read in the FIRAS Covariance matrix and
;                  associated fields
;
;DESCRIPTION: Obsolete UIDL routine which reads in fields of the FIRAS IP/PDS
;             files. Not valid for Pass 4 (=final release) Data Products.
;
;CALLING SEQUENCE: 
;        read_covar - no arguements lists calling sequence
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    param		I/O   	type		description
;    dsname              I      string		full name of file for input
;    covar              [O]     dblarr          output full covariance matrix
;    real_real_covar    [O]     dblarr          self explanatory NON-ZERO
;    imag_imag_covar    [O]     dblarr               parts of the
;    imag_real_covar    [O]     dblarr             covariance matrix
;    real_imag_covar    [O]     dblarr                    "
;    temps_covar        [O]     dblarr                    " (temperatures)
;    real_temps_covar   [O]     dblarr                    "      "
;    imag_temps_covar   [O]     dblarr                    "      "
;    model_label        [O]     string          calib. model soln label
;    cov_label          [O]     string          covariance matrix name
;    chanscan           [O]     string          channel/scan mode
;    nu_zero            [O]     float           freq. of inital data pt (GHz)
;    delta_nu           [O]     float           freq. interval (GHz)
;    num_freq           [O]     long            # freq. points
;    num_ifgs           [O]     long            # contributing ifgs
;    num_coadds         [O]     long            # contributing coadds
;    deg_freedom        [O]     long            degrees of freedom
;    bol_avg            [O]     float           average bolometer temp. (K)
;    cbias_avg          [O]     float           average commmand bias (volts)
;    volt_avg           [O]     float           average readout voltage
;    s0                 [O]     double          DC responsivity (volts/watt)
;    tau                [O]     double          detector time constant (secs)
;    tbol               [O]     double          derived bolometer temp. (K)
;    avg_rcalspec       [O]     dblarr          ave. real parts calibrated 
;                                               spectrum (MJy/sr)
;    avg_icalspec       [O]     dblarr          ave. imag. parts calibrated 
;                                               spectrum (MJy/sr)
;    bin_total          [O]     lonarr          # calibrated spectra in bin
;    wtd_bin_total      [O]     fltarr          weighted total of calibrated
;                                               spectra in bin (MJy/sr)
;    eff_wt             [O]     dblarr          effective weight per IFG
;                                               contributing to bin
;    real_dispersion    [O]     dblarr          total real part, calibrated
;                                               spec. dispersions/bin (MJy/sr)
;    imag_dispersion    [O]     dblarr          total imag. part, calibrated
;                                               spec. dispersions/bin (MJy/sr)
;    temp               [O]     dblarr          temperature and glitch rate
;                                               dispersions totals in bin
;
;#
; SUBROUTINES CALLED: most of the FXB routines for handling the FITS
;                     I/O.
;
; COMMON BLOCKS:  None
;
; LIBRARY CALLS:  None
;
; WARNINGS:
;
; PROGRAMMING NOTES:
;
; MODIFICATION HISTORY: A. Trenholme/J. Newmark May 1994 Original SPR 11759
;
;.TITLE
;Routine READ_COVAR
;-
if n_params() lt 1 then begin
    print, 'Syntax: read_covar, data_set_name, covar=covar,real_real_covar='+$ 
  'real_real_covar,' + $
  'imag_imag_covar=imag_imag_covar,imag_real_covar=imag_real_covar,'+$
  'real_imag_covar=real_imag_covar,temps_covar=temps_covar,'+$
  'real_temps_covar=real_temps_covar,imag_temps_covar=imag_temps_covar,'+$
  'gigahz_scale=gigahz_scale,model_label=model_label,cov_label=cov_label,'+$
  'chanscan=chanscan,nu_zero-nu_zero,delta_nu=delta_nu,num_freq=num_freq,'+$
  'num_ifgs=num_ifgs,num_coadds=num_coadds,deg_freedom=deg_freedom,'+$
  'bol_avg=bol_avg,cbias_avg=cbias_avg,volt_avg=volt_avg,s0=s0,tau=tau,'+$
  'tbol=tbol,avg_rcalspec=avg_rcalspec,avg_icalspec=avg_icalspec,'+$
  'bin_total=bin_total,wtd_bin_total=wtd_bin_total,eff_wt=eff_wt,'+$
  'real_dispersion=real_dispersion,imag_dispersion=imag_dispersion,temp=temp'
    goto, abort
endif
if n_elements(dsname) lt 1 then begin
    print, 'READ_COVAR: You must supply the name of the data file'
    goto, abort
endif
result = is_fits(dsname,fits_exten)
if result(0) ne 1 then begin
    print, 'READ_COVAR: this routine accepts only FITS format but this is ',$
            'not a fits file'
    goto, abort
endif

; get the primary FITS header
header = headfits(dsname)

instrume  = sxpar(header, 'INSTRUME') 
if strtrim(instrume) ne 'FIRAS' then begin
    print, 'READ_COVAR: This is not a FIRAS IP/PDS File'
    goto, abort
endif

get_lun, data_lun
fxbopen,data_lun,DSNAME,1,bin_header
fields = fxpar(bin_header,'TTYPE*',count=fieldcnt)
for i=1,fieldcnt do begin
  fxbread,data_lun,data,i
  if (i le 15) then begin
    z=WHERE(data NE data(0))
    IF z(0) EQ -1 THEN data=data(0) ELSE PRINT,'Problems with ' + $
       fields(i)
  endif
  if i eq 1 then model_label=data
  if i eq 2 then cov_label=data
  if i eq 3 then chanscan=data
  if i eq 4 then nu_zero=data
  if i eq 5 then delta_nu=data
  if i eq 6 then num_freq=data
  if i eq 7 then num_ifgs=data
  if i eq 8 then num_coadds=data
  if i eq 9 then deg_freedom=data
  if i eq 10 then bol_avg=data
  if i eq 11 then cbias_avg=data
  if i eq 12 then volt_avg=data
  if i eq 13 then s0=data
  if i eq 14 then tau=data
  if i eq 15 then tbol=data
  if i eq 16 then avg_rcal=data
  if i eq 17 then avg_ical=data
  if i eq 18 then bin_total=data
  if i eq 19 then wtd_bin_total=data
  if i eq 20 then eff_wt=data
  if i eq 21 then rdisp=data
  if i eq 22 then idisp=data
  if i eq 23 then temp=data
  if i eq 24 then fip_covar=data
endfor        
fxbclose,data_lun
free_lun, data_lun

GIGAHZ_SCALE=NU_ZERO+DELTA_NU*DINDGEN(NUM_FREQ)
NF_1=NUM_FREQ-1
AVG_RCALSPEC=AVG_RCAL(0:NF_1)
AVG_ICALSPEC=AVG_ICAL(0:NF_1)
IMAG_DISPERSION=IDISP(0:NF_1,*)
REAL_DISPERSION=RDISP(0:NF_1,*)

COVAR = DBLARR(368,368)

;  GENERATE AN INDEX TO ACCESS THE REAL-REAL CORRELATION PART OF THE MATRIX,
;  AND ASSIGN THE VALUES FROM THE "FIP_COVAR" FIELD IN THE STRUCTURE

;    REAL-REAL (17766 VALUES) STORED IN RECORDS 1-66 (0-65 IN IDLSPEAK);
;         UPPER HALF OF SYMMETRIC SUBMATRIX STORED AS A SINGLE VECTOR

INDEX1 = LINDGEN(188)
FOR I=1,187 DO INDEX1 = [ INDEX1, (368L*I + LINDGEN(188-I) + I) ]

COVAR( INDEX1 ) = FIP_COVAR( 0:17765 )

;  GENERATE AN INDEX TO ACCESS THE REAL-IMAG CORRELATION PART OF THE MATRIX,
;  AND ASSIGN THE VALUES FROM THE "FIP_COVAR" FIELD IN THE STRUCTURE

;    REAL-IMAG (33840 VALUES) STORED IN RECORDS 67-192 (66-191 IN IDLSPEAK);
;         ENTIRETY OF RECTANGULAR SUBMATRIX STORED AS A SINGLE VECTOR

INDEX2 = LINDGEN(180) + 188
FOR I=1,187 DO INDEX2 = [ INDEX2, (368L*I + LINDGEN(180) + 188) ]

COVAR( INDEX2 ) = FIP_COVAR(17820:51659) 

;  GENERATE AN INDEX TO ACCESS THE IMAG-IMAG CORRELATION PART OF THE MATRIX,
;  AND ASSIGN THE VALUES FROM THE "FIP_COVAR" FIELD IN THE STRUCTURE

;    IMAG-IMAG (16290 VALUES) STORED IN RECORDS 193-253 (192-252 IN IDLSPEAK)
;         UPPER HALF OF SYMMETRIC SUBMATRIX STORED AS A SINGLE VECTOR

INDEX3 = 368L*188 + LINDGEN(180) + 188
FOR I=1,179 DO INDEX3 = [ INDEX3, ( 368L*(188+I) + LINDGEN(180-I) + 188 + I ) ]

COVAR( INDEX3 ) = FIP_COVAR(51840:68129)

;  SINCE THE UPPER-HALF OF THE (SYMMETRIC) COVARIANCE MATRIX IS NOW FILLED, WE
;  GENERATE THE LOWER HALF BY TRANSPOSING THE ARRAY.  ZERO THE ELEMENTS ALONG 
;  THE DIAGONAL OF THE TRANSPOSE MATRIX.

TCOVMAT = TRANSPOSE( COVAR )
TCOVMAT( LINDGEN(368)*(368+1) ) = 0.

;  GENERATE THE MATRIX BY THE ADDITION OF THE UPPER AND LOWER HALF.  
;  COVAR is the 368 by 368 delivered matrix.  It will usually have LOTS of zeros 
;  in it.
COVAR  = (COVAR + TCOVMAT)

TOP_IX= 188+NF_1


; The following guys (whose names are self-evident) are the designated NON-ZERO 
; parts of the covariance matrices.

REAL_REAL_COVAR=COVAR(0:NF_1,0:NF_1)
IMAG_IMAG_COVAR=COVAR(188:TOP_IX,188:TOP_IX)
IMAG_REAL_COVAR=COVAR(0:NF_1,188:TOP_IX)
REAL_IMAG_COVAR=COVAR(188:TOP_IX,0:NF_1)
TEMPS_COVAR=COVAR(180:187,180:187)
REAL_TEMPS_COVAR=COVAR(180:187,0:NF_1)
IMAG_TEMPS_COVAR=COVAR(180:187,188:TOP_IX)

abort:
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


