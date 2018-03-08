function nyquist, fakeit, ngroup, mtm_speed, iunits
;+ NAME/ONE LINE DESCRIPTION OF ROUTINE:
;      NYQUIST returns the Nyquist frequency of a FIRAS spectrum
;
;  DESCRIPTION:
;      IDL function to return the Nyquist frequency for a specified
;      FIRAS instrument state in either wavenumbers or hertz.  This
;      routine incorporates the frequency shift due to the on-orbit
;      sampling rate of the instrument.
;
;  CALLING SEQUENCE:
;      FNYQ = NYQUIST (fakeit, ngroup, mtm_speed, iunits)
;
;  ARGUMENTS (I = input, O = output, [] = optional):
;      FNYQ           O      flt      Nyquist frequency for the
;                                     specified instrument state
;                                     Units:  Icm if IUNITS = 0
;                                              Hz if IUNITS = 1
;      FAKEIT         I      int      MTM scanning/fakeit mode
;                                     FAKEIT = 0 indicates scanning mode
;                                     FAKEIT = 1 indicates fakeit mode
;      NGROUP         I      int      Adds per group for the spectrum
;                                     NGROUP has a range of 1-12
;      MTM_SPEED      I      int      MTM scan speed
;                                     MTM_SPEED = 0 indicates slow speed
;                                     MTM_SPEED = 1 indicates fast speed
;      IUNITS         I      int      Frequency units flag
;                                     IUNITS = 0 indicates wavenumbers
;                                     IUNITS = 1 indicates hertz
;
;  WARNINGS:
;      1.    If given an "impossible" instrument state (e.g., MTM speed
;            neither fast nor slow), the routine returns a value of -1.
;      2.    Fakeit mode data do not have a meaningful conversion to
;            wavenumbers, so hz must be specified (IUNITS = 1) when
;            FAKEIT=1.
;
;  EXAMPLE:
;      To determine the Nyquist frequency for the LLSS scan mode in icm:
;
;            fakeit    = 0
;            ngroup    = 3
;            mtm_speed = 0
;            iunits    = 0
;            fnyq = nyquist (fakeit, ngroup, mtm_speed, iunits)
;
;      Alternatively, one can invoke the function as follows:
;            fnyq = (0, 3, 0, 0)
;
;      In either case, NYQUIST returns a value of 144.981 icm to the
;      variable fnyq.
;#
;  COMMON BLOCKS:
;      None
;
;  PROCEDURE (AND OTHER PROGRAMMING NOTES)
;      Adapted from the Fortran subroutine FUT_NYQUIST.
;
;  PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;      None
;
;  MODIFICATION HISTORY:
;      Written by Gene Eplee, General Sciences Corp., 23 April 1991.
;
; SPR 9616
;
; SPR 10118  --  The FIRAS frequency scale has been recalibrated.  The
;                frequency correction factor, SHIFT, is changed from 1.0053537
;                to 1.0044579.  Gene Eplee, GSC, 14 October 1992.
;
;.TITLE
;Routine NYQUIST
;-

FRINGES  = 20.00E-04            ;grating spacing in cm
OPT_PATH = 3.464102             ;optical path = 4*cos(pi/6)
SUBDIVID = [6.0, 4.0]           ;fringe multiplier for slow/fast MTM speeds
VELOCITY = [0.799, 1.18]        ;MTM speeds in cm/sec
SHIFT    = 1.0044579            ;on-orbit frequency shift

;
;  First check that the minimum necessary scan mode and units 
;  information makes sense.
;
        IF ((NGROUP LT 1) OR (NGROUP GT 16)) THEN BEGIN
           NYQUIST = -1.0
           RETURN,NYQUIST
        ENDIF

        IF ((FAKEIT NE 0) AND (FAKEIT NE 1)) THEN BEGIN
           NYQUIST = -1.0
           RETURN,NYQUIST
        ENDIF

        IF ((IUNITS NE 0) AND (IUNITS NE 1)) THEN BEGIN
           NYQUIST = -1.0
           RETURN,NYQUIST
        ENDIF

;
;  State is OK. Now do the real stuff.
;
        IF (FAKEIT EQ 0) THEN BEGIN
;        ... we're in MTM mode. First make sure we've got a legal speed.
;
           IF ((MTM_SPEED EQ 0) OR (MTM_SPEED EQ 1)) THEN BEGIN
;          ...calculate Nyquist freq in wavenumbers.  Convert to hz if desired.
;
              NYQUIST = SHIFT*SUBDIVID(MTM_SPEED)/FRINGES/OPT_PATH/NGROUP/2.0
              IF (IUNITS EQ 1) THEN NYQUIST = NYQUIST * VELOCITY(MTM_SPEED) 

;          ...bail out due to undefined speed.
           ENDIF  ELSE BEGIN 
              NYQUIST = -1.0
              RETURN,NYQUIST
           ENDELSE
        ENDIF ELSE BEGIN

;
;  ... we're in fakeit mode, so can only supply frequency in Hz.
;
           IF (IUNITS EQ 1) THEN NYQUIST = 256.0/NGROUP ELSE NYQUIST = -1.0
        ENDELSE

        RETURN,NYQUIST
        END
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


