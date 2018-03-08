;+NAME\ONE LINE DESCRIPTION OF ROUTINE: 
;     DIRBE_TOD reads DIRBE data and corrects for order and gains.
;
; DESCRIPTION:
;     This routine reads DIRBE data and converts it from process order
;     to science-detector order. If data is read from NBS_CTOD or 
;     BMD_MA, DIRBE_TOD will correct for auto gain and commandable gain 
;     using pipeline gain correction ratios.  If data is read from 
;     BLI_TOIRSCLD, DIRBE_TOD will decompress it.
;
; CALLING SEQUENCE:
;     DATA = DIRBE_TOD(START,STOP [,MFTIME,DET=detnum,FMT=infmt,
;                      DS=dataset,OMS=mode,AUTO=again,COMM=cgain])
;
; ARGUMENTS (I = input, O = output, [] = optional):
;     DATA     O   fltarr   DIRBE Intensities per band. By default, data 
;                           will be read from BMD_MA and returned for 
;                           all bands and all operating modes in 1x DNs.
;     START    I   str      Start time in zulu or VAX format.
;     STOP     I   str      Stop time in zulu or VAX format.
;     MFTIME  [O]  strarr   Array of MF times in zulu format.
;     DET     [I]  int      Detector number (0-15) = (1A,1B,...10).  If 
;                           specified, only one band will be returned.
;     DS      [I]  str      Dataset to read (Default='BMD_MA').  Valid 
;                           datasets include 'BMD_MA','NBS_CTOD', and 
;                           'BLI_TOIRSCLD'.
;     FMT     [I]  int      Output format of data array. 
;                           If present and set to:
;                            0: DATA = FLTARR(256*MF,16)  (Default)
;                            1: DATA = FLTARR(16,256,MF)
;                           where MF is number of Major Frames read.
;     OMS     [I]  int      Operating mode. If present and set to:
;                           -1: SDM and CAL modes returned (Default)
;                            0: SDM mode data will be returned
;                            1: CAL mode data will be returned
;
;     AUTO    [I]  int      Auto-gain correction switch. If set to:
;                            0: Data will be corrected to 1x auto-gain 
;                               (Default).
;                            1: Data will be corrected to 16x auto-gain.
;                            2: No auto-gain correction, but gain-bit 
;                               will be factored out.
;
;     COMM    [I]  int      Commanded-gain correction switch. If set to:
;                            0: Data corrected for commanded-gain 
;                               (Default).
;                            1: No commanded-gain correction.
;
; WARNINGS:
;     1.  This routine processes the MUX for the first major frame of 
;         data ONLY. It DOES NOT check after that.
;     2.  DIRBE_TOD will not work for mux sequences with redundant 
;         processes.
;     3.  BLI CAL mode data is not corrected for commanded gain.
;     4.  COMM and AUTO switches are ignored when reading from 
;         BLI_TOIRSCLD.
;     5.  The amount of data which can be read at one time is usually
;         limited to approximately three hours, depending on available
;         memory.
;
; EXAMPLE:
;     To get data from BMD_MA for all bands, both CAL and SDM modes,
;     in 1x DNs for the 1st hour of 01-January-1990:
;       data = DIRBE_TOD('900010000','900010100')
;
;     If just CAL mode is desired:
;       data = DIRBE_TOD('900010000','900010100',oms=1)
;
;     If just CAL mode is desired, you want the array divided into 
;     MFs, and you want the MF times:
;       data = DIRBE_TOD('9000100','9000101',mftime,oms=1,fmt=1)
;
;     If just CAL mode is desired, and no commanded gain correction:
;       data = DIRBE_TOD('9000100','9000101',oms=1,comm=1)
;
;     If you want to read BLI data for SDM mode only:
;       data = DIRBE_TOD('9000100','9000101',ds='bli_toirscld',oms=0)
;
;     If you want to read BLI data for band 2A only:
;       data = DIRBE_TOD('9000100','9000101',det=3,ds='bli_toirscld')
;#
; COMMON BLOCKS:
;     None.
;
; PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     Reads the science data, mux sequence, and oms from the specified
;     or default dataset found in CSDR$DIRBE_ARCHIVE.  Gets the 
;     gain-correction ratios from DIRBE_GAINS.  Note: If the system 
;     variable !QUIET is non-zero, no informational messages will be 
;     printed.
;
; PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     Calls routines DIRBE_GAINS and DECOMPRESS_MF.
;
; MODIFICATION HISTORY:
;     Written: Spiesman (GSC) 89
;     Mod: Spiesman 89 - Added !quiet support
;     Mod: Spiesman Jun 91 - Solved Temporary file creation problem.
;          SPR # 8642
;     Mod: Franz Mar 92 - Removed temporary file creation.  Modified 
;          to use pipeline auto and commanded-gain ratios, and 
;          changed default output to 1x DNs.  Added optional output 
;          format and optional dataset specification.  Changed default 
;          dataset to BMD_MA.  Enhanced to allow specification of 
;          operating mode and desired gain corrections.  Added ability 
;          to read and decompress BLI data. Modified to use READ_TOD.
;     Mod: Franz Apr 92 - Changed call to READ_TOD such that only one
;          channel of data will be read if only one detector is 
;          requested.  Changed detector selection to a keyword input.
;          Added MF times as optional output.
;.title
;Routine DIRBE_TOD
;-
function dirbe_tod,start,stop,datimas,det=band,ds=dataset,fmt=outfmt,$
                                  oms=mode,auto=again,comm=cgain

;
; Set default values for unspecified keyword inputs
;
if (n_elements(dataset) eq 0) then begin
     dataset = 'BMD_MA'
     if (!quiet eq 0) then print,'Default dataset, '+dataset+', will be used' 
endif else begin
     dataset = strupcase(strtrim(dataset,2))
     if  (dataset ne 'BMD_MA'        and $
          dataset ne 'NBS_CTOD'      and $
          dataset ne 'BLI_TOIRSCLD') then begin
          if (!quiet eq 0) then  print,'** Invalid dataset name, '+dataset
          return,-1
     endif
endelse
if (n_elements(outfmt) eq 0)  then outfmt = 0
if (n_elements(mode)   eq  0) then mode   =-1
if (n_elements(again)  eq 0)  then again  = 0
if (n_elements(cgain)  eq 0)  then cgain  = 0
if (n_elements(band)   eq 0)  then bflag  = 0  $
                              else bflag  = 1

;
; Get MF times and engineering attributes
;
if (!quiet eq 0) then print,'Reading Engineering Information'
istat = 0
istat = read_tod(dataset,'datimas,dameps,telemetry_format,daoms',$
                 start,stop,datimas,meps,tlm,daoms)
if (not istat) then begin
    print,'READ_TOD returned error status'
    return,-1
endif
nrec = n_elements(daoms)
if (!quiet eq 0) then print,strtrim(string(nrec),2),' MFs in timerange'

;
; Select OMS=mode, if desired, and eliminate MFs with non-science
; telemetry format.  Note: BMD_MA and BLI_TOIRSCLD should already have
; non-science T/M removed.
;
if (!quiet eq 0) then print,'Checking OMS and T/M format'
tlm = shift(tlm,-1) or tlm
if ((mode ge 0) and (mode le 1)) then begin
    select = where((daoms eq mode) and (tlm eq 1))
    if (select(0) eq -1) then begin
        if (!quiet eq 0) then $
          print,'** No OMS=',strtrim(string(mode),2),' data in selected range.'
        return,-1
    endif
endif else begin
    select = where(((daoms eq 0) or (daoms eq 1)) and (tlm eq 1))
    if (select(0) eq -1) then begin
        if (!quiet eq 0) then $
            print,'** No science telemetry in selected range.'
        return,-1
    endif
endelse
nrec = n_elements(select)
if (!quiet eq 0) then $
    print,strtrim(string(nrec),2),' MFs with selected attributes'


;
; Correct MEPS for late field.  Note: If only 1-MF exists after selecting, 
; then there is not enough info in timerange to correct for late field.
;
if (nrec eq 1) then begin
   if (!quiet eq 0) then begin
       print,'** No trailing MF in selected timerange'
       print,'** Unable to get MEPS information.'
   endif
   return,-1
endif
meps = meps(*,select(0)+1)

;
; Calculate detector sequence
;
idet  = bindgen(16)
mm    = meps(2*idet+1)

mux  = bytarr(16)
for i=0,15 do begin
    it   = where(mm eq i)          ; Note:  This will not work for 
    mux(i) = it(0)                   ; mux > mux2 (redundant processes).
endfor

;    
; Reorder from 1a 2a 3a... to 1a 1b 1c...
;
s1=mux(1)
s2=mux(2)
s3=mux(3)
s5=mux(5)
s6=mux(6)
s7=mux(7)

mux(1)=s3
mux(2)=s6
mux(3)=s1
mux(5)=s7
mux(6)=s2
mux(7)=s5

;
; Read data from archive
;
if (!quiet eq 0) then print,'Extracting data from archive'
if (bflag eq 0) then                                    $
    drange = 'dadrbsci2(1:16,1:256)'                    $
else                                                    $
    drange = 'dadrbsci2('+strtrim(string(mux(band)+1),2)+',1:256)'
istat = 0
istat = read_tod(dataset,drange,start,stop,data)
if (not istat) then begin
    print,'READ_TOD returned error status'
    return,-1
endif

;
; Put data in science detector order
;
if (bflag eq 0) then data = data(mux(idet),*,*)

;
; Select MFs with desired attributes
;
if (bflag eq 1) then begin
    temp = data
    data = intarr(256,nrec)
    data(0:255,0:nrec-1) = temp(0,0:255,select)
endif else data = data(*,*,select)
datimas = datimas(select)

if (dataset ne 'BLI_TOIRSCLD') then begin
    ;
    ; Correct to 1x DNs
    ;
    if (!quiet eq 0) then print,'Gain correcting'

    ;
    ; Get gain correction factors and separate-out saturated data
    ;
    x16_x1   = dirbe_gains(0,datimas(0))
    ape      = dirbe_gains(1,datimas(0))
    isat     = where(data eq 32767)

    ;
    ; Factor-out gain bit and remove effects of onboard truncation
    ;
    temp = float(data/2)

    if (again le 1) then begin
        itrunc_h = where(data/2 gt  255) 
        itrunc_l = where(data/2 lt -255) 
        if (itrunc_h(0) ne -1) then temp(itrunc_h) = temp(itrunc_h) - 0.5
        if (itrunc_l(0) ne -1) then temp(itrunc_l) = temp(itrunc_l) + 0.5
    endif

    if (bflag eq 0) then begin
        ;
        ; Loop through detectors and factor-out gains as requested
        ;
        for i=0,15 do begin
            ;
            ; Reduce to 1x auto-gain
            ;
            if (again eq 0) then          $
                temp(i,*,*) = temp(i,*,*) $
                 / ( 1.+ ((data(i,*,*) and 1) ne 1) * (x16_x1(i)-1.) )

            ;
            ; Or increase to 16x auto-gain
            ;
            if (again eq 1) then          $
                temp(i,*,*) = temp(i,*,*) $
                 * ( 1.+ (data(i,*,*) and 1) * (x16_x1(i)-1.) )

            ;
            ; Divide-out commanded gain
            ;
            if (cgain le 0) then temp(i,*,*) = temp(i,*,*)/ape(i)

        endfor

    endif else begin
        ;
        ; For desired detector, factor-out gains as requested
        ;
        if (again eq 0) then          $
                temp(*,*) = temp(*,*) $
                 / ( 1.+ ((data(*,*) and 1) ne 1) * (x16_x1(band)-1.) )

        ;
        ; Or increase to 16x auto-gain
        ;
        if (again eq 1) then          $
            temp(*,*) = temp(*,*) $
            * ( 1.+ (data(*,*) and 1) * (x16_x1(band)-1.) )

        ;
        ; Divide-out commanded gain
        ;
        if (cgain le 0) then temp(*,*) = temp(*,*)/ape(band)

    endelse

    data = temp
    if (isat(0) ne -1) then data(isat) = -16382.      ; Saturation sentinel

endif else begin
    ;
    ; Decompress BLI data
    ;
    if (!quiet eq 0) then print,'Decompressing BLI data'
    data = decompress_mf(data)

endelse

;
; Convert to standard "DIRBE_DATA" format if desired.
;
if (outfmt eq 0) then begin
    if (bflag eq 1) then begin
        temp = fltarr(256*nrec)
        for j = 0L,nrec-1 do temp(j*256:(j+1)*256-1) = data(0:255,j)
        data = temp
    endif else begin
        temp = fltarr(256*nrec,16)
        for i=0,15 do begin
              for j = 0L,nrec-1 do temp(j*256:(j+1)*256-1,i) = data(i,0:255,j)
        endfor
        data = temp
    endelse
endif 

if (bflag eq 1) then begin
    if (!quiet eq 0) then print,'Returning detector ',strtrim(string(band),2)
endif else begin
    if (!quiet eq 0) then print,'Returning all bands'
endelse

return,data
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


