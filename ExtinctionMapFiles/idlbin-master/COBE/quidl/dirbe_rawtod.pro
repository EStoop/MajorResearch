function dirbe_rawtod,start,stop,ds=dataset
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     DIRBE_RAWTOD returns auto-gain corrected data in process order.
;
; DESCRIPTION:
;     IDL function to read the DIRBE time-ordered data and do a 
;     first-order correction for the auto-gain.  The data is returned 
;     in process order at 16x auto-gain and full commanded gain.
;
; CALLING SEQUENCE:
;     DATA = DIRBE_RAWTOD(START,STOP [,DS=dataset])
;
; ARGUMENTS (I = input, O = output, [] = optional):
;     DATA   O    intarr    Array 16 processes x 256 hmnfs x n MFs.
;     START  I    str       Start time in zulu or VAX format.
;     STOP   I    str       Stop time in zulu or VAX format.
;     DS    [I]   str       Dataset to read (Default='BMD_MA'). Valid 
;                           datasets include 'BMD_MA' and 'NBS_CTOD'.
;
; WARNINGS:
;     1. This routine is not accurate for data taken at 1X auto-gain.  
;        It uses a simple factor of 16 to correct the gain bit.  This
;        routine is meant for first-order examination of warm-era data.  
;        For cold-era data, the preferred routine is DIRBE_DATA.
;     2. Data will be returned for all operating modes.
;     3. Data will be returned in the order that it was processed, NOT 
;        in detector order.
;
; COMMON BLOCKS:
;     None.
;
; EXAMPLE:
;     To get data from BMD_MA for the 1st hour of 01-January-1990:
;       data = DIRBE_RAWTOD('900010000','900010100')
;
;     If you want to read NBS_CTOD instead:
;       data = DIRBE_RAWTOD('900010000','900010100',ds='nbs_ctod')
;#
; PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     Reads the science data word from the specified or default dataset 
;     found in CSDR$DIRBE_ARCHIVE.  Checks the gain-bit for each data 
;     sample and multiplies by 16 if the bit is set.
;
; PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     None.
;
; MODIFICATION HISTORY:
;     Written: Spiesman (GSC) 89 (formerly DIRBE_SDM_GAIN).
;     Mod: Franz Mar 92 - Added optional optional dataset specification.  
;          Changed default dataset to BMD_MA.
;.title
;Routine DIRBE_RAWTOD
;-
;
if (n_elements(dataset) eq 0) then begin
     dataset = 'BMD_MA'
     if (!quiet eq 0) then print,'Default dataset, '+dataset+', will be used' 
endif else begin
     dataset = strupcase(strtrim(dataset,2))
     if  (dataset ne 'BMD_MA'        and $
          dataset ne 'NBS_CTOD'      and $
          dataset ne 'BLI_TOIRSCLD') then begin
          if (!quiet eq 0) then  print,'Invalid dataset name, '+dataset
          return,-1
     endif
endelse

;
; Read data from archive
;
print,'Reading data from archive'
istat = 0
istat = read_tod(dataset,'dadrbsci2',start,stop,sdm)
if (not istat) then begin
    print,'READ_TOD returned error status'
    return,-1
endif

;
; Find the number of records read
;
s = size(sdm)
n_rec = s(3)
print,strtrim(string(n_rec),2),' MFs read'

;
; Gain-bit correct
;
print,'Gain-bit correcting'
sdm  = long(sdm)
mask = 1+((sdm and 1)*15)      ;x1 or x16
sdm  = mask*(sdm and (-2))/2

return,sdm
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


