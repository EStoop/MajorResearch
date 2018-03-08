function dirbe_meps,start,stop,ds=dataset
;+NAME:
;     DIRBE_MEPS returns the MUX sequence and prints it out.
;
; DESCRIPTION:
;     IDL function to extract the mux sequence from the DIRBE archive
;     and generate a table of band versus MUX process number.
;
; CALLING SEQUENCE:
;     MUX = MEPS(START,STOP,[DS=dataset])
;
; ARGUMENTS (I = input, O = output, [] = optional):
;     START     I    str      Start time in zulu or VAX format.
;     STOP      I    str      Stop time in zulu or VAX format.
;     DS       [I]   str      Dataset to read (Default='BMD_MA').  Valid 
;                             datasets include 'BMD_MA','NBS_CTOD', and 
;                             'BLI_TOIRSCLD'.
;     MUX       O    bytarr   Array of engineering detector numbers in 
;                             MUX process order.  Mux (0,1,2,3...) = 
;                             Bands (1a,2a,3a,1b...).   
;
; WARNINGS:
;     1.     DIRBE_MEPS will only return data for the first major frame
;            in the requested time range.
;     2.     If an invalid dataset name is specified, -1 is returned.
;
; COMMON BLOCKS:
;     None.
;
; EXAMPLE:
;     To get the mux sequence at the start of 1 Jan 1990.
;         mux = dirbe_meps('900010000','900010001')
;
;     Note: If the BMD_MA file for 1 Jan 1990 was not online, but the 
;           BLI_TOIRSCLD file was online, then you would have to type:
;         mux = dirbe_meps('900010000','900010001',ds='bli_toirscld')
;
;#
; PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     Reads the the field DAMEPS from the specified or default dataset 
;     found in CSDR$DIRBE_ARCHIVE.  Converts DAMEPS to an array of 
;     engineering detector numbers in the mux process order.  If the 
;     system variable !QUIET is non-zero, no informational messages will be 
;     printed.
;
; PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     None.
;
; MODIFICATION HISTORY:
;     Written: Spiesman (GSC) 89
;     Mod: Spiesman 89 - Added !quiet support
;     Revised by BAF 16 Mar 1992 to allow specification of the dataset name.
;
; SPR 9616
;.TITLE
;Routine DIRBE_MEPS
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

;read meps
if (!quiet eq 0) then print,'[Reading meps data]'
istat = 0
istat = read_tod(dataset,'dameps',start,stop,m,maxrec=2)
if (not istat) then begin
    print,'READ_RAWTOD returned error status'
    return,-1
endif

meps=bytarr(16)
;remember that meps applies to FOLLOWING Mf
for i=0,15 do meps(i)=m(2*i+1,1)

print,'Band 1a = process ',where(meps eq 0)
print,'Band 1b = process ',where(meps eq 3)
print,'Band 1c = process ',where(meps eq 6)
print,'Band 2a = process ',where(meps eq 1)
print,'Band 2b = process ',where(meps eq 4)
print,'Band 2c = process ',where(meps eq 7)
print,'Band 3a = process ',where(meps eq 2)
print,'Band 3b = process ',where(meps eq 5)
print,'Band 3c = process ',where(meps eq 8)
print,'Band  4 = process ',where(meps eq 9)
print,'Band  5 = process ',where(meps eq 10)
print,'Band  6 = process ',where(meps eq 11)
print,'Band  7 = process ',where(meps eq 12)
print,'Band  8 = process ',where(meps eq 13)
print,'Band  9 = process ',where(meps eq 14)
print,'Band 10 = process ',where(meps eq 15)
print,'Chop  A = process ',where(meps eq 32)
print,'Chop  B = process ',where(meps eq 33)

return,meps

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


