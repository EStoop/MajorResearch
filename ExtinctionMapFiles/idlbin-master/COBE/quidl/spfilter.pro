function spfilter,inmap,npix,type,badval=bad_data
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    SPFILTER is a general spatial filtering tool.
;
;DESCRIPTION:                                   
;    IDL function which takes as input a 2-D array and replaces each 
;    element with the result of a user-specified function applied to 
;    all data within an NPIX x NPIX window centered on that element.
; 
;CALLING SEQUENCE:
;    outmap = SPFILTER(map,npix,filter_func,badval=badval)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    map           I     arr       2-D array containing the input map
;    npix          I     scl       Window size in pixels (full width 
;                                  and height) 
;    filter_func   I     str       String containing name of filter 
;                                  function.  The function can be any  
;                                  IDL function which takes as an array  
;                                  as input and returns a scalar.  
;                                  Standard examples include 'MIN',
;                                  'MEDIAN','AVG'.  
;    badval       [I]    scl       All data at this value will be ignored 
;                                  (Default=0).
;    outmap        O     arr       2-D array containing the filtered map
;
;WARNINGS:
;    None.
;
;EXAMPLES:
;    To apply a 5x5 minimizing filter to the 2-D array called MAP:
;      filtered_map = SPFILTER(map,5,'min')
;
;    To apply a 10x10 averaging filter to the 2-D array called MAP:
;
;      Create within your !path a function of the form:
;
;        function myavg,data
;         n = n_elements(data)
;         return,total(data)/n
;        end
;
;      Then call SPFILTER as:
;        filtered_map = SPFILTER(map,10,'myavg')
;
;COMMON BLOCKS:
;    None.
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    The window associated with point i,j is defined by the indices 
;    (i-fix(fsize/2):i+fix(fsize/2),j-fix(fsize/2):j+fix(fsize/2)),
;    where any indices less than zero or greater than the max elements
;    per dimension are ignored.  Thus, setting Fsize=4 is equivalent
;    to setting Fsize=5, since both will produce a 5x5 window of the 
;    form (i-2:i+2,j-2:j+2).
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    None.
;
;MODIFICATION HISTORY
;    Written by B.A. Franz, Applied Research Corp.,   Feb 1993.
;    Replaced 'Execute' statement with 'Call_Function', BAF, May 93
;
;------------------------------------------------------------------------------------
;
on_error,2
;
; Check Inputs
;
if (n_params(0) lt 3) then begin
    message,'Not enough arguments passed, returning input map.',/continue
    return,inmap
endif
;
s  = size(inmap)
if (s(0) ne 2) then begin
    message,'Input map must be 2-D array, returning input map.',/continue
    return,inmap
endif
nx = s(1)
ny = s(2)
;
if (npix le 1) then begin
    message,'Filter size is <= 1, returning input map.',/continue
    return,inmap
endif
;
s  = size(type)
ns = n_elements(s)
if (s(ns-2) ne 7) then begin
    message,'Filter type must be a string, returning input map.',/continue
    return,inmap
endif
;
if (n_elements(bad_data) eq 0) then bad_data=0.

;
; Initialize Ouput Array
;
outmap = inmap*0.0 + bad_data

;
; Form map of filtered data values within each subregion
;
w = fix(npix/2)
xrange = intarr(2)
yrange = intarr(2)

for i=0,nx-1 do for j=0,ny-1 do begin
    ;
    ; Select range in x for subregion
    ; 
    xrange(0) = max([   0,i-w])
    xrange(1) = min([nx-1,i+w])
    ;
    ; Select range in y for subregion
    ; 
    yrange(0) = max([   0,j-w])
    yrange(1) = min([ny-1,j+w])
    ;
    ; Form subregion and apply filter to good data
    ; 
    subreg = inmap(xrange(0):xrange(1),yrange(0):yrange(1))
    select = where(subreg ne bad_data)
    if (select(0) ne -1) then                               $
        outmap(i,j) = call_function(type,subreg(select))    $
    else                                                    $
        outmap(i,j) = bad_data
endfor


return,outmap
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


