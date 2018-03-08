      FUNCTION SortIt, XIn, SortVec
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     Sorts a vector in ascending order of the contents of a second vector.
;
;DESCRIPTION:
;     Sorts a vector in ascending order of the contents of a second vector.
;     This function is needed because the IDL SORT function does not do a 
;     double sort correctly. In the example given below, this routine works
;     as expected, whereas the IDL SORT function returns
;             XYX = [ 1, 1, 1, 1, 2, 2, 3, 3, 4, 4 ] 
;             YYX = [ 5, 2, 6, 3, 2, 3, 4, 4, 4, 2 ]
;    
;CALLING SEQUENCE:  
;     SortedVector = SortIt ( XIn, SortVec )
;
;ARGUMENTS (I = input, O = output, [] = optional):
;     XIn           I   Numeric vector   Vector to be sorted.
;     SortVec       I   Numeric vector  Vector specifying the sort order.
;     XOut          O   Numeric vector Sorted vector.
;
;WARNINGS:
;     1. Will only sort in ascending order.
;     2. For 2-d sort, the order of the commands is important -- don't 
;        rearrange the vector you are sorting on until you've finished
;        using it!
;     3. If XIn and SortVec are different lengths, SortIt will sort the
;        first N elements of XIn based on the first N elements of SortVec,
;        where N = number of elements in the shorter of XIn and SortVec.
;
;EXAMPLE: If 
;        X = [0,0,0,0,1,1,1,1,2,2,2,2,5,5,6]
;        I = [1,2,4,5,3,7,6,11,10,9,12,13,15,14,8]
;     then the command 
;        OUT = SORTIT(X,I)
;     gives
;        OUT = [0,0,1,0,0,1,1,6,2,2,1,2,2,5,5]
;
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES): 
;     Repeatedly extracts elements of XIn at MIN(SortVec), then 
;     eliminates those elements from XIn and the minimum values from 
;     SortVec, until SortVec is empty.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     None
;  
;MODIFICATION HISTORY:
;    Written by Jack Saba, STX, Oct 1992.
;    Dave Bazell, Dec. 92  Fixed two endifs and changed example.
;
;.TITLE
; Routine SortIt
;-

; -- --- 1 --- --- 2 ---|--- 3 --- --- 4 --- --- 5 --- --- 6 --- --- 7 |

;...Check input -- must be vectors. If only one element, sort is 
;...trivial.

      IX = N_ELEMENTS ( XIn )
      IF ( IX EQ 0 ) THEN BEGIN
         MESSAGE, 'SortIt: Input vector is empty', /TRACEBACK
      ENDIF ELSE IF ( IX EQ 1 ) THEN BEGIN
         XOut = XIn
         GOTO, AllDone
      ENDIF
      IS = N_ELEMENTS ( SortVec )
      IF ( IS EQ 0 ) THEN BEGIN
         MESSAGE, 'SortIt: Input sorting vector is empty', /TRACEBACK
      ENDIF ELSE IF ( IS EQ 1 ) THEN BEGIN
         XOut = XIn
         GOTO, AllDone
      ENDIF

;...Check lengths and adjust if necessary.

      IF ( IX LT IS ) THEN BEGIN
         LocalSortVec = SortVec(0:IX-1)
         LocalXIn     = XIn
      ENDIF ELSE IF ( IS LT IX ) THEN BEGIN
         LocalSortVec = SortVec
         LocalXIn     = XIn(0:IS-1)
      ENDIF ELSE BEGIN
         LocalSortVec = SortVec
         LocalXIn     = XIn
      ENDELSE

;...Set elements corresponding to minimum values in the sorting vector.

      IS     = WHERE ( LocalSortVec EQ MIN ( LocalSortVec ) )
      IF ( IS(0) EQ -1 ) THEN BEGIN
         MESSAGE, 'SortIt: Sort vector is empty', /TRACEBACK
      ENDIF
      XOut   = LocalXIn ( IS )

;...Eliminate just-selected elements from XIn and SortVec.

      IS     = WHERE ( LocalSortVec NE MIN ( LocalSortVec ) )
      IF ( IS(0) EQ -1 ) THEN BEGIN
         RETURN, XOut
      ENDIF ELSE BEGIN
         Sorter = LocalSortVec ( IS )
         XSort  = LocalXIn ( IS )
      ENDELSE

;...Repeat the above two steps until all elements are accounted for.

      WHILE ( N_ELEMENTS(Sorter) GT 0 ) DO BEGIN

         IS     = WHERE ( Sorter EQ MIN ( Sorter ) )
         XOut   = [ XOut, XSort ( IS ) ]
         IS     = WHERE ( Sorter NE MIN ( Sorter ) )
         IF ( IS(0) EQ -1 ) THEN BEGIN
            RETURN, XOut
         ENDIF ELSE BEGIN
            Sorter = Sorter ( IS )
            XSort  = XSort ( IS )
         ENDELSE

      ENDWHILE

AllDone:

      RETURN, XOut
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


