Pro PIX_TO_XY, Resolution, Pixno, Orientation, Origin, Ix, Iy, Nface
;
;-----------------------------------------------------------------
;   Subroutine to return the location of a quadtree pixel on an
;   unfolded cube.  This routine assumes the following ordering.
;                   F0
;       F4  F3  F2  F1
;                   F5
;   The output x,y coordinates are counted relative to one.
;   The unfolded cube fits on a 1024 x 768 screen without
;   loss of DIRBE resolution.  Variable sized pixels are handled
;   by this routine via the input variable resolution.  The value
;   of resolution must be set to the quadtree index level, .i.e.
;   for DIRBE maps use resolution = 9, FIRAS and DMR use resolution
;   set to 6.
;   
;
;   Jeff Pedelty, ARC, August, 1988
;   modified for sky viewing, January, 1988
;
;------------------------------------------------------------------------
;
; Inputs: 
;       Resolution         6  or  9
;       Pixno              Integer value of Pixel Number
;       Orientation        'L' or 'R' (Left or Right)
;       Origin             'F' or 'U' (Single face or Unfolded)
; Outputs:
;         Ix               X-coordinate of pixel
;         Iy               Y-coordinate of pixel
;         Nface            Face number of pixel
; Modification History
;
; S.R.K. VIDYA SAGAR
;
; MAY 1991
;
; Changes:  Changed the code from FORTRAN to IDL_2 and added more 
;           arguments in the subroutine.
;
;
;--------------------------------------------------------------------
;
ix = 0
iy = 0
nface = 0
face_width = 0
face_area = 0
ipix = 0
jpix = 0
ip = 0
id = 0
tot_pix = 0
tot_width = 0
;
origin = strupcase(origin)
first_time = 'T'
      if (first_time) then begin
         face_width = long(2 ^(resolution - 1))
         face_area  = long(1.0 * (face_width ^ 2.0))
         tot_pix    = long(6.0 * face_area)
         if (origin eq 'U') then begin
             tot_width  = long(4.0 * face_width)
         endif else begin
             tot_width  = long(face_width)
         endelse
         first_time = 'F'
      endif
;
      if (pixno gt tot_pix-1) then begin
         print,' pixno ',pixno, 'tot_pix-1 ',tot_pix-1
         print , 'ERROR: PIXEL NUMBER IS TOO LARGE FOR GIVEN RESOLUTION'
         stop
      endif
;
      nface = long(pixno / face_area)              ; Face number (0 rel)
      ipix  = long(pixno - nface * face_area)      ; Pixel num within face
;
;  Compute the X and Y position of the pixel within its cube face
;  by breaking up the pixel number.  This block of code is lifted 
;  directly from E. Wright's DMR_CENPIX routine, via the code for
;  USA_TPL_INTPIX
;
      jpix = ipix
      ix = 0
      iy = 0
      ip = 1
;
    if (jpix gt 0) then begin
      repeat begin
         id = (jpix) - ((jpix/2)*2)
         jpix = jpix / 2
         ix = id * ip + ix
         id = (jpix) - ((jpix/2)*2)
         jpix = jpix / 2
         iy = id * ip + iy
         ip = 2 * ip
       endrep until (jpix le 0)
    endif
;
;  Now calculate the position of the pixel in the unfolded cube.
;  
 if (origin eq 'U') then begin
      if (nface eq 2) then  ix = ix + face_width
      if (nface eq 3) then  ix = ix + 2 * face_width
      if (nface eq 4) then  ix = ix + 3 * face_width
;
      if (nface lt 5) then begin
         iy = iy + face_width
         if (nface eq 0) then iy = iy + face_width
      endif
endif
;
;   Take care of to flip to astronomical convention (skyward looking)
;      and to 1-relative pixel ordering
;
orientation = strupcase(orientation)
if (orientation eq 'R') then begin
      ix = tot_width - ix - 1
      iy = iy
endif
;
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


