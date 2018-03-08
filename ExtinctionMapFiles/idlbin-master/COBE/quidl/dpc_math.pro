FUNCTION dpc_math,in0,in1,in2,code=code

;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    DPC_MATH performs double precision mathmatical operations.
;
;DESCRIPTION:
;    This function performs a number of mathmatical operations on
;    pseudo-complex IDL arrays.  These arrays are double precision
;    IDL arrays of dimension (*,n) with the real part stored in first
;    n rows and the imaginary part stored in the second n rows.  The
;    user must first convert the IDL complex arrays to the double
;    precision arrays (either using the 'c2d' feature of this routine
;    or manually) before using this routine.
;
;
;CALLING SEQUENCE:
;    out = dpc_math(arr1,[arr2],code=code)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    arr1        I    double,complex      input array 1
;   [arr2]       I    double,complex      input array 2
;    code        I    string              operation code
;    out         O    double,complex      output array
;
;
;EXAMPLES:
;
; Convert from IDL COMPLEXARR to Double precision array
; -----------------------------------------------------
;    cplx = COMPLEX(real,imag)
;
;    d_cplx = dpc_math(cplx,code='c2d')
;
;    d_cplx = [[DOUBLE(real)],[DOUBLE(imag)]]
;
;
;
; Matrix Multiplication
; ---------------------
;    prod = dpc_math(d_cplx1,d_cplx2,code='m')
;
;    d_cplx1 = n1 by 2*n2 double array 
;	(representing an n1 by n2 complex array)
;
;    d_cplx2 = n2 by 2*n3 double array 
;	(representing an n2 by n3 complex array)
;
;    prod    = n1 by 2*n3 double array 
;	(representing an n1 by n3 complex array)
;
;
;
; Matrix Inversion
; ----------------
;    inv = dpc_math(d_cplx,code='i')
;
;    d_cplx = n1 by 2*n1 double array 
;	(representing an n1 by n1 complex array)
;
;    inv    = n1 by 2*n1 double array 
;	(representing an n1 by n1 complex array)
;
;    dpc_math(d_cplx,inv,code='m') = [[Identity],[Zero]]
;
;
;
; Complex Conjugation
; -------------------
;    conjg = dpc_math(d_cplx,code='c')
;
;    d_cplx = [[real],[+imag]],  n1 by 2*n2 double array 
;	(representing an n1 by n2 complex array)
;
;    conjg  = [[real],[-imag]],  n1 by 2*n2 double array 
;	(representing an n1 by n2 complex array)
;
;
;
; Matrix Transposition
; --------------------
;    trans = dpc_math(d_cplx,code='t')
;
;    d_cplx   = [[real],[imag]], n1 by 2*n2 double array 
;	(representing an n1 by n2 complex array)
;
;    trans    = [[transpose(real)],[transpose(imag)]], 
;			n2 by 2*n1 double array 
;			(representing an n2 by n1 complex array)
;
;
;
; Hermitian Conjugation
; ---------------------
;    h_conjg = dpc_math(d_cplx,code='h')
;
;    d_cplx   = [[real],[imag]], n1 by 2*n2 double array 
;	(representing an n1 by n2 complex array)
;
;    h_conjg  = [[transpose(real)],[-transpose(imag)]], 
;			n2 by 2*n1 double array 
;			(representing an n2 by n1 complex array)
;
;
;
; Diagonal Elements
; -----------------
;    diag = dpc_math(d_cplx,code='d')
;
;    d_cplx = n1 by 2*n2 double array 
;	(representing an n1 by n2 complex array)
;
;    diag   = n1 by 2 double array
;	(representing an n1 complex vector)
;
;
;
; Real Part
; ---------
;    real = dpc_math(d_cplx,code='re')
;
;    d_cplx   = [[real],[imag]], n1 by 2*n2 double array 
;	(representing an n1 by n2 complex array)
;
;    real     = [real], n1 by n2 double array
;
;
;
; Imaginary Part
; --------------
;    imag = dpc_math(d_cplx,code='d')
;
;    d_cplx   = [[real],[imag]], n1 by 2*n2 double array 
;	(representing an n1 by n2 complex array)
;
;    imag     = [imag], n1 by n2 double array
;
;
;
; Convert from Double precision array to IDL COMPLEXARR
; -----------------------------------------------------
;    d_cplx = [[real],[imag]]
;
;    cplx = dpc_math(cplx,code='d2c')
;
;    cplx = COMPLEX(FLOAT(real),FLOAT(imag))
;
;
;#
;COMMON BLOCKS:
;    None
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    None
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Jan 94
;    Initial Delivery   SPR XXXXX
;
;
;.TITLE
;Routine DPC_MATH
;-



; Convert from IDL COMPLEXARR to Double precision array
;
; cpx(*,n) --> dbl(*,2n) where first n lines contains real part
;                        and second n lines contains imag part
; ------------------------------------------------------------
IF (STRUPCASE(code) EQ 'C2D') THEN BEGIN

	sz = SIZE(in0)
	IF sz(0) EQ 1 THEN n_row=1 ELSE n_row = sz(2)
	n = 2*n_row
	out = DBLARR(sz(1),n)
	FOR i=0,n_row-1 DO BEGIN
		out(*,i) = DOUBLE(in0(*,i))
		out(*,i+n_row) = DOUBLE(IMAGINARY(in0(*,i)))
	ENDFOR
	GOTO, exit

ENDIF



; Convert from Double precision array to IDL COMPLEXARR
; -----------------------------------------------------
IF (STRUPCASE(code) EQ 'D2C') THEN BEGIN

	sz = SIZE(in0)
	n = sz(2)/2
	out = COMPLEXARR(sz(1),n)
	FOR i=0,n-1 DO out(*,i) = COMPLEX(in0(*,i),in0(*,i+n))
	GOTO, exit

ENDIF


cde = STRUPCASE(STRMID(code,0,1))



; Extract Real Part
; -----------------
IF (STRUPCASE(code) EQ 'RE') THEN BEGIN

	sz = SIZE(in0)

	IF (sz(0) EQ 1) THEN out = in0(0:sz(1)/2-1) ELSE $
			     out = in0(*,0:sz(2)/2-1)
	GOTO, exit

ENDIF



; Extract Imaginary Part
; ----------------------
IF (STRUPCASE(code) EQ 'IM') THEN BEGIN

	sz = SIZE(in0)

	IF (sz(0) EQ 1) THEN out = in0(sz(1)/2:*) ELSE $
			     out = in0(*,sz(2)/2:*)
	GOTO, exit

ENDIF


cde = STRUPCASE(STRMID(code,0,1))

CASE cde OF 


; Matrix Multiplication
; ---------------------
'M':	BEGIN

		sz0 = SIZE(in0)
		sz1 = SIZE(in1)

		IF (sz0(0) EQ 1) THEN BEGIN
			a_real = in0(0:(sz0(1)/2)-1)
			a_imag = in0(sz0(1)/2:*)
		ENDIF ELSE BEGIN
			a_real = in0(*,0:(sz0(2)/2)-1)
			a_imag = in0(*,sz0(2)/2:*)
		ENDELSE

		IF (sz1(0) EQ 1) THEN BEGIN
			b_real = in1(0:(sz1(1)/2)-1)
			b_imag = in1(sz1(1)/2:*)
		ENDIF ELSE BEGIN
			b_real = in1(*,0:(sz1(2)/2)-1)
			b_imag = in1(*,sz1(2)/2:*)
		ENDELSE

		c_real = a_real#b_real - a_imag#b_imag
		c_imag = a_imag#b_real + a_real#b_imag

		out = [[c_real],[c_imag]]


	END




; Diagonal Matrix Multiplication
; ------------------------------
'G':	BEGIN

		sz0 = SIZE(in0)		; regular matrix
		sz1 = SIZE(in1)		; diagonal of diag matrix


		IF (sz0(0) EQ 1) THEN BEGIN
			a_real = in0(0:(sz0(1)/2)-1)
			a_imag = in0(sz0(1)/2:*)
		ENDIF ELSE BEGIN
			a_real = in0(*,0:(sz0(2)/2)-1)
			a_imag = in0(*,sz0(2)/2:*)
		ENDELSE

		b_real = in1(*,0)
		b_imag = in1(*,1)

		c_real = DBLARR(sz0(1),sz0(2)/2)
		c_imag = c_real

		FOR i=0,sz0(2)/2-1 DO $
			c_real(*,i) = a_real(*,i)*b_real(i) - $
				      a_imag(*,i)*b_imag(i)

		FOR i=0,sz0(2)/2-1 DO $
			c_imag(*,i) = a_imag(*,i)*b_real(i) + $
				      a_real(*,i)*b_imag(i)

		out = [[c_real],[c_imag]]


	END



; Transpose
; ---------
'T':	BEGIN

		sz = SIZE(in0)

		a_real = in0(*,0:(sz(2)/2)-1)
		a_imag = in0(*,sz(2)/2:*)

		out = [[TRANSPOSE(a_real)],[TRANSPOSE(a_imag)]]

	END



; Hermetian Conjugate (Transpose + Complex Conjugate)
; ---------------------------------------------------
'H':	BEGIN

		sz = SIZE(in0)

		a_real = in0(*,0:(sz(2)/2)-1)
		a_imag = in0(*,sz(2)/2:*)

		out = [[TRANSPOSE(a_real)],[-TRANSPOSE(a_imag)]]

	END




; Complex Conjugate
; -----------------
'C':	BEGIN

		sz = SIZE(in0)

		a_real = in0(*,0:(sz(2)/2)-1)
		a_imag = in0(*,sz(2)/2:*)

		out = [[a_real],[-a_imag]]

	END




; Get Diagonal Elements
; ---------------------
'D':	BEGIN

		sz = SIZE(in0)

		a_real = in0(*,0:(sz(2)/2)-1)
		a_imag = in0(*,sz(2)/2:*)

		diag_r = a_real(INDGEN(sz(1))*(sz(1)+1))
		diag_i = a_imag(INDGEN(sz(1))*(sz(1)+1))

		out = [[diag_r],[diag_i]]

	END



; Matrix Inversion
; ----------------
'I':	BEGIN

		sz = SIZE(in0)
		n = sz(2)/2

		real = in0(*,0:n-1)
		imag = in0(*,n:*)

		inv_real = INVERT(real,stat)

;		IF (stat NE 0) THEN $
;			MESSAGE,'Bad Inversion Status: '+STRING(stat),/CONT

		q = (imag # inv_real # imag) + real
		b_real = INVERT(q,stat)

;		IF (stat NE 0) THEN $
;			MESSAGE,'Bad Inversion Status:'+STRING(stat),/CONT

		b_imag = - b_real # imag # inv_real

		out = [[b_real],[b_imag]]

	END



; Take Inner Product
; ------------------
'N':	BEGIN

		sz = SIZE(in1)

		o_r = TRANSPOSE(in0(*,0))
		o_i = TRANSPOSE(in0(*,1))

		arr_r = in1(*,0:(sz(2)/2)-1)
		arr_i = in1(*,(sz(2)/2):*)

		vec1_r = in2(*,0)
		vec1_i = in2(*,1)

		rrr = o_r # arr_r # vec1_r
		iir = o_i # arr_i # vec1_r
		rii = o_r # arr_i # vec1_i
		iri = o_i # arr_r # vec1_i

		real = rrr(0)+iir(0)-rii(0)+iri(0)

		rri = o_r # arr_r # vec1_i
		rir = o_r # arr_i # vec1_r
		iii = o_i # arr_i # vec1_i
		irr = o_i # arr_r # vec1_r

		imag = rri(0)+iii(0)+rir(0)-irr(0)

		out = [real,imag]

	END

ENDCASE

exit:

RETURN,out
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


