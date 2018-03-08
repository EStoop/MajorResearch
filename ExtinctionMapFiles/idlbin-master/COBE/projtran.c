/*;DISCLAIMER:
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
;Please send bug reports to CGIS@COBECL.DNET.NASA.GOV.*/

/*******************************************************************************
/* MODULE:	projtran.c

   FUNCTION:	This program generates a projection from a sky cube 
                input.  This input can be of byte, short integer,
		long integer, float, or double data type.  The output
		array will be the same type as the input.  It is called
		by the IDL JPRO procedure.  There is a corresponding 
		IDL function by the same name with the same functionality.
		If the linkimage is not implemented, the IDL version 
		will be compiled and run instead.

   AUTHOR:	J.M. Gales (Applied Research Corp)
   DATE:	7/92, per CCR 619 & SPR 9818 to speed up reprojection.

	MODIFICATION HISTORY:
	DATE	SPR/SER		PROGRAMMER	DESCRIPTION
	----	-------		----------	-----------
     11/9/92    10212           Gales           Renamed from CPROJ
	03/94	11614		Dan Grogan	Changed "export.h" to 
						"export_wrap.h" for IDL 3.0+
       09/95                    Newmark         Upgrade to IDL 4.0


   COMMAND LINE ARGUMENT: input array, i_arr, j_arr, mask, face, min

   INPUT PARAMETERS:
	Name		Type		Description
	----		----		-----------
	input array	B,S,L,F,D arr	unfolded skycube or single face

	i_arr		int arr		corresponding col # of input 
					array for each projection array
					element

	j_arr		int arr		corresponding row # of input 
					array for each projection array
					element

	mask		byte arr	projection mask (255 if outside
					projection, face number within)

	face		int		face number (-1 if skycube)

	min		float		minimum value of input	

   RETURN:		projected output

   FUNCTIONS CALLED:	None

*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include "export_wrap.h"

IDL_VPTR projtran(argc,argv)

int	argc;
IDL_VPTR	argv[];


{

IDL_LONG	i, j, k;
short		row, col, msk_eq, msk_ne, face;
IDL_LONG	inp_dim[2], prj_dim[2];
float		min;
IDL_UCHAR		ma;
char		typ;

/* IDL pointer VARIABLEs */

IDL_UCHAR		*b_input, *b_image;
short		*s_input, *s_image;
IDL_LONG	*l_input, *l_image;
float		*f_input, *f_image;
double		*d_input, *d_image;

short		*iarr, *jarr;
IDL_UCHAR		*mask;

IDL_LONG	*dims;

IDL_VARIABLE	*retdat;


	/* set up pointers to input */

	iarr = (short *) argv[1]->value.arr->data;
			/* get pointer to iarr */

	jarr = (short *) argv[2]->value.arr->data;
			/* get pointer to jarr */

	mask = (IDL_UCHAR *) argv[3]->value.arr->data;
			/* get pointer to mask */

	face = (short) argv[4]->value.i;
			/* get face number */

	min = (float) argv[5]->value.f;
			/* get minumum value */


	dims = (IDL_LONG *) argv[0]->value.arr->dim;
			/* get input dimensions */

	for (i=0; i<=1; ++i)
		inp_dim[i] = *dims++;
			/* load dimensions */

	dims = (IDL_LONG *) argv[1]->value.arr->dim;
			/* get projection dimensions */

	for (i=0; i<=1; ++i)
		prj_dim[i] = *dims++;
			/* load dimensions */


	typ = (char) argv[0]->type;
			/* get data type of input array */


	/* setup pointer to input (skycube/face) array
	   generate temporary output (projected) array */

	switch (typ)
	
	{

	case 1:		/* byte input */

	b_input = (IDL_UCHAR *) argv[0]->value.arr->data;
			/* get pointer to input array */


	b_image = (IDL_UCHAR *)
		   IDL_MakeTempArray(IDL_TYP_BYTE,2,prj_dim, 
		   IDL_BARR_INI_ZERO, &retdat);
			/* allocate space for output array */

	break;


	case 2:		/* short integer input */

	s_input = (short *) argv[0]->value.arr->data;
			/* get pointer to input array */


	s_image = (short *)
		   IDL_MakeTempArray(IDL_TYP_INT,2,prj_dim, 
		   IDL_BARR_INI_ZERO, &retdat);
			/* allocate space for output array */

	break;


	case 3:		/* IDL_LONG integer input */

	l_input = (IDL_LONG *) argv[0]->value.arr->data;
			/* get pointer to input array */


	l_image = (IDL_LONG *)
		   IDL_MakeTempArray(IDL_TYP_LONG,2,prj_dim, 
		   IDL_BARR_INI_ZERO, &retdat);
			/* allocate space for output array */

	break;


	case 4:		/* float input */

	f_input = (float *) argv[0]->value.arr->data;
			/* get pointer to input array */


	f_image = (float *)
		   IDL_MakeTempArray(IDL_TYP_FLOAT,2,prj_dim, 
		   IDL_BARR_INI_ZERO, &retdat);
			/* allocate space for output array */

	break;


	case 5:		/* double input */

	d_input = (double *) argv[0]->value.arr->data;
			/* get pointer to input array */


	d_image = (double *)
		   IDL_MakeTempArray(IDL_TYP_DOUBLE,2,prj_dim, 
		   IDL_BARR_INI_ZERO, &retdat);
			/* allocate space for output array */

	}


	/* generate projected array */


	if (face == -1)		/* sky cube input */
	{

	switch (typ)
	
	{

	case 1:		/* byte input */

	for (i=0; i<prj_dim[1]; i++)
	{
		for (j=0; j<prj_dim[0]; j++)
		{
			k = i*prj_dim[0] + j;	/* offset for lut/proj arrays */
			col = *(iarr + k);	/* sky cube col */
			row = *(jarr + k);	/* sky cube row */
			ma = *(mask + k);	/* mask value */
			msk_ne = (ma != 255) ? 1 : 0;

			*(b_image + k) = *(b_input + row*inp_dim[0] + col) * 
					   msk_ne + (min * (1 - msk_ne));

				/* build image, fill boundary with
				   minimum value */

		}
	}

	break;


	case 2:		/* short integer input */

	for (i=0; i<prj_dim[1]; i++)
	{
		for (j=0; j<prj_dim[0]; j++)
		{
			k = i*prj_dim[0] + j;	/* offset for lut/proj arrays */
			col = *(iarr + k);	/* sky cube col */
			row = *(jarr + k);	/* sky cube row */
			ma = *(mask + k);	/* mask value */
			msk_ne = (ma != 255) ? 1 : 0;

			*(s_image + k) = *(s_input + row*inp_dim[0] + col) * 
					   msk_ne + (min * (1 - msk_ne));

				/* build image, fill boundary with
				   minimum value */

		}
	}


	break;


	case 3:		/* IDL_LONG integer input */


	for (i=0; i<prj_dim[1]; i++)
	{
		for (j=0; j<prj_dim[0]; j++)
		{
			k = i*prj_dim[0] + j;	/* offset for lut/proj arrays */
			col = *(iarr + k);	/* sky cube col */
			row = *(jarr + k);	/* sky cube row */
			ma = *(mask + k);	/* mask value */
			msk_ne = (ma != 255) ? 1 : 0;

			*(l_image + k) = *(l_input + row*inp_dim[0] + col) * 
					   msk_ne + (min * (1 - msk_ne));

				/* build image, fill boundary with
				   minimum value */

		}
	}

	break;


	case 4:		/* float input */

	for (i=0; i<prj_dim[1]; i++)
	{
		for (j=0; j<prj_dim[0]; j++)
		{
			k = i*prj_dim[0] + j;	/* offset for lut/proj arrays */
			col = *(iarr + k);	/* sky cube col */
			row = *(jarr + k);	/* sky cube row */
			ma = *(mask + k);	/* mask value */
			msk_ne = (ma != 255) ? 1 : 0;

			*(f_image + k) = *(f_input + row*inp_dim[0] + col) * 
					   msk_ne + (min * (1 - msk_ne));

				/* build image, fill boundary with
				   minimum value */

		}
	}

	break;


	case 5:		/* double input */

	for (i=0; i<prj_dim[1]; i++)
	{
		for (j=0; j<prj_dim[0]; j++)
		{
			k = i*prj_dim[0] + j;	/* offset for lut/proj arrays */
			col = *(iarr + k);	/* sky cube col */
			row = *(jarr + k);	/* sky cube row */
			ma = *(mask + k);	/* mask value */
			msk_ne = (ma != 255) ? 1 : 0;

			*(d_image + k) = *(d_input + row*inp_dim[0] + col) * 
					   msk_ne + (min * (1 - msk_ne));

				/* build image, fill boundary with
				   minimum value */

		}
	}


	}


	}
	else			/* single face input */
	{

	switch (typ)
	
	{

	case 1:		/* byte input */

	for (i=0; i<prj_dim[1]; i++)
	{
		for (j=0; j<prj_dim[0]; j++)
		{
			k = i*prj_dim[0] + j;	/* offset for lut/proj arrays */
			col = *(iarr + k);	/* face col */
			row = *(jarr + k);	/* face row */
			ma = *(mask + k);	/* mask value */
			msk_eq = (ma == face) ? 1 : 0;

			*(b_image + k) = *(b_input + row*inp_dim[0] + col) * 
					   msk_eq + (min * (1 - msk_eq));

				/* build image, fill boundary with
				   minimum value */

		}
	}

	break;


	case 2:		/* short integer input */

	for (i=0; i<prj_dim[1]; i++)
	{
		for (j=0; j<prj_dim[0]; j++)
		{
			k = i*prj_dim[0] + j;	/* offset for lut/proj arrays */
			col = *(iarr + k);	/* face col */
			row = *(jarr + k);	/* face row */
			ma = *(mask + k);	/* mask value */
			msk_eq = (ma == face) ? 1 : 0;

			*(s_image + k) = *(s_input + row*inp_dim[0] + col) * 
					   msk_eq + (min * (1 - msk_eq));

				/* build image, fill boundary with
				   minimum value */

		}
	}


	break;


	case 3:		/* IDL_LONG integer array */


	for (i=0; i<prj_dim[1]; i++)
	{
		for (j=0; j<prj_dim[0]; j++)
		{
			k = i*prj_dim[0] + j;	/* offset for lut/proj arrays */
			col = *(iarr + k);	/* face col */
			row = *(jarr + k);	/* face row */
			ma = *(mask + k);	/* mask value */
			msk_eq = (ma == face) ? 1 : 0;

			*(l_image + k) = *(l_input + row*inp_dim[0] + col) * 
					   msk_eq + (min * (1 - msk_eq));

				/* build image, fill boundary with
				   minimum value */

		}
	}

	break;


	case 4:		/* float input */

	for (i=0; i<prj_dim[1]; i++)
	{
		for (j=0; j<prj_dim[0]; j++)
		{
			k = i*prj_dim[0] + j;	/* offset for lut/proj arrays */
			col = *(iarr + k);	/* face col */
			row = *(jarr + k);	/* face row */
			ma = *(mask + k);	/* mask value */
			msk_eq = (ma == face) ? 1 : 0;

			*(f_image + k) = *(f_input + row*inp_dim[0] + col) * 
					   msk_eq + (min * (1 - msk_eq));

				/* build image, fill boundary with
				   minimum value */

		}
	}

	break;


	case 5:		/* double input */

	for (i=0; i<prj_dim[1]; i++)
	{
		for (j=0; j<prj_dim[0]; j++)
		{
			k = i*prj_dim[0] + j;	/* offset for lut/proj arrays */
			col = *(iarr + k);	/* face col */
			row = *(jarr + k);	/* face row */
			ma = *(mask + k);	/* mask value */
			msk_eq = (ma == face) ? 1 : 0;

			*(d_image + k) = *(d_input + row*inp_dim[0] + col) * 
					   msk_eq + (min * (1 - msk_eq));

				/* build image, fill boundary with
				   minimum value */

		}
	}


	}


	}
return(retdat);
}
