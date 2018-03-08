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

/****************************************************************************

/* MODULE:	rastr.c

   FUNCTION:	This program rasterizes a pixel and data list.
		Called by IDL PIX2XY routine.  There is an IDL
		procedure by the same name with the same function-
		ality.  If the linkimage is not implemented, the
		IDL routine will be compiled and run instead.

   AUTHOR:	J.M. Gales (Applied Research Corp)
   DATE:	8/92, SPR 10025, Rewrite in "C" for faster BUILDMAP execution.

	MODIFICATION HISTORY:
	DATE	SPR/SER		PROGRAMMER	DESCRIPTION
	----	-------		----------	-----------
     11/9/92    10212           Gales           Renamed from CRASTER
     03/2/93    10643           Gales           Renamed from RASTER
     03/22/93   10729           Gales           Initialize raster array
                                                to bad pixel value
      02/12/93  11468           Newmark         Change to compile under
                                                VAX and Sun C not GNU C
	03/94	11614		Dan Grogan	Changed "export.h" to 
						"export_wrap.h" for IDL 3.0+
       09/95                    Newmark         Upgrade to IDL 4.0


   COMMAND LINE ARGUMENT: pixel array, resolution, face #, sixpack flag,
                          data array, raster output, column array, row array,
			  bad pixel value

   INPUT PARAMETERS:

	Name		Type		Description
	----		----		-----------
	pixel array	IDL_LONG arr	Input pixel array

	resolution	int		skymap resolution

	face number	int		face number (-1 if whole cube)

	sixpack flag	int		0 = unfolded raster
					1 = sixpacked raster

	data array	B,S,L,F,D,C arr	Input data array

	raster output	B,S,L,F,D,C arr	Rasterized image

	column array	int arr		array of sky map column #

	row array	int arr		array of sky map row #

	bad pixel value	flt		bad pixel value

   RETURN:		(void)

   FUNCTIONS CALLED:	None

*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include "export_wrap.h"

void rastr(argc,argv)

int	argc;
IDL_VPTR	argv[];

{

int		i,j,bit,ras_bld,sixpack;
int		res,face,col,row,len,hgt,thrd;
int		offx[6],offy[6];
IDL_LONG	out_dim[2],num_pix_face,cube_side,fpix,one,pix;
IDL_LONG	num_pix,tot_pix,ras_dim[3],off_r,off_p,k;
unsigned int	pow_2[16];
char		dat_typ;
IDL_UCHAR		n_dims;
float		bad_pixval;

/* IDL VARIABLEs */

IDL_UCHAR		*b_dat, *b_raster;
short		*i_dat, *i_raster;
IDL_LONG	*l_dat, *l_raster;
float		*f_dat, *f_raster;
double		*d_dat, *d_raster;
IDL_COMPLEX		*c_dat, *c_raster;


short		*x_out, *y_out;
IDL_LONG	*pixel,*dims;

IDL_VARIABLE	*ret_r, *ret_x, *ret_y;

one=1;
ras_bld= -1;

	pixel = (IDL_LONG *) argv[0]->value.arr->data;
			/* get pointer to pixel numbers */

	num_pix = (IDL_LONG) argv[0]->value.arr->n_elts;
			/* get number of pixels */

	res = (short) argv[1]->value.i;
			/* get resolution */

	face = (short) argv[2]->value.i;
			/* get face value */

	sixpack = (short) argv[3]->value.i;
			/* get sixpack value */

	dat_typ = (char) argv[4]->type;
			/* get data type */

	bad_pixval = (float) argv[8]->value.f;
			/* get bad pixel value */


	n_dims = (IDL_UCHAR) argv[4]->value.arr->n_dim;
	dims = (IDL_LONG *) argv[4]->value.arr->dim;
	thrd = (n_dims == 1) ? 1 : *(dims + n_dims - 1);
		/* get number of data dimensions
		   read dimensions
		   if n_dims=1 then thrd=1 (2D output array)
		   else thrd=2nd dimension of input data 
                   (3D output raster) */

	if ((IDL_LONG) argv[4]->value.arr->n_elts == 1) ras_bld = 0;
		/* if single element in data array then no raster output */

	cube_side = (short) (pow(2.,res-1.));
	num_pix_face = (IDL_LONG) (pow(4.,res-1.));
			/* determine # of pixels on a side & in a face */



	/* unfolded skycube dimensions and offsets */

	if (face == -1 && sixpack == 0)
	{
		tot_pix = 12 * num_pix_face;
		len = 4 * cube_side;
		hgt = 3 * cube_side;
		offx[0]=0; offx[1]=0; offx[2]=1;
		offx[3]=2; offx[4]=3; offx[5]=0;
		offy[0]=2; offy[1]=1; offy[2]=1;
		offy[3]=1; offy[4]=1; offy[5]=0;
	}


	/* sixpacked skycube dimensions and offsets */

	if (face == -1 && sixpack == 1)
	{
		tot_pix = 6 * num_pix_face;
		len = 3 * cube_side;
		hgt = 2 * cube_side;
		offx[0]=0; offx[1]=0; offx[2]=1;
		offx[3]=2; offx[4]=2; offx[5]=1;
		offy[0]=1; offy[1]=0; offy[2]=0;
		offy[3]=0; offy[4]=1; offy[5]=1;
	}


	/* single face dimensions and offsets */

	if (face != -1)
	{
		tot_pix = num_pix_face;
		len = cube_side;
		hgt = cube_side;
		offx[0]=0; offx[1]=0; offx[2]=0;
		offx[3]=0; offx[4]=0; offx[5]=0;
		offy[0]=0; offy[1]=0; offy[2]=0;
		offy[3]=0; offy[4]=0; offy[5]=0;
	}
		
	ras_dim[0] = len;
	ras_dim[1] = hgt;
	ras_dim[2] = thrd;
		/* output raster dimensions */


	for (i=0; i<=15; i++) pow_2[i] = one << i;
		/* build power of 2 array (bit mask) */


	x_out = (short *)
		IDL_MakeTempArray(IDL_TYP_INT,1,argv[0]->value.arr->dim,
	        IDL_BARR_INI_ZERO,&ret_x);

	IDL_VarCopy(ret_x,argv[6]);
		/* build col output array and 
		   associate with 6th VARIABLE  */



	y_out = (short *)
		IDL_MakeTempArray(IDL_TYP_INT,1,argv[0]->value.arr->dim,
	        IDL_BARR_INI_ZERO,&ret_y);

	IDL_VarCopy(ret_y,argv[7]);
		/* build row output array and 
		   associate with 6th VARIABLE  */




	/* decode pixel numbers into face, col, & row numbers */

	for (pix=0; pix<num_pix; pix++)
	{

	face = (int) (*(pixel + pix) / num_pix_face);
	fpix = *(pixel + pix) - num_pix_face*face;
				/* get face and pixel # within face */

	col = 0;		/* zero column # */
	row = 0;		/* zero row # */

		for (bit=0; bit<=res-1; bit++)
		{
		col = col | (int) (pow_2[bit] * (one & fpix));
		fpix = fpix >> 1;
		row = row | (int) (pow_2[bit] * (one & fpix));
		fpix = fpix >> 1;
		}

	*(x_out + pix) = len - (offx[face] * cube_side + col + 1);
	*(y_out + pix) = offy[face] * cube_side + row;
				/* determine row & col of raster
				   and store in output arrays */


	}


	if (ras_bld == -1)	/* if data supplied then build raster */
	{

	switch (dat_typ)
	{

	case 1:		/* byte data */

		b_dat = (IDL_UCHAR *) argv[4]->value.arr->data;

		b_raster = (IDL_UCHAR *)
			   IDL_MakeTempArray(IDL_TYP_BYTE,3,ras_dim, 
			   IDL_BARR_INI_ZERO, &ret_r);
				/* allocate space for raster output */

		for (i=0; i<ras_dim[0]*ras_dim[1]*ras_dim[2]; i++)
			*(b_raster + i) = bad_pixval;
				/* initialize to bad pixel value */

		IDL_VarCopy(ret_r,argv[5]);
				/* associate with 5th VARIABLE */

		for (i=0; i<thrd; i++)	/* for all data sets */
		{
		off_r = i * len * hgt;
			/* raster offset */
		off_p = i * num_pix;
			/* pixel offset */

			for (pix=0; pix<num_pix; pix++)
			{
			k = *(y_out + pix) * len + *(x_out + pix);
				/* offset within raster */

			*(b_raster+off_r+k) = *(b_dat+off_p+pix);
				/* store data into raster */

			}
		}

	break;

	case 2:		/* short integer data */

		i_dat = (short *) argv[4]->value.arr->data;

		i_raster = (short *)
			   IDL_MakeTempArray(IDL_TYP_INT,3,ras_dim, 
			   IDL_BARR_INI_ZERO, &ret_r);
				/* allocate space for raster output */

		for (i=0; i<ras_dim[0]*ras_dim[1]*ras_dim[2]; i++)
			*(i_raster + i) = bad_pixval;
				/* initialize to bad pixel value */

		IDL_VarCopy(ret_r,argv[5]);
				/* associate with 5th VARIABLE */

		for (i=0; i<thrd; i++)	/* for all data sets */
		{
		off_r = i * len * hgt;
			/* raster offset */
		off_p = i * num_pix;
			/* pixel offset */

			for (pix=0; pix<num_pix; pix++)
			{
			k = *(y_out + pix) * len + *(x_out + pix);
			*(i_raster+off_r+k) = *(i_dat+off_p+pix);
			}
		}

	break;

	case 3:		/* IDL_LONG integer data */

		l_dat = (IDL_LONG *) argv[4]->value.arr->data;

		l_raster = (IDL_LONG *)
			   IDL_MakeTempArray(IDL_TYP_LONG,3,ras_dim, 
			   IDL_BARR_INI_ZERO, &ret_r);
				/* allocate space for raster output */

		for (i=0; i<ras_dim[0]*ras_dim[1]*ras_dim[2]; i++)
			*(l_raster + i) = bad_pixval;
				/* initialize to bad pixel value */

		IDL_VarCopy(ret_r,argv[5]);
				/* associate with 5th VARIABLE */

		for (i=0; i<thrd; i++)	/* for all data sets */
		{
		off_r = i * len * hgt;
			/* raster offset */
		off_p = i * num_pix;
			/* pixel offset */

			for (pix=0; pix<num_pix; pix++)
			{
			k = *(y_out + pix) * len + *(x_out + pix);
			*(l_raster+off_r+k) = *(l_dat+off_p+pix);
			}
		}

	break;

	case 4:		/* float data */


		f_dat = (float *) argv[4]->value.arr->data;

		f_raster = (float *)
			   IDL_MakeTempArray(IDL_TYP_FLOAT,3,ras_dim, 
			   IDL_BARR_INI_ZERO, &ret_r);
				/* allocate space for raster output */

		for (i=0; i<ras_dim[0]*ras_dim[1]*ras_dim[2]; i++)
			*(f_raster + i) = bad_pixval;
				/* initialize to bad pixel value */

		IDL_VarCopy(ret_r,argv[5]);
				/* associate with 5th VARIABLE */

		for (i=0; i<thrd; i++)	/* for all data sets */
		{
		off_r = i * len * hgt;
			/* raster offset */
		off_p = i * num_pix;
			/* pixel offset */

			for (pix=0; pix<num_pix; pix++)
			{
			k = *(y_out + pix) * len + *(x_out + pix);
			*(f_raster+off_r+k) = *(f_dat+off_p+pix);
			}
		}

	break;


	case 5:		/* double data */

		d_dat = (double *) argv[4]->value.arr->data;

		d_raster = (double *)
			   IDL_MakeTempArray(IDL_TYP_DOUBLE,3,ras_dim, 
			   IDL_BARR_INI_ZERO, &ret_r);
				/* allocate space for raster output */

		for (i=0; i<ras_dim[0]*ras_dim[1]*ras_dim[2]; i++)
			*(d_raster + i) = bad_pixval;
				/* initialize to bad pixel value */

		IDL_VarCopy(ret_r,argv[5]);
				/* associate with 5th VARIABLE */

		for (i=0; i<thrd; i++)	/* for all data sets */
		{
		off_r = i * len * hgt;
			/* raster offset */
		off_p = i * num_pix;
			/* pixel offset */

			for (pix=0; pix<num_pix; pix++)
			{
			k = *(y_out + pix) * len + *(x_out + pix);
			*(d_raster+off_r+k) = *(d_dat+off_p+pix);
			}
		}

	break;


	case 6:

		c_dat = (IDL_COMPLEX *) argv[4]->value.arr->data;

		c_raster = (IDL_COMPLEX *)
			   IDL_MakeTempArray(IDL_TYP_COMPLEX,3,ras_dim, 
			   IDL_BARR_INI_ZERO, &ret_r);
				/* allocate space for raster output */


		for (i=0; i<ras_dim[0]*ras_dim[1]*ras_dim[2]; i++)
			(*(c_raster + i)).r = bad_pixval;
				/* initialize to bad pixel value */

		IDL_VarCopy(ret_r,argv[5]);
				/* associate with 5th VARIABLE */

		for (i=0; i<thrd; i++)	/* for all data sets */
		{
		off_r = i * len * hgt;
			/* raster offset */
		off_p = i * num_pix;
			/* pixel offset */

			for (pix=0; pix<num_pix; pix++)
			{
			k = *(y_out + pix) * len + *(x_out + pix);
			(*(c_raster+off_r+k)).r = (*(c_dat+off_p+pix)).r;
			(*(c_raster+off_r+k)).i = (*(c_dat+off_p+pix)).i;
			}
		}


	break;


	}

	}
return;
}
