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
;Please send bug reports to CGIS@COBECL.DNET.NASA.GOV. */

/*******************************************************************************
/* MODULE:	pixavg.c

   FUNCTION:	This program averages the multiple observations per
                pixel for DIRBE skymaps.
		Called by IDL BUILDMAP routine.

   AUTHOR:	J.M. Gales (Applied Research Corp)
   DATE:	8/92, SPR 10025, Rewrite in "C" for faster BUILDMAP execution.

	MODIFICATION HISTORY:
	DATE	SPR/SER		PROGRAMMER	DESCRIPTION
	----	-------		----------	-----------
	5/25/93	10962		Gales		Remove sentinal values
						from average
	03/94	11614		Dan Grogan	Changed "export.h" to 
						"export_wrap.h" for IDL 3.0+
       09/95                    Newmark         Upgrade to IDL 4.0


   COMMAND LINE ARGUMENT: input pixel array, input data array, sent_val

   INPUT PARAMETERS:

	Name		Type		Description
	----		----		-----------
	pixel array	IDL_LONG arr	Input (non-condensed) pixel array
					(Assumes input pixels are sorted)

	data array	B,S,L,F,D arr	Input (non-condensed) data array

	sentinal value	float		Sentinal trigger value
					If x < sent_val then considered
					a sentinal value and not 
					averaged in.


   RETURN:		  condensed pixel array, averaged data array
			  (passed back through corresponding input arrays)

   FUNCTIONS CALLED:	None

*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include "export_wrap.h"

void pixavg(argc,argv)

int	argc;
IDL_VPTR	argv[];

{

int		i,j;
IDL_LONG	pix,n_sent;
IDL_LONG	num_pix,num_dat,tot_pix,dim[3],brk,cnt;
float		sent_val;
char		dat_typ;
IDL_UCHAR		n_dims,sec;

/* IDL VARIABLEs */

IDL_UCHAR		*b_dat, *b_newdat, b_dat0;
short		*i_dat, *i_newdat, i_dat0;
IDL_LONG	*l_dat, *l_newdat, l_dat0;
float		*f_dat, *f_newdat, f_dat0;
double		*d_dat, *d_newdat, d_dat0;

short		*n_obs,*n_obs0;
IDL_LONG	*pixel,*dims,*new_pix;

IDL_VARIABLE	*ret_pix, *ret_n, *ret_dat;


	pixel = (IDL_LONG *) argv[0]->value.arr->data;
			/* get pointer to pixel numbers */

	num_pix = (IDL_LONG) argv[0]->value.arr->n_elts;
			/* get number of pixels */

	dat_typ = (char) argv[1]->type;
			/* get data type */

	n_dims = (IDL_UCHAR) argv[1]->value.arr->n_dim;
	dims = (IDL_LONG *) argv[1]->value.arr->dim;
	sec = (n_dims == 1) ? 1 : *(dims + n_dims - 1);
			/* get number of data dimensions
			   read dimensions
			   if n_dims=1 then sec=1 else = 2nd dim */


	sent_val = (float) argv[2]-> value.f;
			/* get sentinal value (upper bound) */


	/* count number of distinct pixels */

	brk = 0;
	cnt = 1;
	for (pix=0; pix<num_pix; pix++)
	{
		if (*(pixel + pix) != *(pixel + brk))
			/* if change in pixel number */
		{
			brk = pix;
				/* set new pixel break value */
			cnt++;
				/* increment count */
		}
	}



	dim[0] = cnt;
	dim[1] = sec;
		/* set dimensions of output arrays */

	new_pix = (IDL_LONG *)
		  IDL_MakeTempArray(IDL_TYP_LONG,1,dim,
	          IDL_BARR_INI_ZERO,&ret_pix);
		/* allocate condensed pixel array */



	n_obs = (short *)
	        IDL_MakeTempArray(IDL_TYP_INT,1,dim,
	        IDL_BARR_INI_ZERO,&ret_n);
		/* allocate # of observations per pixel array */

	n_obs0 = n_obs;
		/* set auxilliary pointer to this array */




	/* condense pixel array and
	   count # of observations of each pixel */

	brk = 0;
	*n_obs0 = 0;
	*new_pix = *pixel;
	for (pix=0; pix<num_pix; pix++)
	{
		if (*(pixel + pix) != *(pixel + brk))
			/* if change in pixel number */
		{
			brk = pix;
			/* set new break point value */

			*(++new_pix) = *(pixel + pix);
			/* increment new pixel ptr
			   store new pixel value */

			*(++n_obs0) = 1;
			/* initialize # of observations */
		}
		else
		{
			(*n_obs0)++;
			/* increment # of observations */
		}
	}


	switch (dat_typ)
	{

	case 1:		/* byte data */

		b_dat = (IDL_UCHAR *) argv[1]->value.arr->data;
				/* get pointer to raw data */

		b_newdat = (IDL_UCHAR *)
			   IDL_MakeTempArray(IDL_TYP_BYTE,2,dim, 
			   IDL_BARR_INI_ZERO, &ret_dat);
				/* allocate space for averaged output */


		for (i=0; i<sec; i++)
		{
		n_obs0 = n_obs;	
		brk = 0;
		n_sent = 0;
		for (pix=0; pix<num_pix; pix++)
		{	/* for all pixels in raw pixel list */

			if (*(pixel + pix) != *(pixel + brk))
			{	/* if change in pixel # */

				brk = pix;
				if (*n_obs0 != n_sent)
				*b_newdat = *b_newdat / (*n_obs0 - n_sent);
				/* divide by number of observations 
				   minus number of sentinals */

				n_obs0++;
				/* increment # of observations  */

				b_dat0 = *(b_dat + i*num_pix + pix);
				/* get next input data */

				b_newdat++;
				/* increment point to accumulated data */

				n_sent = 0;
				/* initialize # of sentinal values */

				if (b_dat0 > sent_val)
					*b_newdat = b_dat0;
				else
					n_sent++;
				/* if data value greater than sentinal
				   value then initialize accumulated
				   data else increment number of 
				   sentinal values. */

			}
			else
			{
				b_dat0 = *(b_dat + i*num_pix + pix);
				if (b_dat0 > sent_val)
					*b_newdat = *b_newdat + b_dat0;
				else
					n_sent++;
				/* if data value greater than sentinal
				   value then add to accumulated
				   data else increment number of 
				   sentinal values. */

			}
		}

			if (*n_obs0 != n_sent)
				*b_newdat = *b_newdat / (*n_obs0 - n_sent);
				/* divide by number of observations 
				   minus number of sentinals */
			b_newdat++;

		}

	break;


	case 2:		/* short integer data */

		i_dat = (short *) argv[1]->value.arr->data;

		i_newdat = (short *)
			   IDL_MakeTempArray(IDL_TYP_INT,2,dim, 
			   IDL_BARR_INI_ZERO, &ret_dat);


		for (i=0; i<sec; i++)
		{
		n_obs0 = n_obs;	
		brk = 0;
		n_sent = 0;
		for (pix=0; pix<num_pix; pix++)
		{		
			if (*(pixel + pix) != *(pixel + brk))
			{	
				brk = pix;
				if (*n_obs0 != n_sent)
				*i_newdat = *i_newdat / (*n_obs0 - n_sent);
				n_obs0++;
				i_dat0 = *(i_dat + i*num_pix + pix);
				i_newdat++;
				n_sent = 0;
				if (i_dat0 > sent_val)
					*i_newdat = i_dat0;
				else
					n_sent++;
			}
			else
			{
				i_dat0 = *(i_dat + i*num_pix + pix);
				if (i_dat0 > sent_val)
					*i_newdat = *i_newdat + i_dat0;
				else
					n_sent++;
			}
		}

			if (*n_obs0 != n_sent)
				*i_newdat = *i_newdat / (*n_obs0 - n_sent);
			i_newdat++;

		}

	break;

	case 3:		/* IDL_LONG integer data */

		l_dat = (IDL_LONG *) argv[1]->value.arr->data;

		l_newdat = (IDL_LONG *)
			   IDL_MakeTempArray(IDL_TYP_LONG,2,dim, 
			   IDL_BARR_INI_ZERO, &ret_dat);


		for (i=0; i<sec; i++)
		{
		n_obs0 = n_obs;	
		brk = 0;
		n_sent = 0;
		for (pix=0; pix<num_pix; pix++)
		{		
			if (*(pixel + pix) != *(pixel + brk))
			{	
				brk = pix;
				if (*n_obs0 != n_sent)
				*l_newdat = *l_newdat / (*n_obs0 - n_sent);
				n_obs0++;
				l_dat0 = *(l_dat + i*num_pix + pix);
				l_newdat++;
				n_sent = 0;
				if (l_dat0 > sent_val)
					*l_newdat = l_dat0;
				else
					n_sent++;
			}
			else
			{
				l_dat0 = *(l_dat + i*num_pix + pix);
				if (l_dat0 > sent_val)
					*l_newdat = *l_newdat + l_dat0;
				else
					n_sent++;
			}
		}

			if (*n_obs0 != n_sent)
				*l_newdat = *l_newdat / (*n_obs0 - n_sent);
			l_newdat++;

		}

	break;

	case 4:		/* float data */

		f_dat = (float *) argv[1]->value.arr->data;

		f_newdat = (float *)
			   IDL_MakeTempArray(IDL_TYP_FLOAT,2,dim, 
			   IDL_BARR_INI_ZERO, &ret_dat);


		for (i=0; i<sec; i++)
		{
		n_obs0 = n_obs;	
		brk = 0;
		n_sent = 0;
		for (pix=0; pix<num_pix; pix++)
		{		
			if (*(pixel + pix) != *(pixel + brk))
			{	
				brk = pix;
				if (*n_obs0 != n_sent)
				*f_newdat = *f_newdat / (*n_obs0 - n_sent);
				n_obs0++;
				f_dat0 = *(f_dat + i*num_pix + pix);
				f_newdat++;
				n_sent = 0;
				if (f_dat0 > sent_val)
					*f_newdat = f_dat0;
				else
					n_sent++;
			}
			else
			{
				f_dat0 = *(f_dat + i*num_pix + pix);
				if (f_dat0 > sent_val)
					*f_newdat = *f_newdat + f_dat0;
				else
					n_sent++;
			}
		}

			if (*n_obs0 != n_sent)
				*f_newdat = *f_newdat / (*n_obs0 - n_sent);
			f_newdat++;

		}


	break;


	case 5:		/* double data */

		d_dat = (double *) argv[1]->value.arr->data;

		d_newdat = (double *)
			   IDL_MakeTempArray(IDL_TYP_DOUBLE,2,dim, 
			   IDL_BARR_INI_ZERO, &ret_dat);


		for (i=0; i<sec; i++)
		{
		n_obs0 = n_obs;	
		brk = 0;
		n_sent = 0;
		for (pix=0; pix<num_pix; pix++)
		{		
			if (*(pixel + pix) != *(pixel + brk))
			{	
				brk = pix;
				if (*n_obs0 != n_sent)
				*d_newdat = *d_newdat / (*n_obs0 - n_sent);
				n_obs0++;
				d_dat0 = *(d_dat + i*num_pix + pix);
				d_newdat++;
				n_sent = 0;
				if (d_dat0 > sent_val)
					*d_newdat = d_dat0;
				else
					n_sent++;
			}
			else
			{
				d_dat0 = *(d_dat + i*num_pix + pix);
				if (d_dat0 > sent_val)
					*d_newdat = *d_newdat + d_dat0;
				else
					n_sent++;
			}
		}

			if (*n_obs0 != n_sent)
				*d_newdat = *d_newdat / (*n_obs0 - n_sent);
			d_newdat++;

		}

	break;

	}

	IDL_VarCopy(ret_pix,argv[0]);
		/* replace old pixel list with condensed list */

	IDL_VarCopy(ret_dat,argv[1]);
		/* replace old data list with averaged data */

	IDL_DELTMP(ret_n);
		/* free n_obs temp array */

return;
}
