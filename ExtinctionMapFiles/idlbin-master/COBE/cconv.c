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

   FUNCTION:	Pixel <--> Unit Vector <--> Lon/Lat conversion routine
		Called by IDL COORCONV routine.  Note: There is a IDL
		function of the same name that is a shell around other
		IDL functions which perform conversions between pixels
		and unit vectors and lon/lat.  If this C module is not
		implememented then the IDL function is compiled and run
		instead.  If implemented, then this module takes prece-
		dence.

   AUTHOR:	J.M. Gales (Applied Research Corp)
   DATE:	9/92, SPR 10025, Rewrite in "C" for faster BUILDMAP execution.

	MODIFICATION HISTORY:
	DATE	SPR/SER		PROGRAMMER	DESCRIPTION
	----	-------		----------	-----------
	12/9/92	  10324		Gales		Add ability to read and
						write double precision
						arrays for uv and ll 
						input.
      02/12/93  11468           Newmark         Change to compile under
                                                VAX and Sun C not GNU C
	03/94	11614		Dan Grogan	Changed "export.h" to 
						"export_wrap.h" for IDL 3.0+
       09/95    12279           Newmark         Upgrade to IDL 4.0

       03/97                    Gales           upgrade for use at res 15.


   COMMAND LINE ARGUMENT: input array, resolution/action code

   INPUT PARAMETERS:

	Name		Type		Description
	----		----		-----------
	input array	LONG arr	Input array
			(pixel)
			float/double arr
			(unit vector)
			float/double arr
			(lon/lat vector)

	resolution	int		skymap resolution/action code
					(For uv<-->ll conversion, set
					 resolution to -1)

   RETURN:		output array (pixel, unit vector, lon/lat)

   FUNCTIONS CALLED:	pix2uv, uv2pix, uv2ll, ll2uv,
				uv2pix_d, uv2ll_d, ll2uv_d

*******************************************************************************/

#include <stdio.h>
#include <math.h>
#include "export_wrap.h"

IDL_VPTR cconv(argc,argv)

int	argc;
IDL_VPTR	argv[];

{

short		i,j,res,in_type;
IDL_LONG	in_dim[2], out_dim[2];
IDL_UCHAR		n_dims,type;

/* IDL pointer IDL_VARIABLEs */

IDL_LONG	*pix,*dims;
float		*uv,*ll;
double		*uv_d,*ll_d;

IDL_VARIABLE	*ret;

/*	printf("New cconv\n");*/

	n_dims = (IDL_UCHAR) argv[0]->value.arr->n_dim;
			/* get number of dimensions of input */

	dims = (IDL_LONG *) argv[0]->value.arr->dim;
	for (i=0; i<n_dims; i++) in_dim[i] = *dims++;
			/* load dimensions */

	type = (IDL_UCHAR) argv[0]->type;
			/* get input VARIABLE data type */

	if (n_dims == 1)
		in_type = 0;
	else
		in_type = (in_dim[1] == 2) ? 1 : 2;
			/* in_type = 0 for pixel input
			   in_type = 1 for lon/lat input
			   in_type = 2 for uv input */


	res = (short) argv[1]->value.i;
			/* get resolution/action code */


	if (in_type == 0)	/* pixel input */
	{
		pix = (IDL_LONG *) argv[0]->value.arr->data;
				/* get pointer to pixel numbers */

		out_dim[0] = in_dim[0];
		out_dim[1] = 3;
				/* uv output n by 3 array */

		uv = (float *)
		       IDL_MakeTempArray(IDL_TYP_FLOAT,2,out_dim, 
		       IDL_BARR_INI_ZERO, &ret);
				/* allocate space for uv output */

		pix2uv(pix,uv,res,in_dim);
				/* pixel to unit vector conversion */
	}



	if (in_type == 1)	/* longitude/latitude input */
	{
		out_dim[0] = in_dim[0];
		out_dim[1] = 3;

		if (type == IDL_TYP_FLOAT) 
		{
			ll = (float *) argv[0]->value.arr->data;
				/* get float pointer to lon/lat input */

			uv = (float *)
		       		IDL_MakeTempArray(IDL_TYP_FLOAT,2,out_dim, 
		       		IDL_BARR_INI_ZERO, &ret);
				/* allocate space for uv output */

			ll2uv(ll,uv,in_dim);
				/* lon/lat to unit vector conversion */
		}

		if (type == IDL_TYP_DOUBLE) 
		{
			ll_d = (double *) argv[0]->value.arr->data;
				/* get double pointer to lon/lat input */

			uv_d = (double *)
		       		IDL_MakeTempArray(IDL_TYP_DOUBLE,2,out_dim, 
		       		IDL_BARR_INI_ZERO, &ret);
				/* allocate space for uv output */

			ll2uv_d(ll_d,uv_d,in_dim);
				/* lon/lat to unit vector conversion */
		}

	}



	if (in_type == 2)	/* unit vector input */
	{

		if (type == IDL_TYP_FLOAT)
			uv = (float *) argv[0]->value.arr->data;
				/* get float pointer to unit vectors */

		if (type == IDL_TYP_DOUBLE)
			uv_d = (double *) argv[0]->value.arr->data;
				/* get double pointer to unit vectors */

		out_dim[0] = in_dim[0];

		if (res == -1)	/* lon/lat output */
		{
			out_dim[1] = 2;

			if (type == IDL_TYP_FLOAT)
			{
				ll = (float *)
			       		IDL_MakeTempArray(IDL_TYP_FLOAT,2,out_dim, 
			       		IDL_BARR_INI_ZERO, &ret);
					/* allocate space for ll output */

				uv2ll(uv,ll,in_dim);
				/* unit vector to lon/lat conversion */
			}

			if (type == IDL_TYP_DOUBLE)
			{
				ll_d = (double *)
			       		IDL_MakeTempArray(IDL_TYP_DOUBLE,2,out_dim, 
			       		IDL_BARR_INI_ZERO, &ret);
					/* allocate space for ll output */

				uv2ll_d(uv_d,ll_d,in_dim);
				/* unit vector to lon/lat conversion */
			}


		}


		if (res > 0)	/* pixel output */
		{

			pix = (IDL_LONG *)
			        IDL_MakeTempArray(IDL_TYP_LONG,1,out_dim, 
			        IDL_BARR_INI_ZERO, &ret);
					/* allocate space for pixel output */

			if (type == IDL_TYP_FLOAT)
			uv2pix(uv,pix,res,in_dim);
				/* unit vector to pixel conversion */

			if (type == IDL_TYP_DOUBLE)
			uv2pix_d(uv_d,pix,res,in_dim);
				/* double unit vector to pixel conversion */
		}

	}


return(ret);		/* return output */
}


pix2uv(pixel,uv,res,dim)	/* pix to uv */

IDL_LONG	*pixel;
float		*uv;
int		res;
int		dim[];

{

int		i, j, k, bit;
IDL_LONG	face,fpix,pix;
float		x,y,xx,yy,eta,xi,half_side,c[3],norm;
double          xtemp,ytemp,xid,etad;
IDL_LONG	num_pix_face,num_pix_side,one;
unsigned int	pow_2[16];
static double	p[29] ;

one=1;
p[0]= 0.0; 
p[1]= -0.27292696;
p[2]= -0.07629969;
p[3]= -0.02819452;
p[4]= -0.22797056;
p[5]= -0.01471565; 
p[6]=0.27058160;
p[7]=0.54852384;
p[8]=0.48051509; 
p[9]= -0.56800938;
p[10]= -0.60441560; 
p[11]= -0.62930065; 
p[12]= -1.74114454;
p[13]=0.30803317; 
p[14]=1.50880086; 
p[15]=0.93412077;
p[16]=0.25795794; 
p[17]=1.71547508; 
p[18]=0.98938102;
p[19]= -0.93678576; 
p[20]= -1.41601920; 
p[21]= -0.63915306;
p[22]=0.02584375; 
p[23]= -0.53022337; 
p[24]= -0.83180469;
p[25]=0.08693841; 
p[26]=0.33887446; 
p[27]=0.52032238; 
p[28]=0.14381585;

	num_pix_side = (IDL_LONG) (pow(2.,res-1.));
	num_pix_face = (IDL_LONG) (pow(4.,res-1.));
		/* determine # of pixels per side & face */

	for (i=0; i<=15; i++) pow_2[i] = one << i;
		/* generate bit mask array */


	half_side = (float) (pow(2.,res-2.));
			/* determine # of pixels per face
			   # pixels per half side */

	for (pix=0; pix<dim[0]; pix++)	/* pixel loop */
	{

	face = *(pixel + pix) / num_pix_face;
	fpix = *(pixel + pix) - num_pix_face*face;
			/* determine face and pixel # within face */

	i = 0;
	j = 0;

	for (bit=0; bit<=res-1; bit++)
	{
		i = i | (int) (pow_2[bit] * (one & fpix));
		fpix = fpix >> 1;
		j = j | (int) (pow_2[bit] * (one & fpix));
		fpix = fpix >> 1;
 	}
		/* decode pixel into column and row */


	x = (i - half_side +0.5) / half_side;
	y = (j - half_side +0.5) / half_side;
			/* get normalized x & y positions */

/*	printf("pixel,face,fpix: %d %d %d\n",*(pixel + pix),face,fpix);
	printf("x, y: %12.7f %12.7f\n",x,y);
*/

	xx = x*x;
	yy = y*y;

	xi=x*(1.+(1.-xx)*( 
	p[1]+xx*(p[2]+xx*(p[4]+xx*(p[7]+xx*(p[11]+xx*(p[16]+xx*p[22]))))) +
	yy*( p[3]+xx*(p[5]+xx*(p[8]+xx*(p[12]+xx*(p[17]+xx*p[23])))) +
	yy*( p[6]+xx*(p[9]+xx*(p[13]+xx*(p[18]+xx*p[24]))) +
	yy*( p[10]+xx*(p[14]+xx*(p[19]+xx*p[25])) +
	yy*( p[15]+xx*(p[20]+xx*p[26]) +
	yy*( p[21]+xx*p[27] + yy*p[28])))) )));

	eta=y*(1.+(1.-yy)*(
	p[1]+yy*(p[2]+yy*(p[4]+yy*(p[7]+yy*(p[11]+yy*(p[16]+yy*p[22]))))) +
	xx*( p[3]+yy*(p[5]+yy*(p[8]+yy*(p[12]+yy*(p[17]+yy*p[23])))) +
	xx*( p[6]+yy*(p[9]+yy*(p[13]+yy*(p[18]+yy*p[24]))) +
	xx*( p[10]+yy*(p[14]+yy*(p[19]+yy*p[25])) +
	xx*( p[15]+yy*(p[20]+yy*p[26]) +
	xx*( p[21]+yy*p[27] + xx*p[28])))) )));

			/* magic quad cube formuli */

/*	printf("xi, eta: %12.7f %12.7f\n",xi,eta);*/

	/* If res 15 then iterate */
	if (res == 15)
	{
		xid = (double) xi;
		etad = (double) eta;

 		incube(xid,etad,&xtemp,&ytemp);
/*	printf("x,y: %12.7f %12.7f\n",xtemp,ytemp);*/


		xid -= (xtemp - x);
		etad -= (ytemp - y);
/*	printf("xi, eta (2): %12.7f %12.7f\n",xid,etad);*/

 		incube(xid,etad,&xtemp,&ytemp);
/*	printf("x,y: %12.7f %12.7f\n",xtemp,ytemp);*/

		xi = xid - (xtemp - x);
		eta = etad - (ytemp - y);
/*	printf("xi, eta (3): %12.7f %12.7f\n",xi,eta);*/
	}


	/* This section taken from XYAXIS.FOR */
	switch(face)	/* determine unnormalized vector */
	{

	case 0 :
	  c[2]=1.;
	  c[0]= -eta;
	  c[1]=xi;
	break;

	case 1 :
	  c[0]=1.;
	  c[2]=eta;
	  c[1]=xi;
	break;

	case 2 :
	  c[1]=1.;
	  c[2]=eta;
	  c[0]= -xi;
	break;

	case 3 :
	  c[0]= -1.;
	  c[2]=eta;
	  c[1]= -xi;
	break;

	case 4 :
	  c[1]= -1.;
	  c[2]=eta;
	  c[0]=xi;
	break;

	case 5 :
	  c[2]= -1.;
	  c[0]=eta;
	  c[1]=xi;
	break;

	}

	norm = pow(c[0]*c[0]+c[1]*c[1]+c[2]*c[2],0.5);
	c[0] = c[0] / norm;
	c[1] = c[1] / norm;
	c[2] = c[2] / norm;
			/* normalize vector */

	*(uv + pix) = c[0];
	*(uv + pix + dim[0]) = c[1];
	*(uv + pix + 2*dim[0]) = c[2];
			/* load into unit vector output array */


	} /* end pixel loop */

return;
}


uv2pix(uv,pixel,res,dim)	/* uv to pix */

float		*uv;
IDL_LONG	*pixel;
short		res;
int		dim[];

{

int		i;
IDL_LONG	pix,pix_num,u2p();
double		c0,c1,c2;
IDL_LONG	num_pix_face,num_pix_side,one;
unsigned int	pow_2[16];

one=1;

	num_pix_side = (IDL_LONG) (pow(2.,res-1.));
	num_pix_face = (IDL_LONG) (pow(4.,res-1.));
		/* determine # of pixels per side & face */

	for (i=0; i<=15; i++) pow_2[i] = one << i;
		/* generate bit mask array */


	for (pix=0; pix<dim[0]; pix++)
	{

	c0 = (double) *(uv + pix);
	c1 = (double) *(uv + pix + dim[0]);
	c2 = (double) *(uv + pix + 2*dim[0]);
			/* read unit vector coordinates */

	pix_num = u2p(c0,c1,c2,res,num_pix_face,num_pix_side,pow_2);

	*pixel++ = pix_num;
			/* load into pixel output array */

	}	/* pix loop */

return;
}


uv2pix_d(uv,pixel,res,dim)	/* uv to pix (double) */

double		*uv;
IDL_LONG	*pixel;
short		res;
int		dim[];

{

int		i;
IDL_LONG	pix,pix_num,u2p();
double		c0,c1,c2;
IDL_LONG	num_pix_face,num_pix_side,one;
unsigned int	pow_2[16];

one=1;

	num_pix_side = (IDL_LONG) (pow(2.,res-1.));
	num_pix_face = (IDL_LONG) (pow(4.,res-1.));
		/* determine # of pixels per side & face */

	for (i=0; i<=15; i++) pow_2[i] = one << i;
		/* generate bit mask array */


	for (pix=0; pix<dim[0]; pix++)
	{

	c0 = *(uv + pix);
	c1 = *(uv + pix + dim[0]);
	c2 = *(uv + pix + 2*dim[0]);
			/* read unit vector coordinates */

	pix_num = u2p(c0,c1,c2,res,num_pix_face,num_pix_side,pow_2);

	*pixel++ = pix_num;
			/* load into pixel output array */

	}	/* pix loop */

return;
}


IDL_LONG u2p(c0,c1,c2,res,num_pix_face,num_pix_side,pow_2)

double		c0,c1,c2;
short		res;
IDL_LONG	num_pix_face,num_pix_side;
unsigned int	pow_2[];

{

int		i, j, k, bit;
IDL_LONG	face,fpix,pix_num;
double		x,y,eta,xi,half_side,norm;
double		abs_yx,abs_zx,abs_zy;

 	/* Determine Face Number */

	abs_yx = (double) fabs(c1/c0);
	abs_zx = (double) fabs(c2/c0);
	abs_zy = (double) fabs(c2/c1);

	if (abs_zx >= 1 && abs_zy >= 1 && c2 >= 0)
	{
		face = 0;
		goto lbl;
	}

	if (abs_zx >= 1 && abs_zy >= 1 && c2 < 0)
	{
		face = 5;
		goto lbl;
	}

	if (abs_zx < 1 && abs_yx < 1 && c0 >= 0)
	{	
		face = 1;
		goto lbl;
	}

	if (abs_zx < 1 && abs_yx < 1 && c0 < 0)
	{
		face = 3;
		goto lbl;
	}

	if (abs_zy < 1 && abs_yx >= 1 && c1 >= 0)
	{
		face = 2;
		goto lbl;
	}

	if (abs_zy < 1 && abs_yx >= 1 && c1 < 0)
		face = 4;



	/* Determine quad cube normalized coordinates */

lbl:	switch(face)
	{

	case 0 :
	  eta = -c0/c2;
	  xi  =  c1/c2;
	break;

	case 1 :
	  eta = c2/c0;
	  xi  = c1/c0;
	break;

	case 2 :
	  eta =  c2/c1;
	  xi  = -c0/c1;
	break;

	case 3 :
	  eta = -c2/c0;
	  xi  =  c1/c0;
	break;

	case 4 :
	  eta = -c2/c1;
	  xi  = -c0/c1;
	break;

	case 5 :
	  eta = -c0/c2;
	  xi  = -c1/c2;
	break;

	}

	incube(xi,eta,&x,&y);

	x=(x+1.)/2;
	y=(y+1.)/2;

	i = (int) (x * num_pix_side);
	j = (int) (y * num_pix_side);
	if (i == num_pix_side) i = i - 1;
	if (j == num_pix_side) j = j - 1;
			/* get column & row numbers
			   knock back one if on boundary */


	fpix = 0;

	for (bit=0; bit<res-1; bit++)
	{
		fpix = fpix | ((pow_2[bit] & i) << bit);
		fpix = fpix | ((pow_2[bit] & j) << (bit+1));
	}
			/* encode column & row into pixel # */

	pix_num = face*num_pix_face + fpix;

return(pix_num);
}





incube(xi,eta,x,y)

double		eta,xi,*x,*y;

{
double		e2,x2,e4,x4,onme2,onmx2,half_side,norm;
static double	gstar,g,m,
		w1,c00,c10,
		c01,c11,c20,
		c02,d0,d1,
		r0;
gstar=1.37484847732;
g= -0.13161671474;
m=0.004869491981;
w1= -0.159596235474;
c00=0.141189631152;
c10=0.0809701286525;
c01= -0.281528535557;
c11=0.15384112876;
c20= -0.178251207466;
c02=0.106959469314;
d0=0.0759196200467;
d1= -0.0217762490699;
r0=0.577350269;


	e2 = eta * eta;
	x2 = xi * xi;
	e4 = e2 * e2;
	x4 = x2 * x2;
	onmx2=1.-x2;
	onme2=1.-e2;


	*x=xi*(gstar+x2*(1.-gstar)+onmx2*(e2*(g+(m-g)*x2
	  +onme2*(c00+c10*x2+c01*e2+c11*x2*e2+c20*x4+c02*e4))
	  +x2*(w1-onmx2*(d0+d1*x2))));


	*y=eta*(gstar+e2*(1.-gstar)+onme2*(x2*(g+(m-g)*e2
	  +onmx2*(c00+c10*e2+c01*x2+c11*e2*x2+c20*e4+c02*x4))
	  +e2*(w1-onme2*(d0+d1*e2))));
		/* Magic formuli to determine x & y face coordinates */

	return;
}


ll2uv(ll,uv,dim)	/* ll to uv */

float		*ll,*uv;
int		dim[];

{

IDL_LONG	cnt;
float		lon,lat,cos_lat,sin_lat,cos_lon,sin_lon;
double		d2r;

	d2r = atan(1.0) / 45;
			/* degrees to radians conv */

	for (cnt=0; cnt<dim[0]; cnt++)
	{
	
		lon = *(ll + cnt) * d2r;
		lat = *(ll + cnt + dim[0]) * d2r;
			/* get lon/lat input */

		cos_lat = cos(lat);
		sin_lat = sin(lat);
		cos_lon = cos(lon);
		sin_lon = sin(lon);

		*(uv + cnt) = cos_lat * cos_lon;
		*(uv + cnt + dim[0]) = cos_lat * sin_lon;
		*(uv + cnt + 2*dim[0]) = sin_lat;
			/* load into unit vector output array */
	}

return;
}


ll2uv_d(ll,uv,dim)	/* ll to uv (double) */

double		*ll,*uv;
int		dim[];

{

IDL_LONG	cnt;
double		lon,lat,cos_lat,sin_lat,cos_lon,sin_lon;
double		d2r;

	d2r = atan(1.0) / 45;
			/* degrees to radians conv */

	for (cnt=0; cnt<dim[0]; cnt++)
	{
	
		lon = *(ll + cnt) * d2r;
		lat = *(ll + cnt + dim[0]) * d2r;
			/* get lon/lat input */

		cos_lat = cos(lat);
		sin_lat = sin(lat);
		cos_lon = cos(lon);
		sin_lon = sin(lon);

		*(uv + cnt) = cos_lat * cos_lon;
		*(uv + cnt + dim[0]) = cos_lat * sin_lon;
		*(uv + cnt + 2*dim[0]) = sin_lat;
			/* load into unit vector output array */
	}

return;
}


uv2ll(uv,ll,dim)

float		*ll,*uv;
int		dim[];

{

IDL_LONG	cnt;
float		c0,c1,c2,lon,lat;
double		r2d;

	r2d = 45 / atan(1.0);
			/* radians to degrees conv */

	for (cnt=0; cnt<dim[0]; cnt++)
	{
		c0 = *(uv + cnt);
		c1 = *(uv + cnt + dim[0]);
		c2 = *(uv + cnt + 2*dim[0]);
			/* read unit vector coordinates */
	
		lon = atan2(c1,c0) * r2d;
		lat = asin(c2) * r2d;

		if (lon < 0) lon = lon + 360;
			/* make lon from 0 to 360 */

		*(ll + cnt) = lon;
		*(ll + cnt + dim[0]) = lat;
			/* load into lon/lat output array */

	}

return;
}


uv2ll_d(uv,ll,dim)

double		*ll,*uv;
int		dim[];

{

IDL_LONG	cnt;
double		c0,c1,c2,lon,lat;
double		r2d;

	r2d = 45 / atan(1.0);
			/* radians to degrees conv */

	for (cnt=0; cnt<dim[0]; cnt++)
	{
		c0 = *(uv + cnt);
		c1 = *(uv + cnt + dim[0]);
		c2 = *(uv + cnt + 2*dim[0]);
			/* read unit vector coordinates */
	
		lon = atan2(c1,c0) * r2d;
		lat = asin(c2) * r2d;

		if (lon < 0) lon = lon + 360;
			/* make lon from 0 to 360 */

		*(ll + cnt) = lon;
		*(ll + cnt + dim[0]) = lat;
			/* load into lon/lat output array */

	}
return;
}
