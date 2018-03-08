/*******************************************************************************
*              
*  Module:       QRDCMPD.C
*
*  Purpose:      QR matrix decomposition.
*
*  Author:       Joel Gales, ARC, 8/26/94, SPR 11893
*
*  Modification History
*  --------------------
*  Trap zero norm input vectors, J.M. Gales, ARC, 10/12/94  SPR# 11932
*        3/27/95 12158           Newmark         Update for v3.5 EXPORT.h
*       09/95                    Newmark         Upgrade to IDL 4.0
*
*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include "export_wrap.h"
#define  BELL  putchar('\07')

void qrdcmp(argc,argv)

int	argc;
IDL_VPTR	argv[];

{

short int	stat,qr0();
IDL_LONG	i,j,k,dim[2],dim_r[2];
char		dat_typ;
IDL_UCHAR		n_dims;


double		*in_arr_d,*q,*r;
float		*in_arr_f;
IDL_LONG	*dims;


IDL_VARIABLE	*ret_q,*ret_r;

	dat_typ = (char) argv[0]->type;
			/* get data type */
	
	if (dat_typ < 4 || dat_typ > 5)
	{
		printf("Input array must be FLOAT or DOUBLE!\n");
		BELL;
		goto bail;
	}	

	dims = (IDL_LONG *) argv[0]->value.arr->dim;
	for (i=0; i<=1; ++i) dim[i] = *dims++;
	if (dim[0] < dim[1])
	{
		printf("# columns must be >= # of rows!\n");
		BELL;
		goto bail;
	}	



	if (dat_typ == 4) in_arr_f = (float *)  argv[0]->value.arr->data;
	if (dat_typ == 5) in_arr_d = (double *) argv[0]->value.arr->data;
		/* get pointer to input array */



	dim_r[0] = dim[1];
	dim_r[1] = dim[1];

	q = (double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,2,dim,
	        	       IDL_BARR_INI_ZERO,&ret_q);
	r = (double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,2,dim_r,
	        	       IDL_BARR_INI_ZERO,&ret_r);

	stat = qr0(dat_typ,in_arr_f,in_arr_d,dim[0],dim[1],q,r);

	IDL_VarCopy(ret_q,argv[1]);
	IDL_VarCopy(ret_r,argv[2]);

bail:

return;
}



short int qr0(dat_typ,af,ad,m,n,mat,r)

double		*ad,*mat,*r;
IDL_LONG	m,n;
float		*af;
char		dat_typ;

{

double		*fac,*hse_vec,*v,*p2,accm,v0;
IDL_LONG	i,j,l;

/*
*       First copy a into mat.
*/
	if (dat_typ == 4)  for (i = 0; i <= m * n; ++i) mat[i] = af[i];
	if (dat_typ == 5)  for (i = 0; i <= m * n; ++i) mat[i] = ad[i];


/* Generate matrices to store right triangular and Householder vectors
   ------------------------------------------------------------------- */
hse_vec = (double *) calloc(m*n,8);
fac = (double *) calloc(n,8);
v   = (double *) calloc(m,8);
p2  = (double *) calloc(m,8);


/* Generate Right Triangular
   ------------------------- */
for (i=0; i<n; ++i)
{

	/* Determine Householder vector
	   ---------------------------- */
	for (accm=0, j=0; j<(m-i); ++j)
	{
		v[j] = mat[i*m + i+j];
		accm += v[j] * v[j];
	}
	accm = pow(accm,0.5);

	if (accm == 0)
	{
		for (j=0; j<(m-i); ++j) p2[j] = 0;
		printf("Row %d all zeros.\n",i);
		BELL;
		goto skip;
	}


	if (v[0] < 0) v0 = v[0] - accm; else v0 = v[0] + accm;

	accm = 1;
	v[0] = 1;
	for (j=1; j<(m-i); ++j)
	{
		v[j] = v[j] / v0;
		accm += v[j] * v[j];
	}
		/* normalize to v(0) = 1 */


	fac[i] = 2 / accm;
	for (j=0; j<(m-i); ++j) p2[j] = fac[i] * v[j];


skip:	for (l=0; l<n; ++l)
	{
		for (accm=0, j=0; j<(m-i); ++j) accm += v[j] * mat[l*m + i+j];
		for (j=0; j<(m-i); ++j) mat[l*m + i+j] -= accm * p2[j];
	}
		/* Apply ith Householder matrix to input matrix */

	for (j=0; j<(m-i); ++j)	hse_vec[i*m + i+j] = v[j];
		/* Store Householder vector for later use */

}	/* Column Loop */


for (i=0; i<n; ++i)  for (j=0; j<n; ++j) r[i*n + j] = mat[i*m + j];
	/* Extract Right Triangular matrix */




/* Generate Orthogonal matrix from product of Householder matrices
   --------------------------------------------------------------- */
for (i=0; i<n*m; ++i)  mat[i] = 0;
for (i=0; i<n; ++i)    mat[i*(m+1)] = 1;


for (l=n-1; l>=0; l--)
{
	for (j=0; j<n; ++j)
	{
		for (accm=0 ,i=0; i<m; ++i) accm += hse_vec[l*m+i] * mat[j*m+i];
		for (i=0; i<m; ++i) mat[j*m+i] -= fac[l] * accm * 
							hse_vec[l*m+i];
	}
}


free(hse_vec);
free(fac);
free(v);
free(p2);


return(1);

}
