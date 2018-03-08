/*******************************************************************************
*              
*  Module:       DBLSVD.C
*
*  Purpose:      Double precision singular value decomposition.
*
*  Author:       Joel Gales, ARC, 8/26/94, SPR 11893
*       SPR 12158  3/27/95  Update for changes in EXPORT.H. J. Newmark
*       09/95                    Newmark         Upgrade to IDL 4.0
*
*  Description:  This is a translation of the Singular Value Decomposition 
*                algorithm of Golub and Reinsch (Numerische Mathematik 14(1970)
*                pp403-470) from Algol 60 to C. The routine takes a single 
*                matrix (A) and computes two additional matrices (U and V) and 
*                a vector W such that A = UWV'  The calling sequence is
*                svd(a,m,n,u,w,v,eps,flags) where a is the original, m by n 
*                matrix; u is the upper resultant; w is the vector of singular 
*                values; v is the lower resultant; eps is a convergence test 
*                constant; and flags tells what to do, values:-
*                                    0 produce only w
*                                    1 produce w and u
*                                    2 produce w and v
*                                    3 produce u, w, and v.
*                NOTE:  M must be greater than n or an error is signaled.
*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include "export_wrap.h"
#define  BELL  putchar('\07')

void dblsvd(argc,argv)

int	argc;
IDL_VPTR	argv[];

{

short int	stat,svd0();
IDL_LONG	i,j,dim[2],dim_v[2],num_elem,flag;
char		dat_typ;
IDL_UCHAR		n_dims;


double		*in_arr_d,*w,*u,*v,eps=1.0e-10;
float		*in_arr_f;
IDL_LONG	*dims;

IDL_VARIABLE	*ret_w,*ret_u,*ret_v;

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
	if (dim[0] > dim[1])
	{
		printf("# columns must be <= # of rows!\n");
		BELL;
		goto bail;
	}	

	if (dat_typ == 4) in_arr_f = (float *)  argv[0]->value.arr->data;
	if (dat_typ == 5) in_arr_d = (double *) argv[0]->value.arr->data;
		/* get pointer to input array */

	dim_v[0] = dim[0];
	dim_v[1] = dim[0];

	w = (double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,1,dim_v,
	        	       IDL_BARR_INI_ZERO,&ret_w);
	u = (double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,2,dim,
	        	       IDL_BARR_INI_ZERO,&ret_u);
	v = (double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,2,dim_v,
	        	       IDL_BARR_INI_ZERO,&ret_v);

	flag = 3;

	stat = svd0(dat_typ,in_arr_f,in_arr_d, 
		    dim[1],dim[0],u,w,v,eps,flag);

	IDL_VarCopy(ret_w,argv[1]);
	IDL_VarCopy(ret_u,argv[2]);
	IDL_VarCopy(ret_v,argv[3]);

bail:

return;
}

#include <math.h>
#define WITHV   1
#define WITHU   2
#define TINY   1.0e-20

short int svd0(dat_typ,af,ad,m,n,u,w,v,eps,flags)
double *ad,                      /* The original matrix [m,n] */
       *u,                      /* The new upper matrix [m,n] */
       *v,                      /* The new lower matrix [m,n] */
       *w,                      /* The vector of singular values [n] */
        eps;                    /* Convergence factor */
IDL_LONG m,                          /* Number of rows */
    n,                          /* Number of columns (n <= m) */
    flags;                      /* Flags controlling what gets computed. */
float	*af;
char	dat_typ;

{
        IDL_LONG i,j,k,l,l1;         /* Mostly loop IDL_VARIABLEs */
        double tol = TINY / eps;    /* tells about machine tolerance */
        double       c,f,g,h,s,x,y,z,*temp; /* temporaries */


	temp = (double *) malloc(8*m);

/*
*       First copy a into u.
*/
	if (dat_typ == 4)  for (i = 0; i <= m * n; ++i) u[i] = af[i];
	if (dat_typ == 5)  for (i = 0; i <= m * n; ++i) u[i] = ad[i];

/*
*       Reduce the u matrix to bidiagonal form with Householder transforms.
*/
        g = (x = 0.0);
        for(i = 0; i < n; ++i) {
                temp[i] = g;
                s = 0.0;
                l = i + 1;
                for (j = i; j < m; j++) {
                        s += u[j * n + i] * u[j * n + i];
                        }
                if (s < tol) {
                        g = 0.0;
                        }
                else {
                        f = u[i * n + i];
                        g = (f < 0.0) ? sqrt(s) : -sqrt(s);
                        h = f * g - s;
                        u[i * n + i] = f - g;
                        for(j = l; j < n; ++j) {
                                s = 0.0;
                                for(k = i; k < m; ++k) {
                                        s += u[k * n + i] * u[k * n + j];
                                        }
                                f = s / h;
                                for(k = i; k < m; ++k) {
                                        u[k * n + j] += f * u[k * n + i];
                                        }
                                }
                        }
                w[i] = g;
                s = 0.0;
                for (j = l; j < n; ++j) {
                        s += u[i * n + j] * u[i * n + j];
                        }
                if (s < tol) {
                        g = 0.0;
                        }
                else {
                        f = u[i * n + i + 1];
                        g = (f < 0.0) ? sqrt(s) : -sqrt(s);
                        h = f * g - s;
                        u[i * n + i + 1] = f - g;
                        for (j = l; j < n; ++j) {
                                temp[j] = u[i * n + j] / h;
                                }
                        for (j = l; j < m; ++j) {
                                s = 0.0;
                                for (k = l; k < n; ++k) {
                                        s += u[j * n + k] * u[i*n+k];
                                        }
                                for (k = l; k < n; ++k) {
                                        u[j*n+k] += s * temp[k];
                                        }
                                }
                        }
                y = fabs(w[i]) + fabs(temp[i]);
                if (y > x) {
                        x = y;
                        }
                }
/*
*       Now accumulate right-hand transforms if we are building v too.
*/
        if (flags & WITHV) for (i = n - 1; i >= 0; --i) {
                if (g != 0.0) {
                        h = u[i * n + i + 1] * g;
                        for (j = l; j < n; ++j) {
                                v[j * n + i] = u[i * n + j] / h;
                                }
                        for (j = l; j < n; ++j) {
                                s = 0.0;
                                for (k = l; k < n; ++k) {
                                        s += u[i * n + k] * v[k * n + j];
                                        }
                                for (k = l; k < n; ++k) {
                                        v[k * n + j] += s * v[k * n + i];
                                        }
                                }
                        }
                for (j = l; j < n; ++j) {
                        v[i * n + j] = (v[j * n + i] = 0.0);
                        }
                v[i * n + i] = 1.0;
                g = temp[i];
                l = i;
                }
/*
*       Now accumulate the left-hand transforms.
*/
        if (flags & WITHU) for (i = n - 1; i >= 0; --i) {
                l = i + 1;
                g = w[i];
                for (j = l; j < n; ++j) {
                        u[i * n + j] = 0.0;
                        }
                if (g != 0.0) {
                        h = u[i * n + i] * g;
                        for (j = l; j < n; ++j) {
                                s = 0.0;
                                for (k = l; k < m; ++k) {
                                        s += u[k * n + i] * u[k * n + j];
                                        }
                                f = s / h;
                                for (k = i; k < m; ++k) {
                                        u[k * n + j] += f * u[k * n + i];
                                        }
                                }
                        for (j = i; j < m; ++j) {
                                u[j * n + i] /= g;
                                }
                        }
                else {
                        for (j = i; j < m; ++j) {
                                u[j * n + i] = 0.0;
                                }
                        }
                u[i * n + i] += 1.0;
                }
/*
*       Now diagonalise the bidiagonal form. BEWARE GOTO's IN THE LOOP!!
*/
        eps = eps * x;
        for (k = n - 1; k >= 0; --k) {
testsplitting:
                for (l = k; l >= 0; --l) {
                        if (fabs(temp[l]) <= eps) 
                                goto testconvergence;
                        if (fabs(w[l - 1]) <= eps)
                                goto cancellation;
                        }
/*
*               Cancellation of temp[l] if l > 0;
*/
cancellation:
                c = 0.0;
                s = 1.0;
                l1 = l - 1;
                for (i = l; i <= k; ++i) {
                        f = s * temp[i];
                        temp[i] *= c;
                        if (fabs(f) <= eps) 
                                goto testconvergence;
                        g = w[i];
                        h = (w[i] = sqrt(f * f + g * g));
                        c = g/h;
                        s = -f/h;
                        if (flags & WITHU) for (j = 0; j < m; ++j) {
                                y = u[j * n + l1];
                                z = u[j * n + i];
                                u[j * n + l1] = y * c + z * s;
                                u[j * n + i] = -y * s + z * c;
                                }
                        }
testconvergence:
                z = w[k];
                if (l == k) goto convergence;
/*
*               Shift from bottom 2x2 minor.
*/
                x = w[l];
                y = w[k - 1];
                g = temp[k - 1];
                h = temp[k];
                f = ((y - z)*(y + z) + (g - h)*(g + h)) / (2 * h * y);
                g = sqrt(f * f + 1);
                f = ((x - z)*(x + z) + h*(y/((f < 0.0)?f-g:f+g) - h)) / x;
/*
*               Next QR transformation.
*/
                c = (s = 1);
                for (i = l + 1; i <= k; ++i) {
                        g = temp[i];
                        y = w[i];
                        h = s * g;
                        g *= c;
                        temp[i - 1] = (z = sqrt(f * f + h * h));
                        c = f / z;
                        s = h/z;
                        f = x * c + g * s;
                        g = -x * s + g * c;
                        h = y * s;
                        y *= c;
                        if (flags & WITHV) for (j = 0; j < n; ++j) {
                                x = v[j * n + i - 1];
                                z = v[j * n + i];
                                v[j * n + i - 1] = x * c + z * s;
                                v[j * n + i] = -x * s + z * c;
                                }
                        w[i - 1] = (z = sqrt(f * f + h * h));
                        c = f / z;
                        s = h / z;
                        f = c * g + s * y;
                        x = -s * g + c * y;
                        if (flags & WITHU) for (j = 0 ; j < m; ++j) {
                                y = u[j * n + i - 1];
                                z = u[j * n + i];
                                u[j * n + i - 1] = y * c + z * s;
                                u[j * n + i] = -y * s + z * c;
                                }
                        }
                temp[l] = 0.0;
                temp[k] = f;
                w[k] = x;
                goto testsplitting;

convergence:
                if (z < 0.0) {
/*
*                       w[k] is made non-negative.
*/
                        w[k] = -z;
                        if (flags & WITHV) for (j = 0; j < n; ++j) {
                                v[j * n + k] = -v[j * n + k];
                                }
                        }
                }
free(temp);

return(1);
        }

