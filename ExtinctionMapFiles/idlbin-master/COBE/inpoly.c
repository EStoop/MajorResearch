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
/* MODULE:	inpoly.c

   FUNCTION:	This program determines whether a given point on a 
		sphere is within a spherical polygon defined by other
		points on the sphere connected by great circle segments.
		It is called by the PXINPOLY IDL procedure.  An IDL 
		version of this function also exists.

   AUTHOR:	J.M. Gales (Applied Research Corp)
   DATE:	6/93

	MODIFICATION HISTORY:
	DATE	SPR/SER		PROGRAMMER	DESCRIPTION
	----	-------		----------	-----------
    07/23/93      11174           Gales         Initial Delivery
	03/94	11614		Dan Grogan	Changed "export.h" to 
						"export_wrap.h" for IDL 3.0+
       09/95                    Newmark         Upgrade to IDL 4.0

   COMMAND LINE ARGUMENT: vertices, edge, plane, index, uv

   INPUT PARAMETERS:
	Name		Type		Description
	----		----		-----------
	vertices	flt arr		polygon vertices in unit vector
			(n_vert,3)		format

	edge		flt arr		edge vectors of triangular
			(3,3*n_tri)	facets of polygon

	plane		flt arr		parameters of plane defined
			(4,n_tri)	by triangular facets
					Ax+By+Cz=D

	index		int arr		indicies of each triangular
			(3,n_tri)	facet in terms of polygon
					vertices indicies

	uv		flt arr		unit vectors of test points
			(3,n_pnts)	on sphere in unit vector format


   RETURN:		flag  [byte array (n_pnts)]
			0 is outside polygon, 1 if inside

   FUNCTIONS CALLED:	None

*******************************************************************************/

/*			Mathematical Interlude
			----------------------


A plane is defined by the equation:

		Ax + By + Cz = D.		(1)

The point given by the intersection of a vector colinear with
a = (u,v,w) and the plane is given by b = (x,y,z) where a x b = 0.
Therefore:

			vz - wy = 0
			wx - uz = 0
			vx - uy = 0	.

If w NE 0 then:

			x = (u/w)*z
			y = (v/w)*z	.

Substituting into (1) and solving for z we get:

			z = D / (A*(u/w) + B*(v/w) + C) .


If w = 0 but v NE 0 then:

			x = (u/v)*y
			z = 0.

Substituting into (1) and solving for y we get:

			y = D / (A*(u/v) + B) .


If w = 0 and v = 0 but u NE 0 then:

			y = 0
			z = 0.

Substituting this into (1) and solving for x we get:

			x = D / A	.


								*/

#include <stdio.h>
#include <math.h>
#include "export_wrap.h"

IDL_VPTR inpoly(argc,argv)

int	argc;
IDL_VPTR	argv[];


{

short		i0, i1, i2;
IDL_LONG	i, j, k, n, n_tri, n_vert, tri, dim[1];
float		uv_pl[3], diff[3][3], a0, a1, a2, y, z, test;
float		u, v, w, A, B, C, D;
float		*diff_adr, *edge_adr;
float		*a, *b;

/* IDL pointer IDL_VARIABLEs */

short		*index;
float		*vert, *edge, *plane, *uv;
IDL_UCHAR		*flag;

IDL_LONG	*dims;

IDL_VARIABLE	*retdat;


	/* set up pointers to input */

	vert = (float *) argv[0]->value.arr->data;
			/* get pointer to vertices */

	edge = (float *) argv[1]->value.arr->data;
			/* get pointer to triangle edge vectors */

	plane = (float *) argv[2]->value.arr->data;
			/* get pointer to triangle plane parameters */

	index = (short *) argv[3]->value.arr->data;
			/* get pointer to triangle vertici index */

	uv = (float *) argv[4]->value.arr->data;
			/* get pointer to test points (in uv format) */


	n_vert = (IDL_LONG) argv[0]->value.arr->n_elts / 3;
			/* get number of vertices */

	n_tri = (IDL_LONG) argv[3]->value.arr->n_elts / 3;
			/* get number of triangles */

	n = (IDL_LONG) argv[4]->value.arr->n_elts / 3;
			/* # of test points */


	dim[0] = n;
	flag = (IDL_UCHAR *)
		IDL_MakeTempArray(IDL_TYP_BYTE,1,dim,IDL_BARR_INI_ZERO, &retdat);
			/* allocate space for output flag array */


	for (i=0; i<n; i++)	/* for each test point */
	{
		u = uv[i];
		v = uv[n+i];
		w = uv[2*n+i];
				/* get uv coordinates */


		if (w != 0)		/* w NE 0 */
		{
			a0 = u/w;
			a1 = v/w;
		}
		else if (v != 0)	/* v NE 0 */
		{
			a0 = u/v;
		}


		for (tri=0; tri<n_tri; tri++)
		{
			i0 = index[3*tri];
			i1 = index[3*tri+1];
			i2 = index[3*tri+2];

			A = plane[4*tri];
			B = plane[4*tri+1];
			C = plane[4*tri+2];
			D = plane[4*tri+3];

			if (w != 0)
			{
				z = D / (A*a0 + B*a1 + C);
				uv_pl[0] = z * a0;
				uv_pl[1] = z * a1;
				uv_pl[2] = z;
			}
			else if (v != 0)
			{
				y = D / (A*a0 + B);
				uv_pl[0] = y * a0;
				uv_pl[1] = y;
				uv_pl[2] = 0;
			}
			else
			{
				uv_pl[0] = D/A;
				uv_pl[1] = 0;
				uv_pl[2] = 0;
			}

		/* get intersection of plane defined by triangular
		   facet and line defined by the sphere center and
		   the test point.
		*/


			diff[0][0] = vert[i0] - uv_pl[0];
			diff[0][1] = vert[i0+n_vert] - uv_pl[1];
			diff[0][2] = vert[i0+2*n_vert] - uv_pl[2];

			diff[1][0] = vert[i1] - uv_pl[0];
			diff[1][1] = vert[i1+n_vert] - uv_pl[1];
			diff[1][2] = vert[i1+2*n_vert] - uv_pl[2];

			diff[2][0] = vert[i2] - uv_pl[0];
			diff[2][1] = vert[i2+n_vert] - uv_pl[1];
			diff[2][2] = vert[i2+2*n_vert] - uv_pl[2];

		/* get difference between this point and vertices
		   of triangular facet
		*/


			for (j=0; j<3; j++)
			{
				a = &edge[9*tri+3*j];
				b = &diff[j][0];
				test = (a[1]*b[2] - a[2]*b[1]) * A +
				       (a[2]*b[0] - a[0]*b[2]) * B +
				       (a[0]*b[1] - a[1]*b[0]) * C;
				if (test < 0) break;

		/* For each of these difference vectors
		   take cross product of edge vector to next vertex
		   and this difference.
		   If in opposite dirction of normal outward from
		   triangular face then outside triangle.
		   Check next facet.
		*/

			}

			if (test > 0)
			{
				flag[i] = 1;
				break;

		/* Inside facet therefore set flag 
		   and go on to next point */

			}

		}	/* tri loop */


	}	/* test point loop */

return(retdat);
}
