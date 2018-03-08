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
/* MODULE:	Adt2zulu.c

   FUNCTION:	This program converts ADT times to ZULU format
                (YYDDDHHMMSSSSS).  It is called by the IDL time
                conversion facility, TIMECONV.  It is passed three arrays, 
		the LONGword array containing the ADT times, the STRING 
		array returning the ZULU times, and the return status array.  

   AUTHOR:	J.M. Gales (Applied Research Corp)
	DATE:	4/92  SPR 9687 This routine, called by IDL time conversion rou-
			       tines, is linked to IDL with the Dataserver (UFC)

	MODIFICATION HISTORY:
	DATE	SPR/SER		PROGRAMMER	DESCRIPTION
	----	-------		----------	-----------
      11/9/92   10212           Gales           Renamed from ADT2ZULU_C
      12/10/92  10337           Gales           Suppress seconds roundup
                                                when sec >= 59.9995 s
      01/05/93  10418           Gales           Change times counter
                                                from short int to IDL_LONG 
      02/12/93  11468           Newmark         Change to compile under
                                                VAX and Sun C not GNU C
	03/94	11614		Dan Grogan	Changed "export.h" to 
						"export_wrap.h" for IDL 3.0+
						Add IDL_STRING* cast in IDL_StrStore()
       09/95                    Newmark         Upgrade to IDL 4.0

   COMMAND LINE ARGUMENT: ADT time, zulu time (returned)
			  return status (returned)

   RETURN:		None

   FUNCTIONS CALLED:	c_subx (internal)
			c_ediv (internal)

*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include "export_wrap.h"

void adt2zulu(argc,argv)

int	argc;
IDL_VPTR	argv[];


{

IDL_ULONG       base_time[2] ;
			/* # of intervals to 01-JAN-1985 */

IDL_ULONG	num_sec_n[2]; 
			/* # of intervals in non-leap year */

IDL_ULONG	num_sec_l[2] ;
			/* # of intervals in leap year */

IDL_ULONG	sec_in_dy[2] ;
			/* # of intervals in day */

IDL_ULONG	sec_in_hr[2] ;
			/* # of intervals in hour */

IDL_ULONG	sec_in_mn[2];
			/* # of intervals in minute */

IDL_ULONG	quot[2], sec_rem[2], sec_rem0[2], sec_in_yr[2];

int		i, temp_len;
int		yr, dy, hr, mn, leap, d_of_y;
IDL_LONG	sec, n_times, k;
char		borrow;
char		time_str[13];
char		temp[6], u_str[17], yr_str[5];

IDL_ULONG	bin_time[2];
char			zulu_time[17];
short			istat;

/* IDL pointer IDL_VARIABLEs */

short			*i_ptr;
IDL_ULONG	*bin_ptr;

IDL_VARIABLE		*retdat;
char			*zulu;
IDL_STRING                  *oupt;

base_time[0] = 154632192l;
base_time[1] = 9266898l;

num_sec_n[0] = 2026291200l;
num_sec_n[1] = 73425l;

num_sec_l[0] = -1557102592l;
num_sec_l[1] = 73626l;

sec_in_dy[0] = 711573504l;
sec_in_dy[1] = 201l;

sec_in_hr[0] = 1640261632l;
sec_in_hr[1] = 8l;

sec_in_mn[0] = 600000000l;
sec_in_mn[1] = 0l;

	bin_ptr = (IDL_ULONG *) argv[0]->value.arr->data;
			/* get pointer to ADT times */

	i_ptr = (short *) argv[2]->value.arr->data;
			/* get pointer to return status array */

	n_times = (IDL_LONG) argv[2]->value.arr->n_elts;
			/* get number of ADT times from status array */

	
	zulu = &zulu_time[0];
			/* assign zulu pointer to 
			   zulu_time work array */

	oupt = (IDL_STRING *) IDL_MakeTempArray(IDL_TYP_STRING, 
			argv[2]->value.arr->n_dim,
			argv[2]->value.arr->dim,
			IDL_BARR_INI_NOP, &retdat);
			/* allocate space for zulu return array
			   using info from status array */

	IDL_VarCopy(retdat,argv[1]);
			/* pass pointer to 2nd argument */


	for (k=0; k<=n_times-1; ++k)
			/* loop through all times */
	{
	istat = 0;
	bin_time[0] = *bin_ptr++;
	bin_time[1] = *bin_ptr++;
			/* initialize status
			   get both words of ADT times */


	c_subx(bin_time,base_time,sec_rem,2,&borrow);
				/* subtract base time */

	if (borrow == 1)
	{
		strcpy(zulu_time,"0000000000000000");
		istat = -1;
		printf("Binary time before 01-JAN-1985");
			/* ADT time before base time  */
			/* set return time to all zeros */
			/* set return status to bad & exit */
	}
	else
	{

		for (yr=1985; yr<=9999;  ++yr)
		{
		leap = 0;
		if (yr % 4   == 0) {leap = -1;}
		if (yr % 100 == 0) {leap = -0;}
		if (yr % 400 == 0) {leap = -1;}

		if (leap == -1)
		{
			sec_in_yr[0] = num_sec_l[0];
			sec_in_yr[1] = num_sec_l[1];
		}
		else
		{
			sec_in_yr[0] = num_sec_n[0];
			sec_in_yr[1] = num_sec_n[1];
		}
			/* decide whether leap year or not
			   get corresponding secs in year */

		c_subx(sec_rem,sec_in_yr,sec_rem0,2,&borrow);
			/* subtract year intervals */

		if (borrow == 0)
		{
			sec_rem[0] = sec_rem0[0];
			sec_rem[1] = sec_rem0[1];
		}
		else
		{
			sprintf(yr_str,"%d",yr);
				/* write yr as STRING */

			goto days;
		}

		} /* endfor */



days:	c_ediv(sec_rem,sec_in_dy,quot,sec_rem);
			/* divide remaining intervals
			   by intervals per day */

	strcpy(time_str,"000000000000");
			/* initialize time STRING with '0' fill */

	d_of_y = quot[0] + 1;
	sprintf(temp,"%d",d_of_y);	/* convert int to STRING */
	temp_len = strlen(temp);	/* get length of STRING */
	for (i=0; i<=temp_len-1; ++i)
	{
		time_str[i+3-temp_len] = temp[i];
	}
			/* copy to time STRING padded with
			   leading zeros */


	c_ediv(sec_rem,sec_in_hr,quot,sec_rem);
			/* divide remaining intervals
			   by intervals per hour */


	hr = quot[0];
	sprintf(temp,"%d",hr);
	temp_len = strlen(temp);
	for (i=0; i<=temp_len-1; ++i)
	{
		time_str[i+5-temp_len] = temp[i];
	}

	c_ediv(sec_rem,sec_in_mn,quot,sec_rem);
		/* divide remaining intervals
		   by intervals per minute */


	mn = quot[0];
	sprintf(temp,"%d",mn);
	temp_len = strlen(temp);
	for (i=0; i<=temp_len-1; ++i)
	{
		time_str[i+7-temp_len] = temp[i];
	}


	if (sec_rem[0] < 599995000l) sec_rem[0] = sec_rem[0] + 5000;
			/* if within a half millisecond of 60 sec
			   then don't add roundup */
	sec = sec_rem[0] / 10000l;
	sprintf(temp,"%d",sec);
	temp_len = strlen(temp);

	for (i=0; i<=temp_len-1; ++i)
	{
		time_str[i+12-temp_len] = temp[i];
	}


	strcpy(u_str,yr_str);		/* copy yr_str to u_str */
	strcat(u_str,time_str);		/* concatonate STRINGs */
	strcpy(zulu_time,u_str);	/* copy u_str to zulu time */

}	/* end else */

	*i_ptr++ = istat;
		/* store status */

	IDL_StrStore((IDL_STRING *)oupt,zulu);
		/* store zulu time in return array */

	oupt++;
		/* increment STRING output pointer */

}	/* end main loop */


return;
}

c_ediv(a,b,quot,rem)

/*
;  NAME:
;    c_ediv
;
;  PURPOSE:                                   
;    Performs quadword integer division
;
;  INPUT:
;    dividend - quadword (2 dim LONGword array)
;    divisor - quadword (2 dim LONGword array)
;
;  OUTPUT:
;    quotient - quadword (2 dim LONGword array)
;    remainder - quadword (2 dim LONGword array)
;
;  SUBROUTINES CALLED:
;    c_subx
;
;  REVISION HISTORY
;    J.M. Gales
;    Feb 92		Initial Release
;-
*/

IDL_ULONG	a[], b[], quot[], rem[];

{

IDL_ULONG	a0[2], b0[2], mask;
int		i, j, bit, bit0;
int		first_1;
char		brw_bit;

bit0=0;

	a0[0] = a[0];
	a0[1] = a[1];
	b0[0] = b[0];
	b0[1] = b[1];

	quot[0] = 0;
	quot[1] = 0;
			/* initialize quotient */


	if ((b[0] == 0) && (b[1] == 0))
	{
		printf("Division by zero");
		rem[0] = a[0];
		rem[1] = a[1];
		return;
	}
			/* division by zero
			   return dividend as remainder */


	if (a[1] < b[1])
	{
		rem[0] = a[0];
		rem[1] = a[1];
		return;
	}
			/* dividend less than divisor
			   return dividend as remainder */


	if ((a[1] == 0) && (b[1] == 0))
	{
		quot[0] = a[0] / b[0];
		rem[0] = a[0] - quot[0]*b[0];
		rem[1] = 0;
		return;
	}
			/* standard LONGword division */


/*	-----------------------
	start quadword division
	-----------------------		*/

	/* Find first nonzero bit in dividend starting from left */

	first_1 = 31;
	mask = 0x80000000;

	while ((a0[1] & mask) != mask)	/*while not found*/
	{
		first_1 = first_1 - 1;	/* decrement bit */
		mask = mask >> 1;	/* build new mask */

	}


/*	Shift divisor to left until first nonzero bits align	*/

	while ((b0[1] & mask) != mask)
	{

		b0[1] = b0[1] << 1;
		b0[1] = b0[1] + ((b0[0] & 0x80000000) >> 31);
		b0[0] = b0[0] << 1;
			/* quadword shift left 1 bit */

		bit0 = bit0 + 1;	/* increment bit counter */
	}


	for (i=bit0; i>=0; --i)
	{
		j = i / 32;		/*get LONGword index*/
		bit = i % 32;		/*get bit number within LONGword*/

		c_subx(a0,b0,rem,2,&brw_bit);
					/*subtract shifted divisor from
					  remaining dividend*/


		if (brw_bit == 0)		/* if a0 > b0 */
		{
			a0[0] = rem[0];
			a0[1] = rem[1];		/*new remaining dividend*/

			quot[j] = quot[j] + (1l << bit);
						/* set bit in quotient */
		}

		b0[0] = b0[0] >> 1;
		b0[0] = b0[0] + ((b0[1] & 1) << 31);
		b0[1] = b0[1] >> 1;
						/* shift 'divisor' right
						   one bit */


	}	/* endfor */

	rem[0] = a0[0];
	rem[1] = a0[1]; 			/* get final remainder */

return;

}

c_subx(a, b, c, n, borrow_bit)

int		n;
IDL_ULONG	a[],b[],c[];
char		*borrow_bit;

{

int		i;
IDL_ULONG	a0[2];
char		borrow[20];

for (i=0; i<=19; ++i) 
  borrow[i]=0;

a0[0] = a[0];
a0[1] = a[1];

for (i=0; i<=n-1; ++i)
{
	if (a0[i] == 0) {borrow[i+1] = borrow[i];}
		/* propogate borrow bit if 0 in minuend */

	a0[i] = a0[i] - borrow[i];
		/* subtract borrow bit for minuend */

	if (a0[i] < b[i]) {borrow[i+1] = 1;}
		/* set borrow bit if subtrahend greater than minuend */

	c[i] = a0[i] - b[i];
		/* calculate difference */

}	/* endfor */

*borrow_bit = borrow[n];

return;
}
