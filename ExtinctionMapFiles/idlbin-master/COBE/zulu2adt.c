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

/****************************************************************************

/* MODULE:	zulu2adt.c

   FUNCTION:	This program converts ZULU times (YYDDDHHMMSSSSS)
		to the VAX ADT time format.  It is called by the IDL time
                conversion facility, TIMECONV.

        AUTHOR:	J.M. Gales (Applied Research Corp)
	DATE:	11/92, SPR 10212

	MODIFICATION HISTORY:
	DATE	SPR/SER		PROGRAMMER	DESCRIPTION
	----	-------		----------	-----------
    01/05/93      10418         J.M. Gales      Change times counter
                                                from short int to
                                                IDL_LONG
      02/12/93  11468           Newmark         Change to compile under
                                                VAX and Sun C not GNU C
	03/94	11614		Dan Grogan	Changed "export.h" to 
						"export_wrap.h" for IDL 3.0+
       09/95                    Newmark         Upgrade to IDL 4.0


   COMMAND LINE ARGUMENT: zulu2adt, zulu_time,adt_time (returned),
				    return status (returned)


   PROCEDURE PARAMETERS:

	Name		Type		Description
	----		----		-----------
	zulu_time	IDL_STRING arr	Zulu ('YYDDDHHMMSSSSS') array

	adt_time	IDL_LONG arr    ADT output long integer array

	return status   int arr		return status array


   RETURN:		None

   FUNCTIONS CALLED:	c_addx (internal)
			c_emul (internal)

*******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "export_wrap.h"

void zulu2adt(argc,argv)

int	argc;
IDL_VPTR	argv[];


{

IDL_ULONG	base_time[2] ;
			/* # of intervals to 01-JAN-1985 */

IDL_ULONG	num_sec_n[2] ;
			/* # of intervals in non-leap year */

IDL_ULONG	num_sec_l[2] ;
			/* # of intervals in leap year */

IDL_ULONG	sec_in_dy[2] ;
			/* # of intervals in day */

IDL_ULONG	sec_in_hr[2] ;
			/* # of intervals in hour */

IDL_ULONG	sec_in_mn[2] ;
			/* # of intervals in minute */

IDL_ULONG	prod[2], sec_in_yr[2];

int		i;
int		yr, dy, hr, mn, leap, d_of_y;
IDL_LONG	sec, n_times, k, dim[1];

IDL_ULONG	bin_time[2];
char			zulu_time[17],u_str[17];
short			istat;

/* IDL pointer VARIABLEs */

IDL_STRING			*zulu_ptr;
IDL_ULONG		*adt_ptr;
short			*i_ptr;
unsigned short		str_len;

IDL_VARIABLE		*retdat;
char			*zulu, *oupt;

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


	zulu_ptr = (IDL_STRING *) argv[0]->value.arr->data;
			/* get pointer to zulu times */

	i_ptr = (short *) argv[2]->value.arr->data;
			/* get pointer to return status array */

	n_times = (IDL_LONG) argv[2]->value.arr->n_elts;
			/* get number of times from status array */

	dim[0] = 2 * n_times;
			/* ADT times are TWO IDL_LONG words each */

	adt_ptr = (IDL_ULONG *) IDL_MakeTempArray(IDL_TYP_LONG, 
			1,dim,IDL_BARR_INI_ZERO, &retdat);
			/* allocate space for adt return array */

	IDL_VarCopy(retdat,argv[1]);
			/* pass pointer to 2nd argument */


	for (k=0; k<=n_times-1; ++k)
			/* loop through all times */
	{

	istat = 0;
	str_len = (unsigned short) zulu_ptr->slen;
			/* get STRING length */

	zulu = IDL_STRING_STR(zulu_ptr);
	for (i=0; i<=str_len; ++i) zulu_time[i] = *zulu++;
			/* get pointer to STRING from descriptor ptr
			   load into zulu_time array */

	if (toupper(zulu_time[str_len-1]) == 'E')
	{
		for (i=str_len-1; i<=15; ++i) zulu_time[i] = '0';
		zulu_time[16] = '\0';
		str_len = 16;
	}
			/* If last character 'E' then pad with zeros
			   place end of STRING delineator at end */

	if (str_len < 16)
	{
		for (i=0; i<str_len; ++i) zulu_time[str_len+1-i] = 
					  zulu_time[str_len-1-i];
		zulu_time[0] = '1';
		zulu_time[1] = '9';
		for (i=str_len+2; i<=15; ++i) zulu_time[i] = '0';
		zulu_time[16] = '\0';
		str_len = 16;
	}
			/* If length < 16 add '19' to beginning
			   pad with zeros */

	zulu_ptr++;
			/* point to next zulu STRING */


	for (i=0; i<=3; ++i) u_str[i] = zulu_time[i];
	u_str[4] = '\0';
	yr = atoi(u_str);
			/* extract year & change to integer */

	for (i=0; i<=2; ++i) u_str[i] = zulu_time[i+4];
	u_str[3] = '\0';
	dy = atoi(u_str);
			/* extract day & change to integer */

	for (i=0; i<=1; ++i) u_str[i] = zulu_time[i+7];
	u_str[2] = '\0';
	hr = atoi(u_str);
			/* extract hour & change to integer */

	for (i=0; i<=1; ++i) u_str[i] = zulu_time[i+9];
	u_str[2] = '\0';
	mn = atoi(u_str);
			/* extract minute & change to integer */

	for (i=0; i<=4; ++i) u_str[i] = zulu_time[i+11];
	u_str[5] = '\0';
	sec = atol(u_str);
			/* extract sec(x10^3)  & change to LONG integer */


	bin_time[0] = base_time[0];
	bin_time[1] = base_time[1];
			/* initialize adt to base time */

	if (yr < 1985)
	{
		istat = -1;
		printf("Binary time before 01-JAN-1985");
			/* ADT time before base time  */
			/* set return time to all zeros */
			/* set return status to bad & exit */
	}
	else
	{

		for (i=1985; i<yr;  ++i)
		{
		leap = 0;
		if (i % 4   == 0) {leap = -1;}
		if (i % 100 == 0) {leap = -0;}
		if (i % 400 == 0) {leap = -1;}

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

		c_addx(bin_time,sec_in_yr,bin_time,2);
			/* add year intervals */

		}  /* end year loop */

	}  /* end year if */


	c_emul(sec_in_dy,dy-1,prod);
		/* multiply # of days times units per day */
	c_addx(bin_time,prod,bin_time,2);
		/* add days */


	c_emul(sec_in_hr,hr,prod);
		/* multiply # of hours times units per hour */
	c_addx(bin_time,prod,bin_time,2);
		/* add hours */


	c_emul(sec_in_mn,mn,prod);
		/* multiply # of minutes times units per minute */
	c_addx(bin_time,prod,bin_time,2);
		/*add minutes */


	prod[0] = sec * 10000l;
		/* convert to 100 nanosecond units */
	prod[1] = 0;
	c_addx(bin_time,prod,bin_time,2);
		/* add seconds */


	*i_ptr++ = istat;
		/* store status */

	*adt_ptr++ = bin_time[0];
	*adt_ptr++ = bin_time[1];
		/* store in adt array and increment adt output pointer */

}	/* end main loop */


return;
}


c_emul(a,b,prod)
/*
;  NAME:
;    c_emul
;
;  PURPOSE:                                   
;    Performs quadword integer multiplication
;
;  INPUT:
;    multiplicand 1 - quadword (2 dim LONGword array)
;    multiplicand 2 - LONGword
;
;  OUTPUT:
;    product - quadword (2 dim LONGword array)
;
;  SUBROUTINES CALLED:
;    c_addx
;
;  REVISION HISTORY
;    J.M. Gales
;    Sep 92		Initial Release
;-
*/
IDL_ULONG	a[], b, prod[];

{

IDL_ULONG	a0[2], mask;
int		i;


	a0[0] = a[0];
	a0[1] = a[1];

	prod[0] = 0;
	prod[1] = 0;
			/* initialize product */


	mask = 1;
	for (i=0; i<=31; ++i)
	{
		if ((b & mask) == mask)	c_addx(prod,a0,prod,2);

		mask = mask << 1;

		a0[1] = a0[1] << 1;
		a0[1] = a0[1] + ((a0[0] & 0x80000000) >> 31);
		a0[0] = a0[0] << 1;
			/* quadword shift left 1 bit */
	}

return;

}


c_addx(a, b, c, n)

int		n;
IDL_ULONG	a[],b[],c[];

{

int		i;
char		carry[20];

for (i=0;i<=19;++i)
  carry[i]=0;

carry[0] = 0;

for (i=0; i<=n-1; ++i)
{
	c[i] = a[i] + b[i] + carry[i];

	if ((c[i] < a[i]) || (c[i] < b[i])) carry[i+1] = 1;

}	/* endfor */

return;
}
