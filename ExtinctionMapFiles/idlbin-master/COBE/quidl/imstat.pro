pro imstat, data, avg, sigma, med, minimum, maximum, npix
;
;+NAME/ONE-LINE DESCRIPTION:
;  IMSTAT returns useful statistics on an image.
;
;DESCRIPTION:
;  IMSTAT calculates a number of useful statistical quantities on an
;  input image.  It will print out the values of avg(image),
;  stdev(image), median(image), min(image), max(image), npixels(image)
;
;CALLING SEQUENCE:
;  IMSTAT(data, avg, stdev, median, min, max, npix)
;
;ARGUMENTS: (I=input, O=output, []=optional)
;  data		I	1-D or 2-D input array
;  [avg]	O	average value of data
;  [stdev]	O	standard deviation of data
;  [median]	O	median of data
;  [min]	O	minimum of data
;  [max]	O	maximum of data
;  [npix]	O	number of pixels in data array
;
;EXAMPLES:
;	imstat,image
;       This prints out statistical values to the screen only.
;
;	imstat,image,mean,stdev,med,minval,maxval,npix
;	This will print out the values to the screen and place them in
;	named variables for later use.
;# 
;SUBROUTINES CALLED:
;  STDEV
;
;COMMON BLOCKS:
;
;LIBRARY CALLS:
; MEDIAN, MIN, MAX
;
;WARNINGS:
;
;PROGRAMMING NOTES:
;
; MODIFICATION HISTORY:
;  Written by Dave Bazell, General Sciences Corporation, Nov. 1993
;
;.TITLE
;Routine IMSTAT
;-
;
ON_ERROR, 2
;
; Get the number of pixels in the array
npix = n_elements(data)
;
; Calculate the mean and standard deviation
sigma = stdev(data,avg)
;
; Calculate the median
med = median(data)
;
; Min and Max
minimum = min(data)
maximum = max(data)
;
; Print it out
print,'      Mean     Stdev       Median       Min         Max         Npix'
print,avg, sigma, med, minimum, maximum, npix, format='(5(e12.4),i8) '

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


