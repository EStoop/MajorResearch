        FUNCTION BITRANGE, inmap, bitlo=bitlo, bithi=bithi, datlo=datlo, $
                                   dathi=dathi, nostretch=nostretch, $
                                   badval=badval
;
;+NAME/ONE-LINE DESCRIPTION:
;   BITRANGE takes an arbitrary range of bits from an array for display
;
; DESCRIPTION:
;   The routine extracts a user-selected range of up to 8 bits from an 
;   assumed 16-bit dynamic range, and creates an output array scaled to 
;   a single byte (0-255).  With a suitable choice of bits (e.g. 
;   chopping out some high-order bits) this operation enhances low-level
;   map features by causing the grey scale or color table to wrap around 
;   repeatedly, creating a "color contour" effect.
;
;   The calling sequence allows the user to specify the bad pixel value
;   (all such pixels will be set to dark) the range of bits, and the 
;   range of data to operate on (by default, the whole good data range). 
;   If the input data is integer and has an intrinsic dynamic range 
;   smaller than 16 bits, the /NOSTRETCH qualifier will cause the 
;   routine to extract bits from the raw data rather than first scaling 
;   it to a 16-bit range. Otherwise, the rescaling is done first.
;
;   All real numeric data types are accepted as input (I*2, I*4, R*4,
;   BYTE), but the output is always in BYTE format.
;
; CALLING SEQUENCE:
;   ARROUT = BITRANGE (arrin, [bitlo=...], [bithi=...], [datlo=...],$
;                             [dathi=...], [badval=...], [/NOSTRETCH])
;
; ARGUMENTS (I=input, O=output, []=optional)
;   ARROUT        O  byte array   Output map or vector, range 0-255
;   arrin         I  any  array   Input map or vector, any real data 
;   [bitlo]      [I] int          Lo-order bit default=8, min=0
;   [bithi]      [I] int          Hi-order bit default=15, max=15
;   [datlo]      [I] int or real  Minimum value of data to consider. 
;                                    Default=minimum value of good data
;   [dathi]      [I] int or real  Maximum value of data to consider. 
;                                    Default=maximum value of good data
;   [badval]     [I] int or real  Bad (i.e. excluded) data value.  Def=0
;   [/NOSTRETCH] [I] qualifier    Do not rescale integer data to 16 bits
;
; WARNINGS:
;    1.  No more than 8 bits can be selected.  Allowable bit range is
;          0-15
;    2.  Complex numbers will give unpredictable results.  Strings and
;          structures are forbidden.
;
; EXAMPLES:
;    Consider a floating point map IN with a bad value of -999 and good
;    data of interest ranging from 0 to 1000.  We can "quantize" the
;    display to 8 grey levels (3 bits) showing only the brightest 
;    features with the call
;
;           OUT = BITRANGE (IN, BADVAL=-999, BITLO=13, BITHI=15)
;
;    If we're only interested in the upper half of the data range but
;    want to emphasize features in the low end of that range, we can
;    expand the dynamic range to eight bits and chop off the highest
;    order bit, i.e.
;
;           OUT = BITRANGE (IN, DATLO=500, BITLO=7, BITHI=14)  
;
;    (Notice that we didn't need to specify BADVAL since the DATLO in
;    this case would exclude it anyway.)  If the input data were integer
;    and we were to give the same call with the /NOSTRETCH qualifier we 
;    would get a different result, because the raw data only extend up 
;    bit 9 (1000 < 2^10).  The previous call would then in effect extract
;    the three high-order bits (7-9).
;#
; COMMON BLOCKS:  none
;
; PROCEDURE:
;    The program creates a buffer array with a 16-bit range.  First, the
;    input is divided by 2^BITLO to chop out the unwanted low end bits.
;    Using modular arithmetic, the integer remainder is divided by the 
;    "width" of the bit slice, 2^(BITHI-BITLO+1).  The remainder is the
;    desired result, which is then scaled to a one-byte range.
;
; LIBRARY CALLS:  none
;                                
; MODIFICATION HISTORY:
;    Written 18 May 1992 by Rich Isaacman, General Sciences Corp.
;
;.TITLE
;Routine BITRANGE
;-
;  First set the deafults and create masks for the good and data data
;  locations.
;
        ON_ERROR, 2                               
        IF (N_ELEMENTS(badval) EQ 0) THEN badval = 0.
        goodmask = WHERE (inmap NE badval)
        badmask = WHERE (inmap EQ badval)
        IF (N_ELEMENTS(bitlo) EQ 0) THEN bitlo = 8
        IF (N_ELEMENTS(bithi) EQ 0) THEN bithi = 15
        IF (N_ELEMENTS(datlo) EQ 0) THEN datlo = MIN(inmap(goodmask))
        IF (N_ELEMENTS(dathi) EQ 0) THEN dathi = MAX(inmap(goodmask))
;
;  Make sure specified bit range makes sense
;
        IF ((bitlo GT bithi) OR (bitlo LT 0) OR (bithi GT 15) OR $
            (bithi-bitlo GT 7)) THEN $
                  MESSAGE, "Invalid bit range specified" 
;
;  Redefine the masks to exclude data outside the desired range
;
        goodmask = WHERE ((inmap NE badval) AND (inmap LE dathi) $
                                            AND (inmap GE datlo))
        badmask = WHERE ((inmap EQ badval) OR (inmap GT dathi) $
                                           OR (inmap LT datlo))
;                 
;  Scale the data to 16 bits (0-65535) if /NOSTRETCH is NOT present.  If 
;  it is present, simply set the data to LONG integer type and force
;  all values to a 16-bit range (-32768 - +32767).  Then crop off the
;  low-order bits by division.
;
        slice = 2^(bithi - bitlo + 1)
        IF (NOT KEYWORD_SET(nostretch)) THEN $
            buff = FLOAT(inmap - datlo)/(dathi - datlo) * (65535./2^bitlo) $
        ELSE buff = LONG(-32768 > inmap < 32767)/2^bitlo
;
;  Replace the bad values, the low-end points, and the too-bright high 
;  end points by the minimum acceptable value.  Then crop off the
;  high-order bits with modular arithmetic.
;
        IF (badmask(0) NE -1) THEN buff(badmask) = MIN(buff(goodmask))
        buff = LONG(buff MOD slice)
        buff = buff - MIN(buff)                ; Shift minimum to zero.
;      
;  Scale to a one-byte range for final result
;
        RETURN, BYTE(buff * 255./MAX(buff))
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


