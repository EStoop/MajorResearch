PRO NICE_HISTO,DATASET,RETURN_HISTO,SCALE,NUM_BINS
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    NICE_HISTO returns a histogram and corresponding x-vector.
;
;DESCRIPTION:
;    Nice_histo is designed to make creation and plotting of histo-
;    grams easier.  The program takes an array of data, scales it 
;    between -1 and 1, forms a histogram with that data with the 
;    (optional) desired number of bins, and also returns a scale 
;    vector of the same length as the histogram running between the 
;    minimum and maximum values in the original array.  The number of 
;    bins defaults to 200 if not specified.  The user may then easily 
;    plot the histogram with the true scaling of the array displayed.
;
;CALLING SEQUENCE:  
;    NICE_HISTO, DATASET, RETURN_HISTO [, SCALE] [,NUM_BINS]
;
;ARGUMENTS: (I = input, O = output, [] = optional)
;     DATASET        I   flt or dbl  The vector containing the 
;                            [arr]   information to be put into 
;                                    histogram form.
;     NUM_BINS      [I]  int         Number of bins to use in the 
;                                    histogram (defaults to 200).
;     RETURN_HISTO   O   int [arr]   The returned histogram
;     SCALE         [O]  dbl [arr]   The returned scaling vector
;
;WARNINGS:
;     None
;
;EXAMPLES:
;     1.  To produce a histogram having 200 bins from vector MYDATA, 
;     returning the histogram vector MY_HISTODATA and scaling vector 
;     SCALE_VEC:
;
;     NICE_HISTO,MYDATA,MY_HISTODATA,SCALE_VEC
;
;     2.  To produce a histogram having 500 bins from vector MYDATA, 
;     returning the histogram vector MY_HISTODATA and scaling vector 
;     SCALE_VEC:
;
;     NICE_HISTO,MYDATA,MY_HISTODATA,SCALE_VEC,500
;
;     In both the above the examples, the user may then plot the 
;     histogram using the scaling vector to show the binned values 
;     using the command
;        PLOT, SCALE_VEC, MY_HISTODATA
;
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     Calculates a bin size for the histogram procedure by dividing 
;     the difference between largest and smallest values in the input
;     array by the desired number of bins.  Calculates histogram with
;     values running between minimum and maximum of array; calculates
;     scale as straight line between minimum and maximum values with 
;     number of elements equal to number of bins.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     NONE
;
;MODIFICATION HISTORY:
;     Alice Trenholme, General Sciences Corp. February, 1991
;.TITLE
;Routine NICE_HISTO
;-
if n_elements(dataset) eq 0 then begin
        print,' '
        print,'Calling sequence: NICE_HISTO,DATASET,RETURN_HISTO,SCALE,NUM_BINS'
        print,'  (Number of bins is optional and defaults to 200)'
        print,'  (Scale is the optional returned scaling vector)'
        print,' '

endif else begin
        if n_elements(num_bins) eq 0 then num_bins = 200

        maxx=max(dataset)
        minn=min(dataset)
        
	;The constant -.1 here ensures that the desired # of bins is made.... 
        binz = double(maxx-minn)/(num_bins-.1)

        RETURN_HISTO = histogram(dataset,binsize=binz,min=minn,max=maxx) 
        nn=n_elements(RETURN_HISTO)
        scale=dindgen(nn)/(nn-1)
        scale=minn+(maxx-minn)*scale
endelse

end
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


