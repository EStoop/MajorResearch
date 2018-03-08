pro gain_calstar

ch1list = '/Users/jkrick/IRAC/iwic210/callist_ch1.txt'
ch2list = '/Users/jkrick/IRAC/iwic210/callist_ch2.txt'

readcol, ch1list, calch1, format="A"
readcol, ch2list, calch2, format = "A"

ch1name = strcompress('/Users/jkrick/IRAC/iwic210/NPM1p67.0636/' + calch1, /remove_all)
ch2name = strcompress('/Users/jkrick/IRAC/iwic210/NPM1p67.0636/' + calch2, /remove_all)

ra = 269.727857
dec = 67.793591


for i = 0, n_elements(ch1name)  - 1 do begin
   fits_read, ch1name(i), coldcaldata, coldcalheader

   ;find the bright star in the frame
   adxy, coldcalheader, ra, dec, x, y
   
   ;run aperture photometry on the star if it is in the field
   if x gt 5 and x lt 250 and y gt 5 and y lt 250 then begin 
      ;doesn't matter which gain I use here, only pertinent for error calculation
      aper, coldcaldata, x, y, flux, fluxerr, sky, skyerr,3.3, [10.], [12,20],[-100,40000], /nan, /exact, /flux, /silent

      ;correct flux for position dependant systematic error
                                ;how to do this efficient;y?  case
                                ;work with 2 variables?  also have
                                ;adifferent table for ch1 and ch2?
                                ;How do I code a table lookup?
   endif

endfor


end
