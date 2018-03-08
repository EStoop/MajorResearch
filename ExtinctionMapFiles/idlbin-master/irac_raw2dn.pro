function irac_raw2dn,data,chan,barrel,fownum
;;
;; Convert IRAC raw fowler sampled data (from dce.fits or mipl.fits files) to floating point DN values.
;; Required inputs:  DATA : the raw data (array or scalar)
;;                   CHAN : the channel number (CHNLNUM in the fits header)
;;                  BARREL: The barrel shift (A0741D00 in the fits header)
;;                  FOWNUM: The fowler number (A0614D00 in the fits header)
;;           
;;  EXAMPLE:  
;;     IDL> dn_values = irac_raw2dn(raw_values,SXPAR(rawhead,'CHNLNUM'),SXPAR(rawhead,'A0741D00'),SXPAR(rawhead,'A0614D00'))
;;
;;     Raw values should have been read in using READFITS, which automatically converts the data to signed integer.
;; 
;; Convert back to unsigned integer   
result = uint(data)       
;;ineg = WHERE(result LT 0,nneg)
;; INSBPOSDOM (PDD Eq. 2.3)
result = 65535 - result
;; CVTI2R4 (16 bit integer to 32 bit real) (PDD EQ after 2.3)
result -= 0.5 * (1. -2.^(-barrel))  
;; Sign wraparound correction (IRACWRAPDET, IRACWRAPCORR)
iwrap = WHERE(result GT 55000.,nwrap)
IF nwrap NE 0 THEN result[iwrap] -= 65536  ;; PDD EQ 2.4, changing the constant from 65535 to 65536
;iwrap = WHERE(result LT -10535.,nwrap)
;IF nwrap NE 0 THEN result[iwrap] += 65536
;; Fowler sampling renormalization (IRACNORM)
IF fownum NE 2^barrel THEN result *=  2.^barrel / fownum  ;; PDD EQ 2.5

RETURN,result
END
