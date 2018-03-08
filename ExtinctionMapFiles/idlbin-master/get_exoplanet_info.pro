FUNCTION get_exoplanet_info_Url_Callback, status, progress, data
 
   ; print the info msgs from the url object
   ;;PRINT, status
 
   ; return 1 to continue, return 0 to cancel
   RETURN, 1
END
 
;-----------------------------------------------------------------
FUNCTION get_exoplanet_info, myurl_host, myurl_path, myurl_query

   CATCH, errorStatus
   IF (errorStatus NE 0) THEN BEGIN
      CATCH, /CANCEL
 
      ; Display the error msg in a dialog and in the IDL output log
      r = DIALOG_MESSAGE(!ERROR_STATE.msg, TITLE='URL Error.  Will just use local file.', $
         /ERROR)
      PRINT, !ERROR_STATE.msg
 
      ; Get the properties that will tell us more about the error.
      oUrl->GetProperty, RESPONSE_CODE=rspCode, $
         RESPONSE_HEADER=rspHdr, RESPONSE_FILENAME=rspFn
      PRINT, 'rspCode = ', rspCode
      PRINT, 'rspHdr= ', rspHdr
      PRINT, 'rspFn= ', rspFn
 
      ; Destroy the url object
      OBJ_DESTROY, oUrl
      RETURN,-1
   ENDIF
 
   ; create a new IDLnetURL object
   oUrl = OBJ_NEW('IDLnetUrl')
 
   ; Specify the callback function
   oUrl->SetProperty, CALLBACK_FUNCTION ='get_exoplanet_info_Url_Callback'
 
   ; Set verbose to 1 to see more info on the transacton
   oUrl->SetProperty, VERBOSE = 0
 
   ; Set the transfer protocol as https
   oUrl->SetProperty, url_scheme = 'https'
 
   ; The Exoplanets.org host
   oUrl->SetProperty, URL_HOST = myurl_host;'exoplanetarchive.ipac.caltech.edu'
 
   ; The FTP server path of the file to download
   oUrl->SetProperty, URL_PATH = myurl_path; 'cgi-bin/nstedAPI/nph-nstedAPI'
   
   oUrl->SetProperty, URL_QUERY = myurl_query; strcompress('table=exoplanets&ra='+string(ra)+'&dec='+string(dec)+'&radius=4%20minutes&select=pl_hostname,st_pmra,st_pmdec',/remove_all)
 
   ; Make a request to the server, output as a string array (each line of output is an element of the array).
    result = oUrl->Get( /STRING_ARRAY )
 
   ; Print the path to the file retrieved from the remote server
;   PRINT, 'filename returned = ', filename
 
   ; Destroy the url object
   OBJ_DESTROY, oUrl
 
 RETURN,result
END

;;http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=exoplanets&ra=227.68531&dec=-61.422365&radius=0.400000%20minutes&select=pl_hostname,st_pmra,st_pmdec

;;http://simbad.u-strasbg.fr/simbad/sim-coo?Coord=227.68531-61.422365&Radius=0.400000&Radius.unit=arcmin&output.format=ASCII&coodisp1=d&list.pmsel=on
