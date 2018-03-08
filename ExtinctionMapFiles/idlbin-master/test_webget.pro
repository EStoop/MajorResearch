FUNCTION exoplanet_light_curve_Url_Callback, status, progress, data
 
   ; print the info msgs from the url object
   PRINT, status
 
   ; return 1 to continue, return 0 to cancel
   RETURN, 1
END
 
pro test_webget

;would like a code that can access current planet characteristics and
;use those in my analysis.

;test = webget('http://exoplanetarchive.ipac.caltech.edu/cgi-bin/DisplayOverview/nph-DisplayOverview?objname=WASP-62+b&type=CONFIRMED_PLANET')
;print, test.text
;ra = 0.0
;dec = 0.9

;a=webget("http://ned.ipac.caltech.edu/cgi-bin/objsearch?search_type=Near+Position+Search&in_csys=Equatorial&in_equinox=J2000.0&lon=0.0&lat=0.0&radius=2.0&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES")

;a=webget('http://www.mpia.de/index.html')
;print, a

;a = webget("exoplanets.org/csv--files/exoplanets.csv")

;la=strsplit(a.text[26],"|",/extract) ;; first entry in the list is 26 lines down, after header
;print, la[1] ;; name is second column



;PointingRA=0.0
; PointingDE=30.0
; QueryURL = strcompress("http://archive.eso.org/dss/dss/image?ra="+$
;                         string(PointingRA)+$
;                          "&dec="+$
;                          string(PointingDE)+$
;                          "&x=10&y=10&Sky-Survey=DSS1&mime-type=download-fits", $
;                          /remove)
; a=webget(QueryURL)
; tvscl,a.Image
; print,a.ImageHead


;test code from Jim

   this_struc = ROUTINE_INFO('test_webget',/SOURCE)   ;;; Get the path to this program
   this_dir = FILE_DIRNAME(this_struc.path[0],/MARK_DIRECTORY)  ;;; Extract directory name
   filename = this_dir + 'exoplanets.csv'
;;; Only download the file once per day:
;;;    

   ; create a new IDLnetURL object
   oUrl = OBJ_NEW('IDLnetUrl')
 
   ; Specify the callback function
   oUrl->SetProperty, CALLBACK_FUNCTION ='exoplanet_light_curve_Url_Callback'
 
   ; Set verbose to 1 to see more info on the transacton
   oUrl->SetProperty, VERBOSE = 0
 
   ; Set the transfer protocol as http
   oUrl->SetProperty, url_scheme = 'http'
 
   ; The Exoplanets.org host
   oUrl->SetProperty, URL_HOST = 'exoplanets.org'
 
   ; The FTP server path of the file to download
   ;http://exoplanets.org/csv-files/exoplanets.csv
   oUrl->SetProperty, URL_PATH = 'csv-files/exoplanets.csv'
 
   ; Make a request to the Exoplanets.org server
   ; Retrieve the csv file and write it
   ; to the current working directory.
 
   filename = oUrl->Get(FILENAME = filename )
 
   ; Print the path to the file retrieved from the remote server
   PRINT, 'filename returned = ', filename
 
   ; Destroy the url object
   OBJ_DESTROY, oUrl
 


end
