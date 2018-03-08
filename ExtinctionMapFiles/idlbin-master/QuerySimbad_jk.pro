Function QuerySimbad_jk,  ra, dec, Found = found, ERRMSG = errmsg, simbad = simbad, $
    Verbose = verbose,  Print = print,Vmag=Vmag,Jmag=Jmag,Hmag=Hmag,Kmag=Kmag,parallax=parallax
;+
; NAME: 
;   QUERYSIMBAD
;
; PURPOSE: 
;   Query the SIMBAD/NED/Vizier astronomical name resolver to obtain coordinates
;
; EXPLANATION: 
;   Uses the IDL SOCKET command to query  the SIMBAD server 
;   over the Web to return a star name given the J2000 coordinates.  
;    
;   For details on the SIMBAD service, see http://simbad.u-strasbg.fr/Simbad 
;
; CALLING SEQUENCE: 
;     QuerySimbad_jk, 280.93945, -32.45678, /VERBOSE
;
; INPUTS: 
;     ra -  Right ascension of the target in J2000.0 in *degrees*, scalar 
;     dec - declination of the target in degrees, scalar
;
; OUTPUTS: 
;     starid - name of the nearest object
;  
; OPTIONAL INPUT KEYWORD:
;     ERRMSG   = If defined and passed, then any error messages will be
;                  returned to the user in this parameter rather than
;                  depending on the MESSAGE routine in IDL.  If no errors are
;                  encountered, then a null string is returned. 
;     /VERBOSE - If set, then the HTTP-GET command is displayed
;     /PRINT - if set, then output coordinates are displayed at the terminal 
;            By default, the coordinates are displayed if no output parameters
;           are supplied to QUERYSIMBAD
;     /SILENT - If set, then don't print warnings if multiple SIMBAD objects
;             correspond to the supplied name.
; OPTIONAL OUTPUT: 
;     id - the primary SIMBAD (or NED) ID of the target, scalar string
;          As of June 2009, a more reliable ID seems to be found when using 
;          CFA (/CFA) server.
;
; OPTIONAL KEYWORD OUTPUTS:
;     found - set to 1 if the translation was successful, or to 0 if the
;           the object name could not be translated by SIMBAD or NED
;     Errmsg - if supplied, then any error messages are returned in this
;            keyword, rather than being printed at the terminal.   May be either
;            a scalar or array.
;     Server - Character indicating which server was actually used to resolve
;           the object, 'S'imbad, 'N'ed or 'V'izier
;     Vmag - supply to receive the SIMBAD V magnitude 
;     Jmag - supply to receive the SIMBAD J magntiude
;     Hmag - supply to receive the SIMBAD H magnitude 
;     Kmag - supply to receive the SIMBAD K magnitude 
;     Parallax - supply to receive the SIMBAD parallax in milliarcseconds
;            
; EXAMPLES:
;     (1) Display the J2000 coordinates for the ultracompact HII region
;         G45.45+0.06 
;
;      IDL> QuerySimbad,'GAL045.45+00.06'
;           ===>19 14 20.77  +11 09  3.6
; PROCEDURES USED:
;       REPSTR(), WEBGET()
; NOTES:
;     The actual  query is made to the Sesame name resolver 
;     ( see http://cdsweb.u-strasbg.fr/doc/sesame.htx ).     The Sesame
;     resolver first searches the Simbad name resolver, then  NED and then
;     Vizier.   
; MODIFICATION HISTORY: 
;     Written by M. Feldt, Heidelberg, Oct 2001   <mfeldt@mpia.de>
;     Minor updates, W. Landsman   August 2002
;     Added option to use NED server, better parsing of SIMBAD names such as 
;          IRAS F10190+5349    W. Landsman  March 2003
;     Turn off extended name search for NED server, fix negative declination
;     with /NED    W. Landsman  April 2003
;     Use Simbad Sesame sever, add /Verbose, /CADC keywords 
;       B. Stecklum, TLS Tautenburg/ W. Landsman, Feb 2007
;    Update NED query to account for new IPAC format, A. Barth  March 2007
;    Update NED query to account for another new IPAC format, A. Barth  
;                                                   July 2007
;     Update message when NED does not find object  W.L.  October 2008
;     Remove CADC keyword, add CFA keyword, warning if more than two
;         matches  W.L. November 2008 
;     Make NED queries through the Sesame server, add Server output 
;          keyword  W.L.  June 2009
;     Don't get primary name if user didn't ask for it  W.L. Aug 2009
;     Added /SILENT keyword W.L. Oct 2009
;     Added /PRINT keyword W.L.   Oct 2011
;     Added ability to get V, J, H, and K magnitudes as well as
;     a parallax - jswift, Jan 2014
;-

;  compile_opt idl2
;  if N_params() LT 1 then begin
;     print,'Syntax - QuerySimbad, name, ra, dec, [ id, ]'
;     print,'                 Found=, /CFA, /NED, ERRMSG=, /VERBOSE]'
;     print,'   Input - object name, scalar string'
;     print,'   Output -  Ra, dec of object (degrees)'
;     return
;  endif
  
  Catch, theError
  IF theError NE 0 THEN BEGIN
     Catch,/CANCEL
     void = cgErrorMsg(/Quiet)
;     RETURN
  ENDIF   
  ;;
  printerr = ~arg_present(errmsg)
  if ~printerr  then errmsg = ''

  ;;to search SIMBAD
  ;;------------------------------------------------------
  if keyword_set(simbad) then begin
     base1 = "http://simbad.u-strasbg.fr/simbad/sim-coo?Coord="
     base2 = "&Radius=2&Radius.unit=arcmin&output.format=ASCII&coodisp1=d"
     QueryURL = strcompress(base1 + string(ra) +  string(dec) + base2,/remove)
  

  if keyword_set(verbose) then print,queryURL
  Result = webget(QueryURL)
  found = 0
  ;;
  if keyword_set(Verbose) then print, Result.Text
  Result=Result.Text

  idx=where(strpos(Result, 'Number of objects') ne -1, nlines)
  if nlines eq 1 then begin
     cnt = strmid(Result[idx], 20)
     ;;but cnt is a string, want it to be an int
     cnt = fix(cnt)
;     print, 'number of returned objects', cnt
  endif
  
  if cnt GE 1 then begin
     if cnt GT 1 then begin 
        message,/INF,strcompress('Warning - '+ string(cnt)+ ' matches found for position '  )
     endif 	
     found=1   

     ;;now pull the name of the first star in the list
     id2 = where(strpos(Result, '1|') ne -1, nlines)
     print, 'nlines 1', nlines
     print, 'Result: ', Result[id2]
     starname = strmid(Result[id2[0]],13,25)
      ;;and remove some blank spaces
     starname = strtrim(starname)

      if keyword_set(Verbose) then print, 'Nearest star: ', starname

   ENDIF ELSE BEGIN 
      errmsg = ['No objects returned by SIMBAD.   The server answered:' , $
                 strjoin(result)]
      if printerr then begin
         message, errmsg[0], /info	
	 message,strjoin(result),/info
      endif	 
  ENDELSE


  endif
  
  ;;to search NEXSCI
  ;;------------------------------------------------------
  
  base1 = "http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=exoplanets&ra="
  base2 = "&dec="
  base3 = "&radius=2%20minutes&select=pl_hostname"
  QueryURL = strcompress(base1 + string(ra) + base2 + string(dec) + base3,/remove_all)
  if keyword_set(verbose) then print,queryURL

  Result = webget(QueryURL)
  found = 0
  ;;Result=Result.Text
  if keyword_set(Verbose) then print, Result
  
  cnt = n_elements(strlen(result.text))
  if cnt ge 2 then begin
     if cnt gt 2 then print, strcompress('Warning - '+ string(cnt)+ ' matches found for position '  )
     firstline = strsplit(result.text[1],',',/extract)
     starname = firstline[0]

  endif else begin
       errmsg = ['No objects returned by SIMBAD.   The server answered:' , $
                 strjoin(result)]
      if printerr then begin
         message, errmsg[0], /info	
	 message,strjoin(result),/info
      endif	 
  ENDELSE

  return, starname
END 
  

    ; Get V mag if present
;;      vi = where(strpos(Result, '%M.V ') ne -1,vcnt)
;;      if vcnt GE 1 then reads,strmid(Result[vi],4),vmag

      ; Get J mag if present
;;      ji = where(strpos(Result, '%M.J ') ne -1,jcnt)
;;      if jcnt GE 1 then reads,strmid(Result[ji],4),jmag

      ; Get H mag if present
;;      hi = where(strpos(Result, '%M.H ') ne -1,hcnt)
;;      if hcnt GE 1 then reads,strmid(Result[hi],4),hmag

      ; Get K mag if present
;;      ki = where(strpos(Result, '%M.K ') ne -1,kcnt)
;;      if kcnt GE 1 then reads,strmid(Result[ki],4),kmag

      ; Get parallax if present
;;      plxi = where(strpos(Result, '%X ') ne -1,plxcnt)
;;      if plxcnt GE 1 then reads,strmid(Result[plxi],2),parallax
  
