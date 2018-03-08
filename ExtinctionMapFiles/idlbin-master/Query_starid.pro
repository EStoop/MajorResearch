Function Query_starid, naxishere, ra_rqst, dec_rqst, Found = found, ERRMSG = errmsg, $
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

  COMMON centroid_block
  
  ;;First search NEXSCI
  ;;------------------------------------------------------
  
  base1 = "http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=exoplanets&ra="
  base2 = "&dec="
  base3 = "&radius="
  base4 = "%20minutes&select=pl_hostname,st_pmra, st_pmdec"
  if naxishere eq 2 then dist = 4
  if naxishere eq 3 then dist = 0.4
  QueryURL = strcompress(base1 + string(ra) + base2 + string(dec) + base3 + string(dist) + base4,/remove_all)
  if keyword_set(verbose) then print,queryURL
  
  Result = webget(QueryURL)
  found = 0
  if keyword_set(Verbose) then print, 'Result', Result

  ;;how many objects were found?
  cnt = n_elements(strlen(result.text))

  ;;found at least one object
  if cnt ge 2 then begin
     if cnt gt 2 then print, strcompress('Warning - '+ string(cnt)+ ' matches found for position '  )
     firstline = strsplit(result.text[1],',',/extract,/preserve_null)
     starname = firstline[0]
     ra = firstline[3]
     dec = firstline[4]
     pmra = firstline[1]  ;mas/yr
     pmdec = firstline[2] ;mas/yr
     print, 'pmra and dec', pmra, pmdec
     ;;assume 2013 as the 'middle' of the warm mission to get
     ;;posistions closer
     pmra = pmra*13/1000./60./60. ;arsec
     pmdec = pmdec*13/1000./60./60. ;arcsec
     ra = ra + pmra
     dec = dec + pmdec
     
     ;;didn't find an object 
  endif else begin
     ;;-----------------------------------------------------------------------------------------
     print, 'No objects returned by NExScI. Searching SIMBAD in sweet spot' ;  The server answered:'+ (result)]
     
     base1 = "http://simbad.u-strasbg.fr/simbad/sim-coo?Coord="
     base2 = "&Radius="
     base3 = "&Radius.unit=arcmin&output.format=ASCII&coodisp1=d&list.pmsel=on"
     ;;if naxis eq '2' then narcmin = 2 else narcmin = 0.4
     narcmin = 0.4
     ;;dec needs to have a sign attached
     if dec lt 0 then sdec = string(dec) else sdec = strcompress('+'+string(dec),/remove)
     QueryURL = strcompress(base1 + string(ra) +  sdec + base2 + string(narcmin) + base3,/remove)
     if keyword_set(verbose) then print, 'query', QueryURL
     Result = webget(QueryURL)
     Result = Result.text
     if keyword_set(Verbose) then print, Result
     
     ;;how many objects did SIMBAD return?
     idx=where(strpos(Result, 'Number of objects') ne -1, nlines)
     if nlines eq 1 then begin  ;got more than one object
        scnt = strmid(Result[idx], 20)
        ;;but cnt is a string, want it to be an int
        scnt = fix(scnt)
        if keyword_set(Verbose) then print, 'SIMBAD number of matches: ', scnt
        
        if scnt GT 1 then begin 
           print,strcompress('Warning - '+ string(cnt)+ ' SIMBAD matches found for position '  )
        endif 	
        found=1   
        
        ;;now pull the name of the first star in the list
        id2 = where(Result.StartsWith('1|'), nlines)
        firstline = Result[id2[0]]
        junk = splitsimbadstring(firstline, starname=starname, ra=ra, dec=dec, pmra=pmra, pmdec=pmdec)
        
        ;;apply the proper motions
        ;;assume 2013 as the 'middle' of the warm mission to get
        ;;posistions closer
        pmra = pmra*13/1000./60./60. ;arsec
        pmdec = pmdec*13/1000./60./60. ;arcsec
        ra = ra + pmra
        dec = dec + pmdec

        if keyword_set(Verbose) then print, 'Nearest star: ', starname      
        
     ENDIF ELSE BEGIN  ;;SIMBAD found either one star, or didn't find anything, now what??
        
        ;;check for single star
        id1=where(strpos(Result, 'Object') ne -1, nlines)
        if nlines eq 1 then begin
           starname = strmid(Result[id1],7,11)
           idcoo = where(Result.StartsWith('Coordinates(ICRS'))
           cooline = Result[idcoo]
           pos = cooline.split(' ')
           ra = pos[1]
           dec = pos[2]
           btest = starname.StartsWith('NVSS')
           if btest gt 0 then GOTO, jumpif
        endif else begin
           jumpif: print, 're-running search because found NVSS target'
           print, ['No objects returned by SIMBAD.   The server answered:' , $
                   strjoin(result)]

           ;;not at the sweet spot, but maybe at the center of the
           ;;array
           if dec_rqst lt 0 then sdec_rqst = string(dec_rqst) else sdec_rqst = strcompress('+'+string(dec_rqst),/remove)
           narcmin = 1.4
           QueryURL = strcompress(base1 + string(ra_rqst) +  sdec_rqst + base2 + string(narcmin) + base3,/remove)
           print, 'query', QueryURL
           Result = webget(QueryURL)
           Result = Result.text
           idx=where(strpos(Result, 'Number of objects') ne -1, nlines)
                         
           if nlines eq 1 then begin ;got more than one object
              scnt = strmid(Result[idx], 20)
              ;;but cnt is a string, want it to be an int
              scnt = fix(scnt)
              found=1   
        
              ;;now pull the name of the first star in the list
              id2 = where(Result.StartsWith('1|') OR Result.StartsWith('1 |'), nlines)
              firstline = Result[id2]
              junk = splitsimbadstring(firstline, starname=starname, ra=ra, dec=dec, pmra=pmra, pmdec=pmdec)

              ;;apply the proper motions
              ;;assume 2013 as the 'middle' of the warm mission to get
              ;;posistions closer
              pmra = pmra*13/1000./60./60. ;arsec
              pmdec = pmdec*13/1000./60./60. ;arcsec
              ra = ra + pmra
              dec = dec + pmdec
              print, 'ra, dec', ra, dec
           endif
           ;;check for single star
           id1=where(strpos(Result, 'Object') ne -1, nlines)
           if nlines eq 1 then begin
              starname = strmid(Result[id1],7,11)
              idcoo = where(Result.StartsWith('Coordinates(ICRS'))
              cooline = Result[idcoo]
              pos = cooline.split(' ')
              ra = pos[1]
              dec = pos[2]
           endif

           
           ;;figure out how to exit the program if the fov is blank
           ;;(ie. subarray off FOV)
           
           none = where(strpos(Result, 'No astronomical') ne -1, nlines)
           if nlines eq 1 then begin
              found = 0
              starname = 'nostar'
           endif
           
              
              ;;
        endelse
        
     ENDELSE
  endelse
  
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
  
;;        starname = strmid(Result[id2[0]],13,25)
;;        ;;and remove some blank spaces
;;        starname = strtrim(starname)
;;        ;;need to keep the found ra and dec
;;        pos = strmid(Result[id2[0]],43,20)
;;        pos = strtrim(pos)
;;        minus = pos.Contains('-')
;;        if minus gt 0 then begin
;;           negdec = pos.Split('\-')
;;           ra = float(negdec[0])
;;           dec = float(negdec[1])
;;           dec = dec*(-1)
;;        endif else begin
;;           posdec = pos.Split('\+')
;;           ra = float(posdec[0])
;;           dec = float(posdec[1])
;;        endelse

;;              starname = strmid(Result[id2[0]],13,25)
;;              ;;and remove some blank spaces
;;              starname = strtrim(starname)
;;                     ;;need to keep the found ra and dec
;;              pos = strmid(Result[id2[0]],43,20)
;;              pos = strtrim(pos)
;;              print, 'pos', pos
;;              minus = pos.Contains('-')
;;              if minus gt 0 then begin
;;                 negdec = pos.Split('\-')
;;                 ra = float(negdec[0])
;;                 dec = float(negdec[1])
;;                 dec = dec*(-1)
;;              endif else begin
;;                 posdec = pos.Split('\+')
;;                 ra = float(posdec[0])
;;                 dec = float(posdec[1])
;;              endelse
