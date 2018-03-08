function vtdef_area,name
; SPR 12060  19-Jan-1995 Variable name change due to IDL3.6. J. Newmark
COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
COMMON uimage_data2,proj_map,proj2_map

   fpnt=FLTARR(1000,2)
   n=0
   forever=1
;
;  Have the user type in the {lon,lat} coordinates of the desired point.
;  ---------------------------------------------------------------------
   PRINT,' '
   PRINT,'You will be prompted for points that define the vertices of an area.'
   PRINT,'You need to define at least 3 points. Enter the Longitude & Latitude'
   PRINT,'in degrees followed by a <CR>. After entering all the points, hit a <CR>.'  
   PRINT,' '
   PRINT,'Please identify the coordinate system you will use.'
   menu = ['Coordinate system:','Ecliptic','Equatorial','Galactic']
   sel = one_line_menu(menu)
   icoord = STRMID('EQG',sel-1,1)
   PRINT,'Next enter the lon/lat pair in degrees (example: '+$
        bold('20. -30.')+').  [<CR> to exit]'
   WHILE(forever EQ 1) DO BEGIN
      pnt = ''
      READ,'lon lat:  ',pnt
      pnt = STRTRIM(pnt,2)
      len = STRLEN(pnt)
      IF(len EQ 0) THEN BEGIN
	 IF(n LT 3) THEN PRINT,'The number of vertices should be >= 3' ELSE forever=0
      ENDIF ELSE BEGIN
         delm = STRPOS(pnt,',')
	 IF(delm EQ -1) THEN delm = STRPOS(pnt,' ')
	 IF(delm EQ -1) THEN BEGIN
            PRINT,'Please use comma or space as a delimiter'
	 ENDIF ELSE BEGIN
	    str1=STRMID(pnt,0,delm)
	    str1=STRTRIM(str1,2)
	    val=validnum(str1)
	    IF(val NE 1) THEN BEGIN
	       PRINT,'The longitude is an invalid number.'
	    ENDIF ELSE BEGIN
	       fpnt(n,0)=FLOAT(str1)
	       IF((fpnt(n,0) LT 0) OR (fpnt(n,0) GT 360)) THEN BEGIN
	          PRINT,'The longitude should be between 0 and 360 degrees'
	       ENDIF ELSE BEGIN
		  str2=STRMID(pnt,delm+1,len-delm)
		  str2=STRTRIM(str2,2)
		  val=validnum(str2)
	          IF(val NE 1) THEN BEGIN
		      PRINT,'The latitude is an invalid number.'
		  ENDIF ELSE BEGIN
		     fpnt(n,1)=FLOAT(str2)
	       	     IF((fpnt(n,1) LT -90) OR (fpnt(n,1) GT 90)) THEN BEGIN
		        PRINT,'The latitude should be between -90 and +90 degrees'
	             ENDIF ELSE BEGIN
	                n=n+1
	             ENDELSE
		  ENDELSE
	       ENDELSE
	    ENDELSE
	 ENDELSE
      ENDELSE
   ENDWHILE
   PRINT,' '
   first_letter = STRMID(name,0,1)
   istatus = EXECUTE('data='+ name + '.data')
   xpnt=FLTARR(n,2)
   pixcoor=FLTARR(n,2)
   xpnt=fpnt(0:n-1,0:1)
   x_out=INTARR(n)
   y_out=INTARR(n)
   CASE first_letter OF
     'M' :  BEGIN 
               res=STRMID(name,3,1)
      	       ocoord='r'+res
      	       pixcoor=coorconv(xpnt,infmt='L',outfmt='P',$
				inco=icoord ,outco=ocoord)
      	       pix2xy,pixcoor,x_out,y_out,res=res
	    END
     'P' :  BEGIN
	       ocoord=icoord
               pixcoor=coorconv(xpnt,infmt='L',outfmt='U',$
				inco=icoord ,outco=ocoord)
	       xyout=INTARR(n,2)
               j = EXECUTE('projection = '+name+'.projection')
	       pr=STRMID(projection,0,1)
               CASE projection OF
                  'AITOFF'            : proj = 'A'
           	  'GLOBAL SINUSOIDAL' : proj = 'S'
           	  'MOLLWEIDE'         : proj = 'M'
            	   ELSE               : PRINT,'Unsupported projection.'
                               
               ENDCASE
               xyout=uv2proj(pixcoor,pr,SIZE(data))
	       x_out=xyout(0:n-1,0)
	       y_out=xyout(0:n-1,1)
	    END
     'F' :  BEGIN
	       res=STRMID(name,4,1)
               ocoord='r'+res
               pixcoor=coorconv(xpnt,infmt='L',outfmt='P',$
				inco=icoord ,outco=ocoord)
               j = EXECUTE('faceno = '+name+'.faceno')
               pix2xy,pixcoor,x_out,y_out,res=res,/face
	    END
      ELSE : print,'VTDEF-AREA,This object cannot be handled by this program-',$
		name
   ENDCASE
   sz=SIZE(data)
   vimage=data(POLYFILLV(x_out(0:n-1),y_out(0:n-1),sz(1),sz(2)))
   RETURN,vimage
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


