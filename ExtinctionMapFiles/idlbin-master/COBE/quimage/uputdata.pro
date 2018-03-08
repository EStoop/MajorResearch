PRO UPUTDATA,data,$
  badpixval=badpixval,$
  bandno=bandno,$
  bandwidth=bandwidth,$
  coordinate_system=coordinate_system,$
  datavar=datavar,$
  error=error,$
  faceno=faceno,$
  frequency=frequency,$
  frequnits=frequnits,$
  freqvar=freqvar,$
  instrume=instrume,$
  name=name,$
  noprompt=noprompt,$
  orient=orient,$
  projection=projection,$
  title=title,$
  weights=weights,$
  units=units
;+NAME/ONE LINE DESCRIPTION:
;    UPUTDATA sets up a UIMAGE data-object.
;
;DESCRIPTION:
;    This routine allows a UIDL user to transfer data from a main-level
;    UIDL variable into the UIMAGE environment.  Ancilliary information
;    about the data may also be passed.
;    One may also use this routine to read in a UIMAGE-type data-object
;    from an IDL saveset.
;
;CALLING SEQUENCE:
;    uputdata,data,[badpixval=...],[bandno=...],[bandwidth=...],$
;      [coordinate_system=...],[faceno=...],[frequency=...],$
;      [instrume=...],[orient=...],[projection=...],[title=...],$
;      [weights=...],[units=...],[/noprompt]
;
;ARGUMENTS (I=input, O=output, []=optional)
;    data      [I]  2-D arr  flt        A 2-D or 3-D array of data,
;                                       or the name of an IDL saveset.
;					if not specified, then the
;					user is prompted for file name
;    badpixval [I]  keyword  flt        Bad pixel value.
;    bandno    [I]  keyword  int        Band number.
;    bandwidth [I]  keyword  int        Band width.
;    coordinat [I]  keyword  str        Coordinate system (usually
;                                       'ECLIPTIC', but could also be
;                                       'EQUATORIAL' or 'GALACTIC').
;    datavar   [I]  keyword  str        Name of variable in a saveset
;                                       that contains the data.
;    error     [O]  keyword  int        0 if no errors, 1 otherwise.
;    faceno    [I]  keyword  int        Face number (for faces only).
;    frequency [I]  keyword  flt        Frequency.
;    frequnits [I]  keyword  str        Frequency units.
;    freqvar   [I]  keyword  str        Name of variable in a saveset
;                                       that contains the frequencies
;                                       for 3-D data.
;    instrume  [I]  keyword  str        Instrument name (either
;                                       'DIRBE', 'DMR', or 'FIRAS').
;    noprompt  [I]  keyword  int        Flag, if set then don't prompt
;                                       for missing keywords.
;    orient    [I]  keyword  str        Orientation ('L' or 'R').
;    projection[I]  keyword  str        Projection (should be 'SKY-CUBE'
;                                       at this level).
;    title     [I]  keyword  str        A descriptive title (preferably
;                                       less than 40 characters).
;    weights   [I]  2-D arr  flt        A 2-D or 3-D array of weights
;    units     [I]  keyword  str        A string which contains both the
;                                       name of the physical quantity
;                                       and its units.
;
;EXAMPLE:
;    uputdata,data,title='DMR: 53A Temperature',instr='DMR',$
;      proj='SKY-CUBE',coor='GALACTIC',ori='R',bad=0.,/noprompt
;#
;COMMON BLOCKS:  uimage_data,uimage_data2
;
;LIBRARY CALLS:  bold, setup_image, setup_object3d, sixunpack, validnum
;
;PROCEDURE:
;    Check if a valid data array was supplied.  Prompt for any missing
;    kewords if NOPROMPT is not set.  Call setup_image.
;
;REVISION HISTORY:
;    Written by John Ewing, ARC, January 1992.
;    
;    Modified:
;	May 1992	P Kryszak-Servin	SPR #9733
;
;       Nov 1992	J Ewing			SPR #10172  allows ints, etc
;	Jan 1993	P Kryszak-Servin	SPR #10429  3-D save sets ok
;       May 1993        J Newmark               SPR #10904  def face=-1
;       Aug 1993        D Ward                              put verbose on
;                                                           restore if var
;                                                           not found
;  SPR 11311 14 Sep 1993  Read in weights in CISS. J. Newmark
;  SPR 11758 17 May 1994  Allow input from UIDL to UIMAGE of 1-D graphs.
;                         J. Newmark
;
;.TITLE
;Routine UPUTDATA
;-
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map

  old_quiet = !QUIET
  !QUIET = 1
  error = 1
  sz = SIZE(data)
  IF(sz(1) EQ 7) THEN BEGIN   ; data is a file name, not a data array
    filename = data
  ENDIF
  IF(sz(0)+sz(1) EQ 0) THEN BEGIN	; no file name or data array 
             				; supplied so get one from user
    PRINT,'You will be asked below to identify an IDL saveset that '+$ 
      'contains a data array.'
    PRINT,' '
    dirname = ''
    READ,'Saveset directory (default = current directory):  ',dirname
    filename = ''
    READ,'Saveset file name (<CR> to exit):  ',filename
    IF(filename EQ '') THEN RETURN
    filename = STRUPCASE(dirname+filename)
  ENDIF
  sz = SIZE(filename)
  IF(sz(1) NE 0) THEN BEGIN   ; filename exists so read from it
                              ; the data array was not provided, but
                              ; will be restored from the save set
    OPENR, lun, filename, ERROR=err, /GET_LUN
    IF(err EQ 0) THEN BEGIN
      FREE_LUN, lun
      IF(NOT KEYWORD_SET(datavar)) THEN BEGIN
        datavar = ''
        READ,'Data variable name (default = "DATA"):  ',datavar
        IF(datavar EQ '') THEN datavar = 'DATA'
      ENDIF
      PRINT,' '
      j = EXECUTE( datavar+' = 0b' )
      RESTORE, filename
      j = EXECUTE('sz = SIZE('+datavar+')')
      IF(sz(0) NE 0) THEN BEGIN
        IF(STRUPCASE(datavar) NE 'DATA') THEN j = EXECUTE( 'data = '+datavar )
      ENDIF ELSE BEGIN
        IF(KEYWORD_SET(noprompt)) THEN BEGIN
          MESSAGE,'DATA was not present in the specified saveset.',/CONT
          RETURN
        ENDIF
        RESTORE, filename, /VERBOSE
        PRINT,''
        PRINT,'The specified variable was not in the saveset.  Which of the'
        PRINT,'variables listed above contains the desired data?'
        datavar = ' '
        PRINT,'' 
        READ, 'Variable name (enter "E" to exit):  ', datavar
        IF(STRUPCASE(datavar) EQ 'E') THEN RETURN
        status = EXECUTE( datavar+' = 0b' )
        PRINT,' '
        RESTORE, filename 
        j = EXECUTE('sz = SIZE('+datavar+')')
        IF(sz(0) EQ 0) THEN BEGIN
          MESSAGE,'There is no variable in the save set with that name.',/CONT
          RETURN
        ENDIF
        status = EXECUTE( 'data = '+datavar )
      ENDELSE
    ENDIF ELSE BEGIN
      MESSAGE,'That file could not be opened.',/CONT
      RETURN
    ENDELSE
  ENDIF
;
;  at this point, data is from passed variable or read from restored save set.
;  ---------------------------------------------------------------------------
  sz = SIZE(data)
; since IDL can handle numeric conversions, and the case of a string
; is already handled, then this code is not really needed.
;  type = sz(sz(0)+1)
;  IF(type EQ 7) THEN BEGIN   ; must not be a string
;    MESSAGE,'Data must not be of type string.',/CONT
;    RETURN
;  ENDIF
  ndim = sz(0)
  dim1 = sz(1)
  dim2 = sz(2)
  IF(ndim NE 3) THEN GOTO, check_dims
;
;  Handle 3-D case.
;  ----------------
  dim3 = sz(3)
  sz = SIZE(frequency)
  IF(sz(0)+sz(1) NE 0) THEN BEGIN
;
;  A value for FREQUENCY was supplied.
;  -----------------------------------
    IF((sz(0) NE 1) OR (sz(2) NE 4)) THEN BEGIN
      MESSAGE,'FREQUENCY/WAVELENGTH must be a 1-D floating point array.',/CONT
      RETURN
    ENDIF
    IF(sz(1) NE dim3) THEN BEGIN
      MESSAGE,'Dimension of frequency/wavelength array is inappropriate.',/CONT
      RETURN
    ENDIF
  ENDIF ELSE BEGIN
;
;  A value for FREQUENCY was not supplied.
;  ---------------------------------------
    sz = SIZE(filename)
    IF(sz(0)+sz(1) NE 0) THEN BEGIN
;
;  If data was read in from a saveset, then ask if there is a 1-D array in
;  the saveset which contains the frequencies, if so, then have the user
;  identify the name of that variable and set FREQUENCY to be equal to that
;  variable.
;  ------------------------------------------------------------------------
      IF(KEYWORD_SET(freqvar)) THEN GOTO, restore_freq
      PRINT,' '
      PRINT,'There was not an array named '+bold('FREQUENCY')+' in the saveset.'
      PRINT,'(A 1-D array of frequencies/wavelengths must be supplied for ' +$
         '3-D data objects).'
      menu = ['Is there a 1-D array of frequencies/wavelengths in the '+$
         'saveset?','Yes','No']
      sel = one_line_menu(menu, init=2)
      IF(sel EQ 2) THEN GOTO,prompt_freq
      freqvar = ' '
      READ, 'Please enter the frequency/wavelength variable name (or "E" to '+$
        'exit):  ', freqvar
      IF(STRUPCASE(freqvar) EQ 'E') THEN RETURN
      PRINT,' '
      PRINT,'Please ignore any "% Identifiers can only..." messages.'
restore_freq:
      status = EXECUTE( freqvar+' = 0b' )
      RESTORE, filename 
      status = EXECUTE('frequency = '+freqvar)
      sz = SIZE(frequency)
      IF(sz(0)+sz(1) EQ 0) THEN BEGIN
        MESSAGE,'There is no variable in the save set with that name.',/CONT
        RETURN
      ENDIF
      IF((sz(0) NE 1) OR (sz(2) NE 4)) THEN BEGIN
        MESSAGE,'FREQUENCY/WAVELENGTH must be a 1-D floating point array.',/CONT
        RETURN
      ENDIF
      IF(sz(1) NE dim3) THEN BEGIN
        MESSAGE,'Dimension of frequency/wavelength array is inappropriate.',$
          /CONT
        RETURN
      ENDIF
      GOTO,ok_freq
    ENDIF ELSE BEGIN
      IF(KEYWORD_SET(noprompt)) THEN BEGIN
        MESSAGE,'FREQUENCY/WAVELENGTH was not set to be a 1-D floating point'+$
         ' array.',/CONT
        RETURN
      ENDIF
      PRINT,bold('FREQUENCY')+' was not passed in.'
      PRINT,'(A 1-D array of frequencies/wavelengths must be supplied for '+$
        '3-D data objects).'
    ENDELSE
prompt_freq:
;
;  Ask if indices should be used for the frequencies.
;  --------------------------------------------------
    menu = ['Use indices as frequencies/wavelengths?','Yes','No']
    sel = one_line_menu(menu, init=2)
    IF(sel EQ 1) THEN BEGIN
      frequency = FINDGEN(dim3)
      frequnits = 'Index value'
      GOTO,ok_freq
    ENDIF
;
;  Ask if frequencies should be defined via a start and delta frequency.
;  ---------------------------------------------------------------------
    menu = ['Enter a start_freq./wave. and a delta_freq./wave.?','Yes','No']
    sel = one_line_menu(menu, init=2)
    IF(sel EQ 1) THEN BEGIN
      str = ' '
      val = 0
      WHILE(val EQ 0) DO BEGIN
        READ,bold('start_frequency/wavelength')+':  ',str
        val = validnum(str)
        IF(val EQ 1) THEN start_freq = FLOAT(str) $
                     ELSE PRINT,'Invalid number, please re-enter.'
      ENDWHILE
      val = 0
      WHILE(val EQ 0) DO BEGIN
        READ,bold('delta_frequency/wavelength')+':  ',str
        val = validnum(str)
        IF(val EQ 1) THEN delta_freq = FLOAT(str) $
                     ELSE PRINT,'Invalid number, please re-enter.'
      ENDWHILE
      frequency = FLTARR(dim3)
      frequency(0) = start_freq
      FOR i=1,dim3-1 DO frequency(i) = frequency(i-1) + delta_freq
      GOTO, ok_freq
    ENDIF
;
;  Ask if frequencies should be manually entered.
;  ----------------------------------------------
    menu = ['Manually enter frequencies/wavelengths or quit?',$
       'Enter freq./wave.','Quit']
    sel = one_line_menu(menu, init=2)
    IF(sel EQ 2) THEN BEGIN
      PRINT,'The 3-D data object was not placed in the UIMAGE data '+$
        'environment.'
      RETURN
    ENDIF
    PRINT,'Please enter the frequency/wavelength value for the indices '+$ 
        'indicated below.'
    frequency = FLTARR(dim3)
    FOR i=0,dim3-1 DO BEGIN
      str = ' '
      val = 0
      WHILE(val EQ 0) DO BEGIN
        READ,bold(STRTRIM(STRING(i+1),2))+':  ',str
        val = validnum(str)
        IF(val EQ 1) THEN frequency(i)=FLOAT(str) $
                     ELSE PRINT,'Invalid number, please re-enter.'
      ENDWHILE
    ENDFOR
ok_freq:
  ENDELSE
check_dims:
;
; check for weights
;
  IF (N_ELEMENTS(weights) NE 0) THEN BEGIN
     szwt=SIZE(weights)
     IF (szwt(0) EQ 2) THEN wgt='Y'
     IF (szwt(0) EQ 3) THEN wgt='3'
     wtitle='Weight of '+title
  ENDIF ELSE wgt=''

  IF((dim2/3)*4 EQ dim1) THEN BEGIN
;
;  MAP6, MAP7, MAP8, or MAP9
;  -------------------------
    type = 'MAP'
    CASE dim1 OF
       128 : res=6
       256 : res=7
       512 : res=8
      1024 : res=9
      ELSE : BEGIN
        MESSAGE,'Sky-cube resolution is outside of valid range [6..9].',/CONT
               RETURN
             END
    ENDCASE
    temp_data = data
    GOTO,good_data
  ENDIF
;
;  Sixpacked MAP6, MAP7, MAP8, or MAP9
;  -----------------------------------
  IF((dim2/2)*3 EQ dim1) THEN BEGIN
    type = 'MAP'
    CASE dim1 OF
        96 : res=6
       192 : res=7
       384 : res=8
       768 : res=9
      ELSE : BEGIN
        MESSAGE,'Sky-cube resolution is outside of valid range [6..9].',/CONT
               RETURN
             END
    ENDCASE

    temp_data = sixunpack(data)
    IF ((wgt EQ '3') OR (wgt EQ 'Y')) THEN weights=sixunpack(weights)
    GOTO,good_data
  ENDIF

  IF(dim1 EQ dim2) THEN BEGIN
;
;  FACE6, FACE7, FACE8, or FACE9
;  -----------------------------
    type = 'FACE'
    CASE dim1 OF
        32 : res=6
        64 : res=7
       128 : res=8
       256 : res=9
      ELSE : BEGIN
        MESSAGE,'Resolution of face is outside of valid range [6..9].',/CONT
               RETURN
             END
    ENDCASE
    temp_data = data
    GOTO,good_data
  ENDIF
  IF (ndim EQ 1) THEN BEGIN
    temp_data = data
    str=' '
    IF (NOT KEYWORD_SET(title)) THEN BEGIN
      READ,'Enter title: ',str
      title=str
    ENDIF
    IF (NOT KEYWORD_SET(frequnits)) THEN BEGIN
      READ,'Enter X axis label:',str
      frequnits=str
    ENDIF
    IF (NOT KEYWORD_SET(units)) THEN BEGIN
      READ,'Enter Y axis label:',str
      units=str
    ENDIF
    noprompt=1
    num=dim1
    nblocks_x=5
    nblocks_y=5
    sz = SIZE(frequency)
    IF(sz(0)+sz(1) NE 0) THEN BEGIN
;
;  A value for FREQUENCY was supplied.
;  -----------------------------------
      IF((sz(0) NE 1) OR (sz(2) NE 4)) THEN BEGIN
        MESSAGE,'FREQUENCY/WAVELENGTH must be a 1-D floating point array.',/CONT
        RETURN
      ENDIF
      IF(sz(1) NE dim1) THEN BEGIN
       MESSAGE,'Dimension of frequency/wavelength array is inappropriate.',/CONT
       RETURN
      ENDIF
     ENDIF ELSE BEGIN
;
;  Ask if indices should be used for the frequencies.
;  --------------------------------------------------
      menu = ['Use indices as frequencies/wavelengths?','Yes','No']
      sel = one_line_menu(menu, init=2)
      IF(sel EQ 1) THEN BEGIN
        frequency = FINDGEN(dim1)
        frequnits = 'Index value'
        GOTO,good_data
      ENDIF
;
;  Ask if frequencies should be defined via a start and delta frequency.
;  ---------------------------------------------------------------------
      menu = ['Enter a start_freq./wave. and a delta_freq./wave.?','Yes','No']
      sel = one_line_menu(menu, init=2)
      IF(sel EQ 1) THEN BEGIN
        str = ' '
        val = 0
        WHILE(val EQ 0) DO BEGIN
          READ,bold('start_frequency/wavelength')+':  ',str
          val = validnum(str)
          IF(val EQ 1) THEN start_freq = FLOAT(str) $
                     ELSE PRINT,'Invalid number, please re-enter.'
        ENDWHILE
        val = 0
        WHILE(val EQ 0) DO BEGIN
          READ,bold('delta_frequency/wavelength')+':  ',str
          val = validnum(str)
          IF(val EQ 1) THEN delta_freq = FLOAT(str) $
                     ELSE PRINT,'Invalid number, please re-enter.'
        ENDWHILE
        frequency = FLTARR(dim1)
        frequency(0) = start_freq
        FOR i=1,dim1-1 DO frequency(i) = frequency(i-1) + delta_freq
        GOTO, good_data
      ENDIF
;
;  Ask if frequencies should be manually entered.
;  ----------------------------------------------
      menu = ['Manually enter frequencies/wavelengths or quit?',$
       'Enter freq./wave.','Quit']
      sel = one_line_menu(menu, init=2)
      IF(sel EQ 2) THEN BEGIN
        PRINT,'The 1-D data object was not placed in the UIMAGE data '+$
          'environment.'
        RETURN
      ENDIF
      PRINT,'Please enter the frequency/wavelength value for the indices '+$ 
          'indicated below.'
      frequency = FLTARR(dim1)
      FOR i=0,dim1-1 DO BEGIN
        str = ' '
        val = 0
        WHILE(val EQ 0) DO BEGIN
          READ,bold(STRTRIM(STRING(i+1),2))+':  ',str
          val = validnum(str)
          IF(val EQ 1) THEN frequency(i)=FLOAT(str) $
                     ELSE PRINT,'Invalid number, please re-enter.'
        ENDWHILE
      ENDFOR
     ENDELSE
     goto,good_data
  ENDIF
  MESSAGE,'Unsupported data-array size.',/CONT
  RETURN
good_data:
;
;  Prompt for missing keywords if NOPROMPT was not set.
;  ----------------------------------------------------
  IF(NOT KEYWORD_SET(noprompt)) THEN BEGIN
    n_fields=12
    field = REPLICATE({fld, name: '', desc: '', type: '', default: '',$
      nlines: 0, capitulize: 0, dimuse: INTARR(2)}, n_fields)
    field(00)={fld,'BADPIXVAL','bad pixel value','float','0.',1,0,[1,1]}
    field(01)={fld,'BANDNO','band number','integer','0',1,0,[1,0]}
    field(02)={fld,'BANDWIDTH','band width','float','0.',1,0,[1,0]}
    field(03)={fld,'COORDINATE_SYSTEM','coordinate system','string',$
      '"ECLIPTIC"',2,1,[1,1]}
    field(04)={fld,'FACENO','cube face number','integer','-1',1,0,[1,1]}
    field(05)={fld,'FREQUENCY','frequency','float','0.',1,0,[1,0]}
    field(06)={fld,'FREQUNITS','frequency description','string',$
      '"Frequency (MHz)"',2,0,[0,1]}
    field(07)={fld,'INSTRUME','instrument name','string','"UNKNOWN"',1,1,[1,1]}
    field(08)={fld,'ORIENT','cube orientation','string','"R"',1,1,[1,1]}
    field(09)={fld,'PROJECTION','projection','string','"SKY-CUBE"',1,1,[1,1]}
    field(10)={fld,'TITLE','image title','string','"Object"',2,0,[1,1]}
    field(11)={fld,'UNITS','data description/units','string',$
      '"Pixel value (unknown units)"',2,0,[1,1]}
;
;  See whether or not values were supplied for all keywords.
;  ---------------------------------------------------------
    any_missing = 0
    FOR i=0,n_fields-1 DO BEGIN
      IF(field(i).dimuse(ndim-2) EQ 1) THEN BEGIN
        j = EXECUTE('sz = SIZE('+field(i).name+')')
        IF((type NE 'MAP') OR (field(i).name NE 'FACENO')) THEN $
          IF(sz(0)+sz(1) EQ 0) THEN any_missing = 1
      ENDIF
    ENDFOR
    IF(any_missing EQ 1) THEN BEGIN
;
;  Print out introductory text.
;  ----------------------------
      PRINT,' '
      PRINT,'You will be prompted below for certain scalar '+$
        'quantities which were'
      PRINT,'not supplied.  You may hit RETURN to use the indicated'+$
        ' default values.'
      PRINT,' '
;
;  Read in values for missing keywords.
;  ------------------------------------
      FOR i=0,n_fields-1 DO BEGIN
;        IF((ndim EQ 2) OR (field(i).iii_d EQ 1)) THEN BEGIN
        IF(field(i).dimuse(ndim-2) EQ 1) THEN BEGIN
          j = EXECUTE('sz = SIZE('+field(i).name+')')
          set = (sz(0)+sz(1)) < 1
          IF((type EQ 'MAP') AND (field(i).name EQ 'FACENO')) THEN set = 1 
          IF(set NE 1) THEN BEGIN
read_str:
            str = ''
            IF(field(i).nlines EQ 1) THEN BEGIN
              READ,bold(field(i).name)+' = '+field(i).desc+' ('+field(i).type+$
                ', default = '+field(i).default+'):  ',str
            ENDIF ELSE BEGIN
              PRINT,bold(field(i).name)+' = '+field(i).desc+' ('+field(i).type+$
                ', default = '+field(i).default+')'
              READ,':  ',str
            ENDELSE
            IF(str NE '') THEN BEGIN
              IF(field(i).type EQ 'string') THEN BEGIN
                IF(field(i).capitulize EQ 1) THEN str = STRUPCASE(str)
                first_char = STRMID(str,0,1)
                IF((first_char EQ '"') OR (first_char EQ "'")) THEN BEGIN
                  str = STRMID(str,1,STRLEN(str)-1)
                  last_char = STRMID(str,STRLEN(str)-1,1)
                  IF((last_char EQ '"') OR (first_char EQ "'")) THEN $
                    str = STRMID(str,0,STRLEN(str)-1)
                ENDIF
                j = EXECUTE(field(i).name+' = str')
              ENDIF ELSE BEGIN
                IF(validnum(str) NE 1) THEN BEGIN
                  PRINT,'Invalid numerical value.  Please re-enter.'
                  GOTO,read_str
                ENDIF
                IF(field(i).type EQ 'integer') THEN val = FIX(str) $
                                               ELSE val = FLOAT(str)
                j = EXECUTE(field(i).name+' = val')
              ENDELSE
            ENDIF ELSE j = EXECUTE(field(i).name+' = '+field(i).default)
          ENDIF ELSE IF(field(i).capitulize EQ 1) THEN $
             j=EXECUTE(field(i).name +'= STRUPCASE('+field(i).name+')')
        ENDIF
      ENDFOR
    ENDIF
  ENDIF

  CASE ndim OF
     3:	BEGIN
           wlink = -1
	   if wgt eq 'Y' then begin
               PRINT,'Object has a 2D weight array associated with it.'
               scale_min=MIN(weights(WHERE(weights NE 0)),MAX=scale_max)
	       wtitle=append_number(wtitle)
               name = setup_image( $
                 data=weights, $
		 hidden=1, $
                 badpixval=0, $
                 bandno = bandno, $
                 frequency=freqency, $
                 bandwidth=bandwidth, $
                 units=units, $
                 faceno = faceno, $
                 title=wtitle, $
                 instrume=instrume, $
                 version=version, $
                 orient = orient, $
                 scale_min=scale_min,$
                 scale_max=scale_max,$
                 projection=projection, $
                 coordinate_system=coordinate_system )
		 j = execute('wlink = '+name+'.window')
	   endif

	   if wgt eq '3' then begin
             PRINT,'Object has a 3D weight array associated with it.'
             title=append_number(title)
             name = setup_object3d( $
               data=temp_data, $
               badpixval=badpixval, $
               frequency=frequency, $
	       frequnits=frequnits, $
               units=units, $
               faceno = faceno, $
               title=title, $
	       linkweight=wlink, $
	       weight3d=weights, $	       
               instrume=instrume, $
               version=version, $
               orient = orient, $
               projection=projection, $
               coordinate_system=coordinate_system)
	   endif

	   if wgt ne '3' then begin
             title=append_number(title)
             name = setup_object3d( $
               data=temp_data, $
               badpixval=badpixval, $
               frequency=frequency, $
	       frequnits=frequnits, $
               units=units, $
               faceno = faceno, $
               title=title, $
	       linkweight=wlink, $
               instrume=instrume, $
               version=version, $
               orient = orient, $
               projection=projection, $
               coordinate_system=coordinate_system)
	   Endif
        END
      2: BEGIN    ; handle 2D images
	     wlink = -1
	     if wgt Eq 'Y' then begin
               PRINT,'Object has a 2D weight array associated with it.'
	       wtitle=append_number(wtitle)
               scale_min=MIN(weights(WHERE(weights NE 0)),MAX=scale_max)
               name = setup_image( $
                 data=weights, $
		 hidden=1, $
                 badpixval=0, $
                 bandno = bandno, $
                 frequency=freqency, $
                 bandwidth=bandwidth, $
                 units=units, $
                 faceno = faceno, $
                 title=wtitle, $
                 instrume=instrume, $
                 version=version, $
                 orient = orient, $
                 scale_min=scale_min,$
                 scale_max=scale_max,$
                 projection=projection, $
                 coordinate_system=coordinate_system)
		 j = execute('wlink = '+name+'.window')
	      endif

              title=append_number(title)
              scale_min=MIN(temp_data,MAX=scale_max)
              name = setup_image( $
                data=temp_data, $
                badpixval=0, $
                bandno = bandno, $
                frequency=freqency, $
                bandwidth=bandwidth, $
                units=units, $
                faceno = faceno, $
                title=title, $
		linkweight=wlink, $
                instrume=instrume, $
                version=version, $
                orient = orient, $
                scale_min=scale_min,$
                scale_max=scale_max,$
                projection=projection, $
                coordinate_system=coordinate_system)
         END
        1:BEGIN
              name = setup_graph(frequency, temp_Data, dim1,$
                title = title, x_title = frequnits, y_title = units, $
                nblocks_x = nblocks_x, nblocks_y = nblocks_y)
              j=EXECUTE('winnum =' + name + '.window')
              IF ((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN WSET,winnum
              plot_graph, name
          END
    ELSE:
  ENDCASE
  !QUIET = old_quiet
  IF(name NE 'ERROR') THEN error = 0
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


