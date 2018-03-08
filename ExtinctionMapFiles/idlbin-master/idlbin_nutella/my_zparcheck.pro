;******************************************************************************
;
;MY_ZPARCHECK checks user parameters to a procedure
;
;       progname  - scalar string name of calling procedure
;       parameter - parameter passed to the routine
;       parnum    - integer parameter number
;       types     - integer scalar or vector of valid types
;                1 - byte        2 - integer   3 - int*4
;                4 - real*4      5 - real*8    6 - complex
;                7 - string      8 - structure 9 - double complex
;       dimens   - integer scalar or vector giving number
;                     of allowed dimensions.
;       message - string message describing the parameter to be printed if an 
;               error is found
;
;******************************************************************************

PRO MY_ZPARCHECK,progname,parameter,parnum,types,dimens,message

IF N_PARAMS() LT 4 THEN BEGIN
	
	PRINT,'Syntax-MY_ZPARCHECK(progname,parameter,parnum,types,dimens,[message])
	
	RETURN

ENDIF

;get type and size of parameter

s=SIZE(parameter)

ndim=s[0]

type=s[ndim+1]

;check if parameter defined.

IF type EQ 0 THEN BEGIN

        err = ' is undefined.'

        GOTO, ABORT 

ENDIF

;check for valid dimensions

valid=WHERE(ndim EQ dimens,Nvalid)

IF Nvalid LT 1 THEN BEGIN

        err='has wrong number of dimensions'

        GOTO,ABORT   

ENDIF

;check for valid type

valid=WHERE(type EQ types,Ngood)

IF ngood LT 1 THEN BEGIN

        err='is an invalid data type'

        GOTO, ABORT   

ENDIF

RETURN

;bad parameter

ABORT:
  mess=' '

IF N_PARAMS() LT 6 THEN message=''

IF message NE '' THEN mess='('+message+')'

PRINT,STRING(7b)+'Parameter'+STRTRIM(parnum,2)+mess,'of routine',STRUPCASE(progname)+' ',erR

sdim=' '

FOR i=0,N_ELEMENTS(dimens)-1 DO BEGIN

	IF dimens[i] EQ 0 THEN sdim=sdim+'scalar' ELSE sdim=sdim+STRING(dimens[i],'(i3)')

END

PRINT,'Valid dimensions are:'+sdim

stype=' '
  
FOR i=0, N_ELEMENTS(types)-1 DO BEGIN
        CASE types[i] OF
                1: stype = stype + ' byte'
                2: stype = stype + ' integer'
                3: stype = stype + ' longword'
                4: stype = stype + ' real*4'
                5: stype = stype + ' real*8'
                6: stype = stype + ' complex'
                7: stype = stype + ' string'
                8: stype = stype + ' structure'
                9: stype = stype + ' dcomplex'
        ENDCASE
ENDFOR
  
PRINT,'Valid types are:' + stype

RETALL  ; my_zparcheck

END
