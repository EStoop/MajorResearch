pro getpro,proc_name
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE
;     GETPRO copies IDL procedure into current default directory.
;
;DESCRIPTION:
;    IDL procedure to extract a procedure from an IDL Library or
;    directory given in the !PATH system variable and place it in the
;    current default directory (presumably to be edited by the user).
;
;CALLING SEQUENCE:
;     GETPRO,PROC_NAME
;
;ARGUMENTS (I = input, O = output, [] = optional):
;     PROC_NAME     I   str        Character string giving the name
;                                  of the IDL procedure or function.
;                                  Do not give an extension.
;
;WARNINGS:
;    User will be unable to obain source code for a native IDL
;    function or procedure, or for a FORTRAN or C routine added with
;    CALL_EXTERNAL.  User must have write privilege to the current
;    directory.
;
;EXAMPLE:
;     To put a copy of the USER library procedure CURVEFIT on the
;       default directory:
;
;       GETPRO,'CURVEFIT'
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    The system variable !PATH is parsed into individual libraries or 
;    directories.   Each library or directory is then searched for the
;    procedure name.   When found, a SPAWN is used to extract or copy
;    the procedure into the user's directory.    If not found in
;    !PATH, then the ROUTINES.HELP file is checked to see if it is an
;    intrinsic IDL procedure.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    zparcheck,fdecomp,getenv,gettok,libname,message,cd
;
;MODIFICATION HISTORY:
;    Written W. Landsman, STX Corp.   June 1990
;    Converted to SUN IDL,  M. Greason, STX, August 1990.
;    Revised to work on both UNIX and VMS, N. Collins, STX, Nov.1990
;    Display the directory where the procedure was found, W. Landsman,
;    Dec 1990
;
; SPR 9616
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;.TITLE
;ROUTINE GETPRO
;-
on_error,2                           ;Return to caller on error
os = !cgis_os                     ;VMS or Unix operating system
if n_params() EQ 0 then begin 	     ;Prompt for procedure name?
	proc_name = ' ' 
      	read,'Enter name of procedure you want a copy of: ',proc_name     
endif else zparcheck,'getpro',proc_name,1,7,0,'Procedure name'
;
fdecomp,proc_name,disk,dir,name      ;Don't want file extensions
name = strtrim(name,2)  
;                                    ;Make sure user has write privileges
openw,lun,name+'.pro',/DELETE,/GET_LUN,ERROR=ERR
if err NE 0 then begin               ;Problem writing a temporary file?
     cd,current=curdir   
     print,!ERR_STRING
     message,curdir + ' has insufficient privilege or file protection violation'
endif
free_lun,lun
;
if os EQ "vms" then begin   
     sep = ','
     dirsep = ''
     copy = 'COPY '
     name = strupcase(name)
endif 
IF (os EQ 'unix') then begin
     sep = ':'
     dirsep = '/'
     copy = 'cp '
endif   
IF (os EQ 'windows') then begin
     sep = ':'
     dirsep = '\'
     copy = 'cp '
endif   
;
if strlowcase(name) EQ 'startup' then begin     
       TEST = getenv("IDL_STARTUP")
       if test NE "" then begin         ;Copy the STARTUP file
           spawn,copy + test + '.pro startup.pro' 
           message,'Created file startup.pro on the current directory',/INF
       ENDIF ELSE message,'The environment IDL_STARTUP has not been defined'
       return
endif
;
temp = !PATH                     ;Get current IDL path of directories
IF os eq "vms" then temp = strupcase(temp)
;
;    Loop over each directory in !PATH until procedure name found
;
while temp NE '' do begin   
   dir = gettok(temp,sep)
;
   if strmid(dir,0,1) EQ '@' then begin          ;Text Library?
     if os ne "vms" then message, $
        '!path contains a invalid VMS directory specification',/cont $
     else begin
        libname = strmid( dir,1,strlen(dir)-1 )         ;Remove the "@" symbol
        spawn,'library/extract='+name+'/out='+name+'.pro '+ $
                               libname,out,count=i
        if i EQ 0 then begin                           ;Success?
            message,name + '.PRO extracted from ' + libname,/INF
            return
        endif
      endelse
;
   endif else begin                              ;Directory
        a = findfile(dir + dirsep + name+'.pro',COUNT=i)
        if I GE 1 then begin                     ;Found by FINDFILE?
;
          if os ne "vms" then spawn,'cp ' +a(0) + ' .' else begin
             semi = strpos(a(0),';')             ;No version number
             if semi gt 0 then fname = strpos(a(0),0,semi)
             spawn,'copy '+disk+dir+name+ '.pro' + ' *'
          endelse
          message,'Procedure '+ NAME + ' copied from directory '+ dir,/INF
          return
        endif
    endelse
endwhile
;
; At this point !PATH has been searched and the procedure has not been found
; Now check if it is an intrinsic IDL procedure or function
;
openr,inunit,filepath('routines.help',SUBDIR='help'),/GET_LUN    ;Open help files
readf,inunit,n                                           ;# OF RECORDS
lv2_topics = strarr(n)                                  
readf,inunit,lv2_topics
lv2_topics = strtrim(strmid(lv2_topics,0,15),2)
test = where(lv2_topics EQ strupcase(name),count)
if count EQ 0 then begin   
     message,'Procedure '+NAME+' not found in the !PATH search string',/CONT
     message,'Check your spelling or search the individual directories',/INF
endif else begin     
     message,'Procedure '+NAME+' is an intrinsic IDL procedure',/CONT
     message,'No source code is available',/INF
endelse
return
;
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


