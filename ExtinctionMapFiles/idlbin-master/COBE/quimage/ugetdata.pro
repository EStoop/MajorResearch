PRO ugetdata,data
;+NAME/ONE LINE DESCRIPTION:
;    UGETDATA copies data from the UIMAGE environment to a UIDL variable.
;
;DESCRIPTION:
;    This routine allows a UIDL user to transfer data from the UIMAGE
;    environment into a main-level variable.
;
;CALLING SEQUENCE:
;    ugetdata,data
;
;ARGUMENTS (I=input, O=output, []=optional)
;    data       O   2-D arr  flt        A variable which will receive
;                                       an array of data.
;
;EXAMPLE:
;    ugetdata,data
;#
;COMMON BLOCKS:  uimage_data, uimage_data2, info_3d
;
;LIBRARY CALLS:  select_object
;
;PROCEDURE:
;    The SELECT_OBJECT function is called to allow the user to select
;    the desired data object.  If a selection is made, then the data
;    associated with that object is returned through the variable DATA.
;    If there were no data objects stored in the common statement, then
;    a message is printed out and DATA will be undefined.
;
;REVISION HISTORY:
;    Written by John Ewing, ARC, January 1992.
;  SPR 10829  Apr  1993  Change info_3d common block.   J. Newmark
;  SPR 11100  Jun 28 93  Enhance choices for 3d objects. J. Newmark
;  SPR 11169  Jul 21 93  Work with zoomed output for bckgrnd fits. J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;
;.TITLE
;Routine UGETDATA
;-
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
   freq3D_2,wt3d_0,wt3d_1,wt3d_2
  COMMON zback,zbgr,zdsrcmap,zbsub
  IF(N_PARAMS(0) NE 1) THEN $
    MESSAGE,'Please designate an output variable on the line of invocation.'
  name=select_object(/map6,/map7,/map8,/map9,/face6,/face7,/face8,$
    /face9,/proj_map,/proj2_map,/graph,/zoomed,/object3d,title='Select object')
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE,'There is no UIMAGE-type data to choose from.',/CONT
    RETURN
  ENDIF
;
;  Handle the different types of data-objects (MAP, FACE, GRAPH, & ZOOMED).
;  ------------------------------------------------------------------------
  first_letter = STRMID(name,0,1)
  CASE first_letter OF
    'M' : j = EXECUTE('data='+name+'.data')
    'F' : j = EXECUTE('data='+name+'.data')
    'P' : j = EXECUTE('data='+name+'.data')
    'G' : BEGIN
            j = EXECUTE('data='+name+'.data')
            j = EXECUTE('num = '+name+'.num')
            j = EXECUTE('scatter = '+name+'.scatter')
            IF(scatter EQ 1) THEN BEGIN
              j = EXECUTE('win_orig1 = '+name+'.win_orig1')
              j = EXECUTE('win_orig2 = '+name+'.win_orig2')
              name1 = get_name(win_orig1)
              name2 = get_name(win_orig2)
              j = EXECUTE('badpixval1 = '+name1+'.badpixval')
              j = EXECUTE('badpixval2 = '+name2+'.badpixval')
              j = EXECUTE('good = WHERE(('+name1+'.data NE badpixval1) AND ('+$
                name2+'.data NE badpixval2))')
              IF(good(0) EQ -1) THEN BEGIN
                MESSAGE,'The two images have no overlapping good pixels.',/CONT
                RETURN
              ENDIF
              j = EXECUTE('data = [['+name1+'.data(good)],['+name2+$
                '.data(good)]]')
            ENDIF ELSE BEGIN
              data = data(0:num-1,*)
            ENDELSE
            PRINT,'A 2-D array is being returned which contains the data'+$
              ' that was used'
            PRINT,'to make the selected graph.  Please index the array'+$
              ' by "(*,0)" to'
            PRINT,'access the X data, and index it by "(*,1)" to access'+$
              ' the Y data.'
          END
    'Z' : BEGIN
            j = EXECUTE('win_orig = '+name+'.win_orig')
            j = EXECUTE('start_x = '+name+'.start_x')
            j = EXECUTE('stop_x = '+name+'.stop_x')
            j = EXECUTE('start_y = '+name+'.start_y')
            j = EXECUTE('stop_y = '+name+'.stop_y')
            start_x = start_x > 0
            start_y = start_y > 0
            bflag=0
            IF (win_orig LT 0) THEN BEGIN
              win_orig=ABS(win_orig)
              bflag=1
            ENDIF
            j=EXECUTE('zoomflag = ' + name + '.zoomflag')
            IF (zoomflag NE 0) THEN bflag=1
            name_orig=get_name(win_orig)
            j = EXECUTE('sz = SIZE('+name_orig+'.data)')
            dim1 = sz(1)
            dim2 = sz(2)
            stop_x = stop_x < (dim1-1)
            stop_y = stop_y < (dim2-1)
            j = EXECUTE('data = '+name_orig+$
              '.data(start_x:stop_x,start_y:stop_y)')
            IF (bflag EQ 1) THEN BEGIN
              IF (zoomflag EQ 1) THEN data=zbgr
              IF (zoomflag EQ 2) THEN data=zdsrcmap
              IF (zoomflag EQ 3) THEN data=zbsub
            ENDIF
          END
    'O' : BEGIN
            index3d = FIX(STRMID(name,9,1))
            str = STRING(index3d,'(i1)')
            menu3=['Select appropriate array','data','freq/wave','weights']
            sel=one_line_menu(menu3,init=1)
            CASE sel OF
              1: j = EXECUTE('data = data3D_'+str)
              2: j = EXECUTE('data = freq3D_'+str)
              3: j = EXECUTE('data = wt3D_'+str)
           ENDCASE
          END
    ELSE:
  ENDCASE
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


