PRO ubefit, name
;+
;  UBEFIT - a UIMAGE-specific routine.  This is UIMAGE's driver routine
;  for the FITMU procedure. This fits spectra to a Bose Einstein
;  spectrum plus a residual blackbody spectrum.
;#
;  Written by J Newmark
;  SPR 10880  29 Apr 93 Creation.
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
    freq3D_2,wt3d_0,wt3d_1,wt3d_2
  COMMON journal,journal_on,luj
  COMMON xscreen,magnify,block_usage,scrsiz
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  menu_title = 'Fit Spectra with a Bose-Ein. plus BB model'
  IF(N_PARAMS(0) EQ 1) THEN GOTO, have_name
;
;  Put up a menu of potential operands.
;  ------------------------------------
show_menu:
  name = select_object(/graph,/object3d,/help,/exit,title=menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no appropriate objects available.', /CONT
    RETURN
  ENDIF
;
;  An operand has been selected - make sure it is valid.
;  -----------------------------------------------------
have_name:
  first_letter=STRMID(name,0,1)
  parent = name
  j = EXECUTE('title = ' + name + '.title')
  j = EXECUTE('badpixval = ' + name + '.badpixval')
  j=EXECUTE('linkweight= ' +name+ '.linkweight')
  wsel=0
;
; Handle Case where selected object is a 3-dim one.
;
 IF (first_letter EQ 'O') THEN BEGIN
   str_index3d=STRMID(name,9,1)
   index3d=FIX(str_index3d)
   dim1=object3d(index3d).dim1
   dim2=object3d(index3d).dim2
   dim3=object3d(index3d).dim3
   j=EXECUTE('freq = freq3d_' +str_index3d)
   j=EXECUTE('inspec = data3d_' + str_index3d)
   j=EXECUTE('havewt= wt3d_' + str_index3d)
   sz=SIZE(havewt)
   IF (sz(0) EQ 3) THEN j=EXECUTE('inweight=wt3d_'+str_index3d) $
   ELSE BEGIN
    PRINT,'There are no weights associated with the chosen object.'
    READ,'Enter an estimate for the Sigma: ',sig
    inweight=inspec*0. + 1./sig^2
   ENDELSE
   bad=WHERE(inspec EQ badpixval)
   inweight(bad)=0.
;
; Now need to extract out only those spectra which are nonzero.
; Place these in a 2D array for ease of use within FITMU.
;
  IF((dim2 / 3) * 4 EQ dim1) THEN BEGIN
    type = 'MAP'
    CASE dim1 OF
       128 : res = 6
       256 : res = 7
       512 : res = 8
      1024 : res = 9
      ELSE : 
    ENDCASE
  ENDIF
  IF(dim1 EQ dim2) THEN BEGIN
    type = 'FACE'
    faceno = object3d(index3d).faceno
    CASE dim1 OF
        32 : res = 6
        64 : res = 7
       128 : res = 8
       256 : res = 9
      ELSE : 
    ENDCASE
  ENDIF
  IF (type EQ 'FACE') THEN BEGIN
    tempix=FINDGEN(FLOAT(dim1)*dim2)
    tempix=LONG(tempix)
    pix2xy,tempix,data=tempix,res=res,raster=outpix,/face
  ENDIF ELSE BEGIN
   tempix=FINDGEN(FLOAT(dim1)*dim2/2)
   tempix=LONG(tempix)
   pix2xy,tempix,data=tempix,res=res,raster=outpix
  ENDELSE
  pixlist=outpix(WHERE(TOTAL(inspec,3) NE badpixval))
  spec=pix2dat(pixel=pixlist,raster=inspec)
  weight=pix2dat(pixel=pixlist,raster=inweight)
;
; Handle case where selected object is a 1-dim graph
;
 ENDIF ELSE BEGIN
  j=EXECUTE('data = ' +name + '.data')
  j=EXECUTE('num = ' +name+ '.num')
  freq=data(0:num-1,0)
  spec=data(0:num-1,1)
  IF( linkweight NE -1) THEN BEGIN
    wsel=1
    namew=get_name(linkweight)
    j = EXECUTE('titlew = ' + namew + '.title')
    PRINT, 'Weights will be taken from ' + bold(titlew)
    j = EXECUTE('badvalw = ' + namew + '.badpixval')
    j = EXECUTE('weight = ' + namew + '.data(0:num-1,1)')
    badw = WHERE(weight EQ badvalw)
    IF(badw(0) NE -1) THEN weight(badw) = 0.
   ENDIF ELSE BEGIN
    PRINT,'There are no weights associated with the chosen object.'
    READ,'Enter an estimate for the Sigma: ',sig
    weight=spec*0. + 1./sig^2
   ENDELSE
   FOR i=0,num-1 DO IF (spec(i) EQ badpixval) THEN weight(i)=0.
 ENDELSE
;        
; Define various input parameters.
; Choose any input arrays used to define fixed parameter values for
; individual spectra.
;
pix_parm=BYTARR(9)
funits=''
get_parms:
IF (first_letter EQ 'O') THEN BEGIN
parmenu = ['Define Input Parameters ',$
          'BE Temperature value array', $
          'BE Chemical Pot. value array', $
          'BB Temperature value array', $
          'Optical Depth value array', $
          'Emissivity value array', $
          'Mask (def=no masking)', $
          'Tolerance (def=1e-4)', $
          'Input Frequency Units (def=icm)',$
          'Input/Output Spectral units (def=eplees)', $
          'Abort','HELP', $
          'Done Defining Inputs']
inparm = UMENU(parmenu, title=0, init=inparm)
CASE inparm OF
  1: BEGIN
    title='Choose object array for BE temp value'
    cmb=select_object(/map6,/map7,/map8,/map9,/face6,/face7,/face8,/face9,$
       /exit,title=title)
    IF(cmb EQ 'EXIT') THEN GOTO, get_parms
    j=EXECUTE('incmb = ' + cmb + '.data')
    cmb_temp=pix2dat(pixel=pixlist,raster=incmb)
    pix_parm(0)=1
    GOTO, get_parms
    END
  2: BEGIN
    title='Choose object array for BE chem pot value'
    cmu=select_object(/map6,/map7,/map8,/map9,/face6,/face7,/face8,/face9,$
       /exit,title=title)
    IF(cmu EQ 'EXIT') THEN GOTO, get_parms
    j=EXECUTE('incmu = ' + cmu + '.data')
    cmb_mu=pix2dat(pixel=pixlist,raster=incmu)
    pix_parm(1)=1
    GOTO, get_parms
    END
  3: BEGIN
    title='Choose object array for BB temp value'
    bb1=select_object(/map6,/map7,/map8,/map9,/face6,/face7,/face8,/face9,$
       /exit,title=title)
    IF(bb1 EQ 'EXIT') THEN GOTO, get_parms
    j=EXECUTE('inbb1 = ' + bb1 + '.data')
    bb_temp=pix2dat(pixel=pixlist,raster=inbb1)
    pix_parm(2)=1
    GOTO, get_parms
    END
  4: BEGIN
    title='Choose object array for OPDEP value'
    opd1=select_object(/map6,/map7,/map8,/map9,/face6,/face7,/face8,/face9,$
       /exit,title=title)
    IF(opd1 EQ 'EXIT') THEN GOTO, get_parms
    j=EXECUTE('inopd = ' + opd1 + '.data')
    opdep=pix2dat(pixel=pixlist,raster=inopd)
    pix_parm(3)=1
    GOTO, get_parms
    END
  5: BEGIN
    title='Choose object array for EMISSIVITY value'
    em1=select_object(/map6,/map7,/map8,/map9,/face6,/face7,/face8,/face9,$
       /exit,title=title)
    IF(em1 EQ 'EXIT') THEN GOTO, get_parms
    j=EXECUTE('inem = ' + em1 + '.data')
    emis=pix2dat(pixel=pixlist,raster=inem)
    pix_parm(4)=1
    GOTO, get_parms
    END
  6: BEGIN
    index=''
    IF (first_letter EQ 'O') THEN stval=dim3 else stval=num
    PRINT,'The range of the frequency/wavelength array is 1 -',FIX(stval)
    READ,'Enter the masked element array values',index
    mask,index,weight,mask
    weight=mask
    pix_parm(5)=1
    GOTO, get_parms
    END
  7: BEGIN
    READ,'Enter Tolerance (Default=1e-4) ',tolerance
    pix_parm(6)=1
    GOTO, get_parms
    END
  8: BEGIN
    fmenu=['Choose Frequency Units (Default=icm)','icm','microns','GHz']
    fchoice=one_line_menu(fmenu,init=1)
    CASE fchoice OF
     1:funits="'icm'"
     2:funits="'microns'"
     3:funits="'GHz'"
    ENDCASE
    pix_parm(7)=1
    GOTO, get_parms
    END
  9: BEGIN
    sunit=''
    smenu=['Choose Spectral Units (Default=eplees)','eplees','wcms','mjy']
    schoice=one_line_menu(smenu,init=1)
    CASE schoice OF
     1:sunit='/eplees'
     2:sunit='/wcms'
     3:sunit='/mjy'
    ENDCASE
    pix_parm(8)=1
    GOTO, get_parms
    END
  10: GOTO, show_menu
  11: BEGIN
      uimage_help, menu_title
      GOTO, get_parms
    END
  12: PRINT,'Initial parameter arrays are defined, proceeding..'
 ENDCASE
ENDIF ELSE BEGIN
parmenu = ['Define Input Parameters ',$
          'Mask (def=no masking)', $
          'Tolerance (def=1e-4)', $
          'Input Frequency Units (def=icm)',$
          'Input/Output Spectral units (def=eplees)', $
          'Abort','HELP', $
          'Done Defining Inputs']
inparm = UMENU(parmenu, title=0, init=inparm)
CASE inparm OF
  1: BEGIN
    index=''
    IF (first_letter EQ 'O') THEN stval=dim3 else stval=num
    PRINT,'The range of the frequency/wavelength array is 1 -',FIX(stval)
    READ,'Enter the masked element array values',index
    mask,index,weight,mask
    weight=mask
    pix_parm(5)=1
    GOTO, get_parms
    END
  2: BEGIN
    READ,'Enter Tolerance (Default=1e-4) ',tolerance
    pix_parm(6)=1
    GOTO, get_parms
    END
  3: BEGIN
    fmenu=['Choose Frequency Units (Default=icm)','icm','microns','GHz']
    fchoice=one_line_menu(fmenu,init=1)
    CASE fchoice OF
     1:funits="'icm'"
     2:funits="'microns'"
     3:funits="'GHz'"
    ENDCASE
    pix_parm(7)=1
    GOTO, get_parms
    END
  4: BEGIN
    sunit=''
    smenu=['Choose Spectral Units (Default=eplees)','eplees','wcms','mjy']
    schoice=one_line_menu(smenu,init=1)
    CASE schoice OF
     1:sunit='/eplees'
     2:sunit='/wcms'
     3:sunit='/mjy'
    ENDCASE
    pix_parm(8)=1
    GOTO, get_parms
    END
  5: GOTO, show_menu
  6: BEGIN
      uimage_help, menu_title
      GOTO, get_parms
    END
  7: PRINT,'Initial parameters are defined, proceeding..'
 ENDCASE
ENDELSE
;
; Define initial parameter values and parameter variability
;
PRINT,'Set singular initial parameters:'
inp = ''
IF ((pix_parm(7) EQ 1) AND (funits EQ "'microns'")) THEN BEGIN
  dflt_str = ['<2.7>','<0>','<20>','<2e-7>','<-1.5>']  
  default = [2.7,0,20,2e-7,-1.5]
ENDIF ELSE BEGIN
  dflt_str = ['<2.7>','<0>','<20>','<2e-7>','<1.5>']
  default = [2.7,0,20,2e-7,1.5]
ENDELSE
parm_title=['Bose Einstein Temperature ', $
	    'Bose Einstein Chemical Potential ', $
            'Residual Black Body Temperature ', $
            'Residual Black Body Optical Depth ', $
            'Residual Black Body Emissivity Spectral Index']
parm0 = FLTARR(5)
FOR i=0,4 DO BEGIN
   IF (pix_parm(i) EQ 1) THEN GOTO,lbl1
   IF ((i EQ 3 OR i EQ 4) AND (parm0(2) EQ 0) $
   AND (pix_parm(2) EQ 0)) $
   THEN GOTO,lbl1
   READ,parm_title(i)+dflt_str(i)+': ',inp
   IF (inp EQ '') THEN parm0(i) = default(i) $
   ELSE parm0(i) = FLOAT(inp)
 lbl1:
ENDFOR
active = BYTARR(5)
FOR i=0,4 DO BEGIN
    IF (pix_parm(i) EQ 1) THEN GOTO,lbl2
    IF ((i GE 2 AND i LE 4) AND (parm0(2) EQ 0) $
    AND (pix_parm(1) EQ 0)) $
    THEN GOTO,lbl2
    READ,'Vary ' + parm_title(i) + '<Yes>? ',inp
    IF (inp EQ '' OR STRUPCASE(STRMID(inp,0,1)) EQ 'Y') $
    THEN active(i) = 1
  lbl2:
ENDFOR
;
; Define which outputs are wanted
;
out_par=INTARR(7)
menu=[' Do you want the Final Fit Parameter Sigmas? ', 'Yes', 'No']
out_par(0)=one_line_menu(menu,init=2)
menu=[' Do you want the Final Fit Residuals? ', 'Yes', 'No']
out_par(1)=one_line_menu(menu,init=2)
menu=[' Do you want the Spectrum w/dust fit subtracted?', 'Yes', 'No']
out_par(2)=one_line_menu(menu,init=2)
menu=[' Do you want the Spectrum w/cmb fit subtracted?', 'Yes', 'No']
out_par(3)=one_line_menu(menu,init=2)
menu=[' Do you want the Final Fit Failure Flag?', 'Yes', 'No']
;out_par(4)=one_line_menu(menu,init=1)
out_par(4)=1
menu=[' Do you want the Fit Chi_squared/df ?', 'Yes', 'No']
out_par(5)=one_line_menu(menu,init=2)
menu=[' Do you want the Number of Iterations for fit?', 'Yes', 'No']
out_par(6)=one_line_menu(menu,init=2)
; 
;set up proper strings for call
;
exstr= 'fit= FITMU(freq,spec,weight,parm0=parm0,active=active'
; Input values
IF (pix_parm(0)) THEN exstr=exstr+',cmb_temp=cmb_temp'
IF (pix_parm(1)) THEN exstr=exstr+',cmb_mu=cmb_mu'
IF (pix_parm(2)) THEN exstr=exstr+',bb_temp=bb_temp'
IF (pix_parm(3)) THEN exstr=exstr+',opdep=opdep'
IF (pix_parm(4)) THEN exstr=exstr+',emis=emis'
IF (pix_parm(6)) THEN exstr=exstr+',tolerance=tolerance'
IF (pix_parm(7)) THEN exstr=exstr+',funits=' +funits
IF (pix_parm(8)) THEN exstr=exstr+',' +sunit
;output values
exstr=exstr+',out_parm=out_val'
IF (out_par(0) EQ 1) THEN exstr=exstr+',parm_sig=parm_sig'
IF (out_par(1) EQ 1) THEN exstr=exstr+',residual=residual'
IF (out_par(2) EQ 1) THEN exstr=exstr+',nodust=nodust'
IF (out_par(3) EQ 1) THEN exstr=exstr+',dust=dust'
IF (out_par(4) EQ 1) THEN exstr=exstr+',fail=fail'
IF (out_par(5) EQ 1) THEN exstr=exstr+',chi2=chi2'
IF (out_par(6) EQ 1) THEN exstr=exstr+',n_iter=n_iter'
exstr=exstr+')'
;
;  Call FITMU
;  ------------
 j = EXECUTE(exstr)
;
;  Store the result in a UIMAGE data-object.
;  -----------------------------------------
;
; Output for graphs
;
IF (first_letter EQ 'G') THEN BEGIN
ntitle='Fit to ' +title
gname= setup_graph(freq,fit,num,title=ntitle,x_title=funits,y_title=sunit)
plot_graph,gname
IF (out_par(1) EQ 1) THEN BEGIN
 ntitle='Residual = spec-fit'
 gname= setup_graph(freq,residual,num,title=ntitle,x_title=funits,y_title=sunit)
 plot_graph,gname
ENDIF
IF (out_par(2) EQ 1) THEN BEGIN
 ntitle='Spectrum w/dust fit subtracted'
 gname= setup_graph(freq,nodust,num,title=ntitle,x_title=funits,y_title=sunit)
 plot_graph,gname
ENDIF
IF (out_par(3) EQ 1) THEN BEGIN
 ntitle='Spectrum w/cmb fit subtracted'
 gname= setup_graph(freq,dust,num,title=ntitle,x_title=funits,y_title=sunit)
 plot_graph,gname
ENDIF
ptext=STRARR(9)
ptext(0)='Output Parameter Values: '
ptext(1)='Bose Einstein Temperature: ' + STRING(out_val(0))
ptext(2)='Bose Einstein Chemical Pot.' + STRING(out_val(1))
ptext(3)='Black Body Temperature: ' + STRING(out_val(2)) 
ptext(4)='Black Body Optical Depth: ' + STRING(out_val(3)) 
ptext(5)='Black Body Emissivity Spectral Index: ' + STRING(out_val(4)) 
IF (out_par(4) EQ 1) THEN BEGIN
 IF (fail(0)) THEN str_fail='1' ELSE str_fail='0'
 ptext(6)='fail (1=failure) = '+ str_fail
ENDIF
IF (out_par(5) EQ 1) THEN ptext(7)='chi2/df = '+STRING(chi2)
IF (out_par(6) EQ 1) THEN ptext(8)='n_iter = '+STRING(n_iter)
IF (out_par(0) EQ 1) THEN BEGIN
 FOR i=1,5 DO ptext(i)=ptext(i) + ' +/-' + STRING(parm_sig(i-1))
ENDIF
FOR i=0,8 DO PRINT, ptext(i)
ENDIF ELSE BEGIN
;
; Output for 3-dim objects
;
; If fit fails then set all values equal to 0.
;
bad=WHERE(fail)
IF (bad(0) NE -1) THEN BEGIN
 fit(bad,*)=0.
 out_val(bad,*)=0.
 IF (out_par(0) EQ 1) THEN parm_sig(bad,*)=0.
 IF (out_par(1) EQ 1) THEN residual(bad,*)=0.
 IF (out_par(2) EQ 1) THEN nodust(bad,*)=0.
 IF (out_par(3) EQ 1) THEN dust(bad,*)=0.
 IF (out_par(5) EQ 1) THEN chi2(bad)=0.
 IF (out_par(6) EQ 1) THEN n_iter(bad)=0.
ENDIF
tempfit=fit
tout=out_val
IF (out_par(0) EQ 1) THEN tempsig=parm_sig 
IF (out_par(1) EQ 1) THEN tempres=residual
IF (type EQ 'FACE') THEN BEGIN
 pix2xy,pixlist,res=res,/face,bad_pixval=badpixval,data=tempfit,raster=fit
 IF (out_par(0) EQ 1) THEN pix2xy,pixlist,res=res,/face,bad_pixval=badpixval,$
  data=tempsig,raster=parm_sig
 pix2xy,pixlist,res=res,/face,bad_pixval=badpixval,data=tout,raster=out_val
 IF (out_par(1) EQ 1) THEN pix2xy,pixlist,res=res,/face,bad_pixval=badpixval,$
  data=tempres,raster=residual
ENDIF ELSE BEGIN
 pix2xy,pixlist,res=res,bad_pixval=badpixval,data=tempfit,raster=fit
 IF (out_par(0) EQ 1) THEN pix2xy,pixlist,res=res,bad_pixval=badpixval,$
  data=tempsig,raster=parm_sig
 pix2xy,pixlist,res=res,bad_pixval=badpixval,data=tout,raster=out_val
 IF (out_par(1) EQ 1) THEN pix2xy,pixlist,res=res,bad_pixval=badpixval,$
  data=tempres,raster=residual
ENDELSE
; Extract a slice to select out parts of out_val.
;
FOR i=0,4 DO BEGIN
 IF (active(i) EQ 1) THEN BEGIN
  data=out_val(*,*,i)
  ntitle=append_number(parm_title(i))
  badpixval=0.
  good = WHERE(data NE badpixval)
  IF(good(0) NE -1) THEN BEGIN
    scale_min = MIN(data(good), MAX = scale_max)
  ENDIF ELSE BEGIN
    scale_min = 0.
    scale_max = 0.
  ENDELSE
  IF (type EQ 'FACE') THEN BEGIN
    name=setup_image(data=data,badpixval=badpixval,scale_min=scale_min,$
     instrume=object3d(index3d).instrume,title=ntitle,faceno=faceno,$
     scale_max=scale_max)
  ENDIF ELSE BEGIN
    name=setup_image(data=data,badpixval=badpixval,scale_min=scale_min,$
     instrume=object3d(index3d).instrume,title=ntitle,scale_max=scale_max)
  ENDELSE
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay,name
  IF (out_par(0) EQ 1) THEN BEGIN
   data=parm_sig(*,*,i)
   ntitle=append_number('Sigma of ' +parm_title(i))
   badpixval=0.
   good = WHERE(data NE badpixval)
   IF(good(0) NE -1) THEN BEGIN
     scale_min = MIN(data(good), MAX = scale_max)
   ENDIF ELSE BEGIN
     scale_min = 0.
     scale_max = 0.
   ENDELSE
   IF (type EQ 'FACE') THEN BEGIN
     name=setup_image(data=data,badpixval=badpixval,scale_min=scale_min,$
      instrume=object3d(index3d).instrume,title=ntitle,faceno=faceno,$
      scale_max=scale_max)
   ENDIF ELSE BEGIN
     name=setup_image(data=data,badpixval=badpixval,scale_min=scale_min,$
      instrume=object3d(index3d).instrume,title=ntitle,scale_max=scale_max)
   ENDELSE
   IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay,name
  ENDIF
 ENDIF
ENDFOR
;
; Create Residual 3D object
;
IF (out_par(1) EQ 1) THEN BEGIN
 menu=[' Create a 3D object for the Final Fit Residuals? ', 'Yes', 'No']
 choice=one_line_menu(menu,init=2)
 IF (choice EQ 1) THEN BEGIN
  ntitle='Final Fit Residuals'
  name=setup_object3d(data=residual,badpixval=badpixval,frequency=freq,$
    title=ntitle,units=object3d(index3d).units,$
    frequnits=object3d(index3d).frequnits,$
    instrume=object3d(index3d).instrume,$
    orient=object3d(index3d).orient,$
    coordinate_system=object3d(index3d).coordinate_system,$
    faceno=object3d(index3d).faceno)
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay,name
 ENDIF
ENDIF
ENDELSE
;
;  If journaling is enabled, then write info to the journal file.
;  --------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:  ' + title
    PRINTF, luj, '  the following statement was executed:'
    PRINTF, luj, '  ' + STRMID(exstr, 0, 76)
    exstr = STRMID(exstr, 76, STRLEN(exstr))
    WHILE(STRLEN(exstr) GE 1) DO BEGIN
      PRINTF, luj, '    ' + STRMID(exstr, 0, 74)
      exstr = STRMID(exstr, 74, STRLEN(exstr))
    ENDWHILE
    IF(wsel EQ 1) THEN BEGIN
      PRINTF, luj, '    weights came from:  ' + titlew
    ENDIF
    PRINTF,luj,'Input Parameter Values: (flag=0 means held fixed)'
    FOR i=0,4 DO PRINTF,luj,parm_title(i),parm0(i),active(i)
    PRINTF,luj,' '
    IF (first_letter EQ 'G') THEN FOR i = 0,8 DO PRINTF, luj, '  ' +ptext(i)
    PRINTF, luj, '----------------------------------------' + $
                 '--------------------------------------'
  ENDIF
  IF(uimage_version EQ 2) THEN GOTO, show_menu
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


