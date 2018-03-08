FUNCTION mkweight, link1,link2
;+NAME/ONE LINE DESCRIPTION:
;    MKWEIGHT adds quadratically 2 weight arrays (properly handling bad pixels).
;
;DESCRIPTION:
;    MKWEIGHT receives 2 window numbers which are pointers to weight
;    arrays, creates a new weight array and outputs a pointer to the
;    new array. This function adds quadratically two weight arrays
;    to create a new weight array and puts this in a UIMAGE object.
;
;CALLING SEQUENCE:
;    newlink = mkweight(link1,link2)
;
;ARGUMENTS (I=input, O=output, []=optional)
;    link1         I   Scalar   A window number of a weight array.
;    link2         I   Scalar   A window number of a weight array.
;
;#
;COMMON BLOCKS:  none.
;
;LIBRARY CALLS:  SETUP_IMAGE
;
;PROCEDURE: The routine recieves as input two window numbers which are
;           pointers to weight arrays. The weight arrays and bad pixel
;           values are read in. A new weight is created assuming that
;           the SIGMAS add quadratically, i.e Snew^2=S1^2+S2^2. This
;           transfers to weights as Wnew=W1*W2/(W1+W2). A new UIMAGE
;           weight object is created and the window number is passed
;           back to the calling routine for use in .LINKWEIGHT field.
;
;REVISION HISTORY:
;    SPR 10829 Creation: J. Newmark, 4 April 1993
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;
;.TITLE
;Routine MKWEIGHT
;-
;
;  Check the validity of the passed parameters.
;  --------------------------------------------
  common uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  common uimage_data2,proj_map,proj2_map
  IF(N_PARAMS() LT 2) THEN BEGIN
    MESSAGE,'Insufficient parameters. (Two Link numbers are expected).',/CONT
    RETURN,-1
  ENDIF
badval_in=FLTARR(2)
weight1=get_name(link1)
weight2=get_name(link2)
IF (STRMID(weight1,0,5) EQ 'GRAPH') THEN graphs = 1 ELSE graphs = 0
j = EXECUTE('badval_in(0) = ' + weight1 + '.badpixval')
j = EXECUTE('badval_in(1) = ' + weight2 + '.badpixval')
IF (NOT graphs) THEN BEGIN
 j = EXECUTE('wt1 = ' + weight1 + '.data')
 j = EXECUTE('wt2 = ' + weight2 + '.data')
 good = WHERE((wt1 NE badval_in(0)) AND (wt2 NE badval_in(1)))
 outbadval=MIN(badval_in)
 outwt=wt1*0. + outbadval
 outwt(good)=(wt1(good)*wt2(good))/(wt1(good)+wt2(good))
;  -----------------------------------------------------------------------
 w = WHERE(outwt(good) EQ outbadval)
 IF(w(0) NE -1) THEN BEGIN
      outbadval = MIN(outwt(good))-1.
      temp = outwt(good)
      outwt = outwt*0. + outbadval
      outwt(good) = temp
 ENDIF
 scale_min= MIN(outwt(good), MAX=scale_max)
 parent=weight1
 title=append_number('New weight')
;
; Set up new UIMAGE weight object.
;
 newweight=setup_image(data=outwt,title=title,badpixval=outbadval,$
    scale_min=scale_min, scale_max=scale_max, hidden=1, parent=parent)
ENDIF ELSE BEGIN
 j = EXECUTE('num =' +weight1 + '.num')
 j = EXECUTE('wt1 = ' + weight1 + '.data(0:num-1,1)')
 j = EXECUTE('freq = ' + weight1 + '.data(0:num-1,0)')
 j = EXECUTE('wt2 = ' + weight2 + '.data(0:num-1,1)')
 good = WHERE((wt1 NE badval_in(0)) AND (wt2 NE badval_in(1)))
 outbadval=MIN(badval_in)
 outwt=wt1*0. + outbadval
 outwt(good)=(wt1(good)*wt2(good))/(wt1(good)+wt2(good))
 good = WHERE(outwt NE outbadval)
 IF (good(0) EQ -1) THEN BEGIN
   PRINT,'Resulting graph consists of all bad values - operation ' + $
     'aborted.'
   RETURN, -1
 ENDIF
 j=EXECUTE('x_title = ' +weight1 + '.x_title')
 j=EXECUTE('y_title = ' +weight1 + '.y_title')
 wtitle=append_number('New weight')
 newweight=setup_graph(freq,outwt,num,title=wtitle,badpixval=outbadval,$
   x_title=x_title,y_title=y_title,hidden=1)
ENDELSE
j = EXECUTE('link3= ' + newweight + '.window')
RETURN, link3
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


