PRO psym_chg, name
;+
;  PSYM_CHG- a UIMAGE-specific routine.  This routine allows
;  an user to select a graph and then to select a new plot
;  symbol at which the graph will be plotted.
;#
;  SPR 11102  Jun 15 93  Creation J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON color_values,c_badpix,c_draw,c_scalemin
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  first = 1
  menu_title = "Change a graph's plot symbol"
  IF(N_PARAMS(0) EQ 1) THEN GOTO, have_name
;
;  Show a menu in which the user can select an available graph.
;  ------------------------------------------------------------
show_menu:
  name = select_object(/graph, /help, /exit, title = menu_title)
  IF(name EQ 'EXIT') THEN BEGIN
    IF(journal_on AND (first NE 1)) THEN $
      PRINTF, luj, '----------------------------------------' + $
                   '--------------------------------------'
    RETURN
  ENDIF
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no graphs currently available.', /CONT
    RETURN
  ENDIF
have_name:
;
;  Let the user select a plot symbol to be used for that graph.
;  ------------------------------------------------------------
  psymbol=0
  READ,'Enter choice of plotting symbols (IDL defaults = 0-7,10)',psymbol
  IF(psymbol LT 0) THEN RETURN
;
;  Re-plot the graph (using the new plot symbol).
;  ----------------------------------------------
  j = EXECUTE(name + '.psymbol = psymbol')
  j = EXECUTE('window = ' + name + '.window')
  IF ((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
   WSET, window
   WSHOW, window
  ENDIF
  plot_graph, name
;
;  If journaling is enabled, then report what was done to the journal file.
;  ------------------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    IF(first) THEN PRINTF, luj, menu_title
    ig = FIX(STRMID(name, 6, STRPOS(name, ')') - 6))
    PRINTF, luj, ' plot symbol ' + STRTRIM(STRING(psymbol),2) + $
      ' was assigned to:  ' + graph(ig).title
  ENDIF
  first = 0
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


