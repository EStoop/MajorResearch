PRO LABELIT, MAINLABEL=Mainlabel, BOTLABEL=Botlabel, YLABEL=Ylabel, $
        SECLABEL=Seclabel, RIGHTLABEL=Rightlabel, NAME_TIME=Name_time
;+ NAME/ONE LINE DESCRIPTION OF ROUTINE
;       LABELIT is a graph labeling procedure for IDL Version II.
;
;Description:  
;     At user option, labels top, bottom, left and right sides and a 
;          top subtitle.  Automatically centers each label.
;     Always prints user ID, date, and time at bottom.
;     If called with no keywords, program will print an informative 
;             message regarding the calling sequence and exit.
;
;CALLING SEQUENCE:
;     LABELIT [, MAINLABEL=mainlabel] [, BOTLABEL=botlabel] $
;           [, YLABEL=ylabel] [, SECLABEL=seclabel] $
;           [, RIGHTLABEL=rightlabel] [,/NAME_TIME]
;
;ARGUMENTS (I  = input, O = output, [] = optional)
;     Mainlabel     [I]  str        Label at top of graph 
;     Botlabel      [I]  str        Label at bottom of graph 
;     Ylabel        [I]  str        Y-axis label (left hand side)
;     Seclabel      [I]  str        Secondary label at top
;     Rightlabel    [I]  str        Y-axis label (right hand side)
;     Name_Time     [I]  key        Indicator to print name and time;
;                                      default is 1 IF any other
;                                      keyword is specified.  Defaults 
;                                      to 0 if no keywords specified.
;
;OUTPUT:  None
;
;EXAMPLES:  
;        1.  Label a graph with 'Original Time Series' as the main 
;        title, 'Bin Number' on the x-axis, 'Counts' on the y-axis, 
;        'FIRAS Data' as the subtitle, and 'TIME_SERIES.PS_PLOT' on 
;        the right hand side of the y-axis:
;
;        LABELIT,Mainlabel = 'ORIGINAL TIME SERIES', Botlabel = $
;                'Bin Number', Ylabel = 'Counts', $
;                Seclabel = 'FIRAS Data', $
;                Rightlabel = 'TIME_SERIES.PS_PLOT'
;
;        2.  Label a graph with 'Original Time Series' as the main
;        title, 'Bin Number' on the x-axis, 'Counts' on the y-axis,
;        no subtitle, and 'TIME_SERIES.PS_PLOT' on the right 
;        hand side of the y-axis:
;
;        LABELIT,Mainlabel = 'ORIGINAL TIME SERIES', Botlabel = $
;                'Bin Number', Ylabel = 'Counts', $
;                Rightlabel = 'TIME_SERIES.PS_PLOT'
;
;        3.  Label a graph with 'Original Time Series' as the main 
;        title, no x-axis label, and 'Counts' on the y-axis
;
;        LABELIT, Mainlabel ='ORIGINAL TIME SERIES', Ylabel ='Counts'
;
;       4.  Print only user name and current time on bottom of graph.
;
;       LABELIT, /Name_Time       or LABELIT, /N
;
;       5.  Print informative message regarding calling sequence.
;
;       LABELIT
;#
;COMMON BLOCKS:
;       None
;
;PROCEDURE:
;       Using normal coordinate system and the XYOUTS procedure,
;       centers all labels given and writes them in varying size
;       script (top label is largest, and so on.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.
;      None
;
;REVISION HISTORY:
;
;       Written by Alice Trenholme 2/4/91
;       Revised by ART 3/17/92 to use keywords
;
;.title
;Routine LABELIT
;-
        toploc = .98
        secloc = .96
        yloc = .05
        rloc = 1
        botloc = .025
        timeloc = .025
        !Fancy = 1.5

        name = getenv('USER')

        keys_there = n_elements(mainlabel)+n_elements(seclabel)+ $
           n_elements(ylabel)+n_elements(rightlabel)+n_elements(botlabel)

        if (n_elements(mainlabel) ne 0) then  $
                xyouts,.5,toploc,$
                '!6' + mainlabel,size = 1.1,/normal,alignment = .5
        if (n_elements(seclabel) ne 0) then  $
                xyouts,.5,secloc,$
                '!3' + seclabel,size = 0.75,/normal,alignment = .5
        if n_elements(ylabel) ne 0 then  $
                xyouts,yloc,.5,        '!3' + ylabel,size = 1.1, $
                orientation = 90,/normal,alignment = .5
        if n_elements(rightlabel) ne 0 then  $
                xyouts,rloc,.5,        '!3' + rightlabel,size = 1.1, $
                orientation = 90,/normal,alignment = .5
        if n_elements(botlabel) ne 0 then  $
                xyouts,.5,botloc,$
                '!3' + botlabel,size = .9, /normal,alignment = .5

        if ((keys_there gt 0) or keyword_set(name_time)) then begin
                 xyouts,0,timeloc,name,size=0.8,/normal,alignment = 0.
                xyouts,1.,timeloc,!stime,size = 0.65,/normal,alignment = 1.

        endif else begin
                print, $
'     LABELIT [, MAINLABEL=mainlabel] [, BOTLABEL=botlabel] $'
                print, $
'           [, YLABEL=ylabel] [, SECLABEL=seclabel] [, RIGHTLABEL=rightlabel] $'
                print, $
'           [,/NAME_TIME]'
                print, ' '
                print, $
'     Mainlabel   -- Optional string --    Label at top of graph '
                print, $
'     Botlabel    -- Optional string --    Label at bottom of graph '
                print, $
'     Ylabel      -- Optional string --    Y-axis label (left hand side)'
                print, $
'     Seclabel    -- Optional string --    Secondary label at top'
                print, $
'     Rightlabel  -- Optional string --    Y-axis label (right hand side)'
                print, $
'     Name_Time   (Specify with /Name_Time to set)'
                print, $
'                                           Indicator to print name and time;'
                print, $
'                                             specify ONLY if you are giving no '
                print, $
'                                             other labels; otherwise name and'
                print, $
'                                             time will be printed automatically.'
        endelse
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


