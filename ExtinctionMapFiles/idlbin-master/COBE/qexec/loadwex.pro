 PRO LoadWEX, WEXFile, Text, Title, SubTopic, SubFile, ErrStat
;+
;   Description
;   -----------
;   LoadWEX reads and parses a WEX format file in order to extract text,
;   file title, related topics to the text, and what files are they
;   located in.  The information is used to build menus and text screens.
;
;   Calling Sequence:
;
;        LoadWEX, WEXFile, Text, Title, SubTopic, SubFile, ErrStat
;
;   where,
;
;        WEXFile  = input string variable containing the WEX file name
;
;        Text     = output array of strings containing text found in the
;                   WEX file.
;
;        Title    = output string variable which passes the title of
;                   the WEX component defined in the WEX file.
;
;        SubTopic = output string array containing a list of menu items
;                   found in the WEX file.
;
;        SubFile  = output string array containing a list of file names
;                   or command associated with the items in SubTopic.
;
;        ErrStat  = output variable returning the I/O error status.
;#
;   Called by
;   ---------
;   DispHelp
;   DispWEX
;   Extract
;   KeyHelp
;   MenuHelp
;   SpecHelp
;
;   Routines Called
;   ---------------
;   InitWEX      WEXName
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  November 1992
;
;   01/14/93  KRT  Added test and error handling for null WEX acquistion...
;                  the answer to SPR 10455.
;
;.Title
;   Procedure LoadWEX
;-
;
 On_Error,2  ; Return to caller if an error occurs
;
;Begin,
;
;   Initialize Variables :
;
    InitWEX
;
    Title     = ''      ;  Title is set to blank if not found
    Record    = ''      ;  String for holding input record
    KeyWord   = ''      ;  String for holding control keywords
;
    CMark     = ';'     ;  Comment marker
    GetComm   = 0       ;  Flag comment input: 1 = input  0 = cull
    StripMode = 0       ;  Mode flag for striping off first character
    KeyStart  = 0       ;  Logical indicating that a keyword mode has started
    KeyMark   = '.'     ;  Keyword marker
    SFNMark   = '==>'   ;  Subfile name symbol
    MrkLen    = StrLen( SFNMark )
;
;   Open WEX Files :
;
    FileName = WEXName( WEXFile )
    Get_LUN, LUN
;
    OpenR, LUN, FileName, Error=ErrStat
;
;   Open Error Exception Handling
;   -----------------------------
;
    If (ErrStat ne 0) then begin
       Free_LUN, LUN
       RETURN
    EndIf
;
;   WEX File Input Loop
;   -------------------
;
    While (not EoF( LUN )) do begin
;
       ReadF, LUN, Record                ; Read record
       ProCheck  = StrMid( Record, 0, 2) ; check for keyword marker
;
       Case ProCheck of
;
          ';+'    :Begin
                      KeyStart  = 1
                      KeyWord   = 'HELP'
                      StripMode = 1
                      GetComm   = 0
                   End
;
          ';#'    :Begin
                      KeyStart  = 1
                      KeyWord   = 'CODING'
                      StripMode = 1
                      GetComm   = 0
                   End
;
          ';-'    :Begin
                      KeyStart  = 1
                      KeyWord   = 'SOURCE'
                      StripMode = 0
                      GetComm   = 1
                   End
       Else:
       EndCase
;
       If (StripMode) then Record = StrMid( Record, 1, StrLen( Record ) )
;
       FirstChr = StrMid( Record, 0, 1 )
;
       If ((GetComm) or (FirstChr ne CMark)) then begin
          If (FirstChr eq KeyMark) then begin
;
;            Keyword set the mode for a certain type of input action:
;
             KeyStart = 1
             KeyWord  = StrUpcase( StrTrim( Record, 2 ) )
             KeyLen   = StrLen( KeyWord ) - 1
             KeyWord  = StrMid( KeyWord, 1, KeyLen )
;
;            Check to see if the keyword is properly set:
;
             If (not ((KeyWord eq      'HELP') or $
                      (KeyWord eq      'TEXT') or $
                      (KeyWord eq    'SOURCE') or $
                      (KeyWord eq  'EXAMPLES') or $
                      (KeyWord eq     'TITLE') or $
                      (KeyWord eq 'SUBTOPICS') or $
                      (KeyWord eq 'HYPERTEXT') or $
                      (KeyWord eq       CMark)) ) then begin
;
                Print, 'The following record in the WEX file has a bad keyword'
                Print, Record
                RETURN         ; Error Exception Handling: Return to Caller!
;
             EndIf
;
          EndIf Else begin
;
;            Actions taken for each mode :
;
             Case KeyWord of

                'HELP'      :Begin ; reading in text
                                If (KeyStart) then begin
                                   Text     = [ Record ]
                                   KeyStart = 0
                                EndIf Else begin
                                   Text     = [ Text, Record ]
                                EndElse
                             End
;
                'TEXT'      :Begin ; reading in text
                                If (KeyStart) then begin
                                   Text     = [ Record ]
                                   KeyStart = 0
                                EndIf Else begin
                                   Text     = [ Text, Record ]
                                EndElse
                             End
;
                'TITLE'     :Begin ; reading in title
                                Title    = StrTrim( Record, 2 )
                                KeyStart = 0
                             End
;
                'SUBTOPICS' :Begin ; reading in subtopics and subfiles
;
;                                If (Record ne '') then begin
;                                  getting subtopic and corresponding subfile,
;
;                                  Parse out the subfile name :
                                   RecLen   = StrLen( Record )
                                   SFNPos   = StrPos( Record, SFNMark )+MrkLen
;
;                                  Parse out the subtopic name :
                                   DescEnd  = SFNPos-MrkLen-1
                                   Descrptn = StrMid( Record, 0, DescEnd )
                                   Descrptn = StrTrim( Descrptn, 2 )
                                   SFName   = StrMid( Record, SFNPos, RecLen )
                                   SFName   = StrTrim( SFName, 2 )
;
                                   If (KeyStart) then begin
;                                     getting the first subtopic and subfile :
;
                                      SubTopic = [ Descrptn ]
                                      SubFile  = [ SFName ]
                                      KeyStart = 0
;
                                   EndIf Else begin
;                                     building subtopic and subfile arrays :
;
                                      SubTopic = [ SubTopic, Descrptn ]
                                      SubFile  = [ SubFile, SFName ]
;
                                   EndElse
;                                EndIf
                             End
             Else:
             EndCase
;
          EndElse
       EndIf
    EndWhile
;
    If (not ErrStat) then begin
       Close, LUN
       Free_LUN, LUN
    EndIf
;
    If ((N_Elements( Subtopic ) eq 0) and (N_Elements( Text ) eq 0)) then begin
        ErrStat = -999
        Print, "WEX File '"+FileName+"' was empty or had no WEX keywords."
    EndIf
;
    Return
 End
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


