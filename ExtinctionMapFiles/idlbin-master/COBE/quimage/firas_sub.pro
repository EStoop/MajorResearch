Pro firas_sub, name, freq, fsize, subscripts

; initialize this to all data
subscripts = ''

Restart:
   Print, ' '
   Print, 'You have selected field '+ name + '.'
   Print, ' '
   Print, 'This field has an array of values for each pixel, so you should'
   Print, 'now select a range of elements or the retrieve the entire field.'
   Print, ' '
   Print, '   Dimensions of '+Name+' :   (1:'+StrTrim(String(FSize(1)),2)+')'
   Print, ' '
   Print, 'For a range, use a colon to separate the upper and lower bounds.'
   Print, 'Type "ALL" if you want the data for the entire field '+Name+'.'
   Print, ' '

      c = 29.97925
      N = FSize(1)-1

      LoGiga = StrN( Freq(0) )
      UpGiga = StrN( Freq(N) )
      LoWave = StrN( Freq(0)/c )
      UpWave = StrN( Freq(N)/c )
;
      Print, 'Specify units with the range using one of the following:'
      Print, ' '
      Print, '                                          Valid Ranges'
      Print, '  IDX  (default)   - indices             (1,'+StrN(Fsize(1))+')'
      Print, '  GHZ              - gigahertz           ('+LoGiga+','+UpGiga+')'
      Print, '  ICM              - wave numbers        ('+LoWave+','+UpWave+')'
      Print, ' '
      Print, 'Examples:   ---> 1:50  IDX     (gets the first 50 indices.)'
      Print, '            ---> 60:80 ICM    (gets all indices within the'
      Print, '                               range of wave numbers 60 and 80.)'
;
   Units   = [ 'IDX', 'GHZ', 'ICM' ]
   N_Units = N_Elements( Units )
   Done    = 0

   While (not Done) do begin
      InputStr = ''
      Read, '---> ', InputStr
      InputStr = STRUPCASE(InputStr)

      If (StrTrim( InputStr, 2 )eq 'ALL') then Goto, NoCheck

      If InputStr ne '' then Done  = 1

      If (StrTrim(InputStr, 2 ) eq 'BACK') then Goto, Back

      On_ioerror, Io_error
      Test = Fix(inputstr)

      RangStr = InputStr
      UnitStr = Units(0)

      For j = 0, N_Units-1 do begin
         StrEnd  = StrLen( RangStr )
         UnitPos = StrPos( InputStr, Units(j) )
         If (UnitPos ge 0) then begin
            RangStr = StrMid( RangStr, 0, UnitPos+1 )
            UnitStr = Units(j)
         EndIf
      EndFor
   EndWhile

   Channels    = (UnitStr eq Units(0))
   GigaHertz   = (UnitStr eq Units(1))
   WaveNumbers = (UnitStr eq Units(2))

   If Channels then begin
;     Check subscripts against rdf info
      Groups = StrParse( InputStr, Delim=',' )

      If (N_Elements( Groups ) ne FSize(0)) then begin
         Print, ' '
         Print, 'Error: Too many or too few sets of subscripts given.'
         Print, '       You must specify ',FSize(0),' ranges.'
         Print, ' '
         Print, 'Press any key to continue...'
         Ans = Get_Kbrd( 1 )
         Goto, Restart
      EndIf

      For i = 1, FSize(0) do begin

         Dims   = StrParse( Groups(i-1), Delim=':' )
         N_Dims = N_Elements( Dims )

         If (N_Dims gt 2) then begin
            Print, ' '
            Print, 'Error: format of ranges not correct. [ '+groups(i-1)+' ]'
            Print, ' '
            Print, 'Press any key to continue...'
            Ans = Get_Kbrd( 1 )
            Goto, Restart
         EndIf
	
         For j = 0, N_Dims-1 do begin
            If (Fix(Dims(j)) gt FSize(i)) then begin
               Print, ' '
               Print, 'Error: subscript is too large. [ ', Dims(j), ' ]'
               Print, '       maxium allowed value is ', FSize(i)
               Print, ' '
               Print, 'Press any key to continue...'
               Ans = Get_Kbrd( 1 )
               Goto, Restart
            EndIf 
            If (Fix(Dims(j)) le  0) then begin
               Print, ' '
               Print, 'Error: subscripts must be greater than 0.'
               Print, ' '
               Print, 'Press any key to continue...'
               Ans = Get_Kbrd( 1 )
               Goto, Restart
            EndIf 
         EndFor
;
;        Drop down the input indices and reconstruct 'subscripts':
         Str0 = StrTrim( String( Fix( Dims(0) )), 2 )

         If (N_Dims eq 1) then begin
            If (i eq 1) then subscripts = Str0 $
                        else subscripts = subscripts + ',' +  Str0
         EndIf Else begin
            Str1 = StrTrim( String( Fix( Dims(1) ) ) )
            If (i eq 1) then subscripts = Str0 + ':' + Str1 $
                        else subscripts = subscripts + ',' +  Str0 + ':' + Str1
         EndElse

      EndFor
   EndIf
;
   If (not Channels) then begin
      Dim = Float( StrParse( RangStr, Delim=':' ) )
      If (WaveNumbers) then Dim = Dim*c  ; Convert to GigaHertz

      If (N_Elements( Dim ) gt 1) then begin

         If (N_Elements( Dim ) gt 2) then begin
            Print, ' '
            Print, 'Error: only specify one upper bound and one lower bound.'
            Print, ' '
            Print, 'Press any key to continue...'
            Ans = Get_Kbrd( 1 )
            Goto, Restart
         EndIf

         If (Fix( Dim(1) ) lt Fix( Dim(0) )) then begin
            Print, ' '
            Print, 'Error: upper bound is smaller than lower bound.'
            Print, ' '
            Print, 'Press any key to continue...'
            Ans = Get_Kbrd( 1 )
            Goto, Restart
         EndIf

         If (Dim(1) lt Freq(0) or Dim(0) gt Freq(N)) then begin
            Print, ' '
            Print, 'Error: in range specification [ '+ RangStr +' ].'
            Print, ' '
            Print, 'Press any key to continue...'
            Ans = Get_Kbrd( 1 )
            Goto, Restart
         EndIf

         Diff    = Min(Abs(Freq - Dim(0)),Start)
         Diff    = Min(Abs(Freq - Dim(1)),Stop)
	 IRange = [Start, Stop]

         subscripts = StrTrim( String( IRange(0)+1 ), 2 ) + ':' +  $
                      StrTrim( String( IRange(1)+1 ), 2 )

         If (WaveNumbers) then begin
            LoFreq = String( Freq(IRange(0))/c, Format='(F7.1)' ) + ' '
            UpFreq = String( Freq(IRange(1))/c, Format='(F7.1)' )  + ' '
         EndIf Else begin
            LoFreq = String( Freq(IRange(0)), Format='(F7.1)' ) + ' '
            UpFreq = String( Freq(IRange(1)), Format='(F7.1)' )  + ' '
         EndElse

         NChans = IRange(1) - IRange(0) + 1

         Print, NChans, ' indices are selected, covering the range:'
         Print, ' '
         Print, 'From ' + LoFreq + UnitStr + ' at index ', IRange(0)+1
         Print, 'To   ' + UpFreq + UnitStr + ' at index ', IRange(1)+1
         Print, ' '

      EndIf Else begin

         Diff    = Min(Abs(Freq - Dim(0)),Closest)

         If (WaveNumbers) then                                            $
            SelFreq = String( Freq(Closest(0))/c, Format='(F7.1)' ) + ' ' $
         Else                                                             $
            SelFreq = String( Freq(Closest(0)), Format='(F7.1)' ) + ' '

         subscripts = StrTrim( String( Closest(0)+1 ), 2 )

         Print, 'The index selected nearest to your specification is :'
         Print, ' '
         Print, SelFreq + UnitStr + ' at index ', Closest(0)+1
         Print, ' '
;
      EndElse
   EndIf
;
;EndIf  ; (FSize(0) gt 0)
;
NoCheck:
;     Everything A.O.K.
      status = 1
      !err   = 1
      !error = 1
      Return

Back :
      subscripts = 'back'
      !ERR   = 88
      !ERROR = 88
      Return

Abort:
;     return to calling program
      Return
;
Io_Error: print,'There was an error in your input, please reenter'
          On_Ioerror, null
          goto, restart

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


