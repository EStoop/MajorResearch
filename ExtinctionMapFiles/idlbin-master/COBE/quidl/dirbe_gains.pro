function dirbe_gains,iflag,time
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    DIRBE_GAINS returns the auto or commanded gain ratios per detector
;
; DESCRIPTION:
;    IDL function to return the auto-gain or commanded gain for a 
;    specified time.
;
; CALLING SEQUENCE:
;     GAIN = DIRBE_GAIN (gtype,ztime)
;
; ARGUMENTS (I = input, O = output, [] = optional)
;     GAIN         O     fltarr    Array of gain ratios in science 
;                                  detector order.
;     GTYPE        I     int       Gain selection flag (0 = auto_gain,
;                                  1 = commanded_gain).
;     ZTIME        I     str       Time in zulu format.
;
; WARNINGS:
;     1.  Auto-gain ratios were derived from optimized cold-mission 
;         data, but they are currently being applied to all data.
;     2.  Gain ratios are not automatically updated when pipeline 
;         ratios change.
;     3.  If out-of-range arguments are supplied, -1 will be returned.
;
; EXAMPLE:
;     To get auto-gain ratios on 01-January-1990:
;       auto = DIRBE_GAINS(0,'900010000')
;     To get commanded-gains on 01-January-1990:
;       comm = DIRBE_GAINS(1,'900010000')
;
; COMMON BLOCKS:
;     None.
;
; PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     Gains are hard-coded within the procedure based on information
;     contained in the BLI baseline scripts, BEX_GCAFT, and BEX_NOM as
;     of 18-March-1992.
;
; PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     None.
;
; MODIFICATION HISTORY:
;     Written by BA Franz, Applied Research Corp.      18 March 1991
;     CCR 585
;.TITLE
;Routine DIRBE_GAINS
;-

if (n_params(0) lt 2) then begin
    print,'Invalid call to get_gain_ratios.'
    return,-1
endif  

;
; Warning: Currently using mux2 x16 ratios for all time ********
;

if ((time(0) ge '89322142400000') and $
    (time(0) lt '89345000000000')) then begin
  ;
  ; Mux1
  ;
  case iflag of
    0:begin
        ;
        ; X16 to X1 ratios
        ;
        gain = [15.78289,15.77305,15.7690, $    ;band 1
                15.88996,15.98606,15.97883, $   ;band 2
                15.7648,16.00026,16.03776, $    ;band 3
                15.95889, $                     ; 4
                15.79972, 15.93593, $           ;5,6
                15.84680, 15.9267,  $           ;7,8
                15.83417, 15.73266 ]            ;9,10
      end

    1:begin
        ;
        ; Commanded gain ratios
        ;
        gain = [1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,$
                1.00,1.00,1.00,1.00,1.00,1.00]

      end

    else:begin
        print,'Invalid call to get_gain_ratios.'
        return,-1
      end

  endcase

endif else if ((time(0) ge '89345000000000') and $
               (time(0) lt '90266001635660')) then begin
  ;
  ; Mux2
  ;
  case iflag of
    0:begin
        ;
        ; X16 to X1 ratios
        ;
        gain = [15.78289,15.77305,15.7690, $    ;band 1
                15.88996,15.98606,15.97883, $   ;band 2
                15.7648,16.00026,16.03776, $    ;band 3
                15.95889, $                     ; 4
                15.79972, 15.93593, $           ;5,6
                15.84680, 15.9267,  $           ;7,8
                15.83417, 15.73266 ]            ;9,10
      end

    1:begin
        ;
        ; Commanded gain ratios
        ;
        gain = [3.01,3.01,3.01,3.01,3.01,3.01,3.01,3.01,3.01,3.01,$
                1.00,1.00,1.00,1.00,3.01,3.01]

      end

    else:begin
        print,'Invalid call to get_gain_ratios.'
        return,-1
      end

  endcase


endif else if ((time(0) ge '90266001635660') and $
               (time(0) lt '90270013115670')) then begin
  ;
  ; Mux3
  ;
  case iflag of
    0:begin
        ;
        ; X16 to X1 ratios
        ;
        gain = [15.78289,15.77305,15.7690, $    ;band 1
                15.88996,15.98606,15.97883, $   ;band 2
                15.7648,16.00026,16.03776, $    ;band 3
                15.95889, $                     ; 4
                15.79972, 15.93593, $           ;5,6
                15.84680, 15.9267,  $           ;7,8
                15.83417, 15.73266 ]            ;9,10
      end

    1:begin
        ;
        ; Commanded gain ratios
        ;
        gain = [3.01,3.01,3.01,3.01,3.01,3.01,3.01,3.01,3.01,3.01,$
                1.00,1.00,1.00,1.00,3.01,3.01]

      end

    else:begin
        print,'Invalid call to get_gain_ratios.'
        return,-1
      end

  endcase


endif else if ((time(0) ge '90270013115670') and $
               (time(0) lt '91021000000000')) then begin
  ;
  ; Mux4
  ;
  case iflag of
    0:begin
        ;
        ; X16 to X1 ratios
        ;
        gain = [15.78289,15.77305,15.7690, $    ;band 1
                15.88996,15.98606,15.97883, $   ;band 2
                15.7648,16.00026,16.03776, $    ;band 3
                15.95889, $                     ; 4
                15.79972, 15.93593, $           ;5,6
                15.84680, 15.9267,  $           ;7,8
                15.83417, 15.73266 ]            ;9,10
      end

    1:begin
        ;
        ; Commanded gain ratios
        ;
        gain = [9.07,9.07,9.07,9.07,9.07,9.07,9.07,9.07,9.07,9.07,$
                1.00,1.00,1.00,1.00,3.01,3.01]

      end

    else:begin
        print,'Invalid call to get_gain_ratios.'
        return,-1
      end

  endcase

endif else if ((time(0) ge '91021000000000') and $
               (time(0) lt '99000000000000')) then begin
  ;
  ; Mux5
  ;
  case iflag of
    0:begin
        ;
        ; X16 to X1 ratios
        ;
        gain = [15.78289,15.77305,15.7690, $    ;band 1
                15.88996,15.98606,15.97883, $   ;band 2
                15.7648,16.00026,16.03776, $    ;band 3
                15.95889, $                     ; 4
                15.79972, 15.93593, $           ;5,6
                15.84680, 15.9267,  $           ;7,8
                15.83417, 15.73266 ]            ;9,10
      end

    1:begin
        ;
        ; Commanded gain ratios
        ;
        gain = [9.07,9.07,9.07,9.07,9.07,9.07,9.07,9.07,9.07,9.07,$
                1.00,1.00,1.00,1.00,3.01,3.01]

      end

    else:begin
        print,'Invalid call to get_gain_ratios.'
        return,-1
      end

  endcase

endif else begin
  print,'Invalid call to get_gain_ratios.'
  return,-1
endelse

return,gain
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


