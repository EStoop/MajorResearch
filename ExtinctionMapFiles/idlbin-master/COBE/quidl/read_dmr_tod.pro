;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    READ_DMR_TOD- can be used to validate the COBE DMR Time-Ordered Data
;
;DESCRIPTION:
;     This code is mainly intended to provide a platform from which
;     user specfifc code can be written. Presently, this code reads in a
;     single days worth of DMR TOD one major_frame at a time. The data can
;     be written to a log file, stored as IDL save sets (1 major frame = 
;     1 saveset) and 1 large saveset created for the calibrated temperatures
;     of the selected major frames. 
;
;     A few functions are included to check for anomalous or out-of-range
;     values
;
;     With simple modifications, this program could be converted into an 
;     analysis tool for any of the interested fields.
;
;CALLING SEQUENCE:
;pro read_dmr_tod, tod_file, verbose, every_frame,saveframe,savetemp,logfile
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    tod_file    I     string              input file
;    verbose    [I]    byte                turn on verbose output
;    everyframe [I]    byte                turn off query for each frame
;    saveframe  [I}    byte                create saveset for each frame
;    savetemp   [I]    byte                create saveset for all temps.
;    logfile    [I]    string              dump output to named file
;
;WARNINGS:
;  1) Except for temperatures stored in IDL saveset (if savetemp flag is
;    set) all quoted values are in uncalibrated Data Units (DUs). Values
;    stored in saveset TEMPS.DAT are calibrated using:
;            temp=(du-mean_baseline)/diode_gain
;  2) If the /saveframe keyword is specified this will create a large
;     (>3 MB) temporary file and then delete it.
;
;EXAMPLE: 
; read_dmr_tod,'cacdm:[dmr_pass3.archive]dca4_tod.p3_8933100',1,0,1,1,'temp.log'
;   where the tod_file = one day of data = dca4_tod.p3_8933100
;   which uses verbose mode, queries for each major frame, creates an
;   IDL save set (name = dmr_tod.frame#), saves all the calibrated temperature
;   measurements (units=mK) in a file (name=temps.dat), and dumps output 
;   to a log file named temp.log
;
;#
; Code written by D. Leisawitz
;                 NASA/GSFC Code 631
;                 Astrophysics Data Facility
;  SPR 11774 6-Jun-1994 Initial delivery and modifications J. Newmark
;
;PROGRAMMING NOTES:
; Because the data were written in native format on a VAX, some byte
; manipulation is required.  For floating point numbers, the IDL
; Astronomy User's Library routine conv_vax_unix is used in the 
; function fconvert for UNIX platforms.
;
;-

; This function converts a byte array into an array of 2- or 4-byte integers
; of specified length output_length.  It operates successfully on an ULTRIX
; platform on the DMR time-ordered data.
function convert, byte_array, output_length,prt,lun
length = n_elements(byte_array)
bytes_per_word = length/output_length
byte_array = reform(byte_array, bytes_per_word, output_length)
case bytes_per_word of
2: output_array = fix(byte_array,0,output_length)
4: output_array = long(byte_array,0,output_length)
else: stat=execute(prt+",'WARNING: trying to make integer but bytes_per_word = ',bytes_per_word")
endcase
return, output_array
end

; This function converts a byte array into an array of floating point 
; numbers of specified length output_length.  If output_length is 1, a
; scalar is returned.  The function operates successfully on an ULTRIX 
; platform on the DMR time-ordered data.
function fconvert, byte_array, output_length,prt,lun
if (!version.arch ne 'vax') then $
  output_array = conv_vax_unix(float(byte_array,0,output_length)) $
else output_array = float(byte_array,0,output_length)
if (output_length eq 1) then output_array = output_array(0)
return, output_array
end

; This function converts an 8-byte array into the VAX ADT (time)
; of a DMR major frame.  ADT is the number of 100 ns intervals elapsed
; since 0 UT on 17 November 1858, not counting leap seconds.
function make_ADT, byte_array,prt,lun
factors=2.D^(8*(indgen(8)))
factors=reform(factors,1,8)
if (n_elements(byte_array) ne 8) then begin
	help, byte_array
	stat=execute(prt+",'WARNING: not 8 bytes available to make ADT'")
endif
ADT = factors # byte_array
ADT = ADT(0)			; remove extraneous dimension to make scalar
return, ADT
end

; For validation purposes, this function checks the ADT for a reasonable
; value
function check_adt, value,prt,lun
if ((value gt 4.1e16) and (value lt 4.3e16)) then begin
	return, 0
endif else begin
      stat=execute(prt+", 'WARNING: ADT = ', value, ' (~4.1 x 10^16 expected)'")
      return, 1
endelse
end

; For validation purposes, this function checks the SPIKE_LEVEL parameter
function check_spike_level, value,prt,lun
if (value eq 50) then begin
	return, 0
endif else begin
      stat=execute(prt+", 'WARNING: SPIKE LEVEL = ', value, ' (50 expected)'")
      return, 1
endelse
end

; For validation purposes, this function checks the WRITE_HIST parameter
function check_write_hist, value,prt,lun
if (value eq 4) then begin
	return, 0
endif else begin
	stat=execute(prt+", 'WARNING: WRITE_HIST = ', value, ' (4 expected)'")
	return, 1
endelse
end

; For validation purposes, this function checks the DSA_FLAG
function check_dsa_flag, value,prt,lun
if (value eq 0) then begin
	return, 0
endif else begin
	stat=execute(prt+", 'WARNING: DSA_FLAG = ', value, ' (0 expected)'")
	return, 1
endelse
end

; For validation purposes, this function checks the FRAME_FLAG
function check_frame_flag, value,prt,lun
if (value eq 0) then begin
	return, 0
endif else begin
	stat=execute(prt+", value, format='(""Note: some bad data in Major Frame; FRAME_FLAG = "", Z)'")
	return, 1
endelse
end

; For validation purposes, this function checks the CELESTIAL_OBJECTS flag
function check_celestial_objects, value,prt,lun
if (value eq 0) then begin
	return, 0
endif else begin
	stat=execute(prt+", 'According to the CELESTIAL_OBJECTS flag ...'")
	if (value mod 2) then stat=execute(prt+",'... the Moon is close'")
	if ((value/2) mod 2) then stat=execute(prt+",'... the Earth is above 31 GHZ cut angle'")
	if ((value/4) mod 2) then stat=execute(prt+",'... the Earth is above 53 GHZ cut angle'")
	if ((value/8) mod 2) then stat=execute(prt+",'... the Earth is above 90 GHZ cut angle'")
	if ((value/16) mod 2) then stat=execute(prt+",'... Jupiter is close'")
	if ((value/32) mod 2) then stat=execute(prt+",'... Mars is close'")
	if ((value/64) mod 2) then stat=execute(prt+",'... Saturn is close'")
	if ((value/128) mod 2) then stat=execute(prt+",'... the Galaxy is close'")
	return, 1
endelse
end

; For validation purposes, this function checks the ATTITUDE_TYPE
function check_attitude_type, value,prt,lun
case value of
0:	begin
		 stat=execute(prt+",'ERROR: No Attitude (ATTITUDE_TYPE = 0)'")
		return, 1
	end
1:	begin
		stat=execute(prt+",'Note: Attitude from SunEarth (ATTITUDE_TYPE = 1)'")
		return, 1
	end
2:	return, 0	; definitive attitude from DIRBE observations of stars
else:	begin
		stat=execute(prt+", 'ERROR: Attitude type ', value, ' undefined'")
		return, 1
	end
endcase
end

; For validation purposes, this function checks the GALAXY_CUT_ANGLE
; parameter
function check_galaxy_cut, value,prt,lun
if (value eq 15) then begin
	return, 0
endif else begin
	stat=execute(prt+", 'WARNING: Galaxy_Cut_Angle = ', value, ' (15 deg expected)'")
	return, 1
endelse
end

; For validation purposes, this function checks the PLANET_CUT_ANGLE
; parameter
function check_planet_cut, value,prt,lun
if (value eq 15) then begin
	return, 0
endif else begin
	stat=execute(prt+", 'WARNING: Planet_Cut_Angle = ', value, ' (15 deg expected)'")
	return, 1
endelse
end

; For validation purposes, this function checks the MOON_CUT_ANGLE
; parameter
function check_moon_cut, value,prt,lun
if (value eq 21) then begin
	return, 0
endif else begin
	stat=execute(prt+", 'WARNING: Moon_Cut_Angle = ', value, ' (21 deg expected)'")
	return, 1
endelse
end

; For validation purposes, this function checks the EARTH_CUT_ANGLE
; parameters
function check_earth_cut, values,prt,lun
if ((values(0) eq 253) and (values(1) eq 255) and (values(2) eq 255)) then begin
	return, 0
endif else begin
	stat=execute(prt+", 'WARNING: Earth_Cut_Angles = ', values")
	stat=execute(prt+", '         (253, 255, 255 deg expected at 31, 53, and 90 GHz)'")
	return, 1
endelse
end

; This function is a generic angle value checker.  For validation purposes,
; it looks for values outside the 0 - 360 degree range expected.
function check_angle, value, param_name,prt,lun
if ((value ge 0.) and (value le 360.)) then begin
	return, 0
endif else begin
	message = param_name + ' out of range; value = '
	stat=execute(prt+", message, value, ' deg (0 to 360 expected)'")
	return, 1
endelse
end

; This function is a generic angle value checker.  For validation purposes,
; it looks for radian values outside the -pi to 3pi/2 range expected.
function check_angle_rad, value, param_name,prt,lun
if ((value ge -!pi) and (value le 3.*!pi/2.)) then begin
	return, 0
endif else begin
	message = param_name + ' out of range; value = '
	stat=execute(prt+", message, value, ' deg (-pi to 3pi/2 expected)'")
	return, 1
endelse
end

; This function reads and checks values in the DMR Time-ordered data one
; major frame (32 s) at a time.  There is an option to skip over major
; frames, but the program reads daily TOD files sequentially (i.e., 
; chronologically).  Many of the comments come directly from the DMR
; Explanatory Supplement.
;
pro read_dmr_tod, tod_file,verbose,every_frame,saveframe,savetemp,logfile
if (n_elements(verbose) eq 0) then verbose=0
if (n_elements(saveframe) eq 0) then saveframe=0
if (n_elements(savetemp) eq 0) then savetemp=0 
if (savetemp) then begin
   temp_31a=fltarr(64,2700)
   temp_31b=fltarr(64,2700)
   temp_53a=fltarr(64,2700)
   temp_53b=fltarr(64,2700)
   temp_90a=fltarr(64,2700)
   temp_90b=fltarr(64,2700)
endif
if (n_elements(every_frame) eq 0) then every_frame=0
lun=0
if (n_elements(logfile) ne 0) then begin
    prt='printf,lun'
    openw,lun,logfile,/get_lun
endif else prt='print'
bytes_per_maj_frame = 2048
maj_frames_per_day = 2700		; number of frames expected in file
maj_frame_index = 0			; a counter
;
; set up diagnostic display
if ((!d.name eq 'X') or (!d.name eq 'WIN')) then begin
 window, 0, xsize = 1000, ysize=300, title = 'DMR TOD file ' + tod_file
 xmin = 0.
 xmax = 64.*maj_frames_per_day
 ymin = 1000.
 ymax = 4000.
 pcolor = [10,110,50,100,90,60]
 pcolor = reform(pcolor,3,2)
 loadct, 15
 plot,[xmin,xmax],[ymin,ymax],xstyle=1,ystyle=1, $
	xrange=[xmin,xmax],yrange=[ymin,ymax],color=30,/nodata
 xyouts, 1.e4, (ymax-300), '31A', color = pcolor(0,0), size = 1.5
 xyouts, 3.e4, (ymax-300), '31B', color = pcolor(0,1), size = 1.5
 xyouts, 5.e4, (ymax-300), '53A', color = pcolor(1,0), size = 1.5
 xyouts, 7.e4, (ymax-300), '53B', color = pcolor(1,1), size = 1.5
 xyouts, 9.e4, (ymax-300), '90A', color = pcolor(2,0), size = 1.5
 xyouts, 11.e4, (ymax-300), '90B', color = pcolor(2,1), size = 1.5
endif
openr, unit, tod_file, /get_lun		; open the TOD file
on_ioerror, cleanup
bytes = bytarr(bytes_per_maj_frame)	; make an array of bytes
;
next_frame:
maj_frame_index = maj_frame_index + 1
if (maj_frame_index gt maj_frames_per_day) then begin
	stat=execute(prt+", 'no more major frames in file'")
	goto, cleanup
endif
readu, unit, bytes			; read major frame data from file
ADT = make_ADT(bytes(0:7))	; number of 100 ns intervals since 0 UT
				; on 17 Nov 1858 (should be ~4 x 10^16)
;
; The COBETRIEVE_HEADER starts at byte 8 and is 22 bytes long
; -----------------------------------------------------------
MAJOR_FRAME_NUMBER = long(bytes(8:11),0)
ORBIT_NUMBER = long(bytes(12:15),0)
stat=execute(prt+",' '")
stat=execute(prt+",'===================================================='")
stat=execute(prt+",MAJOR_FRAME_NUMBER,format='(""Major Frame Number "",i8)'")
stat=execute(prt+", ADT, format='(""Major frame ADT = "", E18.12)'")
ad = check_adt(ADT,prt,lun)			; <<< diagnostic >>>
stat=execute(prt+",ORBIT_NUMBER,format='(""Orbit Number "",i8)'")
stat=execute(prt+", ' '")
SPIKE_LEVEL = bytes(16)			; x RMS X 10 for spike flagging
sl = check_spike_level(SPIKE_LEVEL,prt,lun)	; <<< diagnostic >>>
WRITE_HIST = bytes(17)			; Last Write: 0=DSA,2=DTO,4=DCA
wh = check_write_hist(WRITE_HIST,prt,lun)	; <<< diagnostic >>>
DSA_FLAG = bytes(18)			; 1=Derived From a Simulated Archive
df = check_dsa_flag(DSA_FLAG,prt,lun)		; <<< diagnostic >>>
FRAME_FLAG = fix(bytes(19:20),0)	; Major-Frame Quality Flag
;					  Bits 0 to 15 set for:
						; Not_clear, Gap, Eng, Dump, 
						; Bad_att, Prev_cal, Cal_up, 
						; Cal_dn, All_bad, 
						; 31A -> 90B_bad, spare
ff = check_frame_flag(FRAME_FLAG,prt,lun)	; <<< diagnostic >>>
SEARCH_FLAG = bytes(21)			; logical OR of all (384) minor
					  ; frame science flags            
CELESTIAL_OBJECTS = bytes(22)		; Bits 0 to 7 for:
						; Moon, Earth (31,53,90),
						; Jupiter, Mars, Saturn,
						; Galaxy
co = check_celestial_objects(CELESTIAL_OBJECTS,prt,lun)	; <<< diagnostic >>>
ATTITUDE_TYPE = bytes(23)		; 0=None, 1=SunEarth, 2=Definitive
at = check_attitude_type(ATTITUDE_TYPE,prt,lun)	; <<< diagnostic >>>
Galaxy_Cut_Angle = bytes(24)		; Degrees
gc = check_galaxy_cut(Galaxy_Cut_Angle,prt,lun)	; <<< diagnostic >>>
Planet_Cut_Angle = bytes(25)		; Degrees
pc = check_planet_cut(Planet_Cut_Angle,prt,lun)	; <<< diagnostic >>>
Moon_Cut_Angle = bytes(26)		; Degrees
mc = check_moon_cut(Moon_Cut_Angle,prt,lun)	; <<< diagnostic >>>
Earth_Cut_Angle = bytes(27:29)		; Degrees; array for 31, 53, 90 GHz
ec = check_earth_cut(Earth_Cut_Angle,prt,lun)	; <<< diagnostic >>>
;
; End of COBETRIEVE_HEADER section
; --------------------------------
;
; SPACECRAFT_INFORMATION starts in byte 30 and is 518 bytes long
; --------------------------------------------------------------
ATTITUDE_QUATERNION = fconvert(bytes(30:173),36,prt,lun); 9 quaternions per frame
							; to interpolate
ATTITUDE_QUATERNION = reform(ATTITUDE_QUATERNION,4,9)	; make 4 x 9 Array
ATTITUDE_RESIDUAL = reform(bytes(174:397),4,56)	; 4 x 56 Array of bytes  
						; Residuals for the
							; interpolator
SPACECRAFT_VELOCITY = fconvert(bytes(398:409),3,prt,lun); 3D vector; meters / sec
SPACECRAFT_POSITION = fconvert(bytes(410:421),3,prt,lun); 3D vector; Celestial coord J2000
EARTH_VELOCITY = fconvert(bytes(422:433),3,prt,lun)	; 3D vector; meters / sec
Orbit_Angle = fconvert(bytes(434:437),1,prt,lun)	; Degrees
if (verbose) then begin
      stat=execute(prt+", ' '")
      stat=execute(prt+", 'Spacecraft velocity (m/s): ',SPACECRAFT_VELOCITY")
      stat=execute(prt+", 'Spacecraft position: ', SPACECRAFT_POSITION")
      stat=execute(prt+", 'Earth velocity (m/s): ', EARTH_VELOCITY")
      stat=execute(prt+", 'Orbit angle (deg) = ', Orbit_Angle")
      stat=execute(prt+", ' '")
endif
ca0 = check_angle(Orbit_Angle, 'Orbit Angle',prt,lun)     	; <<< diagnostic >>>
DIPDUATH_Orbit_Var = fconvert(bytes(438:441),1,prt,lun)  	; D.U. 
DP28VRUA_Orbit_Var = fconvert(bytes(442:445),1,prt,lun)  	; D.U.
;
; MPA_BRACKET_THERM_31A  Major frame average:
MPA_BRACKET_THERMISTOR_31 = fconvert(bytes(446:449),1,prt,lun) ; DMT31MAA
;
; ABSOLUTE_PRT_53A and ABSOLUTE_PRT_90A Major frame averages:
CRITICAL_REGION_53 = fconvert(bytes(450:453),1,prt,lun) 	; DAP53MAA
CRITICAL_REGION_90 = fconvert(bytes(454:457),1,prt,lun) 	; DAP90MAA
;
; Fields from spacecraft archive NSB_SCDB_CS
; Torquer bar currents (modulo 64 values in 2-byte arrays):
AMABI = bytes(458:459)				; MMA A BAR CURRENT  
AMAXBI = bytes(460:461)				; MMA AX BAR CURRENT 
AMBBI = bytes(462:463)				; MMA B BAR CURRENT  
AMBXBI = bytes(464:465)				; MMA BX BAR CURRENT 
AMCBI = bytes(466:467)				; MMA C BAR CURRENT  
AMCXBI = bytes(468:469)				; MMA CX BAR CURRENT 
;
; Magnetometer field components (modulo 32 values in 6-byte arrays),
; 4 per major frame in NSB_SCDB_CS.  For each set of six values, the first
; value is the last value from NSB_SCDB_CS for the previous major frame,
; and the sixth value is the first value from NSB_SCDB_CS for the next
; major frame.
ATAMA = bytes(470:475)				; TAM A  FIELD   
ATAMAP = bytes(476:481)				; TAM AP FIELD(PERPENDICULAR)
ATAMAX = bytes(482:487)				; TAM AX FIELD  
ATAMB = bytes(488:493)				; TAM B  FIELD  
ATAMBP = bytes(494:499)				; TAM BP FIELD    
ATAMBX = bytes(500:505)				; TAM BX FIELD      
ATAMC = bytes(506:511)				; TAM C  FIELD    
ATAMCP = bytes(512:517)				; TAM CP FIELD    
ATAMCX = bytes(518:523)				; TAM CX FIELD
;
; Momentum and Reaction Wheels
; (For each set of six reaction wheel values, the first value is the
; last value from NSB_SCDB_CS for the previous major frame, and the
; sixth value is the first value from NSB_SCDB_CS for the next major
; frame)
; The following momentum wheel speeds are in RPM
AMW1WSP = bytes(524)		; ACS MOMENTUM WHEEL 1 SPEED
AMW1SPER = bytes(525)		; ACS MOMENTUM WHEEL 1 SPEED ERROR
AMW2WSP = bytes(526)		; ACS MOMENTUM WHEEL 2 SPEED
AMW2SPER = bytes(527)		; ACS MOMENTUM WHEEL 2 SPEED ERROR
; The following reaction wheel speeds (in RPM) are 6-byte arrays
ARWASP = bytes(528:533)		; ACS REACTION WHEEL A SPEED
ARWBSP = bytes(534:539)		; ACS REACTION WHEEL B SPEED
ARWCSP = bytes(540:545)		; ACS REACTION WHEEL C SPEED
SPARE = bytes(546:547)		; 2-byte array (Spares)
;
; End of SPACECRAFT_INFORMATION section
; -------------------------------------
;
; MOON_DATA starts in byte 548 and is 12 bytes long
; -------------------------------------------------
ANGLE_SPIN_MOON = fconvert(bytes(548:551),1,prt,lun) 	; Degrees
ca1 = check_angle(ANGLE_SPIN_MOON, 'ANGLE_SPIN_MOON',prt,lun); <<< diagnostic >>>
PHASE = fconvert(bytes(552:555),1,prt,lun)		     ; Radians
ca2 = check_angle_rad(PHASE, 'PHASE',prt,lun)	             ; <<< diagnostic >>>
POLARIZATION = fconvert(bytes(556:559),1,prt,lun)	     ; Radians
ca3 = check_angle_rad(POLARIZATION, 'POLARIZATION',prt,lun)  ; <<< diagnostic >>>
;
; End of MOON_DATA section
; ------------------------
;
; SIGNAL starts in byte 560 and is 6 x 232 bytes long (232 bytes / channel)
; -------------------------------------------------------------------------
; extract SIGNAL data from bytes array
signal_bytes = bytes(560:1951)
signal_bytes = reform(signal_bytes, 232, 6)
;
; make arrays to hold SIGNAL data
DIFFERENTIAL_TEMPERATURE = intarr(6,64)
TOTAL_POWER = intarr(6,4)
AVERAGE = replicate({avg, diff_temp:0., absolute_step:0., $
			  lock_in_amp_therm:0., local_osc_therm:0.},6)
DIODE_GAIN = fltarr(6)
MEAN_BASELINE = fltarr(6)
SPLINE_BASELINE = fltarr(6)
FAST_CTD_BASELINE = fltarr(6)
;
; store the data for 6 DMR channels
for i = 0, 5 do begin
	DIFFERENTIAL_TEMPERATURE(i,*) = convert(signal_bytes(0:127,i), 64,prt,lun)
	TOTAL_POWER(i,*) = convert(signal_bytes(128:135,i), 4,prt,lun)
	AVERAGE(i).diff_temp = fconvert(signal_bytes(136:139,i),1,prt,lun)
	AVERAGE(i).absolute_step = fconvert(signal_bytes(140:143,i),1,prt,lun)
	AVERAGE(i).lock_in_amp_therm = fconvert(signal_bytes(144:147,i),1,prt,lun)
	AVERAGE(i).local_osc_therm = fconvert(signal_bytes(148:151,i),1,prt,lun)
	DIODE_GAIN(i) = fconvert(signal_bytes(152:155,i),1,prt,lun)
	MEAN_BASELINE(i) = fconvert(signal_bytes(156:159,i),1,prt,lun)
	SPLINE_BASELINE(i) = fconvert(signal_bytes(160:163,i),1,prt,lun)
	FAST_CTD_BASELINE(i) = fconvert(signal_bytes(164:167,i),1,prt,lun)
endfor
FLAG = signal_bytes(168:231,*)	; (Minor-Frame (Science) Quality Flag)
				; 6 x 64 Byte array corresponding to
				; DIFFERENTIAL_TEMPERATURE array with
				; Bits 0 to 7 set for
				; Bad_data, Moon, Spike, Offscale,
				; Planet, Galaxy, (spare), Bad_att
;
; reconfigure the arrays into 3 frequencies x 2 channels per frequency
DIFFERENTIAL_TEMPERATURE = reform(DIFFERENTIAL_TEMPERATURE, 3, 2, 64)
TOTAL_POWER = reform(TOTAL_POWER, 3, 2, 4)
AVERAGE = reform(AVERAGE, 3, 2)
DIODE_GAIN = reform(DIODE_GAIN, 3, 2)
MEAN_BASELINE = reform(MEAN_BASELINE, 3, 2)
SPLINE_BASELINE = reform(SPLINE_BASELINE, 3, 2)
FAST_CTD_BASELINE = reform(FAST_CTD_BASELINE, 3, 2)
FLAG = reform(FLAG, 64, 3, 2)
;
; End of SIGNAL section
; ---------------------
;
; HOUSEKEEPING_DATA starts in byte 1952 and is 96 bytes long
; ----------------------------------------------------------
; From SUBCOM_99:
DSCRBY = bytes(1952)				; S/C REG BUS VOLT
DNESBI = bytes(1953)				; S/C NON-ES BUS CURR
DP28VRUA = bytes(1954)				; IPDU +28 V MON
DP5VMRUA = bytes(1955)				; IPDU +5 V MON 
DP15VRUA = bytes(1956)				; IPDU +15V MON 
DN15VRUA = bytes(1957)				; IPDU -15V MON 
DIPDUVMA = bytes(1958)				; IPDU BUSV MON 
DIPDUIMA = bytes(1959)				; IPDU BUSI MON 
DIICBIMA = bytes(1960)				; IPDU CONV BUSI
DIPHBVMA = bytes(1961)				; IPDU HEAT BUSV
DIPHBIMA = bytes(1962)				; IPDU HEAT BUSI
DPR10VMA = bytes(1963)				; IPDU PR1V MON 
DPR20VMA = bytes(1964)				; IPDU PR1V MON 
D31MIVMA = bytes(1965)				; 31 MPA/IPA PRVM
D531PRTA = bytes(1966)				; 53 GHz 1K PRT 
D901PRTA = bytes(1967)				; 90 GHz 1K PRT 
DIPDUTHA = bytes(1968)				; IPDU BOX THERM
DIPLMTHA = bytes(1969)				; IPWR LMAC THERM
DHPLMTHA = bytes(1970)				; HPWR LMAC THERM
DPRGTH1A = bytes(1971)				; PRE-REG 1 THERM
DPRGTH2A = bytes(1972)				; PRE-REG 2 THERM
DINTCTHA = bytes(1973)				; IPDU INT CONV T
D53NSVMA = bytes(1974)				; 53 GHz NSIV MON
D90MIVMA = bytes(1975)				; 90 MPA/IFA PRVM
D90NSVMA = bytes(1976)				; 90 GHz NSIV MON
DDEP5VMA = bytes(1977)				; DE +5V IV MON  
D31TCM1A = bytes(1978)				; 31 GHz TCVM #1 
D53TCM1A = bytes(1979)				; 53 GHz TCMV #1 
D90TCM1A = bytes(1980)				; 90 GHz TCMV #1 
DDEUTCMA = bytes(1981)				; DATA ELEC TCVM 
D31NSVMA = bytes(1982)				; 31 GHz NSIV MON
D53MIVMA = bytes(1983)				; 53 MPA/IFA PRVM
DP28VRUB = bytes(1984)				; +28 IPDU IV MON
DP5VMRUB = bytes(1985)				; +5 IPDU IV MON 
DP15VRUB = bytes(1986)				; +15 IPDU IV MON
DN15VRUB = bytes(1987)				; -15 IPDU IV MON
DIPDUVMB = bytes(1988)				; IPDU I BUSV MON
DIPDUIMB = bytes(1989)				; IPDU I BUSI MON
DIICBIMB = bytes(1990)				; IPDU IC BUSI MON
DIPHBVMB = bytes(1991)				; IPDU HBUSV MON  
DIPHBIMB = bytes(1992)				; IPDU HBUSI MON  
DPR10VMB = bytes(1993)				; IPDU PR1 OV MON 
DPR20MVB = bytes(1994)				; IPDU PR2 OV MON 
D31MIVMB = bytes(1995)				; 31 MPA/IFA PRVM 
D531PRTB = bytes(1996)				; 53 GHz 1K PRT  
D901PRTB = bytes(1997)				; 90 GHz 1K PRT  
DIPDUTHB = bytes(1998)				; IPDU BOX THERM 
DIPLMTHB = bytes(1999)				; IPWR LMAC THERM
DHPLMTHB = bytes(2000)				; HTRP LMAC THERM
DPRGTH1B = bytes(2001)				; PRE-REG 1 THERM
DPRGTH2B = bytes(2002)				; PRE-REG 2 THERM
DINTCTHB = bytes(2003)				; IPDU INT CONV T
D53NSVMB = bytes(2004)				; 53 GHz NSIV MON
D90MIVMB = bytes(2005)				; 90 MPA/IFA PRVM
D90NSVMB = bytes(2006)				; 90 GHz NSIV MON
DDEP5VMB = bytes(2007)				; DE +5V INSV MON
D31TCM2B = bytes(2008)				; 31 GHz TCVM #2 
D53TCM2B = bytes(2009)				; 53 GHz TCVM #2 
D90TCM2B = bytes(2010)				; 90 GHz TCVM #2 
DDEUTCMB = bytes(2011)				; DATA ELEC TCVM 
D31NSVMB = bytes(2012)				; 31 GHz NS IVM  
D53MIVMB = bytes(2013)				; 53 MPA/IFA PRVM
DIPDUATH = bytes(2014)				; IPDU-A BOX S/CT
DIPDUBTH = bytes(2015)				; IPDU-B BOX S/CT
DDEUABTH = bytes(2016)				; DEU-A BOX S/CT 
DDEUBBTH = bytes(2017)				; DEU-B BOX S/CT 
DI8NSCA = bytes(2018)				; NS CALIBRATE   
DRELA99A = bytes(2019:2023)			; 5 Byte array; RelayA group
;                 
; Byte 1:
;  Bit0 = DI7DEVRA  Bit1 = D690LIAA  Bit2 = DI590DSA  Bit3 = D453LIAA
;   DE VOLT REG      90 GHz LI AMP    90 GHz DK SW     53 GHz LI AMP
;  Bit4 = DI353DSA  Bit5 = D231LIAA  Bit6 = DI131DSA  Bit7 = DK831CCA
;   53 GHz DK SW     31 GHz LI AMP    31 GHz DK SW     31 GHz CH CNV
; Byte 2:
;  Bit0 = D7DETG5A  Bit1 = D690NSDA  Bit2 = DK590MIA  Bit3 = D453NSDA
;   DE TIM GEN +5    90 GHz NS DIS    90 MPA/IFA       53 GHz NS DIS 
;  Bit4 = DK353MIA  Bit5 = D231NSDA  Bit6 = DK131MIA  Bit7 = D853FMHA 
;   53 MPA/IFA       31 GHz NS DIS    31 MPA/IFA       53 FAIL M HTR  
; Byte 3:
;  Bit0 = D731FMHA  Bit1 = D690WBCA  Bit2 = D590CRCA  Bit3 = D453WBCA
;   31 FAIL M HTR    90 CASE T CNT    90 CRIT T CNT    53 CASE T CNT
;  Bit4 = D353CRCA  Bit5 = DH231CCA  Bit6 = D131CRCA  Bit7 = D1290LOA
;   53 CRIT T CNT    31 CASE T CNT    31 CRIT T CNT    90 GHz LO CNV
; Byte 4:
;  Bit0 = D1153LOA  Bit1 = D10DETCA  Bit2 = D990FMHA  Bit3 = DK12DECA
;   53 GHz LO CNV    DE TEMP CNT      90 FAIL M HTR    DATA ELEC CNV
;  Bit4 = D1131LOA  Bit5 = D1090CCA  Bit6 = DK953CCA  Bit7 = D31MODEA
;   31 GHz LO CNV    90 GHz CH CNV    53 GHz CH CNV    31A INST MODE
; Byte 5:
;  Bit0 = D53MODEA  Bit1 = D90MODEA  Bit2 = D90GAINA  Bit3 = D53GAINA
;   53A INST MODE    90A INST MODE    90 GHz GAIN      53 GHz GAIN
;  Bit4 = D31GAINA  Bit5 = D31BWTHA  Bit6 = DILMACOA  Bit7 = DHLMACOA
;   31 GHz GAIN      31 GHz BANDWD    LMAC OVTP OVR    HPWR OVTP OVR
;                 
D31NSSTA = bytes(2024)				; 31A NS STATUS
D53NSSTA = bytes(2025)				; 53A NS STATUS
D90NSSTA = bytes(2026)				; 90A NS STATUS
D31NSSTB = bytes(2027)				; 31B NS STATUS
D53NSSTB = bytes(2028)				; 53B NS STATUS
D31SCBTH = bytes(2029)				; 31 GHz BOX TH
D90NSSTB = bytes(2030)				; 90B NS STATUS
D53SCBTH = bytes(2031)				; 53 GHz BOX TH
D90SCBTH = bytes(2032)				; 90 GHz BOX TH
DI8NSCB = bytes(2033)				; NS CALIBRATE
DRELA99B = bytes(2034:2038)			; 5 Byte array; RelayB group
;                 
; Byte 1:
;  Bit0 = DI7DEVRB  Bit1 = D690LIAB  Bit2 = DI590DSB  Bit3 = D453LIAB
;   DE VOLT REG      90 GHz LI AMP    90 DICKE SW      53 GHz LI AMP
;  Bit4 = DI353DSB  Bit5 = D231LIAB  Bit6 = DI131DSB  Bit7 = DK831CCB
;   53 DICKE SW      31 GHz LI AMP    31 DICKE SW      31 GHz CH CNV
; Byte 2:
;  Bit0 = D7DETG5B  Bit1 = D690NSDB  Bit2 = DK590MIB  Bit3 = D453NSDB
;   DE TIM GEN +5    90 GHz NS DIS    90 MPA/IFA       53 GHz NS DIS 
;  Bit4 = DK353MIB  Bit5 = D231NSDB  Bit6 = DK131MIB  Bit7 = D853FMHB
;   53 MPA/IFA       31 GHz NS DIS    31 MPA/IFA       53 FAIL M HTR  
; Byte 3:
;  Bit0 = D731FMHB  Bit1 = D690WBCB  Bit2 = D590CRCB  Bit3 = D453WBCB
;   31 FAIL M HTR    90 CASE T CNT    90 CRIT T CNT    53 CASE T CNT
;  Bit4 = D353CRCB  Bit5 = DH231CCB  Bit6 = D131CRCB  Bit7 = D1290LOB
;   53 CRIT T CNT    31 CASE T CNT    31 CRIT T CNT    90 GHz LO CNV
; Byte 4:
;  Bit0 = D1153LOB  Bit1 = D10DETCB  Bit2 = D990FMHB  Bit3 = DK12DECB
;   53 GHz LO CNV    DE TEMP CNT      90 FAIL M HTR    DATA ELEC CNV
;  Bit4 = D1131LOB  Bit5 = D1090CCB  Bit6 = DK953CCB  Bit7 = D31MODEB
;   31 GHz LO CNV    90 GHz CH CNV    53 GHz CH CNV    31B INST MODE
; Byte 5
;  Bit0 = D53MODEB  Bit1 = D90MODEB  Bit2 = D90GAINB  Bit3 = D53GAINB
;   53A INST MODE    90A INST MODE    90 GHz GAIN      53 GHz GAIN
;  Bit4 = D31GAINB  Bit5 = D31BWTHB  Bit6 = DILMACOB  Bit7 = DHLMACOB
;   31 GHz GAIN      31 GHz BANDWD    LMAC OVTP OVR    HPWR OVTP OVR
;                 
Spare = bytes(2039)
;
; From  SUBCOM_98:
DRELA98 = bytes(2040)				; Relay group
;                 
; Bit0 = DIPIC28A        Bit1 = DIPIC28B        Bit2 = DHSC28VA      
;  +28V IPDU INT CONV A   +28V IPDU INT CONV B   +28V HEATER SWITCH A
; Bit3 = DHSC28VB        Bit4 = DIPPR28A        Bit5 = DIPPR28B
;  +28V HEATER SWITCH B   +28V IPDU PRE REG  A   +28V IPDU PRE REG  B
;                 
D31BOXSH = bytes(2041)				; 31 GHz BOX SURV HEATER
D53BOXSH = bytes(2042)				; 53 GHz BOX SURV HEATER
D90BOXSH = bytes(2043)				; 90 GHz BOX SURV HEATER
DDEUABSH = bytes(2044)				; DEU-A BOX SURV HEATER
DDEUBBSH = bytes(2045)				; DEU-B BOX SURV HEATER
DIPDUBSH = bytes(2046)				; IPDU BAY SURV HEATER
DBUSVOLT = bytes(2047)				; DMR BUS VOLTAGE
;
; End of HOUSEKEEPING_DATA section
; --------------------------------
;
; TOTAL LENGTH OF RECORD:   2048 BYTES
; TOTAL NUMBER OF FIELDS:    161
;
; The Validation section follows
;
; Plot the differential temperature arrays (64 bytes each x 3 frequencies
; x 2 channels per frequency)
mx = max(DIFFERENTIAL_TEMPERATURE)
mn = min(DIFFERENTIAL_TEMPERATURE)
tmin = float(maj_frame_index - 1)*64.
t = tmin + findgen(64)
if ((!d.name eq 'X') or (!d.name eq 'WIN')) then begin
  for i=0,2 do begin
	for j=0,1 do begin
		oplot,t,DIFFERENTIAL_TEMPERATURE(i,j,*),color=pcolor(i,j)
	endfor
  endfor
endif
;
; Print information for validation
;
stat=execute(prt+", ' '")
for i=0,2 do begin
case i of
0: freq = '31 GHz '
1: freq = '53 GHz '
2: freq = '90 GHz '
endcase
for j=0,1 do begin
case j of
0: chan = 'A channel'
1: chan = 'B channel'
endcase
heading = '***** ' + freq + chan + ' *****'
stat=execute(prt+", heading")
;
; calculate average of good differential temperatures for comparison with
; telemetered average
dt = reform(DIFFERENTIAL_TEMPERATURE(i,j,*))
good_dt = FLAG(*,i,j)
dt_index = where(good_dt eq 0)		; good data have FLAG = 0
if (dt_index(0) ge 0) then begin
	dt = dt(dt_index)
	n_good = n_elements(dt)
	avg = total(dt)/float(n_good)	; average of good values
;
;	print arrays if at least one datum is not good
	if (verbose) then begin
                stat=execute(prt+",'    Differential Temperature in DUs'")
		stat=execute(prt+",reform(DIFFERENTIAL_TEMPERATURE(i,j,*))")
                stat=execute(prt+",'    Quality Flags'")
		stat=execute(prt+",reform(FLAG(*,i,j))")
	endif
    stat=execute(prt+", 'Average of ',n_good,' good temperature values =',avg")
endif else begin
	stat=execute(prt+", 'NO GOOD VALUES!'")
endelse
stat=execute(prt+", 'Telemetered average = ', AVERAGE(i,j).diff_temp")
if (verbose) then begin
   stat=execute(prt+", 'Average absolute step = ',AVERAGE(i,j).absolute_step")
   stat=execute(prt+", 'Baseline values:', MEAN_BASELINE(i,j),SPLINE_BASELINE(i,j),"+$
    	"FAST_CTD_BASELINE(i,j),' (mean, spline, fast)")
   stat=execute(prt+", 'Diode gain = ', DIODE_GAIN(i,j),' DU per mK'")
   stat=execute(prt+", 'Total power in 8-second intervals: ',reform(TOTAL_POWER(i,j,*))")
   stat=execute(prt+", 'Average Lock-in Amplifier thermistor reading = ',AVERAGE(i,j).lock_in_amp_therm")
   stat=execute(prt+", 'Average Local Oscillator thermistor reading = ',AVERAGE(i,j).local_osc_therm")
   stat=execute(prt+", ' '")
endif
endfor
endfor
if (verbose) then begin
	stat=execute(prt+", 'ATTITUDE_QUATERNION array:'")
	stat=execute(prt+", ATTITUDE_QUATERNION")
	stat=execute(prt+", ' '")
endif
;
; End of Validation section
;
; A program control section follows.  One may stop after the present
; major frame, continue to the next frame, or skip over some frames.
;
if ((!d.name eq 'X') or (!d.name eq 'WIN')) then begin
  if (mn lt ymin or mx gt ymax) then begin	; pause to read diagnostics
	print,'Check diagnostics, then enter 1 to continue'
	read, contin
  endif
endif
if (saveframe) then begin
     framename='dmr_tod.' + strtrim(major_frame_number,2)
     print,'Saving and correcting major frame data into file= ',framename
     print,'   NOTE: a spawned IDL process will be temporarily created.'
     save,/variables,file=framename
     if !version.os eq 'vms' then cmd='delete ' else cmd='rm '
     if !version.os eq 'vms' then suffix="." else suffix=''
     openw,lun2,'tmp.pro',/get_lun
      printf,lun2,'restore,"'+framename
      printf,lun2,'delvar,temp_31a,temp_31b,temp_53a'
      printf,lun2,'delvar,temp_53b,temp_90a,temp_90b'
      printf,lun2,'spawn,"'+cmd+framename+suffix+'"'
      printf,lun2,'save,/variables,file="'+framename
      printf,lun2,'exit'
     free_lun,lun2
     spawn,'idl tmp.pro'
     print,'  '
     print,'Completed save.'
     print,'  '
     print,'  '
endif
if (savetemp) then begin
   for i=0,2 do begin
     for j=0,1 do begin
       differential_temperature(i,j,*)=(differential_temperature(i,j,*)- $
         mean_baseline(i,j))/diode_gain(i,j)
     endfor
   endfor
   temp_31a(*,maj_frame_index)=differential_temperature(0,0,*)
   temp_31b(*,maj_frame_index)=differential_temperature(0,1,*)
   temp_53a(*,maj_frame_index)=differential_temperature(1,0,*)
   temp_53b(*,maj_frame_index)=differential_temperature(1,1,*)
   temp_90a(*,maj_frame_index)=differential_temperature(2,0,*)
   temp_90b(*,maj_frame_index)=differential_temperature(2,1,*)
endif
if (every_frame) then goto, next_frame
print, 'Enter 1 for next frame, 0 to quit, any other number to skip frames'
read, contin
case contin of
1:	goto, next_frame
0:	goto, cleanup
else:	begin
		print, 'Skip how many frames?'
		read, nskip
                skip_bytes=bytarr(bytes_per_maj_frame*float(nskip))
                maj_frame_index=maj_frame_index+fix(nskip)
                print,'Skipping major frames ',MAJOR_FRAME_NUMBER+1,' to ',$
                         MAJOR_FRAME_NUMBER+fix(nskip)
                readu,unit,skip_bytes
		goto, next_frame
	end
endcase
;
; Close the file and end when finished
;
cleanup:
if (savetemp) then begin
 ncol=where(sum(temp_31a,0) ne 0)
 sz=size(ncol)
 t31a=fltarr(64,sz(1))
 if (ncol(0) ne -1) then t31a(*,*)=temp_31a(*,ncol)
 ncol=where(sum(temp_31b,0) ne 0)
 sz=size(ncol)
 t31b=fltarr(64,sz(1))
 if (ncol(0) ne -1) then t31b(*,*)=temp_31b(*,ncol)
 ncol=where(sum(temp_53a,0) ne 0)
 sz=size(ncol)
 t53a=fltarr(64,sz(1))
 if (ncol(0) ne -1) then t53a(*,*)=temp_53a(*,ncol)
 ncol=where(sum(temp_53b,0) ne 0)
 sz=size(ncol)
 t53b=fltarr(64,sz(1))
 if (ncol(0) ne -1) then t53b(*,*)=temp_53b(*,ncol)
 ncol=where(sum(temp_90a,0) ne 0)
 sz=size(ncol)
 t90a=fltarr(64,sz(1))
 if (ncol(0) ne -1) then t90a(*,*)=temp_90a(*,ncol)
 ncol=where(sum(temp_90b,0) ne 0)
 sz=size(ncol)
 t90b=fltarr(64,sz(1))
 if (ncol(0) ne -1) then t90b(*,*)=temp_90b(*,ncol)
 print,'Calculated temperatures are being saved in file= temps.dat'
 save,t31a,t31b,t53a,t53b,t90a,t90b,file='temps.dat'
endif
free_lun, unit
stat=execute(prt+", ' '")
stat=execute(prt+", 'Major frame index = ', maj_frame_index")
if (n_elements(logfile) ne 0) then free_lun,lun
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


