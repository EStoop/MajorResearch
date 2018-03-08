;-----------------------------------------------------------------------------
; IDL routine GRT.PRO
;
; This IDL routine grabs DIRBE GRT data from the archive, converts it from DN's
; to Kelvins, and plots the results.
;
; Note: This routine labels GRT's by the physical location number assigned
;       to them by:
;
;        "Interface Control Drawing, Temperature Sensor Locations, 
;         DIRBE TU" provided by C. Woodhouse for the Phase 3 ITD test. 
;
; So-called "GID" #'s used by the STOL page displays differ from this 
; numbering system; to convert, add 1 to the GID to get the GRTID, i.e. 
; 
;              GRTID # = GID # + 1
;
;  GRTID #:            Location:
;  -------             --------
;    1			DA-1 3-Shooter
;    2			DA-1 6-Shooter
;    3			DA-2
;    4			DA-2
;    5			DA-3, Band 9
;    6			DA-4, Band 10
;    7			IRS Housing
;    8			IRS Baseplate
;    9			Optics Plate Near IRS
;   10			Optics Plate Near DA-3
;   11			Field Stop #1
;   12			Optics Plate  Near Bracket A
;   13			Optics Plate Near Chopper
;   14			Optics Plate, Center Under Baffle Assembly
;   15			Primary
;   16			Calibration Resistor
;
;
; 
;
; INPUTS:		DATA START TIME, STOP TIME
;
;
; Modification History:
;	09-Oct-1992	Kryszak	location of coefficients file was changed
;				as per SPR 10104
;
;-----------------------------------------------------------------------------
cleanplot
openw,2,'t1.dat'
;Initialize String Variables
start=''
stop=''
outdev=''
grlab=''

;-------------------------------------------------------------------------------

;Get GRT Start and Stop Times
read,'Enter GRT Start Time in VAX VMS or Julian Format ',start
print,' ',start
read,'Enter GRT Stop Time in VAX VMS or Julian Format ',stop
print,' ',stop

;Determine Output Device
read,'Enter Output Plot Device (Tek4014=0,VT240=3,Talaris=5)',outdev
;if outdev eq '' then set_plot,'vt240' else set_plot,outdev
if outdev eq '' then set_plot,'REGIS' else set_plot,outdev
if (!d.name eq 'X') then begin 
    device, retain=2
    set_plot,'X'
endif
if (!d.name eq 'WIN') then begin 
    device, retain=2
    set_plot,'WIN'
endif

;-------------------------------------------------------------------------------

;Get GRT Data From Archive
data=read_arcv('nbs_ctod.dagrt',start,stop)

;Determine Number of MJF Requested
mjfnumb=n_elements(data)/16

;-------------------------------------------------------------------------------

;Get Time Data From Archive

mjftime=read_arcv('nbs_ctod.datimas',start,stop)
time=findgen(mjfnumb)*(32./3600.)

;-------------------------------------------------------------------------------

; Set Up To Make Pretty Plots

!type=16
!fancy=1
!linetype=0
erase
!noeras=1
!ytitle='Deg K'
label=''


;-------------------------------------------------------------------------------

;Convert GRT DN's to Degrees Kelvin


;Read In JJP GRT Conversion Coefficients From a Data File
;January 25, 1989    J.J.P.'
;GRT CALIBRATION POINTS IN DEGREES KELVIN
dummy=fltarr(8,16)
coef=fltarr(16,8)
dirspec=getenv('CGIS_DATA')
filespec=dirspec+'grt_coefficients.dat'
openr,1,filespec
readf,1,'$(e24.16)',dummy
close,1
coef=transpose(dummy)




for j=0,15 do begin				; Main Loop Over All GRT's


;Apply JJP Polynomial Fit To GRT DN's

	temp=fltarr(mjfnumb)
	xx=fltarr(mjfnumb)
	lgxx=fltarr(mjfnumb)

       	xx=data(j,*)
	xx=xx(0:mjfnumb-1)
	xx=(xx+2017.0201)/.20275115
       	lgxx=alog10(xx)
       	pf = (((((((coef(j,7)*lgxx)+coef(j,6))*lgxx $
		+coef(j,5))*lgxx+coef(j,4))*lgxx    $
           	+coef(j,3))*lgxx+coef(j,2))*lgxx    $
  		+coef(j,1))*lgxx+coef(j,0)

       	temp=10 ^ pf


;Determine Plot Lables and Store GRT Data in GRT# Variables For Future Use

	if j eq '0' then begin
		grlab='Calibration Resistor'
		grt0=temp
	endif

	if j eq '1' then begin
		grlab='DA-1 3-Shooter'
		grt1=temp
	endif
	
	if j eq '2' then begin
		grlab='DA-1 6-Shooter'
		grt2=temp
	endif

	if j eq '3' then begin
		grlab='DA-2'
		grt3=temp
	endif

	if j eq '4' then begin
		grlab='DA-2'
		grt4=temp
	endif
	
	if j eq '5' then begin
		grlab='DA-3 Band 9'
		grt5=temp
	endif
	
	if j eq '6' then begin
		grlab='DA-3 Band 10'
		grt6=temp
	endif
	
	if j eq '7' then begin
		grlab='IRS Housing'
		grt7=temp
	endif
	
	if j eq '8' then begin
		grlab='IRS Baseplate'
		grt8=temp
	endif

	if j eq '9' then begin
		grlab='Optics Plate Near IRS'
		grt9=temp
	endif
	
	if j eq '10' then begin
		grlab='Optics Plate Near DA-3'
		grt10=temp
	endif
	

	if j eq '11' then begin
		grlab='Field Stop #1'
		grt11=temp
	endif
	
	if j eq '12' then begin
		grlab='Optics Plate Near Bracket A'
		grt12=temp
	endif
	
	if j eq '13' then begin
		grlab='Optics Plate Near Chopper'
		grt13=temp
	endif
	
	if j eq '14' then begin
		grlab='Optics Plate Under Baffle'
		grt14=temp
	endif
	
	if j eq '15' then begin
		grlab='Primary'
		grt15=temp
	endif
	


;-------------------------------------------------------------------------------

;Start Plots

sep=.03						; Plot Spacing
wid=.13						; Plot Width
z0=1.03						; Origin of Top Plot

;````````````````````````````````````````````````````````
if (j le 5) then begin
	set_viewport,.1,.9,z0-((j+1)*sep+(j+1)*wid),z0-((j+1)*sep+j*wid)
;	maxtempscale=8.0
;	set_xy,min(time),max(time),1.0,maxtempscale
	set_xy,min(time),max(time),min(temp),max(temp)
        printf,2,'min(time) ',min(time),'max(time) ',max(time)
        printf,2,'min(temp) ',min(temp),'max(temp) ',max(temp)
	!xtitle=''
	if (j eq 5) then !xtitle='Hours After ' + start
	!psym=3
	plot,time,temp
;	xyouts,.85*max(time),.25*maxtempscale,grlab 
	xyouts,.80*max(time), min(temp)+ .6* (max(temp)-min(temp)), grlab 
endif
if (j eq 5) then begin 
	erase 
endif

;````````````````````````````````````````````````````````
if (j ge 6) and (j le 11) then begin
	set_viewport,.1,.9,z0-((j-5)*sep+(j-5)*wid),z0-((j-5)*sep+(j-6)*wid)
	set_xy,min(time),max(time),min(temp),max(temp)
	!xtitle=''
	if (j eq 11) then !xtitle='Hours After ' + start
	!psym=3
	plot,time,temp
;	xyouts,.85*max(time),.25*maxtempscale,grlab 
	xyouts,.80*max(time), min(temp)+ .6* (max(temp)-min(temp)), grlab 
endif
if (j eq 11) then begin
	erase 
endif

;```````````````````````````````````````````````````````
if (j ge 12) and (j le 15) then begin
	set_viewport,.1,.9,z0-((j-11)*sep+(j-11)*wid),z0-((j-11)*sep+(j-12)*wid)
	set_xy,min(time),max(time),min(temp),max(temp)
	!xtitle=''
	if (j eq 15) then !xtitle='Hours After ' + start
	!psym=3
	plot,time,temp
;	xyouts,.85*max(time),.25*maxtempscale,grlab 
	xyouts,.80*max(time), min(temp)+ .6* (max(temp)-min(temp)), grlab 
endif
if (j eq 15) then begin
	erase 
endif
;-------------------------------------------------------------------------------

endfor


!noeras=0
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


