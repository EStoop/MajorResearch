
PRO fbcd_med,fBCDlis,OUTdir
;;
;;use input BCDlis to preserve NaNs
;
;; NAME:
;      fbcd_med
;
; PURPOSE:
; Subtracts off median of fbcd
;
; CALLING SEQUENCE:
;      fbcd_cleansm, INlist, OUTdir
;      OUTdir must exist before calling program
;
; INPUTS:
;      INLIST = Name of the list of fBCDs
; OUTPUTs:
;      OUTdir = output directory of cleaned fbcds
;               (output file names same)

if N_params() LT 2 then BEGIN
    print,'Syntax - fbcd_clean, INlist, OUTdir'
    return
endif
;;output file for diagnostic stats
outfile=OUTdir+'/outstats'
openw,lun,outfile,/GET_LUN

readcol,fBCDlis,infile,format='a'

FOR i=0,n_elements(infile)-1 DO BEGIN
   do_med,infile[i],OUTdir,medv,stdv 
   printf,lun,format='(a,f12.6,f12.6)',infile[i],medv,stdv
ENDFOR

free_lun, lun

END

pro do_med,IN_file,OUT_dir,med_out,std_out
;
im=readfits(IN_file,header)
im_out=im
cut=2.5
vec1=im[*,*]
aa=where(finite(vec1))
vec=vec1[aa]
var1=moment(vec)
std1=var1(1)^0.5
avg1=var1(0)
med1=median(vec,/even)
;CLEAN ONCE around median to find noise and median level
bb=where((vec lt med1+cut*std1) and (vec gt med1-cut*std1)) 
vec2=vec[bb]
var2=moment(vec2)
std2=var2(1)^0.5
avg2=var2(0)
med2=median(vec2,/even)
;rms=std2, median=med2
up=med2+cut*std2
lo=med2-cut*std2
;print,lo,up,med2,std2
;
med_out=med2
std_out=std2

;output image
im_out=im-med_out
writefits,OUT_dir+'/'+IN_file,im_out,header

END
