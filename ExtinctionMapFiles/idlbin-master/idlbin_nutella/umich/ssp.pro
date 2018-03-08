PRO ssp

close,/all

openr, lun, "/n/Godiva7/jkrick/galaxev/bc03/models/Padova1994/salpeter/test1", /get_lun

myage = fltarr(200)
mycev = fltarr(200)
mycne = fltarr(200)
karr = fltarr(200)
ekarr = fltarr(200)
i = 0
crap = "hello"
WHILE (NOT EOF(lun)) DO BEGIN
 
    READF, lun,z, junk, age, junk, crf, cne, cev, junk,junk,junk,ekcorr,kcorr
    myage[i] = z
    mycev[i] = cev
    mycne[i] = cne
    karr[i] = kcorr
    ekarr[i] =ekcorr
    i = i + 1
ENDWHILE

myage = myage[0:i-1]
mycev=mycev[0:i-1]
karr = karr[0:i-1]
ekarr = ekarr[0:i-1]
mycne = mycne[0:i-1]

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A3888/final/ssp.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plot, myage, mycev, thick = 3, xtitle = "redshift", $
 ytitle = "V-R color"

oplot, myage, mycev - karr, thick = 3
oplot, myage, mycev - ekarr, thick = 3
;oplot, myage, mycev - ekarr + karr, thick = 3
oplot, myage, karr
oplot, myage, mycne, thick = 1, linestyle = 4
;oplot, myage, mycne-karr , thick = 1, linestyle = 4
;oplot, myage, mycne-ekarr, thick = 1, linestyle = 4

close, lun
free_lun, lun
device, /close
set_plot, mydevice




END
;openr, lun, "/n/Godiva7/jkrick/galaxev/bc03/models/Padova1994/salpeter/test2", /get_lun

;myage = fltarr(200)
;mycev = fltarr(200)
;i = 0
;crap = "hello"
;WHILE (NOT EOF(lun)) DO BEGIN
; 
;    READF, lun,z, junk, age, junk, crf, cne, cev, junk,junk,junk,junk,kcorr
;    myage[i] = z
;    mycev[i] = cne
;    i = i + 1
;ENDWHILE
;
;myage = myage[0:i-1]
;mycev=mycev[0:i-1]
;
;oplot, myage , mycev, thick = 3, linestyle = 3
;close, lun
;free_lun, lun
;
;openr, lun, "/n/Godiva7/jkrick/galaxev/bc03/models/Padova1994/salpeter/bc2003_lr_m62_salp_ssp.2color", /get_lun
;
;bcage = fltarr(250)
;bccev = fltarr(250)
;i = 0
;crap = "hello"
;WHILE (NOT EOF(lun)) DO BEGIN
; 
;    READF, lun,logage, junk, junk,junk,color, junk, junk,junk,junk, junk,junk,junk, junk,junk,junk
;    bcage[i] = logage
;    bccev[i] = color
;    i = i + 1
;ENDWHILE
;
;bcage = bcage[0:i-1]
;bccev=bccev[0:i-1]
;;print, bcage, bccev
;plot, bcage, bccev, thick = 3, linestyle = 2, xrange=[7,10.4], xstyle = 1, yrange = [-0.2, 1.4], ystyle =1
;
;close, lun
;free_lun, lun
;openr, lun, "/n/Godiva7/jkrick/galaxev/bc03/models/Padova1994/salpeter/test.1color", /get_lun
;
;bcage = fltarr(250)
;bccev = fltarr(250)
;i = 0
;crap = "hello"
;WHILE (NOT EOF(lun)) DO BEGIN
; 
;    READF, lun,logage, junk, junk,junk,color, junk, junk,junk,junk, junk,junk,junk, junk,junk,bv
;    bcage[i] = logage
;    bccev[i] = bv
;    i = i + 1
;ENDWHILE
;
;bcage = bcage[0:i-1]
;bccev=bccev[0:i-1]
;oplot, bcage , bccev, thick = 3, linestyle = 3
;
;close, lun
;free_lun, lun
;
;
