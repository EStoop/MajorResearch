pro cosmic_slithelp
close, /all

ps_open, filename='/Users/jkrick/palomar/cosmic/slitmasks/mask11.finder.ps',/portrait,/square,/color
colors = GetColor(/load, Start=1)
restore, '/Users/jkrick/idlbin/object.sav'

maskname = '/Users/jkrick/palomar/cosmic/slitmasks/mask11.txt'
;center values
ch=17.D
cm=40.D;9.;9.D;41.D;39.D;39.D;41.D;39.D;40.D
cs=07.994D;4.637;8.758D;17.012D;35.994;6.097D;03.512;15.895D;05.168D;09.490D
dd=69.D
dm=05.D;59.;3.D;03.D;03.D;9.D;59.D;54.D;54.D;50.D
ds=39.29D;06.33D;37.21D;40.40D;36.98D;3.47D;13.18;46.52D;46.27D;52.03D

racenter = ((cs/60. + cm)/60. + ch) * 15.
deccenter = (ds/60. + dm)/60. + dd

;build up a redshift distribution
openw, zlun, '/Users/jkrick/palomar/cosmic/slitmasks/zdist.txt', /get_lun, /append

openw, outlun, '/Users/jkrick/palomar/cosmic/slitmasks/bmasklist12.txt', /get_lun
;readin the values which are in the current mask.
readcol,maskname ,donenum, junk, junk, junk,junk, junk, junk,junk, junk, junk,junk ,format="A"


k = 0
star = 0
true = 0
galaxynum = fltarr(n_elements(donenum))
starnum = fltarr(100)
h = 0.D
m = h
s = h
dd = h
dm = h
ds = h

;read in the master list, only output those objects which do not already have a slit assaigned to them
openr, lun, '/Users/jkrick/palomar/cosmic/slitmasks/bmasklist9.txt', /get_lun

WHILE (NOT EOF(lun)) DO BEGIN
   READF, lun, masternum, h, m, s, dd, dm, ds, mag
   for i = 0, n_elements(donenum) - 1 do begin
      if masternum eq donenum(i) then begin
         
         ;match to real catalog
         ra = ((s /60. + m )/60. + h ) * 15.
         dec = (ds /60. + dm )/60. + dd 

         m=n_elements(ra)
         ir=n_elements(object.ra)
         
         irmatch=fltarr(ir)
         mmatch=fltarr(m)
         irmatch[*]=-999
         mmatch[*]=-999
         dist=irmatch
         dist[*]=0

         dist=sphdist( ra, dec,object.ra,object.dec,/degrees)
         sep=min(dist,ind)
         
          ; galaxy, not star
         if donenum(i) gt 150  then begin
            true = 1
            galaxynum[i] = ind
            printf, zlun, object[ind].zphot
         endif else begin       ;star
            starnum[star] = ind
            star = star + 1
 ;           print, donenum(i), ind, object[ind].ra, object[ind].dec, object[ind].rmaga, mag
         endelse
         
      endif
      
   endfor
   
   if true ne 1 then printf, outlun, format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',masternum, h, m, s, dd, dm, ds, mag
   true = 0
ENDWHILE

starnum = starnum[0:star-1]

;----------------------------------------------------------
; make a finder chart
;make the stars a different color

fits_read, '/Users/jkrick/palomar/lfc/coadd_r.fits', data, head

;display an image around the center given   15' x 6' = 900x360" = 18000x7200 acspix = 5000x2000
;12'x4' = 720x240" = 4000x1330 palomar pix
adxy, head, racenter, deccenter, xcen, ycen

plotimage, xrange=[xcen - 2500 ,xcen+2500],yrange=[ycen-1000, ycen+1000], bytscl(data, min = -0.5, max = 0.5) ,$; bytscl(data, min=1,max=3),$
 /preserve_aspect, /noaxes, ncolors=8

;draw circles around the slitted objects.  

;
adxy, head,object[galaxynum].ra, object[galaxynum].dec,x, y
tvcircle, 20, x, y, /data, thick = 2, color = colors.red

adxy, head,object[starnum].ra, object[starnum].dec,x, y
tvcircle, 20, x, y, /data, thick = 2, color = colors.green

close, outlun
free_lun, outlun
close, lun
free_lun, lun

ps_close, /noprint,/noid

printf, zlun, "last"
close, zlun
free_lun, zlun
end


;           har[j] = h
;            mar[j] = m
;            sar[j] = s
;            ddar[j] = dd
;            dmar[j] = dm
;            dsar[j] = ds
;            j = j + 1
;         endif else begin
;            harr[k] = h
;            marr[k] = m
;            sarr[k] = s
;            ddarr[k] = dd
;            dmarr[k] = dm
;            dsarr[k] = ds
;            k = k + 1
;         endelse


;har = har[0:j-1]
;mar = mar[0:j-1]
;sar = sar[0:j-1]
;ddar = ddar[0:j-1]
;dmar = dmar[0:j-1]
;dsar = dsar[0:j-1]
;
;harr = harr[0:j-1]
;marr = marr[0:j-1]
;sarr = sarr[0:j-1]
;ddarr = ddarr[0:j-1]
;dmarr = dmarr[0:j-1]
;dsarr = dsarr[0:j-1]


;ra = ((sarr/60. + marr)/60. + harr) * 15.
;dec = (dsarr/60. + dmarr)/60. + ddarr
;adxy, head,ra, dec,x, y


;har = fltarr(100)
;mar = har
;sar = har
;ddar = har
;dmar = har
;dsar = har
;j = 0
;harr = fltarr(100)
;marr = harr
;sarr = harr
;ddarr = harr
;dmarr = harr
;dsarr = harr
