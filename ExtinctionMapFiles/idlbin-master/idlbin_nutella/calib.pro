; code to calibrate SDSS field
pro calib,imagename

readcol,'/Users/jkrick/palomar/LFC/SDSS/sdss.csv',run,rerun,camcol,field,obj,type,ra,dec,u,g,r,i,z,Err_u,Err_g,Err_r,Err_i,Err_z,FORMAT="A"

;print, u

image=readfits(imagename,hd)

;print, hd
exptime=sxpar(hd,'EXPTIME')
filter=sxpar(hd,'FILTER')
;image=image/exptime

exptime=1
filter ='i'

case 1 of
 (filter EQ 'g'):m=g
 (filter EQ 'u'):m=u
 (filter EQ 'r'):m=r
 (filter EQ 'i'):m=i
 (filter EQ 'z'):m=z
endcase

adxy,hd,ra,dec,x,y

gcntrd,image,x,y,xc,yc,3,/silent
print, x, xc
aper,image,xc,yc,flux,eflux,sky,skyerr,1,[8],[30,50],[-200,20000],/nan,/exact,/flux,/silent
instmag=-2.5*alog10(flux)
good=where((instmag EQ instmag) AND (m GT 12) and (m lt 20.))

print, "good", good
for n=0, 20, 1 do begin
      print, xc(n), yc(n),flux(n), sky(n)
endfor





result=ladfit(m(good),instmag(good))

outx=findgen(13)
outx=outx+10
outy=outx+result[0]

plot,m(good),instmag(good),psym=3,xrange=[10,20.5],yrange=[-20,0],color='008800'x ,ystyle=1

print,'First Pass:',result[0],result[1]

; now reject outliers
clean=m*result[1]+result[0]
diff=abs(clean-instmag)
good=where((diff LT 1.0)AND(m GT 15.)and (m lt 20.))
result=ladfit(m(good),instmag(good))
outy=outx+result[0]
print,'Outlier rejection:',result[0],result[1]

; now fix slope =1
clean=m+result[0]
diff=abs(clean-instmag)
good=where((diff LT 0.4)AND(m GT 15.)and (m lt 20.))
magoffset=result[0]-median(diff[good])
print,'Final Magnitude Offset:',magoffset
outy=outx +magoffset
plot,outx,outy,xrange=[10,20.5],yrange=[-20,0],color='000088'x,/noerase,ystyle=1
plot,m(good),instmag(good),psym=3,xrange=[10,20.5],yrange=[-20,0],/noerase,ystyle=1

; compute the sigma
clean=m+magoffset
diff=abs(clean-instmag)
good=where((diff LT 1)AND(m GT 12))
diff=clean-instmag 
k=moment(diff(good))
print,'Moments:',k

;stop

end
