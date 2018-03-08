pro denscontour

fits_read, '/Users/jkrick/umich/icl/A4059/fullr.wcs.fits', data, header
readcol,'/Users/jkrick/umich/icl/A4059/member.txt',x,y,a,a,format="A"


xbinsize=1
ybinsize=1
xmax = 1000
ymax= 1000
xmin = 0
ymin = 0

density = hist_2d(x,y,bin1=xbinsize,bin2=ybinsize,max1=xmax,max2=ymas,min1=xmin,min2=ymin)

;density is the result from HIST_2D

maxvalue=max(density,maxpoint)
maxcoord=ARRAY_INDICES(density,maxpoint)
max_x=xmin+maxcoord[0]*xbinsize+xbinsize/2
max_y=ymin+maxcoord[1]*ybinsize+ybinsize/2
print, 'The maximum density is: ', maxvalue
print, 'It occurs at coordinates:'
print, max_x, max_y

tot_histo=total(density)
histoden=histogram(density)
n_histo=n_elements(histoden)
sum=0
k=0

while sum lt 0.9*tot_histo do begin

	sum = sum + histoden[n_histo-1-k] * (maxvalue-k)
	k = k+1

endwhile

maxlevel=histoden[n_histo-1-k] * (maxvalue-k)

contour, density,c_colors=[100,150,200],levels=maxlevel
xyouts, max_x, max_y, '+', /data, color=100

end
