;-------------------------------------------------------------------------------------------
; takes a list of input images and a list of coords, checks for non-NaN
;  coverage in the images and outputs a flag
PRO COVFLAG,imagelist,positionlist

readcol,imagelist,images,format="a"
nfiles=n_elements(images)

readcol,positionlist,cntr,ra,dec,format="a"
npos=n_elements(cntr)

coverage=intarr(npos)
coverage[*]=0

; start looping through the input images
for i=0l,nfiles-1 do begin

print,images[i]
; for speed, just read the header
bcd_hd=headfits(images[i])
naxis1=sxpar(bcd_hd,'NAXIS1')-2
naxis2=sxpar(bcd_hd,'NAXIS2')-2
crval1=sxpar(bcd_hd,'CRVAL1')
crval2=sxpar(bcd_hd,'CRVAL2')

; get the x and y positions of all points relative to current image
adxy,bcd_hd,ra,dec,x,y
x=round(x)
y=round(y)

; check to see if they are on image
onimage=where(((x LT naxis1) AND (x GT 2) AND (y LT naxis2) AND (y GT 2)),count)
print,count

; begin the loop for checking coverage
if count NE 0 then begin

; read the actual image
bcd=readfits(images[i],bcd_hd)

for j=0l,count-1 do begin
state=0
if (bcd(x(onimage(j)),y(onimage(j))) EQ bcd(x(onimage(j)),y(onimage(j)))) then state=state+1
if (bcd(1+x(onimage(j)),y(onimage(j))) EQ bcd(1+x(onimage(j)),y(onimage(j)))) then state=state+1
if (bcd(-1+x(onimage(j)),y(onimage(j))) EQ bcd(-1+x(onimage(j)),y(onimage(j)))) then state=state+1
if (bcd(x(onimage(j)),1+y(onimage(j))) EQ bcd(x(onimage(j)),1+y(onimage(j)))) then state=state+1
if (bcd(x(onimage(j)),-1+y(onimage(j))) EQ bcd(x(onimage(j)),-1+y(onimage(j)))) then state=state+1

if (state GT 2) then coverage(onimage(j))=1


; end the loop over overlapping objects
endfor
; end the branch if there are overlaps on image
   endif

; end the loop over images
   endfor

; dump the output
openw,unit,'out.dat',/get_lun

for i=0l,npos-1 do begin
   printf,unit,cntr(i),ra(i),dec(i),coverage(i)
endfor

free_lun,unit





END
