pro mips70

;readcol, '/Users/jkrick/spitzer/mips/mips70/sources.reg', ra, dec ,format="A"
restore, '/Users/jkrick/idlbin/object.sav'
;a = where(object.mips70flux gt 400)

;plothyperz, a

!p.multi = [0, 4, 4]

ps_open, file = "/Users/jkrick/nep/mips70.acs.ps", /portrait, xsize = 8, ysize = 8,/color
!p.multi = [0, 4, 4]

ra = [265.19069,265.11907,265.03742,265.06217,264.98049,265.24772,265.22879,264.92418,264.90647,264.89065,264.9000,264.84121,265.08829,265.13558,264.77116,264.97061]

dec = [68.967547,68.944999,68.987134,68.98449,69.002995,69.033995,69.035476,69.027445,69.04931,69.072492,69.068861,69.06022,69.0516,69.058437,69.002295,69.100982]

acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits');
adxy, acshead, ra,dec , xcenter, ycenter
size = 20
a = [0,0.25,0.5,0.75,0,0.25,0.5,0.75,0,0.25,0.5,0.75,0,0.25,0.5,0.75]
b = [0,0,0,0,0.25,0.25,0.25,0.25,0.5,0.5,0.5,0.5,0.75,0.75,0.75,0.75]
c = [0.25,0.5,0.75,1.0,0.25,0.5,0.75,1.0,0.25,0.5,0.75,1.0,0.25,0.5,0.75,1.0]
d = [0.25,0.25,0.25,0.25,0.5,0.5,0.5,0.5,0.75,0.75,0.75,0.75,1,1,1,1]
for n = 0,n_elements(ra) - 1 do begin

;output the acs image of the object
   print, "working on n ", n

   acsdata = mrdfits('/Users/jkrick/hst/raw/wholeacs.fits', range=[ycenter(n) -size/0.05, ycenter(n)+size/0.05])
   plotimage, xrange=[xcenter(n) - size/0.05, xcenter(n)+ size/0.05],$
;                 yrange=[ycenter(n) -size/0.05, ycenter(n)+size/0.05], $
              bytscl(acsdata, min = -0.01, max = 0.1),$
              /preserve_aspect, /noaxes, ncolors=60, region = [a(n),b(n),c(n),d(n)], xmargin = [0.20],ymargin=[0,0]
   
   
endfor
ps_close, /noprint, /noid

ps_open, file = "/Users/jkrick/nep/mips70.ps", /portrait, xsize = 8, ysize = 8,/color
!p.multi = [0, 4, 4]

mipshead = headfits('/Users/jkrick/spitzer/mips/mips70/mosaic/combine/mosaic.fits');
adxy, mipshead, ra,dec , xcenter, ycenter

for n = 0,n_elements(ra) - 1 do begin

;output the acs image of the object
   print, "working on n ", n

   mipsdata = mrdfits('/Users/jkrick/spitzer/mips/mips70/mosaic/combine/mosaic.fits', range=[ycenter(n) -4, ycenter(n)+4])
   plotimage, xrange=[xcenter(n) - 4, xcenter(n)+ 4],$
;                 yrange=[ycenter(n) -size/0.05, ycenter(n)+size/0.05], $
              bytscl(mipsdata, min = -0.02, max = 0.5),$
              /preserve_aspect, /noaxes, ncolors=60, region = [a(n),b(n),c(n),d(n)], xmargin = [0.20],ymargin=[0,0]
   
   
endfor


ps_close, /noprint, /noid

end
