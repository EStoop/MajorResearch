pro makefinder
!p.multi = [0,0, 1]
; 2 brightest gals in each of 2 best clusters

ra = [264.68341,264.69266,264.89816,264.89246 ]
dec=[69.039215,69.045792,69.073517,69.068474]

findobject, ra, dec

catnumber = [1491, 1987, 14723, 2003]

;plothyperz, catnumber, '/Users/jkrick/nep/clusters/brightestclustergals.ps'


;now make some finder charts.
;print image around the sources, maybe 2.5' on each side.
ps_open, filename='/Users/jkrick/nep/clusters/brightfinders.ps',/portrait,/square,/color
redcolor = FSC_COLOR("Red", !D.Table_Size-2)

size = 150
acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits');
adxy, acshead,ra,dec , xcenter, ycenter

for n = 0, n_elements(xcenter) - 1 do begin

   acsdata = mrdfits('/Users/jkrick/hst/raw/wholeacs.fits', range=[ycenter(n) -size/0.05, ycenter(n)+size/0.05])
   plotimage, xrange=[xcenter(n) - size/0.05, xcenter(n)+ size/0.05],$
;;                 yrange=[ycenter(n) -size/0.05, ycenter(n)+size/0.05], $
              bytscl(acsdata, min = -0.01, max = 0.01),title='acs',$
              /preserve_aspect, /noaxes, ncolors=60 ;, position=[0,0.5,0.5,1.0]
   tvcircle, 200, xcenter(n), ycenter(n), /data, thick = 2, color = redcolor

endfor


ps_close, /noprint,/noid


end
