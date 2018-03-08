pro acsthumbnails, keeper, filenameplot

ps_open, filename=filenameplot,/portrait,/square,/color, xsize=5, ysize=5

restore, '/Users/jkrick/idlbin/objectnew.sav'

acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits');
adxy, acshead, objectnew[keeper].ra,objectnew[keeper].dec , xcenter, ycenter

size = 5


for n = 0,n_elements(keeper) - 1 do begin

;output the acs image of the objectnew
   print, "working on n ", n

 
   if objectnew[keeper(n)].acsmag gt 0  and ycenter(n) - size/0.05 gt 0 and ycenter(n) + size/0.05 lt 20300 then begin;and objectnew[keeper(n)].acsmag lt 90 
      acsdata = mrdfits('/Users/jkrick/hst/raw/wholeacs.fits', range=[ycenter(n) -size/0.05, ycenter(n)+size/0.05])
      plotimage, xrange=[xcenter(n) - size/0.05, xcenter(n)+ size/0.05],$
;;                 yrange=[ycenter(n) -size/0.05, ycenter(n)+size/0.05], $
                 bytscl(acsdata, min = -0.01, max = 0.1),title='acs',$
                 /preserve_aspect, /noaxes, ncolors=60;, position=[0,0.5,0.5,1.0]
;      xyouts, xcenter(n)- 0.6*size/0.05, -16., strcompress( string(objectnew[keeper(n)].ra) +'   '+ string(objectnew[keeper(n)].dec)),charthick = 3
   endif else begin
      ;need a placeholder
      plot, findgen(10), findgen(10), thick=0, xthick=0,ythick=0, charthick=0, charsize=0, xstyle=4, ystyle=4, psym=3
      xyouts, 4, 4, 'acs', charthick=3
   endelse


   
endfor


;-----------------------------------------------------------------------

ps_close, /noprint, /noid


end
