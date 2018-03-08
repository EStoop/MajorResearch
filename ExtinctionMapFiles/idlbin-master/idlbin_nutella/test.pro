pro test
size = 6
restore, '/Users/jkrick/idlbin/objectnew_swirc.sav'
keeper = 1269
fits_read, '/Users/jkrick/palomar/wirc/wirc_2008/k/coadd.fits',kdata,kheader
adxy,kheader, objectnew[keeper].ra,objectnew[keeper].dec , kxcenter,kycenter

print, 'x,y', kxcenter, kycenter
print, 'about to plot'
if objectnew[keeper].wirckmag gt 0 then begin ;and objectnew[keeper(n)].acsmag lt 90 
   print, 'about to plot'
   plotimage, xrange=[kxcenter - size/0.25,kxcenter+ size/.25],$
              yrange=[kycenter -size/.25,kycenter+size/.25], $
              bytscl(kdata, min =-43, max = 4300),title='k',$ ;
              /preserve_aspect, /noaxes, ncolors=60
;      xyouts, 0, -10., strcompress( string(objectnew[keeper(n)].ra) +'   '+ string(objectnew[keeper(n)].dec)),charthick = 3
endif 


end
