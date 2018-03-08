pro clustercmd

restore, '/Users/jkrick/idlbin/objectnew.sav'
ps_open, file = "/Users/jkrick/nep/clusters/mips24/clustercmd.ps", /portrait, xsize = 4, ysize = 4,/color
!P.charthick=3
!P.thick=3
!X.thick=3
!Y.thick=3

canddec = [69.045052,  69.06851,  69.087766]
candra = [264.68365, 264.89228, 264.83053]

sep = 0.019
openw, outlunred, '/Users/jkrick/nep/clusters/cluster24_24.reg', /get_lun
printf, outlunred, 'fk5'
;printf, outlunred, 'global color=yellow font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source'

for i = 0, n_elements(candra) -1 do begin
   goodall = where(sphdist(objectnew.ra, objectnew.dec, candra[i], canddec[i], /degrees) lt sep and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.rmag gt 0 and objectnew.rmag lt 23 and objectnew.acsmag gt 0 and objectnew.acsmag lt 90 $
and objectnew.mips24mag gt 0 and objectnew.mips24mag lt 90 )

   print, n_elements(goodall)
   for rc=0, n_elements(goodall) -1 do  printf, outlunred, 'circle( ', objectnew[goodall[rc]].ra, objectnew[goodall[rc]].dec, ' 3")'
   
   plot, objectnew[goodall].rmag,  objectnew[goodall].rmag - objectnew[goodall].acsmag ,psym = 2, xrange=[12,30], xstyle=1, yrange=[-2,4]
endfor
close, outlunred
free_lun, outlunred

ps_close, /noprint, /noid


end
