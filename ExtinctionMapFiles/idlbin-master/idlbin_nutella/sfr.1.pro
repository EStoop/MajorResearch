pro sfr

;device, true=24
;device, decomposed=0
colors = GetColor(/load, Start=1)

!p.multi = [0, 0, 1]
ps_open, file = "/Users/jkrick/ZPHOT/sfr.ps", /portrait, xsize = 4, ysize = 4,/color

restore, '/Users/jkrick/idlbin/object.sav'
;restore,'/Users/jkrick/bin/bc03/bc03.sav'
;restore,'/Users/jkrick/bin/bc03/lambda.sav'
restore, '/Users/jkrick/chary01/rrcsedsforevol_w_yunradio.save'

print,"started at "+systime()

;find some galaxies to start with which have all optical, irac, and mips24
n = 0
for i = 0, n_elements(object.ra) - 1 do begin

   if object[i].umaga gt 1 and object[i].umaga lt 30 then begin
      if object[i].gmaga gt 1 and object[i].gmaga lt 30 then begin
         if object[i].rmaga gt 1 and object[i].rmaga lt 30 then begin
            if object[i].imaga gt 1 and object[i].imaga lt 30 then begin
               if object[i].flamjmag gt 1 and object[i].flamjmag lt 30 then begin
                  if object[i].irac1 gt 0 and object[i].irac1 ne 99 then begin
                     if object[i].irac2 gt 0 and object[i].irac2 ne 99 then begin
                        if object[i].irac3 gt 0 and object[i].irac3 ne 99 then begin
                           if object[i].irac4 gt 0 and object[i].irac4 ne 99 then begin
                              if object[i].flux gt 0 and object[i].flux ne 99 then begin
                                 if object[i].zphot lt 0.2 then begin
                                    
                                    print, "got 1: ",object[i].umaga, object[i].gmaga,object[i].rmaga,object[i].imaga,object[i].flamjmag,object[i].irac1,object[i].irac2,object[i].irac3,object[i].irac4,object[i].flux
                                    n = n + 1
                                 endif
                              endif

                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
      endif
   endif


endfor
print, "n objects ", n
print,"Finished at "+systime()

ps_close, /noprint, /noid


end


