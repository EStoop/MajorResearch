

nmax=99999999
rw1=dblarr(nmax)
rw2=rw1
rw3=rw1
rw4=rw1
time=dblarr(nmax)

dn=0.d0
eu=0.d0
t=0.d0

names=['HD 7924', 'HD 7924', 'Hat-P 2', 'RV Transit 16', 'RV Tranist 5',$
'RV Transit 7', 'RV Transit 22','RV Transit 9','HD 209458', 'XO 3', $
       'XO 3', 'RV Tranist 9','r48072704','r48072960']

startt=[1004594183.d0,1041429073.d0,1062885213.d0,1063530855.d0,$
        1043516340.d0,1047527851.d0,1061298497.d0,1062545924.d0,$
        1062437205.d0,1053584269.d0,1052373562.d0,1062587169.d0,$
        1057752123.d0,1058005982.d0]

endt=[1004622741.d0,1041455531.d0,1062918804.d0,1063561093.d0,$
      1043557330.d0,1047563506.d0,1061317983.d0,1062586968.d0,$
      1062467480.d0,1053614416.d0,1052416649.d0,1062628212.d0,$
      1057783740.d0,1058037598.d0]

reqkey=[44605184,46981888,46477312,48815872,46917120,48407296,48823552,$
        48408320,48014336,46468352,46486016,48408064,48072704,$
        48072960]

remarks=['on',' below','on','below','on','on','on','below','on','below',$
         'below','below','NP bubble','NP bubble']

nplots=n_elements(startt)



line=' '
;openr, 1, 'ctab_data.21480.ecsv'
;openr, 1, 'ctab_data.19709.ecsv'
openr, 1, 'ctab_data.17534.ecsv'

n=0L

while not eof(1) do begin
   readf, 1, line
   if(strmid(line,8,1) eq 'T') then begin
      reads, strmid(line,25,3), pid, format='(i3)'
      ibad=strpos(line,'?')
      if(ibad lt 8) then begin
         fpos=strpos(line,',F')
         qs=strmid(line,29,fpos-29)


         reads, qs, dn, eu, t

         time(n)=t

         case pid of
            
            718: rw1(n)=eu
            719: rw2(n)=eu
            720: rw3(n)=eu
            721: begin
               
               rw4(n)=eu
               n=n+1
            end

         endcase

      endif
   endif
endwhile

close, 1

print, 'n: ', n

time=time[0:n-1]
rw1=rw1[0:n-1]
rw2=rw2[0:n-1]
rw3=rw3[0:n-1]
rw4=rw4[0:n-1]

time=time-1.d9
startt=startt-1.d9
endt=endt-1.d9

set_plot, 'PS'

for k=0,nplots-1 do begin

   device, /encapsulated, file=string(reqkey[k],format='(i8)')+'.eps'
   plot, time, rw1,xr=[startt[k]-3600.*12.,endt[k]+3600.*12.], xtitle='SCLK-1.d9', ytitle='RPM', yr=[-2000.,2000], title=names[k]+'  '+remarks[k]
   oplot, time, rw2, linestyle=1
   oplot, time, rw3, linestyle=2
   oplot, time, rw4, linestyle=3

   oplot, [startt[k],startt[k]], [-5000.,5000.]
   oplot, [endt[k],endt[k]], [-5000.,5000.]

   device, /close

endfor


end

