pro wheel_speed_read
  t1 = systime(1)

  close, /all
  ;read in the data
; FIELD2 = P-0718, 719, 720, 721
;FIELD4 = speeds
;FIELD5 = sclk

  ;read_csv doesn't work on this file (maybe too large, maybe formatting)
  filename = '/Users/jkrick/irac_warm/irac_wheel_speeds/ctab_data.17534.ecsv';test.ecsv'
  openr, 1, filename

  n=0L
  line=' '
  nmax=99999999
  rw1=dblarr(nmax)
  rw2=rw1
  rw3=rw1
  rw4=rw1
  time=dblarr(nmax)
  
  dn=0.d0
  eu=0.d0
  t=0.d0

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
;------------------------------------------------------------------------

  ;read in the second set of momentum data
  filename = '/Users/jkrick/irac_warm/irac_wheel_speeds/ctab_data.18123.ecsv';test.ecsv'
  openr, 2, filename

  n=0L
  line=' '
  nmax=99999999
  m1=dblarr(nmax)
  m2=m1
  m3=m1
  m4=m1
  timem=dblarr(nmax)
  
  dn=0.d0
  eu=0.d0
  t=0.d0

  while not eof(2) do begin
     readf, 2, line
     if(strmid(line,8,1) eq 'T') then begin
        reads, strmid(line,25,3), pid, format='(i3)'
        ibad=strpos(line,'?')
        if(ibad lt 8) then begin
           fpos=strpos(line,',F')
           qs=strmid(line,29,fpos-29)
           s = strsplit( qs, ',',/extract)
           timem(n)=s(1)
           
           case pid of
            
              051: m1(n)=s(0)
              052: m2(n)=s(0)
              053: m3(n)=s(0)
              054: begin
                 
                 m4(n)=s(0)
                 n=n+1
              end
              
           endcase
           
        endif
     endif
  endwhile
  
  close, 1
  
  print, 'n: ', n
  
  timem=timem[0:n-1]
  m1=m1[0:n-1]
  m2=m2[0:n-1]
  m3=m3[0:n-1]
  m4=m4[0:n-1]
;------------------------------------------------------------------------

  save, time, rw1, rw2, rw3, rw4, timem, m1, m2, m3, m4, filename = '/Users/jkrick/irac_warm/irac_wheel_speeds/speeds.sav'

 

  
  print, 'time check', systime(1) - t1

end
