pro read_pitch_roll_yaw
; 

infile='/Users/jkrick/external/initial_drift/ctab_data.6378.ecsv'
savefile = '/Users/jkrick/external/initial_drift/ctab_data.6378.ecsv.sav'
nmax=99999999L

pitch=dblarr(nmax)
roll=dblarr(nmax)
yaw=dblarr(nmax)
pitchb=pitch

time=dblarr(nmax)

dn=0.d0
eu=0.d0
t=0.d0
di=0L

line=' '
openr, 1, infile

n=0L

while not eof(1) do begin
   readf, 1, line
   if(strmid(line,8,1) eq 'T') then begin
      reads, strmid(line,24,4), pid, format='(i4)'
      ibad=strpos(line,'?')
      if(ibad lt 8) then begin
         fpos=strpos(line,',U')
         upos=strpos(line,',F')
         if(upos gt 0) then fpos=upos
         qs=strmid(line,29,fpos-29)

         case pid of

            2049: begin
               reads, qs, di, eu,t
               pitch(n)=eu
            end
            2050: begin
               reads, qs, di, eu,t
               yaw(n)=eu
            end
            2847: begin
               reads, qs, di, eu,t
               roll(n)=eu
            end
            2848: begin
               reads, qs, di, eu,t
               pitchb(n)=eu
               time[n]=t
               n=n+1L
            end


         endcase

      endif
   endif
endwhile

close, 1

print, 'n: ', n

roll=roll[0:n-1]
pitch=pitch[0:n-1]
yaw=yaw[0:n-1]
pitchb=pitchb[0:n-1]
time=time[0:n-1]

outfile=infile+'.sav'

save, /variables, filename=outfile



end

