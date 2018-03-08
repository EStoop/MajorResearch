pro find_latent


aorname = ['r38763008','r37908480','r33220096','r29376000']

for j = 0, n_elements(aorname) - 1 do begin
  dirname = strcompress('/Users/jkrick/iracdata/S19INT/final/'+aorname[j]+'/ch1/bcd',/remove_all)
   CD, dirname                   ; change directories to the correct AOR directory
  command  =  "ls  " + dirname + "/*bimsk.fits > /Users/jkrick/irac_warm/latent/flag_latent/int_aug2011/bimsklist.txt"
  print, command
  spawn, command
  
  readcol,'/Users/jkrick/irac_warm/latent/flag_latent/int_aug2011/bimsklist.txt',fitsname, format = 'A', /silent
  
  for i =1, n_elements(fitsname) - 1 do begin
     
     fits_read, fitsname(i), data, header
     a = where(data eq 1024, latentcount)
     if latentcount gt 0 then begin
        print, 'got one',fitsname(i), latentcount
     endif

  endfor
  
endfor ; for each aor

  
end

