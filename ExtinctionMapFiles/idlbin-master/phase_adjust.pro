pro phase_adjust, planetname, chname, apradius
;adjust phases of secondary eclipse observations to not wrap at 0.5

  planetinfo = create_planetinfo()

  for i = 0, n_elements(planetname)  - 1 do begin
     if chname(i) eq '2' then aorname= planetinfo[planetname(i), 'aorname_ch2'] else aorname = planetinfo[planetname(i), 'aorname_ch1'] 
     savename = strcompress( '/Users/jkrick/irac_warm/pcrs_planets/'+planetname(i)+'/'+planetname(i) + '_phot_ch'+chname(i)+'_'+string(apradius(i))+'.sav',/remove_all)

     restore, savename
     for a = 0, n_elements(aorname) - 1 do begin
        phase = planethash[aorname(a),'phase']        
        n = where(phase lt 0, ncount)
        if ncount gt 0 then begin
           phase(n) = phase(n) + 1
           planethash[aorname(a),'phase']   = phase
        endif
        print, 'test start and end phase',  (planethash[aorname(a),'phase'] )[0],  (planethash[aorname(a),'phase'] )[n_elements(phase) - 1]

     endfor
     
     save, planethash, filename=savename
     print, 'saving planethash', savename
     print, '-----------------------------------------'
endfor

end


;phase_adjust, [ 'WASP-13b', 'WASP-15b', 'WASP-16b', 'WASP-16b', 'WASP-38b', 'WASP-62b', 'WASP-62b', 'HAT-P-22','HAT-P-22'], [ '2', '2', '1', '2', '2', '1', '2', '1','2'],  [2.25, 1.75, 2.0, 2.25, 2.25, 2.25, 2.5, 1.75,2.25]
