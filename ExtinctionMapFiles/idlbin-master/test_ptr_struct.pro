pro test_ptr_struct
restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp14/wasp14_phot_ch2.sav'
 aorname = ['r45838592', 'r45840128', 'r45841408', 'r45842176', 'r45842944', 'r45844480', 'r45845248', 'r45846016', 'r45846784', 'r45839104', 'r45840896', 'r45841664', 'r45842432', 'r45843200', 'r45844736', 'r45845504', 'r45846272', 'r45847040', 'r45839616', 'r45841152', 'r45841920', 'r45842688', 'r45844224', 'r45844992', 'r45845760', 'r45846528'] 

for a = 0, n_elements(aorname) - 1 do begin
   print, 'a' , a, n_elements(planethash[aorname(a),'xcen'])
   if a eq 0 then xcen = [planethash[aorname(a),'xcen']]
   if a gt 0 then xcen = [xcen, planethash[aorname(a),'xcen']]
endfor

print, 'all done', n_elements(xcen)

planethash = hash()
ra_ref = 208.92783
dec_ref = -32.159837
aorname = ['r45675264', 'r45675776'] ;ch2
for a = 0, 1 do begin
   planethash[aorname(a)] = HASH('ra',ra_ref, 'dec', dec_ref)
endfor

print, planethash.keys()
print, planethash[aorname(1)].keys()


testhash['aor1'] = HASH('Y',findgen(100),'X', findgen(100))
testhash['aor2'] = HASH('Y',findgen(5),'X', findgen(100))
print,testhash['aor1'].keys()

;print out an individual value
a = where(testhash['aor1', 'Y'] eq 1)
print, (testhash['aor1', 'Y'])(a)
;or
yarr = testhash['aor1', 'Y']
print, yarr(1)

;test binning using histogram
h = histogram(testhash['aor1', 'Y'], OMIN=om, binsize = 5, reverse_indices = ri)
print, 'omin', om, 'nh', n_elements(h)

for j = 0L, n_elements(h) - 1 do begin
   if (ri[j+1] gt ri[j]) then begin
;      print, 'binning together', n_elements((testhash['aor1', 'Y'])[ri[ri[j]:ri[j+1]-1]])
      meanclip, (testhash['aor1', 'Y'])[ri[ri[j]:ri[j+1]-1]], mean, sigma
;      print, mean
   endif
endfor
 
;test for looping over the AORs
st = ['aor1','aor2']
naors = 2

for n = 0, n_elements(st) -1 do begin
   print, 'inside fpr', testhash[st(n), 'Y']
endfor

help, testhash



end
