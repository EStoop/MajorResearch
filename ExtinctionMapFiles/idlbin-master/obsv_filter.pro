function obsv_filter,inlam,innulnu,boxcar=boxcar,lrange=lrange,$
                     acsf435w=acsf435w,acsf606w=acsf606w,acsf775w=acsf775w,$
                     acsf850lp=acsf850lp, $
                     irac1=irac1,irac2=irac2,irac3=irac3,irac4=irac4,$
                     irs16=irs16,mips24=mips24,mips70=mips70,mips160=mips160,$
                     pacs70=pacs70,pacs100=pacs100,pacs160=pacs160,$
                     spire250=spire250,spire350=spire350,spire500=spire500,$
                     scuba850=scuba850
  
; function that will calculate the observed lnu in a given filter
; window. Input is lambda (in microns), nulnu in say Lsun; output is lnu in Lsun/Hz
; Assumes filter curves are in microns and r_lambda

;READ FILTERS
;================================================================
COMMON FILTERCURVE, filters

if keyword_set(acsf435w) then  begin
  flb = filters[18].lamb 
  transp = filters[18].transp
  goto,jumpt
endif
if keyword_set(acsf606w) then  begin
  flb = filters[19].lamb 
  transp = filters[19].transp
  goto,jumpt
endif
if keyword_set(acsf775w) then  begin
  flb = filters[20].lamb 
  transp = filters[20].transp
  goto,jumpt
endif

if keyword_set(acsf850lp) then  begin
  flb = filters[21].lamb 
  transp = filters[21].transp
  goto,jumpt
endif

if keyword_set(irac1) then  begin
  flb = filters[0].lamb 
  transp = filters[0].transp
  goto,jumpt
endif

if keyword_set(irac2) then  begin
  flb = filters[1].lamb 
  transp = filters[1].transp
  goto,jumpt
endif

if keyword_set(irac3) then  begin
  flb = filters[24].lamb 
  transp = filters[24].transp
  goto,jumpt
endif

if keyword_set(irac4) then  begin
  flb = filters[25].lamb 
  transp = filters[25].transp
  goto,jumpt
endif

if keyword_set(irs16) then  begin
  flb = filters[26].lamb 
  transp = filters[26].transp
  goto,jumpt
endif

if keyword_set(mips24) then  begin
  flb = filters[27].lamb 
  transp = filters[27].transp
  goto,jumpt
endif

if keyword_set(mips70) then  begin
  flb = filters[28].lamb 
  transp = filters[28].transp
  goto,jumpt
endif

if keyword_set(pacs70) then  begin
  flb = filters[29].lamb 
  transp = filters[29].transp
  goto,jumpt
endif

if keyword_set(pacs100) then  begin
  flb = filters[30].lamb 
  transp = filters[30].transp
  goto,jumpt
endif

if keyword_set(mips160) then  begin
  flb = filters[31].lamb 
  transp = filters[31].transp
  goto,jumpt
endif

if keyword_set(pacs160) then  begin
  flb = filters[32].lamb 
  transp = filters[32].transp
  goto,jumpt
endif

if keyword_set(spire250) then  begin
  flb = filters[33].lamb 
  transp = filters[33].transp
  goto,jumpt
endif

if keyword_set(spire350) then  begin
  flb = filters[34].lamb 
  transp = filters[34].transp
  goto,jumpt
endif

if keyword_set(spire500) then  begin
  flb = filters[35].lamb 
  transp = filters[35].transp
  goto,jumpt
endif

if keyword_set(scuba850) then  begin
  flb = filters[36].lamb 
  transp = filters[36].transp
  goto,jumpt
endif


;DO INTEGRATION
;================================================================
jumpt: if keyword_set(boxcar) then begin
   flb = lrange[0] + findgen(101)*(lrange[1]-lrange[0])/100.
   transp = fltarr(101)+1.
endif else begin
   lrange = minmax(flb)
endelse

delta_lam = (lrange[1]-lrange[0])/100.
goodrng = where(inlam ge lrange[0] and inlam le lrange[1])

if (goodrng[0] ne -1) then begin
  ltz = where(innulnu lt 0.)
  if (ltz[0] ne -1) then begin
    ;print,'WARNING: SED has negative values.'
    innulnu[ltz] = 0.
  endif
  inlnu = innulnu/(3d14/inlam)

  newinlnu = interpol(inlnu,inlam,flb)
  notval = where((flb lt min(inlam)) or (flb gt max(inlam)))
  if (notval[0] ne -1) then newinlnu[notval] = 0.

  dfreq = 3d14/(flb-0.5*delta_lam)-3d14/(flb+0.5*delta_lam)
 
  totval = total(newinlnu*transp*dfreq)/total(transp*dfreq)

  return,totval[0]
endif else return,0.     ; if filter doesnt see the flux, it is 0

end
