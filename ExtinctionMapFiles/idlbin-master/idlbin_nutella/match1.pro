pro match1

close, /all

numoobjects = 32000
colorarr = fltarr(numoobjects)
magarr = fltarr(numoobjects)
FITS_READ, '/Users/jkrick/palomar/lfc/coadd_r.fits',rdata, rheader

rgalaxy = replicate({object, rxcenter:0D,rycenter:0D,rfluxa:0D,rmaga:0D,rmagerra:0D,rfwhm:0D,risoarea:0D},numoobjects)
ggalaxy = replicate({gobject, gxcenter:0D,gycenter:0D,gfluxa:0D,gmaga:0D,gmagerrb:0D,gfwhm:0D,gisoarea:0D},numoobjects)
irgalaxy = replicate({irobject, irxcenter:0D,irycenter:0D,tmassJ:0D,tmassh:0D,tmassk:0D,irac1:0D,irac2:0D,irac3:0D, irac4:0D},numoobjects)
galaxy= replicate({what, rra:0D, rdec:0D, g_r:0D, rmag:0D, gmag:0D}, numoobjects)

i = 0

openr, lunr, "/Users/jkrick/palomar/LFC/catalog/SExtractor.r.cat", /get_lun
WHILE (NOT EOF(lunr)) DO BEGIN
    READF, lunr, o, xwcsr, ywcsr, xcenterr, ycenterr, fluxautor, magautor, magerrr, fwhmr, isoarear, fluxmaxr

    if (magautor LT 50) and  (fluxmaxr lt 76) then begin
       rgalaxy[i] ={object, xwcsr, ywcsr,fluxautor, magautor, magerrr,fwhmr, isoarear}
       i = i + 1
;       printf, outlun, xcenterr, ycenterr
    ENDIF

 endwhile
rgalaxy = rgalaxy[0:i - 1]
print, "there are ",i," r galaxies"
close, lunr
free_lun, lunr


l = 0
openr, lung, "/Users/jkrick/palomar/LFC/catalog/SExtractor.g.cat", /get_lun
WHILE (NOT EOF(lung)) DO BEGIN
    READF, lung, o, xwcsg, ywcsg, xcenterg, ycenterg, fluxautog, magautog, magerrg, fwhmg, isoareag, fluxmaxg

    if (fluxautog gt 0 )  and (fluxmaxg lt 20) then begin
       ggalaxy[l] ={gobject, xwcsg, ywcsg,fluxautog, magautog, magerrg,fwhmg, isoareag}
;       print, "wprking on galaxy", l, xcenterg, ycenterg, xwcsg, ywcsg

       l = l + 1
    ENDIF

 endwhile
ggalaxy = ggalaxy[0:l - 1]
print, "there are ",l," g galaxies"
close, lung
free_lun, lung

l = 0
openr, luni, "/Users/jkrick/palomar/LFC/catalog/combined_catalog.prt", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, wcsra, wcsdec, xcenter, ycenter, tmassJa,tmassha,tmasska,irac1a,irac2a,irac3a, irac4a

   irgalaxy[l] ={irobject, wcsra, wcsdec, tmassJa,tmassha,tmasska,irac1a,irac2a,irac3a, irac4a}
   l = l + 1
   
 endwhile
irgalaxy = irgalaxy[0:l - 1]
print, "there are ",l," ir galaxies"
close, luni
free_lun, luni


adxy, rheader, irgalaxy.irxcenter, irgalaxy.irycenter,junk1, junk2
openw, junklun, "/Users/jkrick/palomar/lfc/catalog/irgal.txt", /get_lun
for i = 0, n_elements(junk1) - 1 do   printf, junklun, junk1(i), junk2(i)
close, junklun
free_lun, junklun
; read the input data catalogs
;readcol,'gems2.tbl',gcntr,gra,gdec,g1,g2,g3,g4,format="A"
;readcol,'swire.tbl',scntr,sra,sdec,s1,s2,s3,s4,s5,format="A"

; read the GEMS coverage map to define good data
;ch1cov=readfits('musyc/mosaic_ch1_cov.fits',ch1hd)
;adxy,ch1hd,gra,gdec,x,y
;x=round(x)
;y=round(y)
;gcov1=ch1cov[x,y]

; create initial arrays
r=n_elements(rgalaxy.rxcenter)
g=n_elements(irgalaxy.irxcenter)
rmatch=fltarr(r)
gmatch=fltarr(g)
rmatch[*]=-999
gmatch[*]=-999

; define "good" swire data
;sgood=fltarr(s)
;sgood[*]=10000000
;;sgood(where((s1 EQ s1) and (s2 EQ s2) and (s1 GT 5) AND (s2 GT 7)))=1
;sgood[*]=1

; define "good" GEMS data
;ggood=fltarr(g)
;ggood[*]=0
;ggood(where((gcov1 GT 10) and (g1 EQ g1) and (g1 GT 0.5)))=1
;gmatch(where(ggood NE 1))=-99

print,'Matching  r to ir'
print,"Starting at "+systime()
dist=rmatch
dist[*]=0

;for i=0l,g-1 do begin
for i=0l,r-1 do begin

;   dist=sphdist(irgalaxy[i].irxcenter,irgalaxy[i].irycenter,rgalaxy.rxcenter,rgalaxy.rycenter,/degrees)
   dist = sphdist(rgalaxy[i].rxcenter,rgalaxy[i].rycenter,irgalaxy.irxcenter,irgalaxy.irycenter,/degrees)
   sep=min(dist,ind)
;   if (sep LT 0.0005) then gmatch[i]=ind
    if (sep LT 0.0005) then rmatch[i]=ind
  
endfor

print,"Finished at "+systime()

;save,filename='ch1.dat'
;matched=where(gmatch GT 0)
;nonmatched = where(gmatch lt 0)
matched=where(rmatch GT 0)
nonmatched = where(rmatch lt 0)
print, n_elements(matched), " objects matched"
print, n_elements(nonmatched), " objects not matched"
ratio=irgalaxy(matched).irac1/rgalaxy(gmatch[matched]).rfluxa
plot,irgalaxy(matched).irac1,ratio,/xlog,psym=3

openw, outlun, "/Users/jkrick/palomar/lfc/catalog/rnonmatched.txt", /get_lun
openw, outlunr, "/Users/jkrick/palomar/lfc/catalog/rmatched.txt", /get_lun


;adxy, rheader, irgalaxy(nonmatched).irxcenter, irgalaxy(nonmatched).irycenter, irxcen, irycen
;for i = 0, n_elements(irxcen) - 1 do   printf, outlun, irxcen(i), irycen(i)

adxy, rheader, rgalaxy(nonmatched).rxcenter, rgalaxy(nonmatched).rycenter, rxcen, rycen
for i = 0, n_elements(rxcen) - 1 do   printf, outlun, rxcen(i), rycen(i)

;adxy, rheader, irgalaxy(matched).irxcenter,irgalaxy(matched).irycenter, irxcen2, irycen2
;for i = 0, n_elements(irxcen2) - 1 do printf, outlunr, irxcen2(i), irycen2(i)

adxy, rheader, rgalaxy(matched).rxcenter, rgalaxy(matched).rycenter, rxcen2, rycen2
for i = 0, n_elements(rxcen2) - 1 do printf, outlunr, rxcen2(i), rycen2(i)

close, outlunr
free_lun, outlunr
close, outlun
free_lun, outlun

undefine, rgalaxy
undefine, ggalaxy
undefine, galaxy
undefine, irgalaxy

end
