pro get_dewrap_2,fnam,pstr,im_out=im_out,headinfoonly=headinfoonly,makelog=makelog,$
  badlog=badlog,silent=silent,ierr=ierr,header=head
;*********************************************************************************
;reads in IRAC raw FITS does dewrap2  outputs correct floating image and structure
;  containing useful header parameter values
;  NOTE: this version doesn't know how to find the patter value, so if sets pattern=0
;BG 4-2003 from get_dewrap.pro: added capability to return headerinfo only and
;  pass the name of the log file for bad files
;  Also added the integers to the pstr structure for rdmod (=ARDMOD) and fl (=AFLSTAT)
;
;INPUTS
;  fnam  path to input raw IRAC FITS file (required)
;  headinfoonly (optional): set if don't want to read the image and dewrap it
;  makelog (optional) if set will print error messages to badlog (if supplied) or 
;		'bad_files_get_dewrap.log' if badlog argument is not supplied
;  badlog (optional) file name of log to record missing or erroneous header parameters
;  silent (optional) set this if don't want error messages printed to stdout
;
;OUTPUTS
;  pstr  structure with useful header info. entire header as string array
;  im_out (optional) floating image processed by get_dewrap
;  ierr   (optional) =0 if no errors,  =1 if error(s) occurred
;  head (optional) entire header as string array
;***********************************************************************************
do_dewrap=0
if keyword_set(makelog) then lg=1 else lg=0
if lg and arg_present(badlog) then badlog_file=badlog else badlog_file='bad_files_get_dewrap.log'
get_lun,lu_badlog
openw,lu_badlog,badlog_file	;make the log file anyway --delete it at end if unwanted
if keyword_set(silent) then vb=1 else vb=0

;load parameter name lookup table  (will replace this someday)
load_name_to_raw_correspondence_table,pnametab
if keyword_set(headinfoonly) then begin
  do_dewrap=0
  head=headfits(fnam,/silent)
endif else begin
  if arg_present(im_out) then begin
    im=readfits(fnam,head,/silent)
    do_dewrap=1
  endif else begin
    if vb then print,' must supply im_out=argument to return an image'
    stop
  endelse
endelse
;
;  Structure containing useful header stuff and stuff derived from the header
pstr=create_struct(name='pstr','cbad',' ','fnam',' ','rdmod',' ','exptime',0.0,$
	'frame',0.0,'fowl',0,'wait',0,'ps',0,'barrel',0,'weaside',' ','aintbeg',0L,'atimeend',0L,$
	'sclk',0L,'tx',0,'fl',' ','patt',0,'i_rdmod',0,'i_fl',0,'ich',0)
; note that the "rdmod" and "fl" tags are strings: they contain the 2-character hex values
;
make_header_structure,head,hdr	;used internally
tot_succ=0
isuccess=0
ich=sxpar(head,'CHNLNUM',count=count,/silent)	;MIPL channel identifier
if ich lt 0 or ich gt 4 then begin
  isuccess=1
  if vb then print,' ich =',ich,' out of range in file: ',fnam
  if lg then printf,lu_badlog,' ich =',ich,' out of range in file: ',fnam
  ich=5
endif
tot_succ=tot_succ+isuccess
achanid=-1
get_head_param2,'achanid',hdr,pnametab,cval,achanid,lu_badlog,isuccess	;channel number
tot_succ=tot_succ+isuccess
afowlnum='-1'
ifowlnum=0
get_head_param2,'afowlnum',hdr,pnametab,afowlnum,ifowlnum,lu_badlog,isuccess ;Fowler number
tot_succ=tot_succ+isuccess
apedsig=-1
get_head_param2,'apedsig',hdr,pnametab,capedsig,apedsig,lu_badlog,isuccess ;PedSig indicator
tot_succ=tot_succ+isuccess
awaitper='-1'
iwaitper=-1
get_head_param2,'awaitper',hdr,pnametab,awaitper,iwaitper,lu_badlog,isuccess ;num waits
tot_succ=tot_succ+isuccess
aweaside=-1
get_head_param2,'aweaside',hdr,pnametab,cval,aweaside,lu_badlog,isuccess	;WE side
tot_succ=tot_succ+isuccess
abarrel=-1
get_head_param2,'abarrel',hdr,pnametab,cval,abarrel,lu_badlog,isuccess	;Barrel shift
tot_succ=tot_succ+isuccess
atxstat=-1
get_head_param2,'atxstat',hdr,pnametab,cval,atxstat,lu_badlog,isuccess	;TCAL lamps status
tot_succ=tot_succ+isuccess
aflstat=-1
get_head_param2,'aflstat',hdr,pnametab,cval,aflstat,lu_badlog,isuccess	;FCAL lamps status
tot_succ=tot_succ+isuccess
if (aflstat ge 0) and (aflstat lt 256) then $
   cfl=strupcase(string(aflstat,format='(z2.2)')) else cfl='-1'
aintbeg=0L
get_head_param2,'aintbeg',hdr,pnametab,caintbeg,aintbeg,lu_badlog,isuccess ;time of integ begin
tot_succ=tot_succ+isuccess
atimeend=0L
get_head_param2,'atimeend',hdr,pnametab,catimeend,atimeend,lu_badlog,isuccess ;time of integ end
tot_succ=tot_succ+isuccess
;
;result = SXPAR( Hdr, Name, [ Abort, COUNT=, COMMENT =, /NoCONTINUE  ])
sclk=0L
isuccess=0
sclk=sxpar(head,'H0122D00',count=countsclk,/silent)	;LONG value of sclk
if countsclk eq 0 then begin
  if lg then printf,lu_badlog,format='(2a)','H0122D00 sclk not found by sxpar in header'
  tot_succ=tot_succ+1
endif
rdmod=-1
isuccess=0
rdmod=sxpar(head,'A0617D00',count=countrdmod,/silent)	;LONG value of sclk
if countrdmod eq 0 then begin
  if lg then printf,lu_badlog,format='(2a)','A0617D00 READMODE not found by sxpar in header'
  crdmod='-1'
  tot_succ=tot_succ+1
endif else begin
  if (rdmod ge 0) and (rdmod lt 256) then $
    crdmod=strupcase(string(rdmod,format='(z2.2)')) else crdmod='-1'
endelse
if tot_succ gt 0 then begin
  if lg then printf,lu_badlog,format='(2a)',' ERRORS above are for file: ',fnam
  if lg then printf,lu_badlog,' '
endif
errs=tot_succ
;
;assign header items
weaside='U'
if aweaside eq '-32768' then weaside='A'
if aweaside eq '     0' then weaside='B'
;bit0 of RDMOD is 0 for FULL, 1 for SUBARRAY
;other bits indicate PED/SIG or FOW status of each DSP
subar=rdmod mod 2
case subar of
0: begin
  frametime=0.2*(2*ifowlnum+iwaitper)
  exptime=0.2*(ifowlnum+iwaitper)
  end
1: begin
  frametime=0.01*(2*ifowlnum+iwaitper)
  exptime=0.01*(ifowlnum+iwaitper)
  end
else: begin
  frametime=0.0
  exptime=0.0
  end
endcase
cframetime=string(frametime,format='(f7.2)')
cexptime=string(exptime,format='(f7.2)')
;fnamoz=string(fnamo,format='(a45)')
if errs gt 0 then cbad='*' else cbad=' '
pstr.cbad=cbad
pstr.fnam=fnam
pstr.rdmod=crdmod
pstr.fowl=afowlnum
pstr.ps=apedsig
pstr.exptime=exptime
pstr.frame=frametime
pstr.wait=iwaitper
pstr.barrel=abarrel
pstr.weaside=weaside
pstr.aintbeg=aintbeg
pstr.atimeend=atimeend
pstr.sclk=sclk
pstr.tx=atxstat
pstr.fl=cfl
pattern=0
pstr.patt=pattern
pstr.i_rdmod=rdmod
pstr.i_fl=aflstat
pstr.ich=achanid
if do_dewrap then dewrap2,im,achanid,abarrel,afowlnum,apedsig,pattern,im_out
free_lun,lu_badlog
if not(lg) or (not(arg_present(badlog)) and errs eq 0) then spawn,'rm ./bad_files_get_dewrap.log'
if errs gt 0 then ierr=1 else ierr=0
end
