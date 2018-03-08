pro load_name_to_raw_correspondence_table,pnametab
;input lookup table is name_to_raw_correspondence_table.txt
; makes structure with corresponding names and data type

; skip 7 records to get to data in table.
; names are enclosed in single quotes -- use this to parse
;
in_names=strarr(800)
raw_names=in_names
raw_types=in_names
s=' '
; name correspondence structure
;    struct   local pname   raw FITS     data type
pnametab = {pcorres, localname:'', rawname:'', ptype:''}

close,1
openr,1,'/local/d1/s1/idl/pro_atlo/name_to_raw_correspondence_table.txt'
for i=1,7 do begin
readf,1,s,format='(a)'
endfor
ix=0
while not EOF(1) do begin
readf,1,s,format='(a)'
s2=strsplit(s,/extract)	;assumed separated by spaces or tabs
;print,s2
in_names(ix)=s2(0)
raw_names(ix)=s2(1)
raw_types(ix)=s2(2)
ix=ix+1
endwhile
in_names=in_names(0:ix-1)
raw_names=raw_names(0:ix-1)
raw_types=raw_types(0:ix-1)
close,1
pnametab=replicate({pcorres},ix)
pnametab.localname=in_names
pnametab.rawname=raw_names
pnametab.ptype=raw_types
end
