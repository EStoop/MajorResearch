pro mips70_flux_correct

;ran this on jan 12.  don't run again without changing the mips70 photometry

restore, '/Users/jkrick/idlbin/objectnew.sav'

a = where(objectnew.mips70flux gt 0)
objectnew[a].mips70flux = objectnew[a].mips70flux * 1.25

;save, objectnew, filename='/Users/jkrick/idlbin/objectnew.sav'

end
