;******************************************************************************
;The function DIFF_LUMFUNC takes an absolute magnitude, a value of Mstar, and 
;a value of alpha, and returns the value of the differential luminosity 
;function, with a normalization of 1.0.
;******************************************************************************
FUNCTION DIFF_LUMFUNC,xmag,Mstar,alpha

xmag=double(xmag)
Mstar=double(Mstar)
alpha=double(alpha)

firstexponent=double(0.4)*double(mstar-xmag)*double(1.0+alpha)

secondexponent=double(0.4)*double(Mstar-xmag)

blah1=double(10.0^(firstexponent))

blah2=double(exp(double(-1.0*10.0^(secondexponent))))

RETURN,double(0.4*alog(10.0)*blah1*blah2)

END
