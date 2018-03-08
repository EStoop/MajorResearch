FUNCTION conv, x, P

return, fft(devauc(x,P),1)
;return, fft((fft(devauc(x,P), -1)*fft(mymoffat(x), -1)),1)

end
