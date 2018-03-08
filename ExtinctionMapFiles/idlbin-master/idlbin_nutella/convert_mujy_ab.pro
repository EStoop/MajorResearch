function convert_mujy_ab, mujy

ab = 8.926 - 2.5*alog10(mujy*1E-6)

return, ab

end
