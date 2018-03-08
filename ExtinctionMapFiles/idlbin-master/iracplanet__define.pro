function iracplanet::init
  self.xcen = ptr_new(/allocate)
  self.ycen = ptr_new(/allocate)
  return, 1
end

pro iracplanet::set,value
  if n_elements(value) ne 0 then *(self.xcen) = value
  return
end

function iracplanet::get, value
  if n_elements(*(self.xcen)) ne 0 then value = *(self.xcen)
  return, value
end


function iracplanet::cleanup
  if ptr_valid(self.xcen) then ptr_free, self.xcen
  if ptr_valid(self.ycen) then ptr_free, self.ycen
  return,0
end


pro iracplanet::phot
;here is where the photometry would happen


end


pro iracplanet__define

  void = {iracplanet, xcen:ptr_new(), ycen:ptr_new()}

end

