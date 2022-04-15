function func_s34,x,p
   theta = acos(x)/!dtor
   pcir  = p[0]
   s     = p[1]
   se    = p[2]
   th0   = p[3]
   ex    = p[4]
   ;cc    = cos((theta + s*3.13d0*theta*exp(-7d0*theta/180d0))*!dtor)
   ;cc    = cos((theta + s*3.13d0*theta*exp(-se*(theta-th0)/180d0))*!dtor)
   cc    = cos((theta + s*3.13d0*theta*exp(-se*(theta-th0)^ex/180d0))*!dtor)
   y     = pcir*(1d0-cc^2)/(1d0+cc^2)
   return,y
end
