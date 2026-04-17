function func_s12,x,p
   scale = p[0]
   ;y = scale*(1d0-x^2)/(1+x^2)
   ;y = scale*(1d0-x^2)/(1+abs(x)^6d0)
   ;y = scale*(1d0-(x-p[1])^2)
   y = scale*(1d0-x^2)

   return,y
end
