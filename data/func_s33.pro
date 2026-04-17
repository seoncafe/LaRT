function func_s33,x,p
   ;scale = p[0]
   ;y = scale*2d0*x/(1+x^2)
   ;y = 2d0*x/(1+p[0]*x^2)
   y = 2d0*x/(1+p[0]*sqrt(abs(x)))
   return,y
end
