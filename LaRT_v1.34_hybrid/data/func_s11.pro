function func_s11,x,p
   np = n_elements(p)
   if (np eq 1) then begin
      scale = 1d0
      g     = p[0]
   endif else begin
      scale = abs(p[0])
      g     = p[1]
   endelse
   y = scale*0.5d0*(1d0-g^2)/(1d0+g^2-2d0*g*x)^(1.5d0)

   if (np eq 3) then begin
      scale2 = 1d0 - scale
      g2     = p[2]
   endif else if (np ge 4) then begin
      scale2 = abs(p[2])
      g2     = p[3]
   endif
   if (np gt 2) then begin
      y = y + scale2*0.5d0*(1d0-g2^2)/(1d0+g2^2-2d0*g2*x)^(1.5d0)
   endif
   if (np eq 5) then begin
      y = y + p[4]/0.5d0
      ;y = (1d0-p[4])*y + p[4]/0.5d0
   endif

   return,y
end
