;fname = 'mueller_Lyalpha.dat'
fname = 'mueller_matrix.dat'
readcol,fname,wavl,albedo,hgg,nang,skip=1,numline=1
wavl   = wavl[0]
albedo = albedo[0]
hgg    = hgg[0]
nang   = nang[0]
readcol,fname,cost,s11,s12,s33,s34,skip=3

;rayfun = 0.5d0*(1d0-hgg^2)/(1d0+hgg^2-2d0*hgg*cost)^(1.5d0)

ncost = n_elements(cost)
err   = dblarr(ncost)+1d0

!p.multi=[0,2,3]

xr    = [-1,1]
xr2   = [0.5d0,1]
ylog  = 0
quiet = 1
parm0 = [-0.4,0.2]
parms = mpfitfun('func_s12',cost,s12/s11,err,parm0,yfit=yfit,quiet=quiet)
print,'S12   paramter',parms
 plot,cost,s12/s11,xr=xr,/xst
oplot,cost,yfit,color=cgcolor('red')

 plot,cost,s12,xr=xr2,/xst
oplot,cost,yfit*s11,color=cgcolor('red')

;;s33_func = 2d0*cost/(1d0+abs(cost)^2)
;;s33_func = 2d0*cost/(1d0+abs(cost+0.03))
s33_func = 2d0*cost/(1d0+sqrt(abs(cost)))
parm0 = [1.0]
parms = mpfitfun('func_s33',cost,s33/s11,err,parm0,yfit=yfit,quiet=quiet)
print,'S33   paramter',parms
 plot,cost,s33/s11,xr=xr,/xst
oplot,cost,yfit,color=cgcolor('red')
oplot,cost,s33_func,color=cgcolor('green')

 plot,cost,s33,xr=xr2,/xst
oplot,cost,yfit*s11,color=cgcolor('red')

parm0 = [0.39d0,1.0]
parms = mpfitfun('func_s34',cost,s34/s11,err,parm0,yfit=yfit,quiet=quiet)
print,'S34   paramter',parms
 plot,cost,s34/s11,xr=xr,/xst
oplot,cost,yfit,color=cgcolor('red')

 plot,cost,s34,xr=xr2,/xst
oplot,cost,yfit*s11,color=cgcolor('red')

end
