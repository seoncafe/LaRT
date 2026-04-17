;fname = 'mueller_Lyalpha.dat'
fname = 'mueller_matrix.dat'
readcol,fname,wavl,cext,albedo,hgg,nang,skip=1,numline=1
wavl   = wavl[0]
albedo = albedo[0]
hgg    = hgg[0]
nang   = nang[0]
readcol,fname,cost,s11,s12,s33,s34,skip=3

hgfun = 0.5d0*(1d0-hgg^2)/(1d0+hgg^2-2d0*hgg*cost)^(1.5d0)

ncost = n_elements(cost)
err   = dblarr(ncost)+1d0

!p.multi=[0,2,2]

;xr    = [-1,1]
xr    = [0.5d0,1d0]
ylog  = 0
quiet = 1
parm0 = [hgg]
parms = mpfitfun('func_s11',cost,s11,err,parm0,yfit=yfit,quiet=quiet)
print,'one   paramter',parms
 plot,cost,s11,xr=xr,/xst
oplot,cost,yfit,color=cgcolor('red')
oplot,cost,hgfun,color=cgcolor('green'),linestyle=1

 plot,cost,s11,/ylog,xr=[-1,1],/xst
oplot,cost,yfit,color=cgcolor('red')
oplot,cost,hgfun,color=cgcolor('green'),linestyle=1

;parm0 = [1.0d0,hgg]
;parms = mpfitfun('func_s11',cost,s11,err,parm0,yfit=yfit,quiet=quiet)
;print,'two  paramters',parms
; plot,cost,s11,ylog=ylog,xr=xr,/xst
;oplot,cost,yfit,color=cgcolor('red')
;oplot,cost,hgfun,color=cgcolor('green'),linestyle=1

parm0 = [0.7d0,0.8d0,0.4d0]
parms = mpfitfun('func_s11',cost,s11,err,parm0,yfit=yfit,quiet=quiet)
print,'three paramters',parms
 plot,cost,s11,xr=xr,/xst
oplot,cost,yfit,color=cgcolor('red')
oplot,cost,hgfun,color=cgcolor('green'),linestyle=1

 plot,cost,s11,/ylog,xr=[-1,1],/xst
oplot,cost,yfit,color=cgcolor('red')
oplot,cost,hgfun,color=cgcolor('green'),linestyle=1

;parm0 = [0.7d0,0.8d0,0.3d0,0.4d0]
;parms = mpfitfun('func_s11',cost,s11,err,parm0,yfit=yfit,quiet=quiet)
;print,'four paramters',parms
; plot,cost,s11,ylog=ylog,xr=xr,/xst
;oplot,cost,yfit,color=cgcolor('red')
;oplot,cost,hgfun,color=cgcolor('green'),linestyle=1

;parm0 = [0.6d0,0.8d0,0.1d0,0.4d0, 0.2d0]
;parms = mpfitfun('func_s11',cost,s11,err,parm0,yfit=yfit,quiet=quiet)
;print,'five paramters',parms
; plot,cost,s11,/ylog
;oplot,cost,yfit,color=cgcolor('red')
;oplot,cost,hgfun,color=cgcolor('green'),linestyle=1

end
