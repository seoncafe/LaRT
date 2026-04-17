flist = ['tau1e+0','tau1e+1','tau2e+1','tau2e+2']
n = n_elements(flist)

;xr = [1189d0, 1199d0]
yr = [0d0,    3.2d0]
!p.multi=[0,2,2]
for i=0,n-1 do begin
  fname = flist[i]
  a = mrdfits(fname+'_V000.fits.gz',1)
  b = mrdfits(fname+'_V100.fits.gz',1)
  c = mrdfits(fname+'_V200.fits.gz',1)
  d = mrdfits(fname+'_V400.fits.gz',1)
  ;b = mrdfits(fname+'_V050.fits.gz',1)
  ;c = mrdfits(fname+'_V100.fits.gz',1)
     plot, a.lambda,a.Jout,xr=xr,/xst, yr=yr,/yst
    oplot, b.lambda,b.Jout,color=cgcolor('red')
    oplot, c.lambda,c.Jout,color=cgcolor('green')
    oplot, d.lambda,d.Jout,color=cgcolor('yellow')
endfor

end
