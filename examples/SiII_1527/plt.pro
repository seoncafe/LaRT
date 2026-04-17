;flist = ['t1e5tau1e0','t1e5tau2e0','t1e5tau5e0','t1e5tau1e1']
flist = ['t1e5tau1e0','t1e5tau5e0','t1e5tau1e1','t1e5tau2e1']
n = n_elements(flist)

xr = [1526d0, 1534d0]
yr = [0d0, 0.006d0]
!p.multi=[0,2,2]
for i=0,n-1 do begin
  fname = flist[i]
  a = mrdfits(fname+'.fits.gz',1)
  b = mrdfits(fname+'_V100.fits.gz',1)
  c = mrdfits(fname+'_V050.fits.gz',1)
     plot, a.lambda,a.Jout,xr=xr,/xst, yr=yr,/yst
    oplot, b.lambda,b.Jout,color=cgcolor('red')
    oplot, c.lambda,c.Jout,color=cgcolor('green')
endfor

end
