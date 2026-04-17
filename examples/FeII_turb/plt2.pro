flist = ['FeII_UV1','FeII_UV2','FeII_UV3']
n = n_elements(flist)

;xr = [1526d0, 1534d0]
;yr = [0d0, 0.006d0]
yr = [0d0, 2.8d0]
!p.multi=[0,2,2]
for i=0,n-1 do begin
  fname = flist[i]
  a = mrdfits(fname+'.fits.gz',1)
  b = mrdfits(fname+'_V050.fits.gz',1)
  c = mrdfits(fname+'_V100.fits.gz',1)
     plot, a.wavelength,a.Jout,xr=xr,/xst, yr=yr,/yst,xtit=textoidl('\lambda'),ytit='spectrum',tit=fname,psym=10
    oplot, b.wavelength,b.Jout,color=cgcolor('red'),psym=10
    oplot, c.wavelength,c.Jout,color=cgcolor('green'),psym=10
endfor

end
