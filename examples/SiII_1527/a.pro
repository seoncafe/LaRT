;flist = ['t1e5tau1e0','t1e5tau1e0_V050','t1e5tau1e0_V100']
;flist = ['t1e5tau5e0','t1e5tau5e0_V050','t1e5tau5e0_V100']
flist = ['t1e5tau1e1','t1e5tau1e1_V050','t1e5tau1e1_V100']
;flist = ['t1e5tau2e1','t1e5tau2e1_V050','t1e5tau2e1_V100']

!p.multi=[0,1,3]
!p.charsize=3.0
xr = [1526.d0, 1534d0]
for k=0, n_elements(flist)-1 do begin
    fname = flist[k]
    a = mrdfits(fname+'.fits.gz',1,hdr, /sil)
    p = mrdfits(fname+'_obs.fits.gz',0, /sil)
    q = mrdfits(fname+'_obs.fits.gz',1, /sil)
    ;dx = fxpar(hdr,'CD1_1') * !dpi/180d0
    ;dy = fxpar(hdr,'CD2_2') * !dpi/180d0
    ;s  = total(total(p+q,2),2)/(dx*dy)
    ;s = s/median(s)
    s  = total(total(p+q,2),2)
    s  = s/total(s)*total(a.Jout)
     plot,a.lambda,a.Jout,xr=xr,/xst
    oplot,a.lambda,s,color=cgcolor('red'),psym=4
endfor

end
