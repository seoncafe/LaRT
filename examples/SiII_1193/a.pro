;flist = ['tau1e+3_V000','tau1e+3_V200','tau1e+3_V400']
;flist = ['tau5e+2_V000','tau5e+2_V200','tau5e+2_V400']
;flist = ['tau2e+2_V000','tau2e+2_V200','tau2e+2_V400']
;flist = ['tau1e+2_V000','tau1e+2_V200','tau1e+2_V400']
;flist = ['tau5e+1_V000','tau5e+1_V200','tau5e+1_V400']
;flist = ['tau2e+1_V000','tau2e+1_V200','tau2e+1_V400']
flist = ['tau1e+1_V000','tau1e+1_V200','tau1e+1_V400']
;flist = ['tau1e+1_V000','tau1e+1_V050','tau1e+1_V100']

;flist = ['tau1e+2_V000','tau2e+2_V000','tau5e+2_V000']
;flist = ['tau1e+1_V000','tau2e+1_V000','tau5e+1_V000']
;flist = ['tau1e+0_V000','tau2e+0_V000','tau5e+0_V000']

!p.multi=[0,1,3]
!p.charsize=3.0
;xr = [1526.d0, 1534d0]
yr = [0,3]
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
     plot,a.lambda,a.Jout,xr=xr,/xst,yr=yr,/yst,tit=fname
    oplot,a.lambda,s,color=cgcolor('red'),psym=4
    oplot,!x.crange,[1d0,1d0],linestyle=1
endfor

end
