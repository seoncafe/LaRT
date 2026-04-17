;flist = ['FeII_UV1','FeII_UV1_V050','FeII_UV1_V100']
;flist = ['FeII_UV2','FeII_UV2_V050','FeII_UV2_V100']
;flist = ['FeII_UV3','FeII_UV3_V050','FeII_UV3_V100']
;flist = ['FeII_UV1_a','FeII_UV1_b','FeII_UV1_c']
flist = ['FeII_UV1_Pa','FeII_UV1_Pb','FeII_UV1_Pc']

!p.multi=[0,1,3]
!p.charsize=3.0
for k=0, n_elements(flist)-1 do begin
    fname = flist[k]
    a = mrdfits(fname+'.fits.gz',1, /sil)
    p = mrdfits(fname+'_obs.fits.gz',0, /sil)
    q = mrdfits(fname+'_obs.fits.gz',1, /sil)
    s = total(total(p+q,2),2)
    s = s/median(s)
     plot,a.lambda,a.Jout,/xst
    oplot,a.lambda,s,color=cgcolor('red')
endfor

end
