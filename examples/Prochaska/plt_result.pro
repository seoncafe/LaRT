pro plt_result, fname, mode=mode

if n_elements(fname) ne 1 then fname = 'fiducial'
if n_elements(mode)  ne 1 then mode = 0

if mode eq 0 then begin
   a  = mrdfits(fname+'.fits.gz',1,/sil)
   J0 = mean(a.Jin)
endif else begin
   a        = mrdfits(fname+'.fits.gz',1,/sil)
   Jin_map  = mrdfits(fname+'_obs.fits.gz',2,/sil)
   Jout_map = mrdfits(fname+'_obs.fits.gz',0,/sil) + mrdfits(fname+'_obs.fits.gz',1,/sil)
   npix     = n_elements(Jin_map[0,*,*])
   a.Jin    = total(total(Jin_map, 2),2) /npix
   a.Jout   = total(total(Jout_map,2),2) /npix
   J0       = mean(a.Jin)
endelse

;xr = [2786d0,2812d0]
yr = [0.0d0, 2.7d0]

 plot,a.lambda,a.Jout/J0,xr=xr,/xst,yr=yr,/yst,psym=10
oplot,a.lambda,a.Jin /J0,color=cgcolor('green'),psym=10

oplot,[1d0,1d0]*2803.531d0,!y.crange,linestyle=1,color=cgcolor('green')
oplot,[1d0,1d0]*2796.352d0,!y.crange,linestyle=1,color=cgcolor('green')

readcol,'prochaska_MgII.txt',x,y,format='d,d'
oplot,x,y,color=cgcolor('red')

end
