!p.multi=[0,1,2]
;mode = 0
a = mrdfits('FeII_UV1_a.fits.gz',1)
b = mrdfits('FeII_UV1_b.fits.gz',1)
c = mrdfits('FeII_UV1_c.fits.gz',1)
plot,a.lambda,a.Jout,/xst
oplot,b.lambda,b.Jout,color=cgcolor('red')
oplot,c.lambda,c.Jout,color=cgcolor('green')

Pa = mrdfits('FeII_UV1_Pa.fits.gz',1)
Pb = mrdfits('FeII_UV1_Pb.fits.gz',1)
Pc = mrdfits('FeII_UV1_Pc.fits.gz',1)
 plot,Pa.lambda,Pa.Jout,/xst
oplot,Pb.lambda,Pb.Jout,color=cgcolor('red')
oplot,Pc.lambda,Pc.Jout,color=cgcolor('green')

end
