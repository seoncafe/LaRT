fname = 'mueller_Lyalpha.dat'
readcol,fname,wavl,cext,albedo,hgg,nang,skip=1,numline=1
wavl   = wavl[0]
albedo = albedo[0]
hgg    = hgg[0]
nang   = nang[0]
readcol,fname,cost,s11,s12,s33,s34,skip=3

xr    = [-1,1]
xr2   = [-1,1]

if !d.name eq 'PS' then begin
   set_thick,3.0
   !p.font = 1
   !x.margin=[3,3]
   !y.margin=[2,2]
   margin = 0.07
   psname = 'mueller_Lya.ps'
   device,/color,bit=8,set_font='Times',xsize=19,ysize=13,yoffset=2,filename=psname,/tt_font
   !p.charsize = 1.6
   charsize = 1.6
endif

!p.multi=[0,2,2]
!x.title=textoidl('cos\theta')

;=======
 plot,cost,s11,xr=xr,/xst,ytit=textoidl('S_{11}')
 plot,cost,s12,xr=xr2,/xst,yr=[-0.5,0],ytit=textoidl('S_{12}')
 plot,cost,s33,xr=xr,/xst,ytit=textoidl('S_{33}')
 plot,cost,s34,xr=xr2,/xst,ytit=textoidl('S_{34}')
if !d.name eq 'PS' then begin
 device,/close
 set_plot,'x'
endif

end
