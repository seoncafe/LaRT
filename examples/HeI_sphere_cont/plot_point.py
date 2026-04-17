#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from   matplotlib.ticker import AutoMinorLocator, MaxNLocator, LogLocator
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams.update({'font.size':14})

ps_file = 'fig_point_sphere.pdf'
fname = "Nscatt_sphere_point.txt"
temp, tau0, nscatt = np.loadtxt(fname, comments="#", unpack=True)
#lines = np.loadtxt(fname, comments="#")

for i in np.arange(temp.size):
   print("%5.0e %5.0e %9.6e" % (temp[i], tau0[i], nscatt[i]))

w1 = np.where(temp == 1e1)[0]
w4 = np.where(temp == 1e4)[0]

x1 = tau0[w1]
y1 = nscatt[w1]
x4 = tau0[w4]
y4 = nscatt[w4]

width = 0.9
fig, ax = plt.subplots(1,1,figsize=(6,6))
xmin = 0.6
xmax = 2e9
#xmin = 1.0
#xmax = 1e9
ymin = xmin
ymax = xmax
p1   = np.log10(xmin)
p2   = np.log10(xmax)
q1   = np.log10(ymin)
q2   = np.log10(ymax)
dp   = (p2-p1)*0.035
dq   = (q2-q1)*0.05

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
s = 0.9579
ax.plot([xmin,xmax],[s*ymin,s*ymax],'k',linewidth=width,linestyle='--')
ax.plot(x1,y1,'bo',linewidth=width)
ax.plot(x4,y4,'r.',linewidth=width)
ax.set_xscale("log")
ax.set_yscale("log")
ax.xaxis.set_major_locator(LogLocator(numticks=10))
ax.xaxis.set_minor_locator(LogLocator(subs=np.arange(2,10.)*0.1,numticks=10))
ax.yaxis.set_major_locator(LogLocator(numticks=10))
ax.yaxis.set_minor_locator(LogLocator(subs=np.arange(2,10.)*0.1,numticks=10))
ax.set_xlabel(r'$\tau_0$')
ax.set_ylabel('Number of Scatterings')

label  = r'point source, sphere'
label2 = r'$T_{\rm K} = 10$ K'
label3 = r'$T_{\rm K} = 10^4$ K'
ax.text(np.power(10.0,p1+dp),    np.power(10.0,q2-dq),    label, verticalalignment='top',horizontalalignment='left',fontsize=12)
ax.text(np.power(10.0,p1+dp*2.8),np.power(10.0,q2-dq*2.3),label2,verticalalignment='top',horizontalalignment='left',color='b',fontsize=12)
ax.text(np.power(10.0,p1+dp*2.8),np.power(10.0,q2-dq*3.3),label3,verticalalignment='top',horizontalalignment='left',color='r',fontsize=12)
aa = np.array([0.0]) + p1+dp*1.5
bb = np.array([1.0]) * (q2 - dq*2.6)
ax.plot(np.power(10.0,aa), np.power(10.0,bb),'bo',linewidth=width)
aa = np.array([0.0]) + p1+dp*1.5
bb = np.array([1.0]) * (q2 - dq*3.6)
ax.plot(np.power(10.0,aa), np.power(10.0,bb),'r.',linewidth=width)

aa = np.array([0.0]) + p1+dp*1.5
bb = np.array([1.0]) * (q2 - dq*4.6)
label4 = r'$N_{\rm scatt} = 0.958\tau_0$ ($\tau_0\rightarrow\infty$)'
ax.text(np.power(10.0,aa),np.power(10.0,bb),label4,verticalalignment='top',horizontalalignment='left',fontsize=12)

#plt.tight_layout(pad=0.5)
#plt.savefig(ps_file, papertype='a3')
plt.savefig(ps_file)
plt.show()
