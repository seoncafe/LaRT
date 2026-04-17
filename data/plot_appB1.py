#/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from   scipy.optimize import curve_fit
from   matplotlib.ticker import AutoMinorLocator, LinearLocator, MultipleLocator

pdf_file = 'fig_appB1.pdf'
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams.update({'font.size':8})
plt.rcParams.update({'font.size':9})

def hg_func(cost,g):
    y = 0.5*(1.0-g**2)/np.power(1.0+g**2-2.0*cost*g,3.0/2.0)
    return y

def s11_func(cost,a1,g1,g2):
    y = a1*hg_func(cost,g1) + (1.0-a1)*hg_func(cost,g2)
    return y

def s12_func(cost,p,th0,scal):
    th     = np.arccos(cost)
    y      = -p*(1.0-cost**2)/(1.0+(scal*np.cos(th-th0))**2)
    #y      = -p*(1.0-cost**2)/(1.0+scal*(np.cos(th-th0))**2)
    return y

def s33_func(cost,scal,a0):
    y1 = 2.0*cost/(1.0+cost**2)
    y2 = cost + a0*(cost**2 - 1.0)
    y  = scal*y1 + (1.0-scal)*y2
    return y

def s34_func(cost,pcir,s,scal):
    th = np.arccos(cost)
    cc = np.cos(th + s*th*np.exp(-scal*th/np.pi))
    #y  = -pcir*(1.0-cc**2)/(1.0+cc**2)
    y  = pcir*(1.0-cc**2)/(1.0+cc**2)
    return y

#fig, ax = plt.subplots(3,4,figsize=(12,9))
#fig, ax = plt.subplots(3,4,figsize=(8,4.5))
fig, ax = plt.subplots(3,4,figsize=(8,5))
label = [r'$S_{11}$',r'$S_{12}/S_{11}$',r'$S_{33}/S_{11}$',r'$S_{34}/S_{11}$']
dtype = ['MW','LMC','SMC']

dir   = '../data/'
flist = ['mueller_Lyalpha.dat','mueller_Lyalpha_LMC.dat','mueller_Lyalpha_SMC.dat']
n     = len(flist)
for i in range(n):
   fname = dir+flist[i]
   cost,S11,S12,S33,S34 = np.loadtxt(fname,usecols=(0,1,2,3,4),skiprows=3,unpack=True)
   #S34 = -S34
   file  = open(fname,"r")
   for k in range(2): a = file.readline()
   file.close()
   g_lya = np.float_(a.split()[3])
   print("%-25s: %7.4f" % (fname,g_lya))

   popt11, pcov11 = curve_fit(s11_func, cost, S11, bounds=([0.0,0.5,0.0],[1.0,1.0,0.7]))
   popt12, pcov12 = curve_fit(s12_func, cost, S12/S11)
   popt33, pcov33 = curve_fit(s33_func, cost, S33/S11)
   popt34, pcov34 = curve_fit(s34_func, cost, S34/S11)
   print(popt11)
   print(popt12)
   print(popt33)
   print(popt34)

   #--- 11
   ax[i][0].semilogy(cost,S11,color='k')
   ax[i][0].semilogy(cost,s11_func(cost,*popt11),color='r')
   ax[i][0].semilogy(cost,hg_func(cost,g_lya),color='g')
   #th = np.arccos(cost)*180.0/np.pi
   #ax[0].semilogy(th,S11,color='k')
   #ax[0].semilogy(th,s11_func(cost,*popt),color='r')
   #ax[0].semilogy(th,hg_func(cost,g_lya),color='g')
   #--- 12
   ax[i][1].plot(cost,S12/S11,color='k')
   ax[i][1].plot(cost,s12_func(cost,*popt12),color='r')
   #--- 33
   ax[i][2].plot(cost,S33/S11,color='k')
   ax[i][2].plot(cost,s33_func(cost,*popt33),color='r')
   #--- 34
   ax[i][3].plot(cost,S34/S11,color='k')
   ax[i][3].plot(cost,s34_func(cost,*popt34),color='r')

   ax[i][0].set_ylim(3e-2,3e1)
   ax[i][1].set_ylim(-0.6,0.0)
   ax[i][2].set_ylim(-1.0,1.0)
   #ax[i][3].set_ylim(-0.4,0.0)
   ax[i][3].set_ylim(0.0,0.4)
   ax[i][2].yaxis.set_major_locator(MultipleLocator(0.5))
   ax[i][2].yaxis.set_minor_locator(MultipleLocator(0.1))
   ax[i][3].yaxis.set_major_locator(MultipleLocator(0.1))
   ax[i][3].yaxis.set_minor_locator(MultipleLocator(0.02))

   for k in range(4): 
      ax[i][k].set_xlim(-1,1)
      ax[i][k].set_xlabel(r'$\cos\theta$')
      ax[i][k].xaxis.set_minor_locator(MultipleLocator(0.1))
      ax[i][k].set_ylabel(label[k])
      p1, p2 = ax[i][k].get_xlim()
      #dp     = (p2-p1)*0.035*3.0
      dp     = (p2-p1)*0.035*2.0
      if k==0:
         q1, q2 = np.log10(ax[i][k].get_ylim())
         dq     = (q2-q1)*0.05
         ax[i][k].text(p1+dp, np.power(10.0,q2-dq), dtype[i], color='k', verticalalignment='top',horizontalalignment='left', size=7)
      else:
         q1, q2 = ax[i][k].get_ylim()
         dq     = (q2-q1)*0.05
         ax[i][k].text(p1+dp, q2-dq, dtype[i], color='k', verticalalignment='top',horizontalalignment='left', size=7)

#plt.tight_layout()
plt.tight_layout(h_pad=0.5,w_pad=0.5)
plt.savefig(pdf_file)
plt.show()
