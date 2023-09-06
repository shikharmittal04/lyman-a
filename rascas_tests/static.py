#Use this to analyse the output from a static, uniform and homogeneous medium.
#For a monochromatic source at the centre the results may be compared with analytical
#solution of Dijkstra et al (2006).

import argparse
import sys
from matplotlib import pyplot as plt
from scipy.io import FortranFile as ff
import numpy as np
import photons as p

rtilde=0.49
L=1e8           #box size (NOT the sphere radius) in cm
b=1.29e5        #average thermal speed in cm/s
nhi=3.477e11    #number density in cm^-3

NH=nhi*rtilde*L #column density
T=1e4*(b/12.9e5)**2     #temperature
tau_0=8.3e6*(NH/2e20)*(T/2e4)**(-1/2) #central optical depth

grid = p.from_file('/home/shikhar/work/tau1e7.dat')

x = np.arange(-80,80,0.01)
a       = 4.71e-4 / np.sqrt(T/1e4)
J_x    = x**2 * np.sqrt(np.pi) / (np.sqrt(96)*a*tau_0) / (1+np.cosh(np.sqrt(2*np.pi**3/27)*np.abs(x**3)/a/tau_0))
# Dijkstra 2006, formulae C-17 
# analytic solution from Dijkstra+06
# Normalised to 1/4π

print("\nTemperature (K) =","{0:.2e}".format(T))
#print("Number density (cm^-3) =","{0:.2E}".format(nH))
print("Line centre OD, tau_0 =","{0:.2e}".format(tau_0))
#print("v_th (cm/s) =", "{0:.2e}".format(v),'\n')


delta_nu = p.nu_0 * b / p.clight

xx = (grid.nu - p.nu_0) / delta_nu

nbin=100
hist,bn=np.histogram(xx, bins=nbin, density=True)
hist = hist/(4*np.pi)

#------------------------------------------------------------------

rec=p.from_file('/home/shikhar/work/tau1e7_recoil.dat')
x_rec = (rec.nu - p.nu_0) / delta_nu

histrec,bnrec=np.histogram(x_rec, bins=nbin, density=True)
histrec = histrec/(4*np.pi)

#------------------------------------------------------------------

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig,ax=plt.subplots(figsize=(8.3,7.3))
fig.subplots_adjust(left=0.11, bottom=0.05, right=0.89, top=1.0)
ax.plot((bn[1:]+bn[:-1])/2, 1000*hist,'b',label='No recoil')
ax.plot((bnrec[1:]+bnrec[:-1])/2, 1000*histrec,color='limegreen',ls='dashed',label='With recoil')
ax.plot(x,1000*J_x,color='red',ls=':',lw=2,label='Analytical')


ax.set_xlabel('$x$', fontsize=20)
ax.set_ylabel(r'$10^3\times J(x)$', fontsize=20)
#ax.set_title(r"$T={Temp}, n_{{\mathrm{{H}}}}={numden}, \tau_0={OD}$".format(Temp=T, numden=nH, OD=tau_0),fontsize=18)
#ax.set_title(r'$T=10^4\,$K, $\tau_{0}=1.63\times10^{6}$, no.of photons $=10^5$.',fontsize=18,pad=10)
ax.minorticks_on()
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.set_ylim([-2e-3,1.8])
ax.set_xlim([-80,80])
ax.tick_params(axis='both', which='major', length=5, width=1, labelsize=20,direction='in')
ax.tick_params(axis='both', which='minor', length=3, width=1, direction='in')
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
ax.legend(fontsize=18,frameon=False)
plt.savefig('tau1e7.pdf')
