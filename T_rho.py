#Shikhar Mittal
#Read your RAMSES hydrodynamic output and make T vs rho scatter plot.
#Also, fit an ideal gas law straight line to it.

import yt
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as scop

def Tfit(rho,gam,c):
	return 10**c*rho**(gam-1)

#Enter the path to RAMSES snapshot
ds=yt.load('/mnt/exports/data/mittal/outram/hydro-cool3/output_00015/')

ad=ds.all_data()

den=ad['gas','density']
temp=ad['gas','temperature']


def den_temp(density, temperature):
	sort=np.argsort(density)
	den_sort=density[sort]
	temp_sort=temperature[sort]
	popt,pcov=scop.curve_fit(Tfit,den_sort,temp_sort)

	T=Tfit(den_sort,*popt)                    
	print('\nParameters ',popt)
	print('Error ',np.sqrt(np.diag(pcov)),'\n')
	return den_sort, T

d,T = den_temp(den,temp)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig,ax=plt.subplots(figsize=(8.3,7.5),dpi=200)
fig.subplots_adjust(left=0.12, bottom=0.1, right=0.88, top=0.9)
ax.scatter(np.log10(den),np.log10(temp),color='b',s=2)
ax.plot(np.log10(d),np.log10(T),'k--',label=r'$T\propto\rho^{\gamma-1}$')
ax.set_xlabel(r'$\log_{10}\rho$', fontsize=20)
ax.set_ylabel(r'$\log_{10}T$', fontsize=20)
ax.minorticks_on()
ax.set_rasterized(True)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(axis='both', which='major', length=5, width=1, labelsize=20,direction='in')
ax.tick_params(axis='both', which='minor', length=3, width=1,direction='in')
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
ax.legend(fontsize=18,frameon=False)
plt.savefig('Tvsrho.pdf')
print('Shikhar: done plotting!')
