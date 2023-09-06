#Shikhar Mittal
#Use this code to analyse and plot the total accumulated optical depth (OD) by the photons on their journey.

import numpy as np
from matplotlib import pyplot as plt
from scipy.io import FortranFile as ff
import scipy.integrate as scint


dir = '/mnt/exports/data/mittal/outras/hydro-cool3/' #Folder where you have stored your RASCAS output file

f = ff(dir+'tau.dat')		#Optical depth file (unformatted binary file)
[nphot] = f.read_ints()
print('No.of photons =',nphot)
tau = f.read_reals()
f.close()
x = np.log10(tau)

histo,bn=np.histogram(x, bins=200, density=True)
ltau = (bn[1:]+bn[:-1])/2

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig,ax=plt.subplots(figsize=(8.3,7.5))
fig.subplots_adjust(left=0.12, bottom=0.07, right=0.88, top=0.97)
#ax.hist(x,bins=200,density=True,color='blue',histtype='step')#,log=True)
ax.plot(ltau, histo,'b')
ax.set_xlabel(r'$\log_{10}\tau$', fontsize=20)
ax.minorticks_on()
ax.set_ylim(bottom=0.0)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(axis='both', which='major', length=5, width=1, labelsize=20,direction='in')
ax.tick_params(axis='both', which='minor', length=3, width=1, direction='in')
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plt.savefig(dir+'OD_distribution.pdf')
#plt.show()


#The following lines of code are to find what % of MC photons had a total OD less than 0.1 

#ind = np.where(ltau<-1)
#i = ind[0][-1]
#print('percent of photons with total OD less than 0.1 ',100*scint.trapz(histo[0:i],ltau[0:i]))
