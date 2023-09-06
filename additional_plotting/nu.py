#Shikhar Mittal
#Use this code to analyse and plot the how do the output frequencies compare with input frequencies.
#You can also check the SED of the input frequencies.

import numpy as np
from matplotlib import pyplot as plt
from scipy.io import FortranFile as ff
import scipy.integrate as scint

lam_alp = 1215.67e-10	#Lya wavelength in m
lam_bet = 1025.72e-10	#Ly-beta wavelength in m
c = 2.998e8		#Speed of light in m/s
nu_alp = c/lam_alp	#Lya frequency in Hz
nu_bet = c/lam_bet	#Ly-beta frequency in Hz

dir = '/mnt/exports/data/mittal/outras/hydro-cool3/lowmc/'  #Folder where you have stored your RASCAS output file

#-----------------------------------------------------------
print('Reading the input frequencies...')
f1 = ff(dir+'ic_p1e8.dat')
[nphot] = f1.read_ints()
tot_flux = f1.read_reals()
iseed = f1.read_ints()
id = f1.read_ints()
nu_emit = f1.read_reals()
f1.close()
print('No.of photons =',nphot,'\n')

#-----------------------------------------------------------

print('Reading the output frequencies...')
f2 = ff(dir+'nu_tau_p1e8.dat')
[nphot] = f2.read_ints()
nu_exit = f2.read_reals()
f2.close()
print('Done')
#----------------------------------------------------------
rs=nu_emit/nu_exit

max_rs=np.max(rs)
id=np.where(rs==max_rs)
print('\nLyman-a frequency (Hz) =',nu_alp)
print('Input frequency (Hz) =',nu_emit[id])
print('Exit frequency (Hz) =',nu_exit[id])
print('Maximum redshift (1+z) =',max_rs)
print('Minimum redshift (1+z) =',np.min(rs))

'''
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig,ax=plt.subplots(figsize=(8.3,7.5))
fig.subplots_adjust(left=0.12, bottom=0.07, right=0.88, top=0.97)
ax.semilogx(rs,color='b',lw=0.5)
ax.axhline(y=1,color='r',ls='--')
ax.set_xlabel('Photon ID', fontsize=20)
ax.set_ylabel('$1+z$', fontsize=20)
ax.minorticks_on()
ax.set_xlim([1,nphot])
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(axis='both', which='major', length=5, width=1, labelsize=20,direction='in')
ax.tick_params(axis='both', which='minor', length=3, width=1, direction='in')
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
ax.set_rasterized(True)
plt.savefig(dir+'z_distribution_new.pdf')
'''
#----------------------------------------------------------
#using the following lines of code you can check if the input frequencies have the SED you expect it to be.

'''
a=-0.86
Nab=6520
x = np.linspace(nu_alp,nu_bet,1000)
y = (a+1)*Nab*x**a/(nu_bet**(a+1)-nu_alp**(a+1))

histo,bn=np.histogram(nu_emit, bins=200, density=True)
nu = (bn[1:]+bn[:-1])/2
histo = histo*Nab

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig,ax=plt.subplots(figsize=(8.3,7.5))
fig.subplots_adjust(left=0.12, bottom=0.05, right=0.88, top=0.98)
ax.loglog(nu, histo,'b',label='Code')
ax.loglog(x, y,'r--',label='Analytical')
ax.legend(fontsize=18,frameon=False)
ax.set_xlabel(r'$\nu\,$(Hz)', fontsize=20)
ax.minorticks_on()
ax.set_xlim([nu_alp,nu_bet])
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(axis='both', which='major', length=5, width=1, labelsize=20,direction='in')
ax.tick_params(axis='both', which='minor', length=3, width=1, direction='in')
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
ax.set_rasterized(True)
plt.savefig(dir+'SED.pdf')
print('\nDone.')
'''
