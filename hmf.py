#Use this code to plot to find the HMF and compare it with analytical form. Can be run on both, hydro as well as DMO.
#I have used PS74 for analytical. To run this code:
#python3 hmf.py -f /mnt/exports/data/mittal/outram/dmo3/output_*/

import yt
import numpy as np
import matplotlib.pyplot as plt
import sys
import pathlib
import argparse

from yt.extensions.astro_analysis.halo_analysis import HaloCatalog
from colossus.cosmology import cosmology
from colossus.lss import mass_function
from colossus.lss import peaks
my_cosmo = {'flat': True, 'H0': 67.4, 'Om0': 0.315, 'Ob0': 0.049, 'sigma8': 0.811, 'ns': 0.965,'relspecies': False,'Tcmb0': 2.725}
cosmo = cosmology.setCosmology('my_cosmo', my_cosmo)

h=0.674


parser = argparse.ArgumentParser()
parser.add_argument("-f","--file", help="enter path of RAMSES output folder")
args = parser.parse_args()
info = args.file

hcc=input('Halo catalogue created? (1 or 0): ')

ds = yt.load(info)
z=ds.current_redshift
l=ds.length_unit.in_units('Mpccm')

V=(l.d*h)**3 #Volume in (cMpc/h)^3

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#The following is the halo catalogue path. I store it in RAMSES snapshots directory.
fn=info[:-13]+'info_000'+info[-3:-1]+'/info_000'+info[-3:-1]+'.0.h5'

if hcc=='0':
	hc = HaloCatalog(data_ds=ds, finder_method="hop",output_dir=info[:-13],
                         finder_kwargs={"threshold": 178, "ptype":"all"})
	hc.create()
	dsh = yt.load(fn)
elif hcc=='1':
	dsh = yt.load(fn)
else:
	print('Invalid option. Bye!')
	sys.exit


print('\nShikhar: loading haloes from', fn,'\n')
adh=dsh.all_data()	#This loads all the halos at once

M=adh['halos','particle_mass'].in_units('Msun') #Accessing masses of halos in M_sun units.
M=np.sort(M)	#Arrange them in ascending order.

hist,bins=np.histogram(np.log10(M), bins=50)
widths = np.diff(bins)
hist = hist/widths/V

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig,ax=plt.subplots(figsize=(8.3,7.5))
fig.subplots_adjust(left=0.1, bottom=0.12, right=0.9,top=0.92)
ax.semilogy((bins[1:]+bins[:-1])/2, hist,'r',label='simulation')
ax.semilogy(np.log10(M),
                np.log(10)*mass_function.massFunction(M*h, z, q_in='M',
                q_out='dndlnM', mdef='fof', model = 'press74'),'b',label='analytical')

ax.set_xlabel(r'$\log_{10} (M/\mathrm{M}_{\odot})$',fontsize=20)
ax.set_ylabel(r'$\frac{\mathrm{d}n}{\mathrm{d}\log_{10} M}\,(\mathrm{cMpc}\,h^{-1})^{-3}$',fontsize=20)
ax.tick_params(axis='both', which='major', length=5, width=1, labelsize=20,direction='in')
ax.tick_params(axis='both', which='minor', length=3, width=1, direction='in')
ax.legend(fontsize=18,frameon=False)
ax.minorticks_on()
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plotname='/hmf_000'+info[-3:-1]+'.pdf'
print('\nPlot saved as',plotname)
plt.savefig(plotname)
