#In this I plot the DM density slice/projection from a RAMSES simulation. Can be used on both types of simulation.
#python3 dm_analysis.py -f /mnt/exports/data/mittal/outram/dmo3/output_00053/ -typ 'slc' -lm 8 -ah 1 -dhfe 1

import numpy as np
import matplotlib.pyplot as plt
import yt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from yt.extensions.astro_analysis.halo_analysis import HaloCatalog
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f","--file", help="enter path of text file")
parser.add_argument("-typ","--plotype", help="Slice ('slc') or Projection ('proj')",default='slc')
parser.add_argument("-lm","--levelmin",help="Base level of simulation",default=8)
parser.add_argument("-ah","--haloes",help="Add haloes '1' or '0'",default=0)
parser.add_argument("-dhfe","--halcre", help="Does halo file exist?",default=0)


args = parser.parse_args()

info = args.file
typ=args.plotype
lev=int(args.levelmin)
add_haloes=args.haloes
hcc=args.halcre


ds = yt.load(info)
hp = info[:-13]+'info_000'+info[-3:-1]+'/info_000'+info[-3:-1]+'.0.h5'


h = ds.hubble_constant
z = ds.current_redshift

l=ds.length_unit.in_units('Mpccm').d	#length in cMpc
L=ds.length_unit.in_units('Mpccm').d*h  #length in cMpc/h

cg = ds.arbitrary_grid(left_edge=ds.domain_left_edge, right_edge=ds.domain_right_edge, dims=[2**lev, 2**lev, 2**lev])
den=cg['deposit','all_density']

G=6.67e-11
Om_m = ds.omega_matter
Mpc2km=3.086e19
Ho=100*h
d_crit=3*(Ho/Mpc2km)**2/(8*np.pi*G)
d_DM=d_crit*Om_m*(1+z)**3
print("\nShikhar: mean DM density is ", 1000*np.mean(den.d),'kg/m^3')
print('The analytical estimate is',d_DM,'kg/m^3\n')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig,ax=plt.subplots(figsize=(8.3,7.5))
fig.subplots_adjust(left=0.12, bottom=0.05, right=0.88, top=0.95)

if typ=='slc':
	den=den.to('Msun/Mpccm**3').d
	c=ax.imshow(den[:,:,2**(lev-1)].T, cmap='viridis',aspect='equal',origin='lower',
				extent=[-L/2,L/2,-L/2,L/2],norm=colors.LogNorm())
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	cb=plt.colorbar(c,cax=cax)
	cb.set_label(r'$\rho\,\left(\mathrm{M}_{\odot}\,\mathrm{cMpc}^{-3}\right)$',fontsize=20)

elif typ=='proj':
	den=den.to('Msun/Mpccm**3')
	den=np.mean(den.d,axis=2)
	c=ax.imshow(den.T, cmap='viridis',aspect='equal',origin='lower',
				extent=[-L/2,L/2,-L/2,L/2],norm=colors.LogNorm())
	if add_haloes=='1':
		if hcc=='0':
			hc = HaloCatalog(data_ds=ds, finder_method="hop",output_dir=info[:-13],
			 finder_kwargs={"threshold": 160, "ptype":"all"})
			hc.create()
			dsh = yt.load(hp)
		elif hcc=='1':
			dsh = yt.load(hp)	
		else:
			print('Shikhar: Invalid halo option. Bye!')
			sys.exit()
		adh=dsh.all_data()
		mass=adh['halos','particle_mass'].d
		posi=adh['halos','particle_position'].in_units('Mpccm').d*h-L/2
		m=min(mass)
		ax.scatter(posi[:,0],posi[:,1], facecolors='none',edgecolors='white',
 s=20*(np.log10(mass/m))**4,linewidth=0.5)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	cb=plt.colorbar(c,cax=cax)
	cb.set_label(r'Projected density$\,\left(\mathrm{M}_{\odot}\,\mathrm{cMpc}^{-3}\right)$',fontsize=20)
else:
	print('Invalid plot type!')
	sys.exit()

ax.set_xlabel("$x\,$(cMpc$h^{-1}$)",fontsize=20)
ax.set_ylabel("$y\,$(cMpc$h^{-1}$)",fontsize=20)
ax.tick_params(axis='both', which='major', length=5, width=1, labelsize=20,direction='in')
ax.tick_params(axis='both', which='minor', length=3, width=1, direction='in')
ax.minorticks_on()
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

out_dir='/mnt/exports/data/mittal/outram/dmo3/'
if typ=='slc':
	plt.savefig(out_dir+"dm_000"+info[-3:-1]+'_Slice_z_density.pdf')
elif typ=='proj':
	plt.savefig(out_dir+"dm_000"+info[-3:-1]+'_Projection_z_density.pdf')
