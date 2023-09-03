#Use this code to make density, temperature or neutral hydrogen fraction slice/projection plots.
#Example command: python3 gas_analysis2.py -f /mnt/exports/data/mittal/outram/hydro-cool3/output_00015/ -pq 'density' -typ 'slc' -lm 8 -ah 1 -dhfe 0

import yt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse
import sys
from yt.extensions.astro_analysis.halo_analysis import HaloCatalog

direc='hydro-cool2'

parser = argparse.ArgumentParser()
parser.add_argument("-f","--fln", help="enter path of text file")
parser.add_argument("-pq","--phyquan", help="'density'/'temperature'/'xHI'",default='density')
parser.add_argument("-typ","--plotype", help="Slice ('slc') or Projection ('proj')",default='slc')
parser.add_argument("-lm","--levelmin",help="Base level of simulation",default=8)
parser.add_argument("-ah","--haloes",help="Add haloes '1' or '0'",default=0)
parser.add_argument("-dhfe","--halcre", help="Does halo file exist?",default=0)

args = parser.parse_args()
info = args.fln
fld=args.phyquan
typ=args.plotype
lev=int(args.levelmin)
add_haloes=args.haloes
hcc=args.halcre

ds = yt.load(info);
cg= ds.covering_grid(level=0, left_edge=[0, 0.0, 0.0], dims=[2**lev, 2**lev, 2**lev])
if fld!='xHI':
	quan=cg['gas',fld]
h = ds.hubble_constant

'''
meanfld=np.mean(ad['gas',fld])
if fld=='density':
	G=6.67e-11
	Omb=0.049
	Mpc2km=3.086e19
	H0=67.4
	z=ds.current_redshift
	dcrit=3*(H0/Mpc2km)**2/(8*np.pi*G)
	dgas=dcrit*Omb*(1+z)**3
	print("\nShikhar:Mean value of",fld,"is,", meanfld.to('kg/m**3'))
	print('The analytical estimate is,',dgas,'\n')
elif fld=='temperature':
	z=ds.current_redshift
	Tb=0.0187*(1+z)**2
	print("\nShikhar:Mean value of",fld,"is,", meanfld)
	print('The analytical estimate based on adiabatic evolution is,',Tb,'\n')
else:
	print("Shikhar: I don't know the analytical estimate for this quantity!\n")
'''

L=ds.length_unit.in_units('Mpccm').d*h  #length in cMpc/h
l=ds.length_unit.in_units('Mpccm').d	#length in cMpc

#The following is the halo catalogue path. I store it in RAMSES snapshots directory.
hp=info[:-13]+'info_000'+info[-3:-1]+'/info_000'+info[-3:-1]+'.0.h5'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fs=20

fig,ax=plt.subplots(figsize=(8.3,7.2))
fig.subplots_adjust(left=0.11, bottom=0.05, right=0.89, top=1.0)

if typ=='slc':
	if fld=='density':
		print('\nMean value', np.mean(quan),'\n')
		quan=quan.to('kg/m**3').d
		quan=quan/np.mean(quan)
		c=ax.imshow(quan[:,:,2**(lev-1)].T, cmap='viridis',aspect='equal',origin='lower',
				extent=[-L/2,L/2,-L/2,L/2],norm=colors.LogNorm())
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		cb=plt.colorbar(c,cax=cax)
		cb.set_label(r'$1+\delta$',fontsize=fs)
		#cb.set_label(r'$\rho\,\left(\mathrm{kg}\,\mathrm{m}^{-3}\right)$',fontsize=fs)
	elif fld=='temperature':
		print('\nMean value', np.mean(quan),'\n')
		quan=quan.d	
		c=ax.imshow(quan[:,:,2**(lev-1)].T, cmap='plasma',aspect='equal',origin='lower',
				extent=[-L/2,L/2,-L/2,L/2])
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		cb=plt.colorbar(c,cax=cax)
		cb.set_label(r'$T_{\mathrm{k}}\,(\mathrm{K})$',fontsize=fs)
	elif fld=='xHI':
		nHI = (cg['gas','HI_number_density']).d
		nHII= (cg['gas','HII_number_density']).d
		xHI = nHI/(nHI+nHII)
		print('\nMean value', np.mean(xHI),'\n')
		logxHI = 1e4*np.log10(xHI)
		c=ax.imshow(logxHI[:,:,2**(lev-1)].T, cmap='jet',aspect='equal',origin='lower',
				extent=[-L/2,L/2,-L/2,L/2])#,norm=colors.LogNorm())
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		cb=plt.colorbar(c,cax=cax)
		cb.set_label(r'$\log_{10}x_\mathrm{H\,{\Large\textsc{\lowercase{i}}}}$',fontsize=fs)

elif typ=='proj':
	print('\nShikhar: Full depth is being projected. If you want partial depth go back to the code.\n')
	dz=cg['gas','dz'].in_units('Mpccm').d
	if fld=='density':
		quan=quan.to('kg/m**3').d
		quan=quan/np.mean(quan)
		quan=np.mean(quan,axis=2)
		c=ax.imshow(quan.T, cmap='viridis',aspect='equal',origin='lower',
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
		#cb.set_label(r'Projected density$\,\left(\mathrm{M}_{\odot}\,\mathrm{cMpc}^{-3}\right)$',fontsize=fs)
		cb.set_label(r'$1+\delta$',fontsize=fs)
	elif fld=='temperature':
		quan=np.mean(quan.d,axis=2)
		c=ax.imshow(quan.T, cmap='plasma',aspect='equal',origin='lower',
				extent=[-L/2,L/2,-L/2,L/2])#,norm=colors.LogNorm()
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		cb=plt.colorbar(c,cax=cax)
		cb.set_label(r'$T_{\mathrm{k}}\,(\mathrm{K})$',fontsize=fs)
	elif fld=='xHI':
		nHI = (cg['gas','HI_number_density']).d
		nHII= (cg['gas','HII_number_density']).d
		xHI = nHI/(nHI+nHII)
		logxHI = 1e4*np.log10(xHI)
		proj_logxHI= np.sum(logxHI*dz,axis=2)/l

		c=ax.imshow(proj_logxHI.T, cmap='jet',aspect='equal',origin='lower',
				extent=[-L/2,L/2,-L/2,L/2])
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		cb=plt.colorbar(c,cax=cax)
		#cb.set_label(r'Projected $\log_{10}x_\mathrm{H\,{\Large\textsc{\lowercase{i}}}}\,(\mathrm{cMpc})$',fontsize=20)
else:
	print('Invalid plot type!')
	sys.exit()

cb.ax.tick_params(labelsize=fs)
ax.set_xlabel("$x\,$(cMpc$h^{-1}$)",fontsize=fs)
ax.set_ylabel("$y\,$(cMpc$h^{-1}$)",fontsize=fs)
ax.tick_params(axis='both', which='major', length=5, width=1, labelsize=fs,direction='in')
ax.tick_params(axis='both', which='minor', length=3, width=1, direction='in')
ax.minorticks_on()
#ax.set_yticks([-10,-5,0,5,10])
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')


if typ=='slc':
	if fld=='density':
		plt.savefig(direc+"/gas_000"+info[-3:-1]+'_Slice_z_density.pdf')
		print("Figure saved as "+direc+"/gas_000"+info[-3:-1]+'_Slice_z_density.pdf')
	elif fld=='temperature':
		plt.savefig(direc+"/gas_000"+info[-3:-1]+'_Slice_z_temperature.pdf')
		print("Figure saved as "+direc+"/gas_000"+info[-3:-1]+'_Slice_z_temperature.pdf')
	elif fld=='xHI':
		plt.savefig(direc+"/gas_000"+info[-3:-1]+'_Slice_z_xHI.pdf')
		print("Figure saved as "+direc+"/gas_000"+info[-3:-1]+'_Slice_z_xHI.pdf')

elif typ=='proj':
	if fld=='density':
		plt.savefig(direc+"/gas_000"+info[-3:-1]+'_Projection_z_density.pdf')
		print("Figure saved as "+direc+"/gas_000"+info[-3:-1]+'_Projection_z_density.pdf')
	elif fld=='temperature':
		plt.savefig(direc+"/gas_000"+info[-3:-1]+'_Projection_z_temperature.pdf')
		print(direc+"/gas_000"+info[-3:-1]+'_Projection_z_temperature.pdf')
	elif fld=='xHI':
		plt.savefig(direc+"/gas_000"+info[-3:-1]+'_Projection_z_xHI.pdf')
		print(direc+"/gas_000"+info[-3:-1]+'_Projection_z_xHI.pdf')
