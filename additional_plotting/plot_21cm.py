#Shikhar Mittal
#Use this code to make single panel plots (either projections or slices) of Lya coupling or 21-cm signal.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

L_cMpc_by_h = 10	#Box size in cMpc/h
lev = 8			#levelmin (no.of cells accross each axis in log_2 units)
dir = '/mnt/exports/data/mittal/outras/delta/'	#folder where you have stored your outputs.

print('Reading files ...\n')
xa_new = np.load(dir+'xa.npy')
#T21_new = np.load(dir+'T21.npy')
#xa_old = np.load(dir+'xa_old.npy')
#T21_old = np.load(dir+'T21_1e8.npy')

print('Plotting ...')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def slices(name,data):
	fig,ax=plt.subplots(figsize=(8.5,7.2),dpi=200)
	fig.subplots_adjust(left=0.11, bottom=0.05, right=0.89, top=1.0)
	if name=='T21' or name=='T21_old':
		c=ax.imshow(data[:,:,2**(lev-1)].T, cmap='brg',aspect='equal',origin='lower',extent=[-L_cMpc_by_h/2,L_cMpc_by_h/2,-L_cMpc_by_h/2,L_cMpc_by_h/2], vmin=-1100, vmax=200)
	elif name=='palpha' or name=='xalpha' or 'xalpha_old':
		c=ax.imshow(data[:,:,2**(lev-1)].T, cmap='cool',aspect='equal',origin='lower',extent=[-L_cMpc_by_h/2,L_cMpc_by_h/2,-L_cMpc_by_h/2,L_cMpc_by_h/2], norm=colors.LogNorm())#vmin=1e-3,vmax=1500
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	
	if name=='T21' or name=='T21_old':
		cb=plt.colorbar(c,cax=cax)
		cb.set_label(r'$\Delta T_{\mathrm{b}}\,(\mathrm{mK})$',fontsize=16)
	elif name=='palpha':
		cb=plt.colorbar(c,cax=cax, extend='min')
		cb.set_label(r'$P_{\alpha}\,(\mathrm{s}^{-1})$',fontsize=16)
	elif name=='xalpha_old':
		cb=plt.colorbar(c,cax=cax)
		cb.set_label(r'$x_{\alpha}$',fontsize=16)
	elif name=='xalpha':
		cb=plt.colorbar(c,cax=cax, extend='min')
		cb.set_label(r'$x_{\alpha}$',fontsize=16)
	cb.ax.tick_params(labelsize=16)
	ax.set_xlabel("$x\,$(cMpc$h^{-1}$)",fontsize=20)
	ax.set_ylabel("$y\,$(cMpc$h^{-1}$)",fontsize=20)
	ax.tick_params(axis='both', which='major', length=5, width=1, labelsize=20,direction='in')
	ax.tick_params(axis='both', which='minor', length=3, width=1,direction='in')
	#ax.set_yticks([-10,-5,0,5,10])
	ax.minorticks_on()
	ax.yaxis.set_ticks_position('both')
	ax.xaxis.set_ticks_position('both')
	ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
	plt.savefig('/mnt/exports/data/mittal/outras/'+simid+'Slc_'+name+'.png')
	#plt.savefig('/mnt/exports/data/mittal/outras/'+simid+'Slc_'+name+simid[-6:-1]+'.png')
	#print('Figure saved as /mnt/exports/data/mittal/outras/'+simid+'Slc_'+name+simid[-6:-1]+'.png')

def proj(name,data):
	fig,ax=plt.subplots(figsize=(8.5,7.2),dpi=200)
	fig.subplots_adjust(left=0.11, bottom=0.05, right=0.89, top=1.0)
	data = np.mean(data,axis=2)
	if name=='T21_new' or name=='T21_old':
		c=ax.imshow(data.T/1000, cmap='brg',aspect='equal',origin='lower',extent=[-L_cMpc_by_h/2,L_cMpc_by_h/2,-L_cMpc_by_h/2,L_cMpc_by_h/2],vmin=-0.35,vmax=0.012)
	elif name=='palpha' or name=='xalpha' or 'xalpha_old':
		c=ax.imshow(data.T, cmap='cool',aspect='equal',origin='lower',extent=[-L_cMpc_by_h/2,L_cMpc_by_h/2,-L_cMpc_by_h/2,L_cMpc_by_h/2])#, vmin=0,vmax=2)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	
	if name=='T21_new' or name=='T21_old':
		cb=plt.colorbar(c,cax=cax)
		cb.set_label(r'$T_{21}\,(\mathrm{K})$',fontsize=20)
	elif name=='palpha':
		cb=plt.colorbar(c,cax=cax)
		cb.set_label(r'$P_{\alpha}\,(\mathrm{s^{-1}})$',fontsize=20)
	elif name=='xalpha_old':
		cb=plt.colorbar(c,cax=cax)
		cb.set_label(r'$x_{\alpha}$',fontsize=20)
	elif name=='xalpha':
		cb=plt.colorbar(c,cax=cax)
		cb.set_label(r'$x_{\alpha}$',fontsize=20)
	cb.ax.tick_params(labelsize=20)
	ax.set_xlabel("$x\,$(cMpc$h^{-1}$)",fontsize=20)
	ax.set_ylabel("$y\,$(cMpc$h^{-1}$)",fontsize=20)
	ax.tick_params(axis='both', which='major', length=5, width=1, labelsize=20,direction='in')
	ax.tick_params(axis='both', which='minor', length=3, width=1,direction='in')
	ax.minorticks_on()
	ax.yaxis.set_ticks_position('both')
	ax.xaxis.set_ticks_position('both')
	ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
	plt.savefig(dir+'Proj_'+name+'.png')
	print('Figure saved as '+dir+'Proj_'+name+'.png')

#slices('T21_new',T21_new)
#slices('T21_old',T21_old)
#slices('xalpha',xa_new)
#slices('xalpha_old',xa_old)

#proj('T21_new',T21_new)
#proj('T21_old',T21_old)
proj('xalpha',xa_new)
#proj('xalpha_old',xa_old)
print('\n')
