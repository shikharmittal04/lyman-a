#Shikhar Mittal
#Use this code to make 3 by 2 panel plot of Lya coupling difference as in figure 6.

#norm = colors.SymLogNorm(linthresh=1.0,vmin=-1,vmax=1e4)
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import DivergingNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

L_cMpc_by_h = 10	#Box size in cMpc/h
lev = 8			#levelmin (no.of cells accross each axis in log_2 units)
N=6			#No.of panels.

data = np.zeros((N+1,2**lev,2**lev,2**lev))
proj = np.zeros((N,2**lev,2**lev))
dir = '/mnt/exports/data/mittal/outras/'

txt = ['A','B','C','D','E','F']

print('Reading files ...')
data[0,:,:,:] = np.load(dir+'10cMpc_256/xa.npy')
data[1,:,:,:] = np.load(dir+'compare/none/xa.npy')
data[2,:,:,:] = np.load(dir+'compare/bulk/xa.npy')
#data[3,:,:,:] = np.load(dir+'compare/all_but_tem/xa.npy')
data[3,:,:,:] = np.load(dir+'compare/all_but_bulk/xa.npy')
data[4,:,:,:] = np.load(dir+'compare/all_but_nhi/xa.npy')
data[5,:,:,:] = np.load(dir+'compare/all_but_rec/xa.npy')
data[6,:,:,:] = np.load(dir+'compare/all_but_dir/xa.npy')


print('\nPrinting statistics ...')
for i in range(N):
	diff = (data[i+1,:,:,:]-data[0,:,:,:])/data[0,:,:,:]
	mean_diff = np.mean(diff)
	RMS = np.sqrt(np.mean(diff**2))
	stan_devi = np.sqrt(np.sum((diff-mean_diff)**2)/(8**lev-1))
	print('\n',txt[i])
	print('Mean diff =',mean_diff)
	print('Standard deviation =',stan_devi)
	print('RMS =',RMS)
	
	proj[i,:,:] = np.mean(diff,axis=2)
	print('Max diff =',np.max(proj[i,:,:]))
	print('Min diff =',np.min(proj[i,:,:]),'\n')

'''
mp = 1 - np.max(slc) / (np.max(slc) + abs(np.min(slc)))
orig_cmap = matplotlib.cm.bwr
shifted_cmap = shiftedColorMap(orig_cmap, midpoint=mp, name='shifted')
divnorm=colors.TwoSlopeNorm(vmin=-1, vcenter=0., vmax=1e4)
'''

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fs=13


print('\nPlotting ...')
fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True,sharey=True, figsize=(9,11))
fig.subplots_adjust(left=0.04, bottom=0.05, right=0.96, top=0.99, wspace=-0.24, hspace=0.01)

#fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True,sharey=True, figsize=(10.85,9.5))
#fig.subplots_adjust(left=0.04, bottom=0.06, right=0.96, top=0.99, wspace=-0.15, hspace=0.01)

for i in range(3):
	for j in range(2):
		print('i, j, slice no.:', i,j,2*i+j)
		c=axes[i,j].imshow(proj[2*i+j,:,:].T, cmap='bwr',aspect='equal',origin='lower',extent=[-L_cMpc_by_h/2,L_cMpc_by_h/2,-L_cMpc_by_h/2,L_cMpc_by_h/2],norm=DivergingNorm(vmin=-1,vcenter=0,vmax=30))
		axes[i,j].tick_params(axis='both', which='major', length=5, width=1, labelsize=fs,direction='in')
		axes[i,j].tick_params(axis='both', which='minor', length=3, width=1,direction='in')
		axes[i,j].minorticks_on()
		axes[i,j].yaxis.set_ticks_position('both')
		axes[i,j].xaxis.set_ticks_position('both')
		axes[i,j].text(-4,4,txt[2*i+j], fontsize=1+fs,weight='bold')

divider = make_axes_locatable(axes[0,1])
cax = divider.append_axes("right", size="5%", pad=0.05)
#cax,kw = mpl.colorbar.make_axes([axis for axis in ax.flat], size="5%", pad=0.05)
cb = plt.colorbar(c,cax=cax)#, extend='max')
cb.set_label(r'$\Delta x_{\alpha}$',fontsize=fs)
cb.ax.tick_params(labelsize=fs)
#cb.set_ticks([-30,-20,-10,0,10,20, 30])

divider = make_axes_locatable(axes[1,1])
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(c,cax=cax)
cb.set_label(r'$\Delta x_{\alpha}$',fontsize=fs)
cb.ax.tick_params(labelsize=fs)
#cb.set_ticks([-30,-20,-10,0,10,20])


divider = make_axes_locatable(axes[2,1])
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(c,cax=cax)
cb.set_label(r'$\Delta x_{\alpha}$',fontsize=fs)
cb.ax.tick_params(labelsize=fs)
#cb.set_ticks([-30,-20,-10,0,10,20])


axes[2,0].set_xticks([-4, -2, 0, 2, 4])
axes[2,1].set_xticks([-4, -2, 0, 2, 4])

axes[2,0].set_xlabel("$x\,$(cMpc$h^{-1}$)",fontsize=fs)
axes[2,1].set_xlabel("$x\,$(cMpc$h^{-1}$)",fontsize=fs)
axes[0,0].set_ylabel("$y\,$(cMpc$h^{-1}$)",fontsize=fs)
axes[1,0].set_ylabel("$y\,$(cMpc$h^{-1}$)",fontsize=fs)
axes[2,0].set_ylabel("$y\,$(cMpc$h^{-1}$)",fontsize=fs)

plt.savefig(dir+'compare/compare_xa.pdf')
print('Figure saved as'+dir+'compare/compare_xa.pdf')
