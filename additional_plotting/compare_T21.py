#Shikhar Mittal
#Use this code to make 3 by 2 panel plot of 21-cm signal difference as in figure 7.

#norm = colors.SymLogNorm(linthresh=1.0,vmin=-1,vmax=1e4)
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

L_cMpc_by_h = 10	#Box size in cMpc/h
lev = 8			#levelmin (no.of cells accross each axis in log_2 units)
N=6			#No.of panels.

data = np.zeros((N+1,2**lev,2**lev,2**lev))
proj = np.zeros((N,2**lev,2**lev))
dir = '/mnt/exports/data/mittal/outras/'

print('Reading files ...')
data[0,:,:,:] = np.load(dir+'10cMpc_256/T21.npy')
data[1,:,:,:] = np.load(dir+'compare/none/T21.npy')
data[2,:,:,:] = np.load(dir+'compare/bulk/T21.npy')
#data[3,:,:,:] = np.load(dir+'compare/all_but_tem/T21.npy')
data[3,:,:,:] = np.load(dir+'compare/all_but_bulk/T21.npy')
data[4,:,:,:] = np.load(dir+'compare/all_but_nhi/T21.npy')
data[5,:,:,:] = np.load(dir+'compare/all_but_rec/T21.npy')
data[6,:,:,:] = np.load(dir+'compare/all_but_dir/T21.npy')

data = 1e-3*np.abs(data)
txt = ['A','B','C','D','E','F']
print('\nPrinting statistics ...')
for i in range(N):
	diff = data[i+1,:,:,:]-data[0,:,:,:]
	print('\n',txt[i])
	mean_diff = np.mean(diff)
	RMS = np.sqrt(np.mean(diff**2))
	stan_devi = np.sqrt(np.sum((diff-mean_diff)**2)/(8**lev-1))
	print('Mean diff =',mean_diff)
	print('Standard deviation =',stan_devi)
	print('RMS =',RMS)
	
	proj[i,:,:] = np.mean(diff,axis=2)
	print('Max diff =',np.max(proj[i,:,:]))
	print('Min diff =',np.min(proj[i,:,:]),'\n')


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
		c=axes[i,j].imshow(proj[2*i+j,:,:].T, cmap='bwr',aspect='equal',origin='lower',extent=[-L_cMpc_by_h/2,L_cMpc_by_h/2,-L_cMpc_by_h/2,L_cMpc_by_h/2],vmin=-0.2,vmax=0.2)
		axes[i,j].tick_params(axis='both', which='major', length=5, width=1, labelsize=fs,direction='in')
		axes[i,j].tick_params(axis='both', which='minor', length=3, width=1,direction='in')
		axes[i,j].minorticks_on()
		axes[i,j].yaxis.set_ticks_position('both')
		axes[i,j].xaxis.set_ticks_position('both')
		axes[i,j].text(-4,4,txt[2*i+j], fontsize=1+fs,weight='bold')

divider = make_axes_locatable(axes[0,1])
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(c,cax=cax)#, extend='max')
cb.set_label(r'$\Delta |T_{21}|\,$(K)',fontsize=fs)
cb.ax.tick_params(labelsize=fs)
cb.set_ticks([-0.2,-0.15,-0.10,-0.05,0,0.05,0.1,0.15,0.2])

divider = make_axes_locatable(axes[1,1])
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(c,cax=cax)
cb.set_label(r'$\Delta |T_{21}|\,$(K)',fontsize=fs)
cb.ax.tick_params(labelsize=fs)
cb.set_ticks([-0.2,-0.15,-0.10,-0.05,0,0.05,0.1,0.15])


divider = make_axes_locatable(axes[2,1])
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(c,cax=cax)
cb.set_label(r'$\Delta |T_{21}|\,$(K)',fontsize=fs)
cb.ax.tick_params(labelsize=fs)
cb.set_ticks([-0.2,-0.15,-0.10,-0.05,0,0.05,0.1,0.15])


axes[2,0].set_xticks([-4, -2, 0, 2, 4])
axes[2,1].set_xticks([-4, -2, 0, 2, 4])

axes[2,0].set_xlabel("$x\,$(cMpc$h^{-1}$)",fontsize=fs)
axes[2,1].set_xlabel("$x\,$(cMpc$h^{-1}$)",fontsize=fs)
axes[0,0].set_ylabel("$y\,$(cMpc$h^{-1}$)",fontsize=fs)
axes[1,0].set_ylabel("$y\,$(cMpc$h^{-1}$)",fontsize=fs)
axes[2,0].set_ylabel("$y\,$(cMpc$h^{-1}$)",fontsize=fs)

plt.savefig(dir+'compare/compare_T21.pdf')
print('Figure saved as'+dir+'compare/compare_T21.pdf')
