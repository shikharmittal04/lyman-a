#Shikhar Mittal
#Use this code to make individual difference maps such as figure 4 and B1 from the paper.
#You can also use this to make a distribution (a histogram) of the difference.

#norm = colors.SymLogNorm(linthresh=1.0,vmin=-1,vmax=1e4)
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

L_cMpc_by_h = 10
lev = 8

otpt_dir='/mnt/exports/data/mittal/outras/'

print('Reading files ...\n')
data0 = np.load(otpt_dir+'10cMpc_256/T21.npy')
#data1 = np.load(otpt_dir+'highmc/xa.npy')
#data1 = np.load(otpt_dir+'10cMpc_256/bynumber/xa.npy')
data1 = np.load(otpt_dir+'compare/none/T21.npy')
#data1 = np.load(otpt_dir+'previous_work/bulk/xa_bulk.npy')
#data1 = np.load(otpt_dir+'previous_work/bulk_nhi/xa_bulk_nhi.npy')
#data1 = np.load(otpt_dir+'previous_work/bulk_nhi_T/xa_bulk_nhi_T.npy')
#data1 = np.load(otpt_dir+'previous_work/bulk_nhi_T_dir/xa_bulk_nhi_T_dir.npy')
#data1 = np.load(otpt_dir+'previous_work/bulk_nhi_T_rec/xa_bulk_nhi_T_rec.npy')
#data1 = np.load(otpt_dir+'previous_work/bulk_T_dir_rec/xa.npy')
#data1 = np.load(otpt_dir+'previous_work/nhi_T_dir_rec/xa.npy')
#data1 = np.load(otpt_dir+'previous_work/bulk_nhi_dir_rec/xa.npy')
#data1 = np.load(otpt_dir+'previous_work/T/xa.npy')
#data1 = np.load(otpt_dir+'previous_work/Lorentz/xa.npy')

diff = np.abs(data1)-np.abs(data0)#)/data1

mean_diff = np.mean(diff)
RMS = np.sqrt(np.mean(diff**2))
stan_devi = np.sqrt(np.sum((diff-mean_diff)**2)/(8**lev-1))
print('Mean differance =',mean_diff)
print('Standard deviation =',stan_devi)
print('RMS =',RMS)

print('Max overall diff =',np.max(diff))
print('Min overall diff =',np.min(diff),'\n')

proj_diff = np.mean(diff,axis=2)
print('Max proj diff =',np.max(proj_diff))
print('Min proj diff =',np.min(proj_diff),'\n')

#slc = diff[:,:,2**(lev-1)]
#print('Max slice diff =',np.max(slc))
#print('Min slice diff =',np.min(slc),'\n')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


print('Plotting ...')
fig, axes = plt.subplots(figsize=(8.5,7.2),dpi=200)
fig.subplots_adjust(left=0.11, bottom=0.05, right=0.89, top=1.0)

c=axes.imshow((proj_diff.T), cmap='bwr',aspect='equal',origin='lower',extent=[-L_cMpc_by_h/2,L_cMpc_by_h/2,-L_cMpc_by_h/2,L_cMpc_by_h/2],vmin=-70,vmax=70)#norm=colors.LogNorm(vmin=1e-4)
#c=axes.imshow(proj_diff.T, cmap='bwr',aspect='equal',origin='lower',extent=[-L_cMpc_by_h/2,L_cMpc_by_h/2,-L_cMpc_by_h/2,L_cMpc_by_h/2],vmin=-30,vmax=30)
axes.tick_params(axis='both', which='major', length=5, width=1, labelsize=20,direction='in')
axes.tick_params(axis='both', which='minor', length=3, width=1,direction='in')
axes.minorticks_on()
axes.yaxis.set_ticks_position('both')
axes.xaxis.set_ticks_position('both')
#axes.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

divider = make_axes_locatable(axes)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(c,cax=cax)#, extend='min')
#cb.set_label(r'$10^3\times\Delta x_{\alpha}/x_{\alpha}$',fontsize=20)
#cb.set_label(r'$\Delta x_{\alpha}$',fontsize=20)
cb.set_label(r'$\Delta |T_{21}|\,$(mK)',fontsize=20)
cb.ax.tick_params(labelsize=20)

axes.set_xlabel("$x\,$(cMpc$h^{-1}$)",fontsize=20)
axes.set_xlabel("$x\,$(cMpc$h^{-1}$)",fontsize=20)
axes.set_ylabel("$y\,$(cMpc$h^{-1}$)",fontsize=20)
axes.set_ylabel("$y\,$(cMpc$h^{-1}$)",fontsize=20)

plt.savefig(otpt_dir+'compare/T21_reis.pdf')
#plt.savefig(otpt_dir+'highmc/1e8vs2e8.pdf')
print('Done.')

#----------------------------------------------------------------------------------------------
#For making difference distribution

'''
arr = rel_diff.flatten()
axes.hist(arr,bins=300,density=True,color='blue',histtype='step',log=True)
axes.set_xlabel(r'$\Delta x_{\alpha}/x_{\alpha}$', fontsize=20)
axes.tick_params(axis='both', which='major', length=5, width=1, labelsize=20,direction='in')
axes.tick_params(axis='both', which='minor', length=3, width=1,direction='in')
axes.minorticks_on()
axes.yaxis.set_ticks_position('both')
axes.xaxis.set_ticks_position('both')
axes.set_aspect(1.0/axes.get_data_ratio(), adjustable='box')
plt.savefig('rel_err_dist.png')
'''
