#Shikhar Mittal
#Use this code to make 2 by 2 panel plot as in figure 3 and 5 from the paper.

#norm=colors.LogNorm(vmin=0,vmax=50)
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Change the folllowing settings for different box size.
#dir, L_cMpc_by_h, lev, set_{x/y}ticks, vmin, vmax

L_cMpc_by_h = 10	#Box size in cMpc/h

flnm= '10cMpc_512'
dir_new= '/mnt/exports/data/mittal/outras/'+flnm


print('Reading files...')
xa_mittal = np.load(dir_new+'/xa.npy')
xa_old = np.load(dir_new+'/xa_old.npy')
T21_mittal = np.load(dir_new+'/T21.npy')
T21_old = np.load(dir_new+'/T21_old.npy')

print('Taking projected mean ...')
proj_xa_mittal = np.mean(xa_mittal,axis=2)
proj_xa_old = np.mean(xa_old,axis=2)
proj_T21_mittal = np.mean(T21_mittal,axis=2)
proj_T21_old = np.mean(T21_old,axis=2)

print('proj_xa_mittal min, max = ',np.min(proj_xa_mittal),np.max(proj_xa_mittal))
print('proj_xa_old min, max = ',np.min(proj_xa_old),np.max(proj_xa_old))
print('proj_T21_mittal min, max = ',np.min(proj_T21_mittal),np.max(proj_T21_mittal))
print('proj_T21_old min, max = ',np.min(proj_T21_old),np.max(proj_T21_old))

print('Done.\n')
fs=13

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True,sharey=True, figsize=(10.63,9.5))
fig.subplots_adjust(left=0.04, bottom=0.06, right=0.96, top=0.97, wspace=-0.15, hspace=0.01)

c00=axes[0,0].imshow(proj_xa_old.T, cmap='cool',aspect='equal',origin='lower',extent=[-L_cMpc_by_h/2,L_cMpc_by_h/2,-L_cMpc_by_h/2,L_cMpc_by_h/2], vmin=4,vmax=40)
axes[0,0].set_title('No multiple scatterings',fontsize=fs)

c01=axes[0,1].imshow(proj_xa_mittal.T, cmap='cool',aspect='equal',origin='lower',extent=[-L_cMpc_by_h/2,L_cMpc_by_h/2,-L_cMpc_by_h/2,L_cMpc_by_h/2], vmin=4,vmax=40)
axes[0,1].set_title('With multiple scatterings',fontsize=fs)

divider = make_axes_locatable(axes[0,1])
cax = divider.append_axes("right", size="5%", pad=0.05)
cb=plt.colorbar(c01,cax=cax)
cb.set_label(r'$x_{\alpha}$',fontsize=fs)
cb.ax.tick_params(labelsize=fs)
#cb.set_ticks([1e-2,0.1,1,10,100])

c10=axes[1,0].imshow(1e-3*proj_T21_old.T, cmap='brg',aspect='equal',origin='lower',extent=[-L_cMpc_by_h/2,L_cMpc_by_h/2,-L_cMpc_by_h/2,L_cMpc_by_h/2], vmin=-0.35, vmax=0.0)

c11=axes[1,1].imshow(1e-3*proj_T21_mittal.T, cmap='brg',aspect='equal',origin='lower',extent=[-L_cMpc_by_h/2,L_cMpc_by_h/2,-L_cMpc_by_h/2,L_cMpc_by_h/2], vmin=-0.35, vmax=0.0)

divider = make_axes_locatable(axes[1,1])
cax = divider.append_axes("right", size="5%", pad=0.05)
cb=plt.colorbar(c11,cax=cax)
cb.set_label(r'$T_{21}\,(\mathrm{K})$',fontsize=fs)
cb.ax.tick_params(labelsize=fs)


for i in range(2):
	for j in range(2):
		axes[i,j].tick_params(axis='both', which='major', length=5, width=1, labelsize=fs,direction='in')
		axes[i,j].tick_params(axis='both', which='minor', length=3, width=1,direction='in')
		axes[i,j].minorticks_on()
		axes[i,j].yaxis.set_ticks_position('both')
		axes[i,j].xaxis.set_ticks_position('both')
		axes[i,j].set_aspect(1.0/axes[i,j].get_data_ratio(), adjustable='box')

axes[1,0].set_xlabel("$x\,$(cMpc$h^{-1}$)",fontsize=fs)
axes[1,1].set_xlabel("$x\,$(cMpc$h^{-1}$)",fontsize=fs)


axes[1,1].set_xticks([-4, -2, 0, 2, 4])
axes[1,0].set_xticks([-4, -2, 0, 2, 4])
axes[0,0].set_yticks([-4, -2, 0, 2, 4])
axes[1,0].set_yticks([-4, -2, 0, 2, 4])

'''
axes[1,1].set_xticks([-8, -4, 0, 4, 8])
axes[1,0].set_xticks([-8, -4, 0, 4, 8])
axes[0,0].set_yticks([-8, -4, 0, 4, 8])
axes[1,0].set_yticks([-8, -4, 0, 4, 8])
'''
axes[0,0].set_ylabel("$y\,$(cMpc$h^{-1}$)",fontsize=fs)
axes[1,0].set_ylabel("$y\,$(cMpc$h^{-1}$)",fontsize=fs)

plt.savefig(dir_new+'/'+flnm+'.pdf')
print('Figure saved as '+dir_new+'/'+flnm+'.pdf')

