#Use this code to make DM projection or slice plots. This can be run on both hydro or DMO simulation.
#Eg. python3 dm_analysis.py -f /home/shikhar/research/P5/dmo2/output_00003/ -typ 'proj' -m True -ah True
import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.extensions.astro_analysis.halo_analysis import HaloCatalog
import sys
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-f","--file", help="enter path of text file")
parser.add_argument("-typ","--plotype", help="Slice ('slc') or Projection ('proj')",default='slc')


args = parser.parse_args()
info = args.file
typ=args.plotype

add_haloes=input('Add haloes on the plot? (1 or 0): ')

ds = yt.load(info)
fn=info[:-13]+'dm_000'+info[-3:-1]+'.h5'


h=0.674
z=ds.current_redshift
l=ds.length_unit

Mpchinv=l.d/(3.0857e24)*(1+z)*h		#Box size in cMpc.h^-1


plt.rc('text', usetex=True)
plt.rc('font', family='serif')


L=2**6
cg = ds.arbitrary_grid(left_edge=ds.domain_left_edge, right_edge=ds.domain_right_edge, dims=[L,L,L])
cg.save_as_dataset(fields=[("deposit", "all_density")],filename=fn)
ds_grid = yt.load(fn)

hp=info[:-13]+'info_000'+info[-3:-1]+'/info_000'+info[-3:-1]+'.0.h5'

if typ=='slc':
	p = yt.SlicePlot(ds_grid, 'z', ("grid", "all_density"),width=(Mpchinv/h,"Mpccm"))
	p.set_unit(("grid", "all_density"), "Msun/Mpccm**3")
	p.set_colorbar_label(label=r'DM density $\left(\mathrm{M}_{\odot}\,\mathrm{cMpc}^{-3}\right)$',field=("grid", "all_density"))
elif typ=='proj':
	p = yt.ProjectionPlot(ds_grid, 'z', ("grid", "all_density"),width=(Mpchinv/h,"Mpccm"))
	p.set_unit(("grid", "all_density"), "Msun/Mpccm**2")
	p.set_colorbar_label(label=r'DM column density $\left(\mathrm{M}_{\odot}\,\mathrm{cMpc}^{-2}\right)$',field=("grid", "all_density"))
	if add_haloes=='1':
		hcc=input('Halo catalogue created? (1 or 0): ')
		if hcc=='0':
			hc = HaloCatalog(data_ds=ds, finder_method="hop",output_dir=info[:-13],finder_kwargs={"threshold": 178, "ptype":"all"})
			hc.create()
			dsh = yt.load(hp)
		elif hcc=='1':
			dsh = yt.load(hp)	
		else:
			print('Invalid option. Bye!')
			sys.exit()
		p.annotate_halos(dsh)
else:
	print('Invalid plot type!')
	sys.exit()

p.annotate_timestamp(corner="upper_left", redshift=True, draw_inset_box=True, time_format='$t = {time:.1f}\,${units}',redshift_format='$z = {redshift:.2f}$')

p.set_xlabel("$x\,$(cMpc)")
p.set_ylabel("$y\,$(cMpc)")
p.set_cmap('all', 'viridis')
p.save(name='dmo'+info[-15]+"/dm_000"+info[-3:-1],suffix='pdf')
