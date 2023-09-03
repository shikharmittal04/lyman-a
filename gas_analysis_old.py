#This uses mostly yt functionality to plot. Doesn't give much flixibility in plotting. Use instead gas_analysis.py
#Example command: python3 gas_analysis_old.py -f /mnt/exports/data/mittal/ramses/output_00001/ -pq 'density' -typ 'slc'

#import mpi4py
import yt
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from yt.extensions.astro_analysis.halo_analysis import HaloCatalog
#yt.enable_parallelism()

parser = argparse.ArgumentParser()
parser.add_argument("-f","--fln", help="enter path of text file")
parser.add_argument("-pq","--phyquan", help="'density'/'pressure'/'temperature'/'velocity_magnitude'",default='density')
parser.add_argument("-typ","--plotype", help="Slice ('slc') or Projection ('proj')",default='slc')
parser.add_argument("-ah","--haloes",help="Add haloes '1' or '0'",default=0)
parser.add_argument("-dhfe","--halcre", help="Does halo file exist?",default=0)
parser.add_argument("-str","--stars", help="Add stars?",default=0)

args = parser.parse_args()
info = args.fln
fld=args.phyquan
typ=args.plotype
add_haloes=int(args.haloes)
hcc=int(args.halcre)
st=int(args.stars)

if typ=='slc':
	if fld=='density':
		unit='Msun/Mpccm**3'
	elif fld=='pressure':
		unit='Pa'
	elif fld=='temperature':
		unit='K'
	elif fld=='velocity_magnitude':
		unit='km/s'
	else:
		print('No such field set, yet!')
		sys.exit()
elif typ=='proj':
	if fld=='density':
		unit='Msun/Mpccm**2'
	elif fld=='pressure':
		unit='Pa*m'
	elif fld=='temperature':
		unit='K*Mpccm'
	elif fld=='velocity_magnitude':
		unit='Mpccm*km/s'
	else:
		print('Invalid field!')
		sys.exit()
else:
	print('Invalid plot type!')
	sys.exit()	


h=0.674

ds = yt.load(info);
ad = ds.all_data();

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

L=ds.length_unit.in_units('Mpccm').d*h #length in cMpc/h

hp=info[:-13]+'info_000'+info[-3:-1]+'/info_000'+info[-3:-1]+'.0.h5'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

if typ=='slc':
	slc = yt.SlicePlot(ds, "z", ("gas", fld),axes_unit='Mpccm')
	if fld=='velocity_magnitude':
		slc.annotate_velocity(factor=16)
	elif fld=='density':
		slc.set_colorbar_label(label=r'Gas density $\left(\mathrm{M}_{\odot}\,\mathrm{cMpc}^{-3}\right)$',field=("gas", "density"))
	slc.set_unit(("gas", fld), unit)	
elif typ=='proj':
	box = ds.box(left_edge=[0, 0, 0], right_edge=[1, 1, 1])
	slc = yt.ProjectionPlot(ds, "z", ("gas", fld),axis_unit= 'Mpccm', data_source=box)
	print('\nShikhar: Full depth is being projected. If you want partial depth go back to the code.\n')
	if fld=='velocity_magnitude':
		slc.annotate_velocity(factor=16)
	elif fld=='density':
		slc.set_colorbar_label(label=r'Gas column density $\left(\mathrm{M}_{\odot}\,\mathrm{cMpc}^{-2}\right)$',field=("gas", "density"))
	elif fld=='temperature':
		slc.set_colorbar_label(label=r'Projected temperature $\left(\mathrm{K}\cdot\mathrm{cMpc}\right)$',field=("gas", "temperature"))
		if add_haloes==1:
			print('\nShikhar: You have chosen to add haloes.\n')
			if hcc==0:
				print('\nShikhar: Now finding haloes.\n')
				hc = HaloCatalog(data_ds=ds, finder_method="hop",output_dir=info[:-13], finder_kwargs={"threshold": 160, "ptype":"all"})
				hc.create()
				dsh = yt.load(hp)
			elif hcc==1:
				print('\nShikhar: Since halo catalogue already exists we just load it here.\n')
				dsh = yt.load(hp)	
			else:
				print('Invalid option. Bye!')
				sys.exit()
			slc.annotate_halos(dsh)
			if st==1:
				print('Shikhar: You have chosen to mark the stars.\n')
				slc.annotate_particles((L, "Mpccm"),p_size=50.0, col='r', marker='o', ptype='star')
	slc.set_unit(("gas", fld), unit)
else:
	print('Invalid plot type!')
	sys.exit()

		
#weight_field=("gas", "density")
#slc.set_colorbar_label(label='Pa$\cdot$m',field=("gas", field))
#slc.set_log(("gas", field), False)
slc.annotate_timestamp(corner="upper_left", redshift=True, draw_inset_box=True, time_format='$t = {time:.1f}\,${units}',redshift_format='$z = {redshift:.2f}$')	#To add the timestamp
#slc.annotate_scale(corner="upper_right", unit='Mpccm',scale_text_format='${scale}\,${units}')
slc.set_cmap(field=("gas", fld), cmap="viridis")
slc.set_xlabel("$x\,$(cMpc)")
slc.set_ylabel('$y\,$(cMpc)')
slc.save(name="hydro-cool3/gas_000"+info[-3:-1],suffix='pdf')

