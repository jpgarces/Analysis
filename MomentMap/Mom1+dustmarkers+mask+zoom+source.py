## Import statements
from astropy.io import fits
from astropy.coordinates import SkyCoord
import aplpy
from astropy import units as u  
from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import numpy as np 
 
############### SETUP
## Filenames
filename='/home/juan/Desktop/Research/Data2.0/Serpens/Serpens_Feather_c18o.fits'
filename_mom1 = '/home/juan/Desktop/Research/Data2.0/moment_1_masked_New.fits' 
filename_cont='/home/juan/Desktop/Research/Data2.0/dustContinuum.fits'
 
####
source=45
source2=46
filenametxt='/home/juan/Desktop/Research/Data/sourcecoordinates.txt'
ra,dec = np.loadtxt(filenametxt,unpack=True,usecols=(0,1),dtype='S')
coord=SkyCoord(ra,dec,frame='icrs',unit=(u.hourangle,u.deg)) 

fig = plt.figure('Source '+str(source)+' and '+str(source2),figsize=(15,15)) 

cube1=SpectralCube.read(filename)
cubec=SpectralCube.read(filename_cont)
moment=fits.open(filename_mom1)
moment_1=moment[0].data
#moment_1=cube1.moment(order=1)
  
### Setup some information about your source
### EDIT THE INFORMATION HERE
c = SkyCoord(coord[source-1].ra,coord[source-1].dec, frame='icrs',unit=(u.hourangle,u.deg))
RAcoord=c.ra
Deccoord=c.dec
vc=8.0*u.km/u.s 
 
############### MOMENT 1 MAP
f2 = aplpy.FITSFigure(filename_mom1,figure=fig)
mom1_vmin = ((vc.value-0.8)*1000.)
mom1_vmax = ((vc.value+0.8)*1000.)
f2.show_colorscale(vmin=mom1_vmin,vmax=mom1_vmax,cmap='jet')
 
################ EXTRAS
## you can recenter the plot to your chosen coordinate
xx=np.abs(cube1.header['NAXIS1']*cube1.header['cdelt1'])#np.abs is absolute value. 'Naxis1'*'cdelt1" will give length of x axis in 
yy=np.abs(cube1.header['NAXIS2']*cube1.header['cdelt2'])
f2.recenter(RAcoord+0.0*u.arcsec,Deccoord+0.0*u.arcsec,width=xx*0.03,height=yy*0.03)

##Show marker on each source
filenametxt='/home/juan/Desktop/Research/Data/sourcecoordinates.txt'
ra,dec = np.loadtxt(filenametxt,unpack=True,usecols=(0,1),dtype='S')
coord=SkyCoord(ra,dec,frame='icrs',unit=(u.hourangle,u.deg))
for i in range(len(ra)):
    f2.show_markers(coord[i].ra.value,coord[i].dec.value,edgecolor='red',marker='o',linewidths=2,s=30)
### Show the colorbar
f2.add_colorbar()

### Hack for colorbar labels in km/s

cb_inter=200.*u.m/u.s
labels=(np.arange(mom1_vmin,mom1_vmax+cb_inter.value,cb_inter.value))*(u.m/u.s)
f2.colorbar.set_labels(False)
#f2.colorbar.set_ticks(labels.value)
f2.colorbar.set_width(0.1)
f2.colorbar.set_pad(-0.1)
f2.colorbar.set_axis_label_pad(30)
f2.colorbar.set_axis_label_text('km s$^{-1}$')
f2.colorbar.set_axis_label_font(size='14')
f2.hide_yaxis_label()
f2.hide_ytick_labels()
xlim_mom=f2._ax1.get_xlim()
ylim_mom=f2._ax1.get_ylim()
hack_labels=np.arange(0,1.01,1./(labels.size-1))
for ll in range(labels.size): 
	labelfloat = round(labels[ll].value/1000.,2)-vc.value
	f2.add_label(1.05,hack_labels[ll],str(labelfloat),relative=True,color='black',size='14')

### Add the beamsize 
bmaj=1.043*u.arcsec.to(u.deg)
bmin=0.7046*u.arcsec.to(u.deg)
bpa=76.299
f2.add_beam(major=bmaj,minor=bmin,angle=bpa,pad=0.2)

##Save figure
#f2.save('Source '+str(source)+' to'+str(source2)+'.png')

fig.show()
