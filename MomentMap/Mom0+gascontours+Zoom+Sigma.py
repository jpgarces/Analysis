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
filename='/home/jp/Desktop/Research/Data2.0/Serpens/Serpens_concat_c18o_cube.fits'
#filename_mom0 = 'moment_0.fits' 
source=16
source2=18
filenametxt='/home/jp/Desktop/Research/Data/sourcecoordinates.txt'
ra,dec = np.loadtxt(filenametxt,unpack=True,usecols=(0,1),dtype='S')
coord=SkyCoord(ra,dec,frame='icrs',unit=(u.hourangle,u.deg)) 

cube0=SpectralCube.read(filename)
moment_0=cube0.moment(order=0)
 
fig = plt.figure('Source '+str(source)+' and '+str(source2),figsize=(15,15)) 
 
### Setup some information about your source
### EDIT THE INFORMATION HERE
c = SkyCoord(coord[source-1].ra,coord[source-1].dec, frame='icrs',unit=(u.hourangle,u.deg)) 
RAcoord=c.ra
Deccoord=c.dec
vc=8.0*u.km/u.s 
############### MOMENT 0 MAP
#f2 = aplpy.FITSFigure(filename_mom0,dimensions=[0,1],figure=fig)
f2 = aplpy.FITSFigure(moment_0.hdu, figure=fig)
mom0_vmin = moment_0[~np.isnan(moment_0)].min().value*0.6
mom0_vmax = moment_0[~np.isnan(moment_0)].max().value*0.6
f2.show_grayscale(vmin=mom0_vmin,vmax=mom0_vmax,invert=True)
 
################ EXTRAS
## you can recenter the plot to your chosen coordinate
xx=np.abs(cube0.header['NAXIS1']*cube0.header['cdelt1'])
yy=np.abs(cube0.header['NAXIS2']*cube0.header['cdelt2'])


##Show label on each source
f2.add_label(coord[source-1].ra.value,coord[source-1].dec.value,source,size=14,zorder='9',color='brown')
f2.add_label(coord[source2-1].ra.value,coord[source2-1].dec.value,source2,size=14,zorder='9',color='brown')

##Shw marker on each source
f2.show_markers(coord[source-1].ra.value,coord[source-1].dec.value,edgecolor='green',marker='o',linewidths=2,s=90,zorder='8')
f2.show_markers(coord[source2-1].ra.value,coord[source2-1].dec.value,edgecolor='green',marker='o',linewidths=2,s=90,zorder='8')

## If you want to show blue-shifted or red-shifted contours
velwhole = cube0.spectral_axis.to(u.km/u.s)
omitvel=0.4*u.km/u.s 
minimumvelocity = min(velwhole)
maximumvelocity = vc-omitvel
blucube = cube0.spectral_slab(min(velwhole),(vc-omitvel))
redcube = cube0.spectral_slab((vc+omitvel),max(velwhole))
 
print('BLUE cube v$_{LSR}=$%.1f to %.1f %s'%(min(velwhole).value,(vc-omitvel).value,velwhole.unit))
print('RED cube v$_{LSR}=$%.1f to %.1f %s'%((vc+omitvel).value,max(velwhole).value,velwhole.unit))

moment_0_blu = blucube.moment(order=0)
moment_0_red = redcube.moment(order=0)

##Info for rms (sigma)
mom0_rmsb=moment_0_blu[689:719,939:977].std().value

mom0_rmsr=moment_0_red[2259:2345,1886:1976].std().value

f2.show_contour(moment_0_blu.hdu,colors='blue',levels=np.linspace(3.,18.,12)*mom0_rmsb,linewidths=2.0)

f2.show_contour(moment_0_red.hdu,colors='red',levels=np.linspace(3.,18.,12)*mom0_rmsr,linewidths=2.0)

f2.recenter(RAcoord+0.0*u.arcsec,Deccoord+0.0*u.arcsec,width=xx*0.042,height=yy*0.032)
### Add the beamsize 
### Add the beamsize 
bmaj=1.043*u.arcsec.to(u.deg)
bmin=0.7046*u.arcsec.to(u.deg)
bpa=76.299
f2.add_beam(major=bmaj,minor=bmin,angle=bpa,pad=0.2)

#f2.save('Source '+str(source)+' to'+str(source2)+'.png')
fig.show()

