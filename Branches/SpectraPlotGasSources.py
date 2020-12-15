##Plot Single Pixel intensities over diffrent velocities (channells)
#Import
from astropy.coordinates import SkyCoord
import numpy as np
from astropy import wcs
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import units as u  

##Info for rms (sigma)
filename_gas=fits.open('/home/jp/Desktop/Research/Data2.0/Serpens/Serpens_concat_c18o_cube.fits')
#filename_gas=fits.open('/home/jp/Desktop/Research/Data2.0/Serpens/Serpens_Feather_c18o.fits')
mom0_rms=np.mean(filename_gas[0].data[:,689:719,939:977].std())
levs_3=3*mom0_rms

#Create Pixlist
#From arcsec to deg
filenametxt='/home/jp/Desktop/Research/Data/sourcecoordinates.txt'
ra,dec = np.loadtxt(filenametxt,unpack=True,usecols=(0,1),dtype='S')
coord=SkyCoord(ra,dec,frame='icrs',unit=(u.hourangle,u.deg))

#From deg to pix
filename='/home/jp/Desktop/Research/Data2.0/moment_0_New.fits'
Pixlist=[]
hdulist = fits.open(filename)
w=wcs.WCS(hdulist[0].header)
for i in range(len(coord)):
	arr=np.array([[coord.ra[i].value,coord.dec[i].value]])
	#to make this line work use moment0 map.
	xpix,ypix=w.all_world2pix(arr,0)[0]
	Pixlist.append((xpix,ypix))

#Quantities
source_num=17
a=source_num-1##--> index of the coordinate you want to use in Pixlist
chan=72
xpix=int(round(Pixlist[a][0]))
ypix=int(round(Pixlist[a][1]))

xpix2=1603
ypix2=1216

#open file
hdu=filename_gas

#extract header info
vo=hdu[0].header['CRVAL3']/1000
deltv=hdu[0].header['CDELT3']/1000

#create arrays to plot
l1array=vo+np.arange(chan)*deltv
l2array=hdu[0].data[:,ypix,xpix]
l3array=np.mean(hdu[0].data[:,ypix-2:ypix+2,xpix-2:xpix+2],axis=(1,2))
l4array=np.mean(hdu[0].data[:,ypix2-4:ypix2+4,xpix2-4:xpix2+4],axis=(1,2))
##Creating horizontal line 3*rms
l=[]
for i in range(len(l1array)):
	l.append(levs_3)
#Plot
plt.figure(ra[0]+','+dec[0])
plt.title('Source '+str(source_num))
plt.xlabel('Velocity [km/s]')
plt.ylabel('Intensity [Jy/beam]')
#plt.plot(l1array,l2array,drawstyle='steps-mid')
plt.plot(l1array,l3array,drawstyle='steps-mid')
plt.plot(l1array,l4array,drawstyle='steps-mid',color='red',label='branch')
#plt.plot(l1array,l,color='gray',marker='.')
#plt.plot(l1array,l3smoothed,color='red',label='Smoothed Spectra ('+str(smoothedchan)+' channels)')
plt.plot(l1array,l,color='gray',linestyle='--',label='3$\sigma$')
plt.legend()
#plt.plot(l1array,ls,color='green',linestyle='--',label='3$\sigma$ Smoothed')
plt.show()
