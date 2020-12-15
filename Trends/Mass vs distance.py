import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u  
from astropy.coordinates import SkyCoord

filenametxt='/home/jp/Research/Data/PeakFluxMassIntGasClass.txt'
ra,dec,contmass,sclass= np.loadtxt(filenametxt,unpack=True,usecols=(1,2,5,7),dtype='S')

lcontmass=[]
for i in range(len(contmass)):
	lcontmass.append(float(contmass[i]))
coord=SkyCoord(ra,dec,frame='icrs',unit=(u.hourangle,u.deg))
c1=SkyCoord(ra=277.5148313*u.degree, dec=-2.048732216*u.degree, distance=436*u.pc, frame='icrs')
ldist=[]
for i in range(len(coord)):
	c2 = SkyCoord(ra=coord[i].ra, dec=coord[i].dec, distance=436*u.pc, frame='icrs')
	dist=c1.separation_3d(c2)
	ldist.append(dist.value)
lcontmassN=[]
lcontmass0_I=[]
lcontmassF=[]
ldistN=[]
ldist0_I=[]
ldistF=[]
for i in range(len(contmass)):
	if sclass[i]=='N':
		lcontmassN.append(float(contmass[i]))
		ldistN.append(float(ldist[i]))
	elif sclass[i]=='F':
		lcontmassF.append(float(contmass[i]))
		ldistF.append(float(ldist[i]))
	else:
		lcontmass0_I.append(float(contmass[i]))
		ldist0_I.append(float(ldist[i]))

figname='Continuum Mass vs. distance'
fig=plt.figure('Continuum Mass vs. distance',figsize=(18,12))
plt.title('Continuum Mass vs. distance')
plt.xlabel('Distance [parsec]')
plt.ylabel('Mass [Mo]')
plt.yticks(np.arange(min(lcontmass),max(lcontmass)+0.05,0.05))
plt.xticks(np.arange(min(ldist),max(ldist)+0.025,0.025))
plt.scatter(ldistF,lcontmassF,color='blue',label='F')
plt.scatter(ldistN,lcontmassN,color='gray',label='None')
plt.scatter(ldist0_I,lcontmass0_I,color='red',label='Class 0-I')
plt.rc('axes', labelsize=18, titlesize=18)
plt.rc('xtick', labelsize=13) 
plt.rc('ytick', labelsize=13) 
plt.legend()
fig.savefig(figname+'.png')
plt.show()
plt.close(fig)
