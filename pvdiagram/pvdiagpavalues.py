# PYTHON 2.7 VERSION!
#TO DO:
#enhance subplot. Add angle values, axis (maybe), lines to orient. Vary "amp".
#do subplots for mass.
#PACKAGES:
import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
#plt.ion()
#plt.ioff()
from astropy.io import fits 
import numpy as np  
import aplpy ##JP## I get a WARNING here
from spectral_cube import SpectralCube 
from matplotlib import lines
 
from astropy import wcs
from astropy import coordinates
from astropy import units as u
from astropy.coordinates import SkyCoord
 
import os
import sys
 
from pvextractor.gui import PVSlicer 
from pvextractor import extract_pv_slice
from pvextractor.geometry import Path
import matplotlib.ticker as ticker

#PREVIOUS INFO SETUP:
writechan = False ## CHOOSE WHETHER TO WRITE IMAGE TO A FILE (True), OR DISPLAY DIRECTLY (False). #MOD
source=45 #MOD
vc=8.0*u.km/u.s #central velocity #MOD
#File:
filename='/home/jp/Research/Data/gasC18O_FullMapNew.corrected.fits' #MOD
outname='pv_testS'+str(source)+'.eps'
## Set up your map
cube0 = SpectralCube.read(filename)
ch0='' #MOD
chf='' #MOD
if ch0 != '' and chf !='': 
  cube = cube0.subcube(zlo=ch0,zhi=chf)
else: 
  cube=cube0
  ch0 = 0
  chf = cube.header['NAXIS3']

#Create figure
fig3 = plt.figure(figsize=(9,6)) #MOD
fig3.clf() 

#SOURCE COORDS:
filenametxt='/home/jp/Research/Data/sourcecoordinates.txt'
ra,dec = np.loadtxt(filenametxt,unpack=True,usecols=(0,1),dtype='S')
coord=SkyCoord(ra,dec,frame='icrs',unit=(u.hourangle,u.deg)) 
c = SkyCoord(coord[source-1].ra,coord[source-1].dec, frame='icrs',unit=(u.hourangle,u.deg)) 


### OPTIONAL:
def standard_setup(sp): 
  sp.set_frame_color('black')
  sp.set_tick_labels_font(size='14')
  sp.set_axis_labels_font(size='16')
  sp.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss') 
  sp.set_tick_color('black')
  #sp.set_system_latex(True) #use latex math symbols in plots
def pv_setup(sp):
  sp.set_frame_color('black')
  sp.set_tick_labels_font(size='14')
  sp.set_axis_labels_font(size='16')
  #sp.set_tick_labels_format(xformat='mm:ss',yformat='dd:mm:ss') 
  sp.set_tick_color('black')
  #sp.set_system_latex(True)


#LOCATING SOURCE:
RAcoord=c.ra
Deccoord=c.dec
xcoord=[0,cube.shape[2]] 
ycoord=[0,cube.shape[1]]
x_w,y_w=[RAcoord.value],[Deccoord.value] 		
  

dchan=1 #0.5 to shorten interval
amp=5

filenametxt='/home/jp/Research/Data/sourcecoordinates.txt'
ra,dec = np.loadtxt(filenametxt,unpack=True,usecols=(0,1),dtype='S')
coord=SkyCoord(ra,dec,frame='icrs',unit=(u.hourangle,u.deg)) 
c = SkyCoord(coord[source-1].ra,coord[source-1].dec, frame='icrs',unit=(u.hourangle,u.deg)) 
#SOURCE PV INFO: THERE IS PROBABLY A MUCH SHORTER AND GENERAL WAY TO DO THIS
file2dfit='/home/jp/Research/Data/2DFit/'+str(source)+'/2DGssfit'+str(source)+'Log.txt'
file1=open(file2dfit,'r')
i=1
for line in file1:
	if i==27:
		line1=line
	if i==29:
		line2=line
	i+=1
file1.close()
for j in range(len(line1)):
	if line1[j]==':':
		k=j+1
		majax1=''
		while line1[k]!='+':
			majax1+=line1[k]
			k+=1
for j in range(len(line2)):
	if line2[j]==':':
		k=j+1
		angle1=''
		while line2[k]!='+':
			angle1+=line2[k]
			k+=1
dchan=1 #0.5 to shorten interval
amp=5
majorax=amp*float(majax1.strip())

 
################################
#### PANEL (b): PV PLOT
################################
#PLOT:
panels=[[0.005,0.5,0.326,0.48],[0.336,0.5,0.326,0.48],[0.667,0.5,0.326,0.48],[0.005,0.01,0.326,0.48],[0.336,0.01,0.326,0.48],[0.667,0.01,0.326,0.48]]
angles=[0.0,30.0,60.0,float(angle1.strip()),110.0,140.0]

for i in range(6):
	angle=angles[i]
	majorax=amp*float(majax1.strip())

	#PV-CUT:
	halflength = majorax/2.*u.arcsec  #MOD
	pa = angle*u.deg #Deg: degrees away from vertical.
	avnumb=6 #number of pixels to average perpendicular to pv-cut #MOD
	inter = 1*u.arcsec #interval you want the arcsec to be marked on #MOD

# Calculate the x,y coordinates on either end of the PV cut
	upperx = RAcoord+halflength.to(u.deg)*np.sin(pa.to(u.rad))
	lowerx = RAcoord-halflength.to(u.deg)*np.sin(pa.to(u.rad))
	uppery = Deccoord+halflength.to(u.deg)*np.cos(pa.to(u.rad))
	lowery = Deccoord-halflength.to(u.deg)*np.cos(pa.to(u.rad))
	endpoints_deg = np.array([[upperx.value,uppery.value],[lowerx.value,lowery.value]]) 

	################################
	#### MAKE A PV CUBE
	################################
	pvcube = cube.subcube(xlo=xcoord[0],xhi=xcoord[1],ylo=ycoord[0],yhi=ycoord[1],zlo=ch0,zhi=chf+1) 
	hdu = pvcube
	w = wcs.WCS(hdu.header)
	endpoints_pix = w.sub([wcs.WCSSUB_CELESTIAL]).wcs_world2pix([[x,y] for x,y in endpoints_deg], 0) 
	endpoints_list=endpoints_pix.tolist() 
	path=Path(endpoints_list,width=avnumb) 
	pv = extract_pv_slice(hdu,path,wcs=w) 
	pv_vel=pvcube.spectral_axis 
	minvel = pv_vel.min().value
	maxvel = pv_vel.max().value
	pvstd = pv.data[120:160,90:120].std() #choose a low signal region #MOD
	pvstd2 = pv.data[620:680,90:120].std() #MOD
	pvrms = np.nanmean([pvstd,pvstd2])
	F4 = aplpy.FITSFigure(pv,dimensions=[1,0],figure=fig3,subplot=(panels[i])) 
	#pv_setup(F4)
	F4.show_grayscale(aspect='auto',pmin=.5,pmax=99.75,stretch='linear',invert=True)
	#F4.show_contour(colors='orange',levels=5,dimensions=[1,0]) # does not work as
	plt.axis('off')
	F4.add_label(0.15,0.9,str(pa.value)+' deg',relative=True,color='red',size='large')
	#F4.show_contour(colors='cyan',levels=10,linewidths=1.5,zorder='8') # does not work as intended
	#UNCOMENT THE TEXT BETWEEN ## AND ## TO SHOW VERTICAL AND HORIZONTAL LINES THAT HELP LOCATING STAR POSITION AND CENTRAL VELOCITY
	##
# Set up a second Axis on top of pv plot, to show annotations (if necessary) ##?## invisible axes just to make annotations
	res_elem = 0.15*u.arcsec ##?## Where do I get this?? #resolution element, or beam size #MOD
	x1_res=res_elem.to(u.deg).value 
	x2_res=x1_res+res_elem.to(u.deg).value 
	y1_res=pvcube.spectral_axis.min().value
	y2_res=pvcube.spectral_axis.max().value 
	xlim=F4._ax1.get_xlim() 
	ylim=F4._ax1.get_ylim()
	xlim_lo = minvel 
	xlim_hi = maxvel
	ylim_lo=(ylim[0])*pv.header['CDELT1']+pv.header['CRVAL1']-pv.header['CDELT1']
	ylim_hi=(ylim[1])*pv.header['CDELT1']+pv.header['CRVAL1']-pv.header['CDELT1']
	ax2 = fig3.add_axes(panels[i]) 
	ax2.set_xlim(xlim_lo,xlim_hi)
	ax2.set_ylim(ylim_lo,ylim_hi)
	 
	ax2.set_axis_off()
	 
	# VERTICAL (dashed) LINE TO SHOW v_c in PV Plot
	hx1 = 0.
	hx2 = 2.*halflength.to(u.deg).value
	hy1 = vc.to(u.m/u.s).value
	ax2.vlines(hy1,hx1,hx2,colors='green',linestyles='dashed',label='$v_c$',zorder='3',lw=2) ##?## ##JP## I am getting vertical line shifted #MOD
	ax2.annotate('$v_c$',xy=(hy1,x1_res*2.0),xycoords='data',color='green',size='16')
	 
	# HORIZONTAL (dashed) LINE TO SHOW star position in PV Plot. Separate position to left and right from the source along the position line. 
	hx3 = xlim_lo
	hx4 = xlim_hi
	hy3 = halflength.to(u.deg).value
	#### NOTE, You have to move these around a bit depending on your target. #MOD
	ax2.hlines(hy3,hx3,hx4,colors='magenta',linestyles='dashed',label='protostar',zorder='3',lw=2)
	#ax2.annotate('protostar',xy=(xlim_hi,hy3-2*res_elem.to(u.deg).value),xycoords='data',color='magenta',size='20',verticalalignment='top')
	#ax2.annotate('(b) CO, \n PV diag.',xy=(xlim_hi,hy3+8*res_elem.to(u.deg).value),xycoords='data',color='Black',size='20',verticalalignment='top')
	##
	#MODIFY AXES:
	'''F4.refresh() 
	print([x.get_text()[1:-1] for x in F4._ax1.get_yticklabels()]) 
	xlim=F4._ax1.get_xlim() 
	ylim=F4._ax1.get_ylim()
	xlim_lo = minvel 
	xlim_hi = maxvel
	ylim_lo=(ylim[0])*pv.header['CDELT1']+pv.header['CRVAL1']-pv.header['CDELT1']
	ylim_hi=(ylim[1])*pv.header['CDELT1']+pv.header['CRVAL1']-pv.header['CDELT1']
	npix = round(halflength.to(u.deg) / (cube.header['cdelt2']*u.deg))
	inter_px = inter.to(u.deg).value/(ylim_hi-ylim_lo)*(ylim[1]-ylim[0])
	ntick = (int((ylim_hi-ylim_lo)/inter.to(u.deg).value)+2) 
	print([x.get_text() for x in F4._ax1.get_yticklabels()])
	F4.refresh() 
	F4._ax1.yaxis.set_major_locator(ticker.MultipleLocator(inter_px)) 
	F4._ax1.yaxis.set_minor_locator(ticker.MultipleLocator(inter_px/inter.value)) 
	newticks = np.arange(0,inter.value*ntick,inter.value)-halflength.value 
	F4._ax1.set_yticklabels(newticks)
	F4.refresh()
	print([x.get_text()[:-2] for x in F4._ax1.get_yticklabels()])
	F4._ax1.set_xticklabels([str(float(x.get_text()[:])/1000.*dchan)[:-2] for x in F4._ax1.get_xticklabels()])
	F4._ax1.set_yticklabels([x.get_text()[:] for x in F4._ax1.get_yticklabels()])##NOTE## Here edit the indexing of x.get_text()[:-2]. #with [:-2] we achieve to do '12.0 -> '12'
	F4._ax1.set_ylabel('Offset (arcsec)') 
	F4._ax1.set_xlabel("$V_{LSR} (\mathrm{km\ s}^{-1})$")
	print([x.get_text() for x in F4._ax1.get_yticklabels()])
	F4.refresh()
	F4._ax2.yaxis.set_major_locator(ticker.MultipleLocator(inter_px)) 
	F4._ax2.yaxis.set_minor_locator(ticker.MultipleLocator(inter_px/inter.value))
	F4._ax1.yaxis.labelpad = -9 #Edit to modify spacing #MOD'''

'''
################################
## FINAL ANNOTATIONS
################################
# HORIZONTAL LINE TO SHOW resolution in PV Plot. 
res_elem = 0.15*u.arcsec ##?## Where do I get this?? #resolution element, or beam size #MOD
x1_res=res_elem.to(u.deg).value 
x2_res=x1_res+res_elem.to(u.deg).value 
y1_res=pvcube.spectral_axis.min().value
y2_res=pvcube.spectral_axis.max().value 
 

 
# Set up a second Axis on top of pv plot, to show annotations (if necessary) ##?## invisible axes just to make annotations
ax2 = fig3.add_axes(rightpanel) 
ax2.set_xlim(xlim_lo,xlim_hi)
ax2.set_ylim(ylim_lo,ylim_hi)
 
ax2.set_axis_off()
 
# VERTICAL (dashed) LINE TO SHOW v_c in PV Plot
hx1 = 0.
hx2 = 2.*halflength.to(u.deg).value
hy1 = vc.to(u.m/u.s).value
ax2.vlines(hy1,hx1,hx2,colors='green',linestyles='dashed',label='$v_c$',zorder='3',lw=2) ##?## ##JP## I am getting vertical line shifted #MOD
ax2.annotate('$v_c$',xy=(hy1,x1_res*2.0),xycoords='data',color='green',size='16')
 
# HORIZONTAL (dashed) LINE TO SHOW star position in PV Plot. Separate position to left and right from the source along the position line. 
hx3 = xlim_lo
hx4 = xlim_hi
hy3 = halflength.to(u.deg).value
#### NOTE, You have to move these around a bit depending on your target. #MOD
ax2.hlines(hy3,hx3,hx4,colors='magenta',linestyles='dashed',label='protostar',zorder='3',lw=2)
ax2.annotate('protostar',xy=(xlim_hi,hy3-2*res_elem.to(u.deg).value),xycoords='data',color='magenta',size='20',verticalalignment='top')
ax2.annotate('(b) CO, \n PV diag.',xy=(xlim_hi,hy3+8*res_elem.to(u.deg).value),xycoords='data',color='Black',size='20',verticalalignment='top')'''

##SHOW/SAVE FIGURE:
#plt.show()
if writechan: 
  plt.savefig(outname)
if not writechan:
  print 'print to screen'
  fig3.show() 
#plt.close(fig3)


